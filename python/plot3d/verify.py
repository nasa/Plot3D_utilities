"""Verification of connectivity and periodicity using permutation matrices.

The 8 pre-computed permutation matrices encode every possible face orientation.
Instead of re-extracting points in different traversal orders, we:

1. Extract both faces as canonical 2D grids (ascending index order).
2. Apply ``PERMUTATION_MATRICES[perm_idx]`` to face B's grid.
3. Compare point-by-point within tolerance.

Bit encoding: ``perm_idx = u_reversed | (v_reversed << 1) | (swapped << 2)``

- 0-3 (in-plane): same constant axis, direction flips only.
- 4-7 (cross-plane): different constant axes, loop order changes.

JSON Export Convention
---------------------
When exporting face matches with lb/ub diagonal bounds:

- **In-plane** (perm 0-3): ``permutation_index`` is set to **-1** because the
  traversal direction is fully encoded in block2's lb/ub ordering.
- **Cross-plane** (perm 4-7): ``permutation_index`` contains the actual index
  (4-7) because lb/ub alone cannot represent an axis swap.

The ``permutation_matrix`` field always contains the actual 2x2 matrix regardless
of whether ``permutation_index`` is -1.
"""

from .block import Block
from .blockfunctions import reduce_blocks, rotate_block, compute_min_gcd, scale_face_bounds, constant_axis
from .connectivity import PERMUTATION_MATRICES
from .periodicity import create_rotation_matrix
from typing import List, Tuple
from copy import deepcopy
from math import radians
import numpy as np


# ---------------------------------------------------------------------------
# Core helpers: extract, permute, compare
# ---------------------------------------------------------------------------

def get_bounds(face: dict) -> Tuple[tuple, tuple]:
    """Extract ascending (lo, hi) bounds from either lo/hi or lb/ub format.

    Handles both canonical (lo/hi) and diagonal (lb/ub) JSON formats
    transparently. Always returns ascending bounds.

    Args:
        face: dict with either 'lo'/'hi' or 'lb'/'ub' keys.

    Returns:
        (lo, hi) tuples with lo[d] <= hi[d] for all d.
    """
    if 'lo' in face:
        return tuple(face['lo']), tuple(face['hi'])
    elif 'lb' in face:
        c1, c2 = face['lb'], face['ub']
        return (tuple(min(a, b) for a, b in zip(c1, c2)),
                tuple(max(a, b) for a, b in zip(c1, c2)))
    else:
        raise KeyError("Face dict must have 'lo'/'hi' or 'lb'/'ub' keys")


def extract_canonical_grid(block: Block, lb: list, ub: list) -> Tuple[np.ndarray, int, int]:
    """Extract face as a canonical 2D grid (nu, nv, 3) in ascending index order.

    Finds the constant axis, then extracts points with the first varying axis
    as the outer loop (u) and the second as the inner loop (v), both ascending.

    Args:
        block: Block to extract from.
        lb: Lower diagonal corner [i, j, k].
        ub: Upper diagonal corner [i, j, k].

    Returns:
        (grid, nu, nv) where grid has shape (nu, nv, 3).

    Raises:
        ValueError: If no constant axis is found.
    """
    lo = [min(lb[d], ub[d]) for d in range(3)]
    hi = [max(lb[d], ub[d]) for d in range(3)]

    const_dim = constant_axis(lo, hi)
    if const_dim < 0:
        raise ValueError(f"No constant axis found: lo={lo}, hi={hi}")

    vary = [d for d in range(3) if d != const_dim]
    d0, d1 = vary
    nu = hi[d0] - lo[d0] + 1
    nv = hi[d1] - lo[d1] + 1

    grid = np.empty((nu, nv, 3))
    idx = [0, 0, 0]
    idx[const_dim] = lo[const_dim]
    for u in range(nu):
        idx[d0] = lo[d0] + u
        for v in range(nv):
            idx[d1] = lo[d1] + v
            grid[u, v] = [block.X[idx[0], idx[1], idx[2]],
                          block.Y[idx[0], idx[1], idx[2]],
                          block.Z[idx[0], idx[1], idx[2]]]
    return grid, nu, nv


def apply_permutation(grid: np.ndarray, perm_idx: int) -> np.ndarray:
    """Apply a pre-computed permutation matrix to a 2D face grid.

    Uses bit operations on ``perm_idx`` (0-7) to flip and/or transpose the grid.
    The permutation matrix is looked up from ``PERMUTATION_MATRICES``, not recalculated.

    Bit encoding: ``perm_idx = u_reversed | (v_reversed << 1) | (swapped << 2)``

    Args:
        grid: Face grid with shape (nu, nv, 3).
        perm_idx: Permutation index 0-7.

    Returns:
        Permuted grid with shape (out_nu, out_nv, 3).
    """
    g = grid
    if perm_idx & 1:
        g = g[::-1, :, :]    # flip u
    if perm_idx & 2:
        g = g[:, ::-1, :]    # flip v
    if perm_idx & 4:
        g = g.transpose(1, 0, 2)  # swap u, v
    return np.ascontiguousarray(g)


def verify_match(pts_a: np.ndarray, pts_b: np.ndarray, tol: float) -> bool:
    """Compare two point arrays within tolerance.

    Args:
        pts_a: Array of shape (N, 3) or (nu, nv, 3).
        pts_b: Array of same shape as pts_a.
        tol: Euclidean distance tolerance.

    Returns:
        True if all corresponding points are within tolerance.
    """
    if pts_a.shape != pts_b.shape:
        return False
    diffs = np.linalg.norm(pts_a - pts_b, axis=-1)
    return float(diffs.max()) < tol


def try_all_permutations(grid_a: np.ndarray, grid_b: np.ndarray,
                         tol: float) -> int:
    """Try all 8 permutation matrices on grid_b to find one that matches grid_a.

    For each permutation index 0..8:
    1. Apply the permutation to grid_b via :func:`apply_permutation`.
    2. Check output shape matches grid_a's shape.
    3. Compare point-by-point via :func:`verify_match`.

    Args:
        grid_a: Canonical grid of face A, shape (nu_a, nv_a, 3).
        grid_b: Canonical grid of face B, shape (nu_b, nv_b, 3).
        tol: Euclidean distance tolerance.

    Returns:
        Permutation index (0-7) on success, -1 if no permutation works.
    """
    nu_a, nv_a = grid_a.shape[:2]
    for perm_idx in range(8):
        g = apply_permutation(grid_b, perm_idx)
        if g.shape[0] != nu_a or g.shape[1] != nv_a:
            continue
        if verify_match(grid_a, g, tol):
            return perm_idx
    return -1


def verify_partial_match(grid_a: np.ndarray, grid_b_permuted: np.ndarray,
                         tol: float) -> Tuple[int, int]:
    """Count how many points of face B (small, after permutation) match face A (large).

    Face A is the large face, face B is the small face. After applying the
    permutation to face B, check how many of B's transformed points exist
    within face A (within tolerance). If all of B's points match, the larger
    face A should be split.

    Args:
        grid_a: Canonical grid of face A (large), shape (nu_a, nv_a, 3).
        grid_b_permuted: Permuted grid of face B (small), shape (nu_b, nv_b, 3).
        tol: Euclidean distance tolerance.

    Returns:
        (match_count, total_b_points).
    """
    pts_a = grid_a.reshape(-1, 3)
    pts_b = grid_b_permuted.reshape(-1, 3)
    tol2 = tol * tol
    count = 0
    for b in pts_b:
        diffs = np.sum((pts_a - b) ** 2, axis=1)
        if np.any(diffs <= tol2):
            count += 1
    return count, len(pts_b)


def determine_plane(lb1, ub1, lb2, ub2) -> str:
    """Determine if two faces are in-plane or cross-plane.

    Args:
        lb1, ub1: Diagonal corners of face 1.
        lb2, ub2: Diagonal corners of face 2.

    Returns:
        "in-plane" if both faces share the same constant axis, "cross-plane" otherwise.
    """
    lo1 = [min(lb1[d], ub1[d]) for d in range(3)]
    hi1 = [max(lb1[d], ub1[d]) for d in range(3)]
    lo2 = [min(lb2[d], ub2[d]) for d in range(3)]
    hi2 = [max(lb2[d], ub2[d]) for d in range(3)]
    const_a = constant_axis(lo1, hi1)
    const_b = constant_axis(lo2, hi2)
    return "in-plane" if const_a == const_b else "cross-plane"


def _try_stored_then_all_perms(grid_a: np.ndarray, grid_b: np.ndarray,
                                stored_perm, fm: dict, tol: float) -> int:
    """Try stored permutation first, then brute-force all 8.

    Args:
        grid_a: Canonical grid of face A.
        grid_b: Canonical grid of face B.
        stored_perm: Stored permutation index (int or None).
        fm: Face match dict (used only for lb/ub to determine plane).
        tol: Euclidean distance tolerance.

    Returns:
        Permutation index (0-7) on success, -1 if none works.
    """
    nu_a, nv_a = grid_a.shape[:2]
    if stored_perm is not None and stored_perm >= 0:
        g = apply_permutation(grid_b, stored_perm)
        if g.shape[:2] == (nu_a, nv_a) and verify_match(grid_a, g, tol):
            return stored_perm
    return try_all_permutations(grid_a, grid_b, tol)


# ---------------------------------------------------------------------------
# verify_connectivity
# ---------------------------------------------------------------------------

def verify_connectivity(blocks: List[Block], face_matches: list, tol: float = 1E-6):
    """Verify connectivity face matches using permutation matrices.

    For each face match:
    1. GCD-reduce blocks and scale indices.
    2. Extract both faces as canonical 2D grids.
    3. Try stored permutation_index first (if available).
    4. Fall back to trying all 8 permutations.
    5. On success, store the correct permutation_index in the orientation.

    Args:
        blocks: List of all blocks (original full-resolution).
        face_matches: List of face_match dicts.
        tol: Euclidean distance tolerance.

    Returns:
        (verified, mismatched) lists of face_match dicts.
    """
    gcd_to_use = compute_min_gcd(blocks)
    reduced_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    scaled_matches = deepcopy(face_matches)
    # Normalize to lb/ub format before scaling
    for fm in scaled_matches:
        for side in ['block1', 'block2']:
            lo, hi = get_bounds(fm[side])
            fm[side]['lb'] = list(lo)
            fm[side]['ub'] = list(hi)
    scale_face_bounds(scaled_matches, gcd_to_use, divide=True)

    verified = []
    mismatched = []

    for idx, fm in enumerate(scaled_matches):
        b1 = fm['block1']
        b2 = fm['block2']

        try:
            grid_a, nu_a, nv_a = extract_canonical_grid(
                reduced_blocks[b1['block_index']], b1['lb'], b1['ub'])
            grid_b, nu_b, nv_b = extract_canonical_grid(
                reduced_blocks[b2['block_index']], b2['lb'], b2['ub'])
        except (ValueError, IndexError):
            mismatched.append(face_matches[idx])
            continue

        stored_perm = fm.get('orientation', {}).get('permutation_index')
        perm_idx = _try_stored_then_all_perms(grid_a, grid_b, stored_perm, fm, tol)
        if perm_idx >= 0:
            corrected = deepcopy(face_matches[idx])
            plane = determine_plane(b1['lb'], b1['ub'], b2['lb'], b2['ub'])
            # In-plane: direction encoded in lb/ub → export -1.
            # Cross-plane: lb/ub can't encode axis swap → export actual index.
            export_perm = -1 if plane == 'in-plane' else perm_idx
            corrected['orientation'] = {
                'permutation_index': export_perm,
                'plane': plane,
                'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
            }
            verified.append(corrected)
        else:
            orig = face_matches[idx]
            lo1, hi1 = get_bounds(orig['block1'])
            lo2, hi2 = get_bounds(orig['block2'])
            print(f"verify_connectivity: MISMATCH at index {idx}")
            print(f"  block {orig['block1']['block_index']}: lo={lo1} hi={hi1}")
            print(f"  block {orig['block2']['block_index']}: lo={lo2} hi={hi2}")
            mismatched.append(face_matches[idx])

    return verified, mismatched


# ---------------------------------------------------------------------------
# verify_periodicity
# ---------------------------------------------------------------------------

def verify_periodicity(blocks: List[Block], face_matches: list, theta: float,
                       rotation_axis: str = 'x', tol: float = 1E-6):
    """Verify periodic face matches using permutation matrices with rotation.

    For each face match, rotates block1 by +/- theta and uses the same
    canonical grid + permutation approach as :func:`verify_connectivity`.

    Args:
        blocks: List of all blocks (original full-resolution).
        face_matches: List of face_match dicts from periodicity.
        theta: Rotation angle in degrees.
        rotation_axis: Axis of rotation 'x', 'y', or 'z'.
        tol: Euclidean distance tolerance.

    Returns:
        (verified, mismatched) lists of face_match dicts.
    """
    gcd_to_use = compute_min_gcd(blocks)
    reduced_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    rotation_matrix_pos = create_rotation_matrix(radians(theta), rotation_axis)
    rotation_matrix_neg = create_rotation_matrix(radians(-theta), rotation_axis)

    rotated_blocks_pos = [rotate_block(b, rotation_matrix_pos) for b in reduced_blocks]
    rotated_blocks_neg = [rotate_block(b, rotation_matrix_neg) for b in reduced_blocks]

    scaled_matches = deepcopy(face_matches)
    # Normalize to lb/ub format before scaling
    for fm in scaled_matches:
        for side in ['block1', 'block2']:
            lo, hi = get_bounds(fm[side])
            fm[side]['lb'] = list(lo)
            fm[side]['ub'] = list(hi)
    scale_face_bounds(scaled_matches, gcd_to_use, divide=True)

    verified = []
    mismatched_list = []

    for idx, fm in enumerate(scaled_matches):
        b1 = fm['block1']
        b2 = fm['block2']
        b1_idx = b1['block_index']
        b2_idx = b2['block_index']

        try:
            grid_b, nu_b, nv_b = extract_canonical_grid(
                reduced_blocks[b2_idx], b2['lb'], b2['ub'])
        except (ValueError, IndexError):
            mismatched_list.append(face_matches[idx])
            continue

        found = False

        for rotated_blocks in [rotated_blocks_pos, rotated_blocks_neg]:
            if found:
                break

            try:
                grid_a, nu_a, nv_a = extract_canonical_grid(
                    rotated_blocks[b1_idx], b1['lb'], b1['ub'])
            except (ValueError, IndexError):
                continue

            stored_perm = fm.get('orientation', {}).get('permutation_index')
            perm_idx = _try_stored_then_all_perms(grid_a, grid_b, stored_perm, fm, tol)
            if perm_idx >= 0:
                corrected = deepcopy(face_matches[idx])
                plane = determine_plane(b1['lb'], b1['ub'], b2['lb'], b2['ub'])
                export_perm = -1 if plane == 'in-plane' else perm_idx
                corrected['orientation'] = {
                    'permutation_index': export_perm,
                    'plane': plane,
                    'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
                }
                verified.append(corrected)
                found = True
                break

        if not found:
            orig = face_matches[idx]
            lo1, hi1 = get_bounds(orig['block1'])
            lo2, hi2 = get_bounds(orig['block2'])
            print(f"verify_periodicity: MISMATCH at index {idx}")
            print(f"  block {orig['block1']['block_index']}: lo={lo1} hi={hi1}")
            print(f"  block {orig['block2']['block_index']}: lo={lo2} hi={hi2}")
            mismatched_list.append(face_matches[idx])

    return verified, mismatched_list
