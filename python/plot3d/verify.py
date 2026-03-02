from .block import Block
from .blockfunctions import reduce_blocks, rotate_block
from .periodicity import create_rotation_matrix
from typing import List, Tuple
from copy import deepcopy
from math import radians
import math
import numpy as np


# ---------------------------------------------------------------------------
# Helpers for directed traversal
# ---------------------------------------------------------------------------

def _directed(start: int, end: int) -> range:
    """Inclusive range stepping +1 or -1."""
    return range(start, end + 1) if start <= end else range(start, end - 1, -1)


def _face_dims(lb: list, ub: list) -> Tuple[int, int, int]:
    """Return (di, dj, dk) face dimensions from diagonal corners."""
    return (abs(ub[0] - lb[0]) + 1,
            abs(ub[1] - lb[1]) + 1,
            abs(ub[2] - lb[2]) + 1)


def _extract_face_points(block: Block, lb: list, ub: list, swap_loops: bool = False) -> np.ndarray:
    """Extract all face points in directed traversal order (lb -> ub).

    Returns an (N, 3) array where point n from face A must match
    point n from face B — the diagonal convention preserves the
    node-to-node mapping between blocks.

    If swap_loops is True, the two varying axes swap their loop nesting
    order (outer becomes inner and vice versa). This handles the cross-plane
    case where lb/ub cannot encode loop order.
    """
    ri = list(_directed(lb[0], ub[0]))
    rj = list(_directed(lb[1], ub[1]))
    rk = list(_directed(lb[2], ub[2]))

    if not swap_loops:
        pts = []
        for i in ri:
            for j in rj:
                for k in rk:
                    pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
        return np.array(pts)

    # Swap the two varying axes' loop order
    const_dim = None
    for d in range(3):
        if lb[d] == ub[d]:
            const_dim = d
            break

    if const_dim is None:
        # No constant axis — fall back to normal order
        pts = []
        for i in ri:
            for j in rj:
                for k in rk:
                    pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
        return np.array(pts)

    pts = []
    if const_dim == 0:
        # I-constant: normal is j-outer,k-inner → swapped is k-outer,j-inner
        for k in rk:
            for j in rj:
                pts.append([block.X[ri[0], j, k], block.Y[ri[0], j, k], block.Z[ri[0], j, k]])
    elif const_dim == 1:
        # J-constant: normal is i-outer,k-inner → swapped is k-outer,i-inner
        for k in rk:
            for i in ri:
                pts.append([block.X[i, rj[0], k], block.Y[i, rj[0], k], block.Z[i, rj[0], k]])
    else:
        # K-constant: normal is i-outer,j-inner → swapped is j-outer,i-inner
        for j in rj:
            for i in ri:
                pts.append([block.X[i, j, rk[0]], block.Y[i, j, rk[0]], block.Z[i, j, rk[0]]])
    return np.array(pts)


def _try_all_permutations(block1: Block, b1_lb: list, b1_ub: list,
                          block2: Block, b2_lb: list, b2_ub: list,
                          tol: float) -> Tuple[bool, list, list, int]:
    """Try all 8 permutations of block2's face against block1 using index-based grid manipulation.

    Extracts block2's face points **once** in canonical (ascending) order, then
    applies each of the 8 permutations via numpy index manipulation rather than
    re-extracting points for each permutation.

    Bit encoding: ``perm_idx = u_reversed | (v_reversed << 1) | (swapped << 2)``

    Returns:
        (matched, corrected_lb, corrected_ub, perm_idx) — perm_idx in 0..7
        if matched, -1 otherwise.
    """
    dims1 = _face_dims(b1_lb, b1_ub)
    dims2 = _face_dims(b2_lb, b2_ub)
    n1 = dims1[0] * dims1[1] * dims1[2]
    n2 = dims2[0] * dims2[1] * dims2[2]

    if n1 != n2:
        return False, b2_lb, b2_ub, -1

    pts1 = _extract_face_points(block1, b1_lb, b1_ub)

    # Normalize to ascending ranges
    lo = [min(b2_lb[d], b2_ub[d]) for d in range(3)]
    hi = [max(b2_lb[d], b2_ub[d]) for d in range(3)]

    # Find constant axis
    const_dim = -1
    for d in range(3):
        if lo[d] == hi[d]:
            const_dim = d
            break
    if const_dim < 0:
        return False, b2_lb, b2_ub, -1

    varying = [d for d in range(3) if d != const_dim]
    d0, d1 = varying  # "u" and "v" axes
    nu = hi[d0] - lo[d0] + 1
    nv = hi[d1] - lo[d1] + 1

    # Extract face2 points once in canonical ascending order: u outer, v inner
    grid = np.empty((nu, nv, 3))
    for u in range(nu):
        for v in range(nv):
            idx = [0, 0, 0]
            idx[const_dim] = lo[const_dim]
            idx[d0] = lo[d0] + u
            idx[d1] = lo[d1] + v
            grid[u, v] = [block2.X[idx[0], idx[1], idx[2]],
                          block2.Y[idx[0], idx[1], idx[2]],
                          block2.Z[idx[0], idx[1], idx[2]]]

    # Try each of the 8 permutations
    for perm_idx in range(8):
        u_rev = bool(perm_idx & 1)
        v_rev = bool(perm_idx & 2)
        swap = bool(perm_idx & 4)

        g = grid
        if u_rev:
            g = g[::-1, :, :]
        if v_rev:
            g = g[:, ::-1, :]
        if swap:
            g = g.transpose(1, 0, 2)

        pts2 = g.reshape(-1, 3)
        if pts1.shape != pts2.shape:
            continue
        diffs = np.linalg.norm(pts1 - pts2, axis=1)
        if diffs.max() < tol:
            # Reconstruct corrected lb/ub
            new_lo = list(lo)
            new_hi = list(lo)
            new_lo[const_dim] = lo[const_dim]
            new_hi[const_dim] = lo[const_dim]
            if swap:
                new_lo[d0] = hi[d1] if u_rev else lo[d1]
                new_hi[d0] = lo[d1] if u_rev else hi[d1]
                new_lo[d1] = hi[d0] if v_rev else lo[d0]
                new_hi[d1] = lo[d0] if v_rev else hi[d0]
            else:
                new_lo[d0] = hi[d0] if u_rev else lo[d0]
                new_hi[d0] = lo[d0] if u_rev else hi[d0]
                new_lo[d1] = hi[d1] if v_rev else lo[d1]
                new_hi[d1] = lo[d1] if v_rev else hi[d1]
            return True, new_lo, new_hi, perm_idx

    return False, b2_lb, b2_ub, -1


# ---------------------------------------------------------------------------
# verify_connectivity — full directed point-by-point traversal
# ---------------------------------------------------------------------------

def verify_connectivity(blocks: List[Block], face_matches: list, tol: float = 1E-6):
    """Verifies face_matches using full directed point-by-point traversal.

    For each face_match:
      1. Checks that both faces have the same number of points.
      2. Extracts all points from each face in directed traversal order
         (lb -> ub), stepping +1 or -1 per dimension.
      3. Compares every point: point n from face1 must equal point n
         from face2. If they don't match, tries all 8 permutations
         (4 direct + 4 transposed) of block2's two varying axes.
      4. If a permutation matches, corrects the face_match's block2
         lb/ub and adds it to the verified list.

    Uses GCD reduction for efficient coordinate lookups.

    Args:
        blocks (List[Block]): List of all blocks (original full-resolution)
        face_matches (list): List of face_match dicts from connectivity or periodicity
        tol (float, optional): Euclidean distance tolerance. Defaults to 1E-6.

    Returns:
        (list): verified face_matches whose diagonals are confirmed or corrected
        (list): cross_plane face_matches where points didn't match with standard traversal
    """
    # Compute GCD and reduce blocks
    gcd_array = [math.gcd(b.IMAX - 1, math.gcd(b.JMAX - 1, b.KMAX - 1)) for b in blocks]
    gcd_to_use = min(gcd_array)
    reduced_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # Scale down face_matches indices by GCD
    scaled_matches = deepcopy(face_matches)
    for fm in scaled_matches:
        for side in ['block1', 'block2']:
            fm[side]['lb'] = [x // gcd_to_use for x in fm[side]['lb']]
            fm[side]['ub'] = [x // gcd_to_use for x in fm[side]['ub']]

    verified = []
    cross_plane = []

    for idx in range(len(scaled_matches)):
        fm = scaled_matches[idx]
        b1 = fm['block1']
        b2 = fm['block2']
        block1 = reduced_blocks[b1['block_index']]
        block2 = reduced_blocks[b2['block_index']]

        # Step 1: Dimension check
        dims1 = _face_dims(b1['lb'], b1['ub'])
        dims2 = _face_dims(b2['lb'], b2['ub'])
        n1 = dims1[0] * dims1[1] * dims1[2]
        n2 = dims2[0] * dims2[1] * dims2[2]

        if n1 != n2:
            orig = face_matches[idx]
            print(f"verify_connectivity: DIMENSION MISMATCH at index {idx}")
            print(f"  block {orig['block1']['block_index']}: "
                  f"lb={orig['block1']['lb']} ub={orig['block1']['ub']} dims={dims1}")
            print(f"  block {orig['block2']['block_index']}: "
                  f"lb={orig['block2']['lb']} ub={orig['block2']['ub']} dims={dims2}")
            cross_plane.append(face_matches[idx])
            continue

        # Step 2: Extract face1 points (held constant)
        pts1 = _extract_face_points(block1, b1['lb'], b1['ub'])

        # Step 3: Check stored diagonal first
        pts2 = _extract_face_points(block2, b2['lb'], b2['ub'])
        diffs = np.linalg.norm(pts1 - pts2, axis=1)

        if diffs.max() < tol:
            verified.append(face_matches[idx])
            continue

        # Step 4: Try all permutations of block2's direction
        matched, corr_lb, corr_ub, _perm_idx = _try_all_permutations(
            block1, b1['lb'], b1['ub'],
            block2, b2['lb'], b2['ub'],
            tol
        )

        if matched:
            corrected = deepcopy(face_matches[idx])
            corrected['block2']['lb'] = [x * gcd_to_use for x in corr_lb]
            corrected['block2']['ub'] = [x * gcd_to_use for x in corr_ub]
            verified.append(corrected)
            if b1['block_index'] == b2['block_index']:
                print(f"verify_connectivity: Self-match corrected for block index {b1['block_index']}")
        else:
            orig = face_matches[idx]
            b1_orig = orig['block1']
            b2_orig = orig['block2']
            worst = float(diffs.max())
            n_bad = int(np.sum(diffs > tol))
            print(f"verify_connectivity: POINT MISMATCH at index {idx}")
            print(f"  block {b1_orig['block_index']}: lb={b1_orig['lb']} ub={b1_orig['ub']}")
            print(f"  block {b2_orig['block_index']}: lb={b2_orig['lb']} ub={b2_orig['ub']}")
            print(f"  total points: {len(pts1)}, mismatched: {n_bad}, max dist: {worst:.6e}")
            cross_plane.append(face_matches[idx])

    return verified, cross_plane


# ---------------------------------------------------------------------------
# verify_periodicity — full directed point-by-point traversal with rotation
# ---------------------------------------------------------------------------

def verify_periodicity(blocks: List[Block], face_matches: list, theta: float,
                       rotation_axis: str = 'x', tol: float = 1E-6):
    """Verifies periodic face_matches using full directed point-by-point traversal.

    For each face_match, rotates block1 by +theta (and if needed -theta) and
    checks that every rotated point matches the corresponding point on block2
    in directed traversal order. Tries all direction permutations of block2
    if the stored diagonal doesn't match.

    Uses GCD reduction for efficient coordinate lookups.

    Args:
        blocks (List[Block]): List of all blocks (original full-resolution)
        face_matches (list): List of face_match dicts from periodicity
        theta (float): Rotation angle in degrees
        rotation_axis (str, optional): Axis of rotation 'x', 'y', or 'z'. Defaults to 'x'.
        tol (float, optional): Euclidean distance tolerance. Defaults to 1E-6.

    Returns:
        (list): verified face_matches whose diagonals are confirmed or corrected
        (list): mismatched face_matches where no permutation matched
    """
    # Compute GCD and reduce blocks
    gcd_array = [math.gcd(b.IMAX - 1, math.gcd(b.JMAX - 1, b.KMAX - 1)) for b in blocks]
    gcd_to_use = min(gcd_array)
    reduced_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # Build rotation matrices for +theta and -theta
    rotation_matrix_pos = create_rotation_matrix(radians(theta), rotation_axis)
    rotation_matrix_neg = create_rotation_matrix(radians(-theta), rotation_axis)

    # Pre-rotate all reduced blocks in both directions
    rotated_blocks_pos = [rotate_block(b, rotation_matrix_pos) for b in reduced_blocks]
    rotated_blocks_neg = [rotate_block(b, rotation_matrix_neg) for b in reduced_blocks]

    # Scale down face_matches indices by GCD
    scaled_matches = deepcopy(face_matches)
    for fm in scaled_matches:
        for side in ['block1', 'block2']:
            fm[side]['lb'] = [x // gcd_to_use for x in fm[side]['lb']]
            fm[side]['ub'] = [x // gcd_to_use for x in fm[side]['ub']]

    verified = []
    mismatched = []

    for idx in range(len(scaled_matches)):
        fm = scaled_matches[idx]
        b1 = fm['block1']
        b2 = fm['block2']
        b2_idx = b2['block_index']
        block2 = reduced_blocks[b2_idx]

        # Dimension check
        dims1 = _face_dims(b1['lb'], b1['ub'])
        dims2 = _face_dims(b2['lb'], b2['ub'])
        n1 = dims1[0] * dims1[1] * dims1[2]
        n2 = dims2[0] * dims2[1] * dims2[2]

        if n1 != n2:
            orig = face_matches[idx]
            print(f"verify_periodicity: DIMENSION MISMATCH at index {idx}")
            print(f"  block {orig['block1']['block_index']}: "
                  f"lb={orig['block1']['lb']} ub={orig['block1']['ub']} dims={dims1}")
            print(f"  block {orig['block2']['block_index']}: "
                  f"lb={orig['block2']['lb']} ub={orig['block2']['ub']} dims={dims2}")
            mismatched.append(face_matches[idx])
            continue

        found = False

        # Try +theta rotation first, then -theta
        for rotated_blocks in [rotated_blocks_pos, rotated_blocks_neg]:
            if found:
                break

            block1_rotated = rotated_blocks[b1['block_index']]

            # Check stored diagonal first (full point-by-point)
            pts1 = _extract_face_points(block1_rotated, b1['lb'], b1['ub'])
            pts2 = _extract_face_points(block2, b2['lb'], b2['ub'])
            diffs = np.linalg.norm(pts1 - pts2, axis=1)

            if diffs.max() < tol:
                verified.append(face_matches[idx])
                found = True
                break

            # Try all permutations of block2's direction
            matched, corr_lb, corr_ub, _perm_idx = _try_all_permutations(
                block1_rotated, b1['lb'], b1['ub'],
                block2, b2['lb'], b2['ub'],
                tol
            )

            if matched:
                corrected = deepcopy(face_matches[idx])
                corrected['block2']['lb'] = [x * gcd_to_use for x in corr_lb]
                corrected['block2']['ub'] = [x * gcd_to_use for x in corr_ub]
                verified.append(corrected)
                found = True
                break

        if not found:
            orig = face_matches[idx]
            b1_orig = orig['block1']
            b2_orig = orig['block2']
            print(f"verify_periodicity: MISMATCH at index {idx}")
            print(f"  block {b1_orig['block_index']}: lb={b1_orig['lb']} ub={b1_orig['ub']}")
            print(f"  block {b2_orig['block_index']}: lb={b2_orig['lb']} ub={b2_orig['ub']}")
            mismatched.append(face_matches[idx])

    return verified, mismatched
