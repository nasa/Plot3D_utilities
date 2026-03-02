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


def _extract_face_points(block: Block, lb: list, ub: list) -> np.ndarray:
    """Extract all face points in directed traversal order (lb -> ub).

    Returns an (N, 3) array where point n from face A must match
    point n from face B — the diagonal convention preserves the
    node-to-node mapping between blocks.
    """
    pts = []
    for i in _directed(lb[0], ub[0]):
        for j in _directed(lb[1], ub[1]):
            for k in _directed(lb[2], ub[2]):
                pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
    return np.array(pts)


def _generate_permutations(lb: list, ub: list) -> List[Tuple[list, list]]:
    """Generate all 8 traversal permutations for a face (4 direct + 4 transposed).

    Determines which axis is constant, then for the two varying axes
    generates all 4 direction combinations (increase/decrease per axis).
    Then swaps (transposes) the two varying axes and generates 4 more,
    giving 8 total permutations.
    Returns list of (new_lb, new_ub) tuples.
    """
    perms = []
    for dim in range(3):
        if lb[dim] == ub[dim]:
            varying = [d for d in range(3) if d != dim]
            d0, d1 = varying
            vals = [[lb[d0], ub[d0]], [lb[d1], ub[d1]]]
            # 4 direct permutations
            for s0 in [0, 1]:  # 0=forward, 1=reversed for d0
                for s1 in [0, 1]:  # 0=forward, 1=reversed for d1
                    new_lb = list(lb)
                    new_ub = list(ub)
                    new_lb[d0] = vals[0][s0]
                    new_ub[d0] = vals[0][1 - s0]
                    new_lb[d1] = vals[1][s1]
                    new_ub[d1] = vals[1][1 - s1]
                    perms.append((new_lb, new_ub))
            # 4 transposed permutations (swap which axis values go to d0 vs d1)
            for s0 in [0, 1]:
                for s1 in [0, 1]:
                    new_lb = list(lb)
                    new_ub = list(ub)
                    new_lb[d0] = vals[1][s0]   # d1's values → d0's slot
                    new_ub[d0] = vals[1][1 - s0]
                    new_lb[d1] = vals[0][s1]   # d0's values → d1's slot
                    new_ub[d1] = vals[0][1 - s1]
                    perms.append((new_lb, new_ub))
            break
    return perms


def _try_all_permutations(block1: Block, b1_lb: list, b1_ub: list,
                          block2: Block, b2_lb: list, b2_ub: list,
                          tol: float) -> Tuple[bool, list, list]:
    """Try all direction permutations of block2's face against block1's face.

    Holds block1's traversal fixed (lb->ub). For block2, tries all 8
    permutations (4 direct + 4 transposed) of the two varying axes.

    Returns:
        (matched, corrected_lb, corrected_ub) — if matched is True,
        corrected_lb/ub give the block2 diagonal that makes all points align.
    """
    dims1 = _face_dims(b1_lb, b1_ub)
    dims2 = _face_dims(b2_lb, b2_ub)
    n1 = dims1[0] * dims1[1] * dims1[2]
    n2 = dims2[0] * dims2[1] * dims2[2]

    if n1 != n2:
        return False, b2_lb, b2_ub

    pts1 = _extract_face_points(block1, b1_lb, b1_ub)

    for perm_lb, perm_ub in _generate_permutations(b2_lb, b2_ub):
        pts2 = _extract_face_points(block2, perm_lb, perm_ub)
        if pts1.shape != pts2.shape:
            continue
        diffs = np.linalg.norm(pts1 - pts2, axis=1)
        if diffs.max() < tol:
            return True, perm_lb, perm_ub

    return False, b2_lb, b2_ub


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
        (list): mismatched face_matches where no permutation matched
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
    mismatched = []

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
            mismatched.append(face_matches[idx])
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
        matched, corr_lb, corr_ub = _try_all_permutations(
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
            mismatched.append(face_matches[idx])

    return verified, mismatched


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
            matched, corr_lb, corr_ub = _try_all_permutations(
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
