from .block import Block
from .blockfunctions import reduce_blocks, rotate_block
from .periodicity import create_rotation_matrix
from typing import List
from copy import deepcopy
from math import radians
import math


def verify_connectivity(blocks: List[Block], face_matches: list, tol: float = 1E-6):
    """Verifies that the diagonal corners of face_matches are spatially consistent.

    For each face_match, checks that block1's lower corner coordinates match
    block2's lower corner coordinates (and similarly for upper corners) within
    the specified tolerance. If the stored diagonal doesn't match, tries all
    permutations of block2's face corners. If a valid permutation is found,
    the face_match is corrected and added to the verified list.

    Uses GCD reduction (same as connectivity_fast) for efficient coordinate lookups.

    Args:
        blocks (List[Block]): List of all blocks (original full-resolution)
        face_matches (list): List of face_match dicts from connectivity or periodicity
        tol (float, optional): Euclidean distance tolerance. Defaults to 1E-6.

    Returns:
        (list): verified face_matches whose diagonals are confirmed or corrected
        (list): mismatched face_matches where no corner permutation matched
    """
    # Compute GCD and reduce blocks (same pattern as connectivity_fast)
    gcd_array = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
    gcd_to_use = min(gcd_array)
    reduced_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # Scale down face_matches indices by GCD
    scaled_matches = deepcopy(face_matches)
    for fm in scaled_matches:
        for side in ['block1', 'block2']:
            fm[side]['lb'] = [x // gcd_to_use for x in fm[side]['lb']]
            fm[side]['ub'] = [x // gcd_to_use for x in fm[side]['ub']]

    verified = list()
    mismatched = list()

    for idx in range(len(scaled_matches)):
        fm = scaled_matches[idx]
        b1 = fm['block1']
        b2 = fm['block2']
        b1_idx = b1['block_index']
        b2_idx = b2['block_index']
        block1 = reduced_blocks[b1_idx]
        block2 = reduced_blocks[b2_idx]

        # Block1 diagonal coordinates
        x1_l = block1.X[b1['lb'][0], b1['lb'][1], b1['lb'][2]]
        y1_l = block1.Y[b1['lb'][0], b1['lb'][1], b1['lb'][2]]
        z1_l = block1.Z[b1['lb'][0], b1['lb'][1], b1['lb'][2]]

        x1_u = block1.X[b1['ub'][0], b1['ub'][1], b1['ub'][2]]
        y1_u = block1.Y[b1['ub'][0], b1['ub'][1], b1['ub'][2]]
        z1_u = block1.Z[b1['ub'][0], b1['ub'][1], b1['ub'][2]]

        # Enumerate unique corners of block2's face
        I2 = [b2['lb'][0], b2['ub'][0]]
        J2 = [b2['lb'][1], b2['ub'][1]]
        K2 = [b2['lb'][2], b2['ub'][2]]

        unique_corners = list()
        seen = set()
        for i in I2:
            for j in J2:
                for k in K2:
                    key = (i, j, k)
                    if key not in seen:
                        seen.add(key)
                        unique_corners.append(key)

        # Check stored diagonal first
        x2_l = block2.X[b2['lb'][0], b2['lb'][1], b2['lb'][2]]
        y2_l = block2.Y[b2['lb'][0], b2['lb'][1], b2['lb'][2]]
        z2_l = block2.Z[b2['lb'][0], b2['lb'][1], b2['lb'][2]]

        x2_u = block2.X[b2['ub'][0], b2['ub'][1], b2['ub'][2]]
        y2_u = block2.Y[b2['ub'][0], b2['ub'][1], b2['ub'][2]]
        z2_u = block2.Z[b2['ub'][0], b2['ub'][1], b2['ub'][2]]

        dx = x2_l - x1_l; dy = y2_l - y1_l; dz = z2_l - z1_l
        d_lower = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx = x2_u - x1_u; dy = y2_u - y1_u; dz = z2_u - z1_u
        d_upper = math.sqrt(dx*dx + dy*dy + dz*dz)

        if d_lower < tol and d_upper < tol:
            verified.append(face_matches[idx])
            continue

        # Try all permutations of block2's corners
        found = False
        best_d_lower = d_lower
        best_d_upper = d_upper

        for corner_lower in unique_corners:
            for corner_upper in unique_corners:
                if corner_lower == corner_upper:
                    continue

                il, jl, kl = corner_lower
                iu, ju, ku = corner_upper

                x2_l = block2.X[il, jl, kl]
                y2_l = block2.Y[il, jl, kl]
                z2_l = block2.Z[il, jl, kl]

                x2_u = block2.X[iu, ju, ku]
                y2_u = block2.Y[iu, ju, ku]
                z2_u = block2.Z[iu, ju, ku]

                dx = x2_l - x1_l; dy = y2_l - y1_l; dz = z2_l - z1_l
                dl = math.sqrt(dx*dx + dy*dy + dz*dz)
                dx = x2_u - x1_u; dy = y2_u - y1_u; dz = z2_u - z1_u
                du = math.sqrt(dx*dx + dy*dy + dz*dz)

                if dl < best_d_lower:
                    best_d_lower = dl
                if du < best_d_upper:
                    best_d_upper = du

                if dl < tol and du < tol:
                    corrected = deepcopy(face_matches[idx])
                    corrected['block2']['lb'] = [il * gcd_to_use, jl * gcd_to_use, kl * gcd_to_use]
                    corrected['block2']['ub'] = [iu * gcd_to_use, ju * gcd_to_use, ku * gcd_to_use]
                    verified.append(corrected)
                    if b1_idx == b2_idx:
                        print("verify_connectivity: Self-match corrected for block index {0}".format(b1_idx))
                    found = True
                    break
            if found:
                break

        if not found:
            orig = face_matches[idx]
            b1_orig = orig['block1']
            b2_orig = orig['block2']
            print(f"verify_connectivity: MISMATCH at face_match index {idx}")
            print(f"  block1 (block_index={b1_orig['block_index']}): "
                  f"lower=({b1_orig['lb'][0]},{b1_orig['lb'][1]},{b1_orig['lb'][2]}) "
                  f"upper=({b1_orig['ub'][0]},{b1_orig['ub'][1]},{b1_orig['ub'][2]})")
            print(f"  block2 (block_index={b2_orig['block_index']}): "
                  f"lower=({b2_orig['lb'][0]},{b2_orig['lb'][1]},{b2_orig['lb'][2]}) "
                  f"upper=({b2_orig['ub'][0]},{b2_orig['ub'][1]},{b2_orig['ub'][2]})")
            print(f"  block1 lower xyz = ({x1_l:.6e}, {y1_l:.6e}, {z1_l:.6e})")
            print(f"  block1 upper xyz = ({x1_u:.6e}, {y1_u:.6e}, {z1_u:.6e})")
            print(f"  Closest block2 corner dist to block1 lower: {best_d_lower:.6e}")
            print(f"  Closest block2 corner dist to block1 upper: {best_d_upper:.6e}")
            mismatched.append(face_matches[idx])

    return verified, mismatched


def verify_periodicity(blocks: List[Block], face_matches: list, theta: float, rotation_axis: str = 'x', tol: float = 1E-6):
    """Verifies that the diagonal corners of periodic face_matches are spatially
    consistent after rotating block1 by ±theta about the given axis.

    For each face_match, rotates block1 by +theta (and if needed -theta) and checks
    that the rotated lower/upper corner coordinates match block2's lower/upper
    corners within tolerance. If the stored diagonal doesn't match, tries all
    permutations of block2's face corners. If a valid permutation is found,
    the face_match is corrected and added to the verified list.

    Uses GCD reduction (same as connectivity_fast / periodicity_fast) for
    efficient coordinate lookups.

    Args:
        blocks (List[Block]): List of all blocks (original full-resolution)
        face_matches (list): List of face_match dicts from periodicity or rotated_periodicity
        theta (float): Rotation angle in degrees
        rotation_axis (str, optional): Axis of rotation 'x', 'y', or 'z'. Defaults to 'x'.
        tol (float, optional): Euclidean distance tolerance. Defaults to 1E-6.

    Returns:
        (list): verified face_matches whose diagonals are confirmed or corrected
        (list): mismatched face_matches where no corner permutation matched
    """
    # Compute GCD and reduce blocks (same pattern as connectivity_fast)
    gcd_array = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
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

    verified = list()
    mismatched = list()

    for idx in range(len(scaled_matches)):
        fm = scaled_matches[idx]
        b1 = fm['block1']
        b2 = fm['block2']
        b1_idx = b1['block_index']
        b2_idx = b2['block_index']
        block2 = reduced_blocks[b2_idx]

        # Enumerate unique corners of block2's face
        I2 = [b2['lb'][0], b2['ub'][0]]
        J2 = [b2['lb'][1], b2['ub'][1]]
        K2 = [b2['lb'][2], b2['ub'][2]]

        unique_corners = list()
        seen = set()
        for i in I2:
            for j in J2:
                for k in K2:
                    key = (i, j, k)
                    if key not in seen:
                        seen.add(key)
                        unique_corners.append(key)

        found = False
        best_d_lower = float('inf')
        best_d_upper = float('inf')

        # Try +theta rotation first, then -theta
        for rotated_blocks in [rotated_blocks_pos, rotated_blocks_neg]:
            if found:
                break

            block1_rotated = rotated_blocks[b1_idx]

            # Block1 rotated diagonal coordinates
            x1_l = block1_rotated.X[b1['lb'][0], b1['lb'][1], b1['lb'][2]]
            y1_l = block1_rotated.Y[b1['lb'][0], b1['lb'][1], b1['lb'][2]]
            z1_l = block1_rotated.Z[b1['lb'][0], b1['lb'][1], b1['lb'][2]]

            x1_u = block1_rotated.X[b1['ub'][0], b1['ub'][1], b1['ub'][2]]
            y1_u = block1_rotated.Y[b1['ub'][0], b1['ub'][1], b1['ub'][2]]
            z1_u = block1_rotated.Z[b1['ub'][0], b1['ub'][1], b1['ub'][2]]

            # Check stored diagonal first
            x2_l = block2.X[b2['lb'][0], b2['lb'][1], b2['lb'][2]]
            y2_l = block2.Y[b2['lb'][0], b2['lb'][1], b2['lb'][2]]
            z2_l = block2.Z[b2['lb'][0], b2['lb'][1], b2['lb'][2]]

            x2_u = block2.X[b2['ub'][0], b2['ub'][1], b2['ub'][2]]
            y2_u = block2.Y[b2['ub'][0], b2['ub'][1], b2['ub'][2]]
            z2_u = block2.Z[b2['ub'][0], b2['ub'][1], b2['ub'][2]]

            dx = x2_l - x1_l; dy = y2_l - y1_l; dz = z2_l - z1_l
            d_lower = math.sqrt(dx*dx + dy*dy + dz*dz)
            dx = x2_u - x1_u; dy = y2_u - y1_u; dz = z2_u - z1_u
            d_upper = math.sqrt(dx*dx + dy*dy + dz*dz)

            if d_lower < best_d_lower:
                best_d_lower = d_lower
            if d_upper < best_d_upper:
                best_d_upper = d_upper

            if d_lower < tol and d_upper < tol:
                verified.append(face_matches[idx])
                found = True
                break

            # Try all permutations of block2's corners
            for corner_lower in unique_corners:
                for corner_upper in unique_corners:
                    if corner_lower == corner_upper:
                        continue

                    il, jl, kl = corner_lower
                    iu, ju, ku = corner_upper

                    x2_l = block2.X[il, jl, kl]
                    y2_l = block2.Y[il, jl, kl]
                    z2_l = block2.Z[il, jl, kl]

                    x2_u = block2.X[iu, ju, ku]
                    y2_u = block2.Y[iu, ju, ku]
                    z2_u = block2.Z[iu, ju, ku]

                    dx = x2_l - x1_l; dy = y2_l - y1_l; dz = z2_l - z1_l
                    dl = math.sqrt(dx*dx + dy*dy + dz*dz)
                    dx = x2_u - x1_u; dy = y2_u - y1_u; dz = z2_u - z1_u
                    du = math.sqrt(dx*dx + dy*dy + dz*dz)

                    if dl < best_d_lower:
                        best_d_lower = dl
                    if du < best_d_upper:
                        best_d_upper = du

                    if dl < tol and du < tol:
                        corrected = deepcopy(face_matches[idx])
                        corrected['block2']['lb'] = [il * gcd_to_use, jl * gcd_to_use, kl * gcd_to_use]
                        corrected['block2']['ub'] = [iu * gcd_to_use, ju * gcd_to_use, ku * gcd_to_use]
                        verified.append(corrected)
                        found = True
                        break
                if found:
                    break

        if not found:
            orig = face_matches[idx]
            b1_orig = orig['block1']
            b2_orig = orig['block2']
            print(f"verify_periodicity: MISMATCH at face_match index {idx}")
            print(f"  block1 (block_index={b1_orig['block_index']}): "
                  f"lower=({b1_orig['lb'][0]},{b1_orig['lb'][1]},{b1_orig['lb'][2]}) "
                  f"upper=({b1_orig['ub'][0]},{b1_orig['ub'][1]},{b1_orig['ub'][2]})")
            print(f"  block2 (block_index={b2_orig['block_index']}): "
                  f"lower=({b2_orig['lb'][0]},{b2_orig['lb'][1]},{b2_orig['lb'][2]}) "
                  f"upper=({b2_orig['ub'][0]},{b2_orig['ub'][1]},{b2_orig['ub'][2]})")
            print(f"  Closest rotated block1 corner dist to block2 lower: {best_d_lower:.6e}")
            print(f"  Closest rotated block1 corner dist to block2 upper: {best_d_upper:.6e}")
            mismatched.append(face_matches[idx])

    return verified, mismatched
