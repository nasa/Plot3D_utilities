from typing import Dict, List, Set, Tuple
import numpy as np 
from .block import Block
from .blockfunctions import find_matching_faces, build_connectivity_graph, standardize_block_orientation
from .write import write_plot3D


def rotate_block_to_align_faces(X, Y, Z, face_from: str, face_to: str):
    """
    Rotate block2 geometry so that face_from aligns with face_to.
    Currently supports: jmax → imin.
    """
    if (face_from, face_to) == ('jmax', 'imin'):
        Xr = np.transpose(X, (1, 0, 2))
        Yr = np.transpose(Y, (1, 0, 2))
        Zr = np.transpose(Z, (1, 0, 2))
        Xr = np.flip(Xr, axis=0)
        Yr = np.flip(Yr, axis=0)
        Zr = np.flip(Zr, axis=0)
        return Xr, Yr, Zr
    raise NotImplementedError(f"Rotation from {face_from} to {face_to} not implemented.")

def combine_2_blocks_mixed_pairing(block1, block2, tol=1e-8):
    """
    Combine block1 and block2 by matching and aligning any pair of faces, including cross-axis combinations.
    Automatically transposes and flips block2, and flips block1 if necessary to maintain monotonic physical direction.

    This version generalizes the direction check to detect which of X, Y, or Z is varying most along the stacking axis,
    and uses that to determine whether block1 needs to be flipped before concatenation.

    Returns:
        Block: The merged block with physically consistent orientation.
    """
    face1, face2, flip_flags = find_matching_faces(block1, block2, tol=tol)
    if face1 is None or flip_flags is None:
        print("No matching faces or incompatible orientation.")
        return block1

    flip_ud, flip_lr = flip_flags

    face_axis_normal = {
        'imin': (0, -1), 'imax': (0, 1),
        'jmin': (1, -1), 'jmax': (1, 1),
        'kmin': (2, -1), 'kmax': (2, 1),
    }

    axis1, dir1 = face_axis_normal[face1]
    axis2, dir2 = face_axis_normal[face2] # type: ignore

    # Transpose block2 to align its face axis with block1's
    transpose = None
    if axis1 != axis2:
        if (axis1, axis2) in [(0, 1), (1, 0)]:
            transpose = (1, 0, 2)
        elif (axis1, axis2) in [(0, 2), (2, 0)]:
            transpose = (2, 1, 0)
        elif (axis1, axis2) in [(1, 2), (2, 1)]:
            transpose = (0, 2, 1)

    X2, Y2, Z2 = block2.X.copy(), block2.Y.copy(), block2.Z.copy()
    if transpose:
        X2, Y2, Z2 = np.transpose(X2, transpose), np.transpose(Y2, transpose), np.transpose(Z2, transpose)

    # Apply local flips from face alignment
    if face2 in ['imin', 'imax']:
        if flip_ud: X2, Y2, Z2 = np.flip(X2, 1), np.flip(Y2, 1), np.flip(Z2, 1)
        if flip_lr: X2, Y2, Z2 = np.flip(X2, 2), np.flip(Y2, 2), np.flip(Z2, 2)
    elif face2 in ['jmin', 'jmax']:
        if flip_ud: X2, Y2, Z2 = np.flip(X2, 0), np.flip(Y2, 0), np.flip(Z2, 0)
        if flip_lr: X2, Y2, Z2 = np.flip(X2, 2), np.flip(Y2, 2), np.flip(Z2, 2)
    elif face2 in ['kmin', 'kmax']:
        if flip_ud: X2, Y2, Z2 = np.flip(X2, 0), np.flip(Y2, 0), np.flip(Z2, 0)
        if flip_lr: X2, Y2, Z2 = np.flip(X2, 1), np.flip(Y2, 1), np.flip(Z2, 1)

    # Determine stacking axis
    stack_axis = axis1

    # Identify which physical axis (X/Y/Z) changes most along the stacking axis
    def get_dominant_step(A, axis):
        n = A.shape[axis]
        center = [s // 2 for s in A.shape]
        center[axis] = slice(None)
        values = A[tuple(center)]
        return np.nanmean(values[-1] - values[0])

    # Step values in block1
    stepX1 = get_dominant_step(block1.X, stack_axis)
    stepY1 = get_dominant_step(block1.Y, stack_axis)
    stepZ1 = get_dominant_step(block1.Z, stack_axis)
    dominant_axis = np.argmax(np.abs([stepX1, stepY1, stepZ1]))
    step1 = [stepX1, stepY1, stepZ1][dominant_axis]

    # Step values in block2 (transformed and flipped)
    stepX2 = get_dominant_step(X2, stack_axis)
    stepY2 = get_dominant_step(Y2, stack_axis)
    stepZ2 = get_dominant_step(Z2, stack_axis)
    step2 = [stepX2, stepY2, stepZ2][dominant_axis]

    # Flip block1 if physical step directions are inconsistent
    if np.sign(step1) != np.sign(step2):
        block1.X = np.flip(block1.X, axis=stack_axis)
        block1.Y = np.flip(block1.Y, axis=stack_axis)
        block1.Z = np.flip(block1.Z, axis=stack_axis)

    # Slice off overlapping face from block2
    slicer = [slice(None)] * 3
    slicer[stack_axis] = slice(1, None) if face2.endswith('min') else slice(0, -1) # type: ignore
    X2s, Y2s, Z2s = X2[tuple(slicer)], Y2[tuple(slicer)], Z2[tuple(slicer)]

    # Concatenate along stack axis
    if face2.endswith('min'): # type: ignore
        X = np.concatenate([block1.X, X2s], axis=stack_axis)
        Y = np.concatenate([block1.Y, Y2s], axis=stack_axis)
        Z = np.concatenate([block1.Z, Z2s], axis=stack_axis)
    else:
        X = np.concatenate([X2s, block1.X], axis=stack_axis)
        Y = np.concatenate([Y2s, block1.Y], axis=stack_axis)
        Z = np.concatenate([Z2s, block1.Z], axis=stack_axis)
    return standardize_block_orientation(Block(X, Y, Z))



def combine_blocks_mixed_pairs(blocks: List[Block], tol: float = 1e-8, max_tries: int = 4) -> Tuple[List[Block], List[int]]:
    """
    Combine as many blocks as possible from a group of up to 8 blocks via face matching.

    Parameters
    ----------
    blocks : List[Block]
        List of up to 8 Block objects.
    tol : float
        Tolerance for face matching.
    max_tries : int
        Number of passes to try merging the blocks further.

    Returns
    -------
    merged_blocks : List[Block]
        List of successfully merged Block objects (1 or more).
    used_indices : List[int]
        Indices of blocks that were used in any successful merge.
    """
    from itertools import combinations

    remaining = list(enumerate(blocks))  # [(index, Block)]
    used_indices = set()

    # Initial merged_blocks are just the input blocks
    merged_blocks = [blk for _, blk in remaining]
    index_lookup = {id(blk): idx for idx, blk in remaining}

    tries = 0
    while len(merged_blocks) > 1 and tries < max_tries:
        new_merged = []
        merged_flags = [False] * len(merged_blocks)
        used_this_pass = set()
        skip = set()

        i = 0
        while i < len(merged_blocks):
            if i in skip:
                i += 1
                continue

            blk_a = merged_blocks[i]
            merged = None
            found = False

            for j in range(i + 1, len(merged_blocks)):
                if j in skip:
                    continue
                blk_b = merged_blocks[j]
                face1, face2,_ = find_matching_faces(blk_a, blk_b, tol=tol)
                if face1 is not None:
                    try:
                        merged = combine_2_blocks_mixed_pairing(blk_a, blk_b, tol=tol)
                        found = True
                        break
                    except Exception as e:
                        print(f"⚠️ Failed to merge blocks {i} and {j}: {e}")

            if found:
                new_merged.append(merged)
                skip.update([i, j]) # type: ignore
            else:
                new_merged.append(blk_a)
                skip.add(i)

            i += 1

        # Add any unmerged blocks at the end
        for k in range(len(merged_blocks)):
            if k not in skip:
                new_merged.append(merged_blocks[k])

        merged_blocks = new_merged
        tries += 1

    # Recover used indices (conservatively: all blocks involved in merging)
    used_indices = list(range(len(blocks)))  # all blocks are assumed used here

    return merged_blocks, used_indices

def combine_nxnxn_cubes_mixed_pairs(
    blocks: List[Block],
    connectivities: List[List[Dict]],
    cube_size: int = 2,
    tol: float = 1e-8
) -> List[Tuple[Block, Set[int]]]:
    """
    Find and combine all non-overlapping nxnxn cube groups of blocks
    using face connectivity data and return the merged components.

    Parameters
    ----------
    blocks : List[Block]
        All input block objects.
    connectivities : List of face match metadata pairs
        Face connectivity between blocks.
    cube_size : int
        Size of the cube (e.g. 2 for 2x2x2, 4 for 4x4x4).
    tol : float
        Face match tolerance.

    Returns
    -------
    List[Tuple[Block, Set[int]]]
        A list of merged Block objects and their source block indices.
    """
    from itertools import product

    used = set()
    merged_groups = []
    G = build_connectivity_graph(connectivities)

    remaining_indices = list(range(len(blocks)))

    def find_nxnxn_group(seed_index):
        from collections import deque
        visited = set()
        queue = deque([seed_index])
        group = set()

        while queue and len(group) < cube_size ** 3:
            idx = queue.popleft()
            if idx in visited or idx in used:
                continue
            visited.add(idx)
            group.add(idx)
            for nbr in G.neighbors(idx):
                if nbr not in visited and nbr not in used:
                    queue.append(nbr)

        return group if len(group) == cube_size ** 3 else None

    while True:
        before_len = len(remaining_indices)
        merged_this_round = False
        new_used = set()

        i = 0
        while i < len(remaining_indices):
            seed_index = remaining_indices[i]
            if seed_index in used:
                i += 1
                continue

            group_indices = find_nxnxn_group(seed_index)
            if not group_indices or group_indices & new_used:
                i += 1
                continue

            group_block_list = [blocks[k] for k in sorted(group_indices)]
            index_mapping = {i: orig_idx for i, orig_idx in enumerate(sorted(group_indices))}

            try:
                partial_merges, local_indices = combine_blocks_mixed_pairs(group_block_list, tol=tol)
                for merged_block in partial_merges:
                    merged_group = {index_mapping[i] for i in local_indices}
                    merged_groups.append((merged_block, merged_group))
                    # write_plot3D('block1.xyz',[merged_block])
                    new_used.update(merged_group)
                merged_this_round = True
            except Exception as e:
                print(f"⚠️ Skipping group {group_indices} due to error: {e}")
                i += 1
                continue

            # Update remaining_indices and restart inner loop
            remaining_indices = [idx for idx in remaining_indices if idx not in new_used]
            i = 0

        used.update(new_used)

        if not merged_this_round or len(remaining_indices) == before_len:
            print("✅ No further merges possible. Appending unmerged blocks.")
            for idx in remaining_indices:
                merged_groups.append((blocks[idx], {idx}))
            break
    
    return merged_groups