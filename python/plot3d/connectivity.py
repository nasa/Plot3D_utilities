from .block import Block
from .blockfunctions import reduce_blocks
from .face import Face
from .facefunctions import create_face_from_diagonals, split_face, get_outer_faces
import math
from itertools import product, combinations
from tqdm import trange
import numpy as np
import pandas as pd
from typing import List, Tuple
import math
from .point_match import point_match
from copy import deepcopy


def _directed(start: int, end: int) -> range:
    """Inclusive range stepping +1 or -1."""
    return range(start, end + 1) if start <= end else range(start, end - 1, -1)


def _extract_face_points(block: Block, lb: list, ub: list) -> np.ndarray:
    """Extract face points in directed lb->ub traversal order. Returns (N,3) array."""
    pts = []
    for i in _directed(lb[0], ub[0]):
        for j in _directed(lb[1], ub[1]):
            for k in _directed(lb[2], ub[2]):
                pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
    return np.array(pts)


def _extract_face_indices(lb: list, ub: list) -> List[Tuple[int, int, int]]:
    """Extract face indices in directed lb->ub traversal order. Returns list of (i,j,k)."""
    indices = []
    for i in _directed(lb[0], ub[0]):
        for j in _directed(lb[1], ub[1]):
            for k in _directed(lb[2], ub[2]):
                indices.append((i, j, k))
    return indices


def _face_point_count(lb: list, ub: list) -> int:
    """Total number of points on a face defined by lb/ub."""
    return ((abs(ub[0] - lb[0]) + 1) *
            (abs(ub[1] - lb[1]) + 1) *
            (abs(ub[2] - lb[2]) + 1))


def _generate_face2_permutations(lb: list, ub: list) -> List[Tuple[list, list]]:
    """Generate all traversal permutations for face2.

    Determines which axis is constant, then for the two varying axes
    generates all 4 direction combinations (increase/decrease per axis).
    Returns list of (new_lb, new_ub) tuples.
    """
    perms = []
    for dim in range(3):  # check which dimension is constant
        if lb[dim] == ub[dim]:
            # This dimension is constant; the other two vary
            varying = [d for d in range(3) if d != dim]
            d0, d1 = varying
            vals = [[lb[d0], ub[d0]], [lb[d1], ub[d1]]]
            # 4 direction combos: (fwd/rev for each varying axis)
            for s0 in [0, 1]:  # 0=forward, 1=reversed
                for s1 in [0, 1]:
                    new_lb = list(lb)
                    new_ub = list(ub)
                    new_lb[d0] = vals[0][s0]
                    new_ub[d0] = vals[0][1 - s0]
                    new_lb[d1] = vals[1][s1]
                    new_ub[d1] = vals[1][1 - s1]
                    perms.append((new_lb, new_ub))
            break  # only one constant axis
    return perms


def _constant_axis(lb: list, ub: list) -> int:
    """Return which axis is constant (0=I, 1=J, 2=K), or -1 if none."""
    for d in range(3):
        if lb[d] == ub[d]:
            return d
    return -1


def _varying_dims(lb: list, ub: list) -> Tuple[int, int]:
    """Return (n_outer, n_inner) for the two varying axes of a face.

    The triple-nested I→J→K loop with one constant axis produces:
      I-const: outer=J, inner=K
      J-const: outer=I, inner=K
      K-const: outer=I, inner=J
    """
    ca = _constant_axis(lb, ub)
    dims = [abs(ub[d] - lb[d]) + 1 for d in range(3)]
    if ca == 0:
        return dims[1], dims[2]   # J outer, K inner
    elif ca == 1:
        return dims[0], dims[2]   # I outer, K inner
    else:
        return dims[0], dims[1]   # I outer, J inner


def _compute_orientation(df: pd.DataFrame, lb1: list, ub1: list) -> List[int]:
    """Compute the orientation vector from a face match DataFrame.

    The orientation maps face1's axes to face2's axes:
        orientation[d] = 1-indexed face2 axis corresponding to face1 axis d.
        (1=I, 2=J, 3=K)

    Direction information is NOT in the orientation — it is already encoded
    in the lb/ub values of each face.

    The mapping is derived from the DataFrame by observing which face2 axis
    changes when a specific face1 axis changes.
    """
    ca1 = _constant_axis(lb1, ub1)
    varying1 = [d for d in range(3) if d != ca1]
    inner1 = varying1[1]  # innermost in I→J→K nesting
    outer1 = varying1[0]  # outermost in I→J→K nesting
    inner_size = abs(ub1[inner1] - lb1[inner1]) + 1

    orientation = [0, 0, 0]
    mapped = set()

    ref_f2 = [int(df.iloc[0]['i2']), int(df.iloc[0]['j2']), int(df.iloc[0]['k2'])]

    # Row 1: face1's innermost varying axis changed by 1
    if len(df) > 1:
        row1_f2 = [int(df.iloc[1]['i2']), int(df.iloc[1]['j2']), int(df.iloc[1]['k2'])]
        for d2 in range(3):
            if row1_f2[d2] != ref_f2[d2]:
                orientation[inner1] = d2 + 1
                mapped.add(d2)
                break

    # Row inner_size: face1's outermost varying axis changed by 1
    if inner_size < len(df):
        row_n_f2 = [int(df.iloc[inner_size]['i2']), int(df.iloc[inner_size]['j2']), int(df.iloc[inner_size]['k2'])]
        for d2 in range(3):
            if row_n_f2[d2] != ref_f2[d2] and d2 not in mapped:
                orientation[outer1] = d2 + 1
                mapped.add(d2)
                break

    # Remaining face2 axis → face1's constant axis
    remaining = [d for d in range(3) if d not in mapped]
    if len(remaining) == 1:
        orientation[ca1] = remaining[0] + 1
    elif len(remaining) == 2:
        # Outer axis has size 1 (line face) — both remaining are unmapped.
        # The face2 constant axis maps to face1's constant axis.
        ca2_lb = [int(df.iloc[0]['i2']), int(df.iloc[0]['j2']), int(df.iloc[0]['k2'])]
        ca2_ub = [int(df.iloc[-1]['i2']), int(df.iloc[-1]['j2']), int(df.iloc[-1]['k2'])]
        for d2 in remaining:
            if ca2_lb[d2] == ca2_ub[d2]:
                orientation[ca1] = d2 + 1
                orientation[outer1] = [r for r in remaining if r != d2][0] + 1
                break
        else:
            # Fallback: assign arbitrarily
            orientation[ca1] = remaining[0] + 1
            orientation[outer1] = remaining[1] + 1

    return orientation


def _extract_2d_grid(block: Block, lb: list, ub: list) -> np.ndarray:
    """Extract face points as 2D array (n_outer, n_inner, 3).

    I-const: (nJ, nK, 3), J-const: (nI, nK, 3), K-const: (nI, nJ, 3).
    """
    n_outer, n_inner = _varying_dims(lb, ub)
    pts = _extract_face_points(block, lb, ub)
    return pts.reshape(n_outer, n_inner, 3)


def _extract_2d_indices(lb: list, ub: list) -> np.ndarray:
    """Extract face indices as 2D array (n_outer, n_inner, 3).

    Each element is [i, j, k].
    """
    n_outer, n_inner = _varying_dims(lb, ub)
    indices = _extract_face_indices(lb, ub)
    return np.array(indices).reshape(n_outer, n_inner, 3)


def _find_corner_on_face(block: Block, face_lb: list, face_ub: list,
                         x: float, y: float, z: float, tol: float):
    """Find the [i,j,k] index on a face where point (x,y,z) exists.

    Returns [i,j,k] or None if not found within tolerance.
    """
    X2 = select_multi_dimensional(block.X,
            (min(face_lb[0], face_ub[0]), max(face_lb[0], face_ub[0])),
            (min(face_lb[1], face_ub[1]), max(face_lb[1], face_ub[1])),
            (min(face_lb[2], face_ub[2]), max(face_lb[2], face_ub[2])))
    Y2 = select_multi_dimensional(block.Y,
            (min(face_lb[0], face_ub[0]), max(face_lb[0], face_ub[0])),
            (min(face_lb[1], face_ub[1]), max(face_lb[1], face_ub[1])),
            (min(face_lb[2], face_ub[2]), max(face_lb[2], face_ub[2])))
    Z2 = select_multi_dimensional(block.Z,
            (min(face_lb[0], face_ub[0]), max(face_lb[0], face_ub[0])),
            (min(face_lb[1], face_ub[1]), max(face_lb[1], face_ub[1])),
            (min(face_lb[2], face_ub[2]), max(face_lb[2], face_ub[2])))

    result = point_match(x, y, z, X2, Y2, Z2, tol)
    if sum(result) == -2:
        return None

    # point_match returns 2D indices into the sliced array.
    # Convert back to 3D block indices.
    p2, q2 = int(result[0]), int(result[1])
    ca = _constant_axis(face_lb, face_ub)
    i_min = min(face_lb[0], face_ub[0])
    j_min = min(face_lb[1], face_ub[1])
    k_min = min(face_lb[2], face_ub[2])

    if ca == 0:  # I-const: X2 is 2D (J, K)
        return [face_lb[0], p2 + j_min, q2 + k_min]
    elif ca == 1:  # J-const: X2 is 2D (I, K)
        return [p2 + i_min, face_lb[1], q2 + k_min]
    else:  # K-const: X2 is 2D (I, J)
        return [p2 + i_min, q2 + j_min, face_lb[2]]


def _try_permutations_with_transpose(block1: Block, lb1: list, ub1: list,
                                      block2: Block, lb2: list, ub2: list,
                                      tol: float) -> Tuple[bool, pd.DataFrame]:
    """Try all 8 permutations (4 direct + 4 transposed) to match face1 to face2.

    Returns:
        (matched, df) where df has columns i1,j1,k1,i2,j2,k2 if matched.
        Returns (False, empty DataFrame) if no permutation works.
    """
    pts1 = _extract_face_points(block1, lb1, ub1)
    idx1 = _extract_face_indices(lb1, ub1)
    n_outer1, n_inner1 = _varying_dims(lb1, ub1)

    for perm_lb, perm_ub in _generate_face2_permutations(lb2, ub2):
        # --- Direct comparison ---
        pts2 = _extract_face_points(block2, perm_lb, perm_ub)
        if pts1.shape == pts2.shape:
            diffs = np.linalg.norm(pts1 - pts2, axis=1)
            if diffs.max() < tol:
                idx2 = _extract_face_indices(perm_lb, perm_ub)
                match_location = [{'i1': idx1[n][0], 'j1': idx1[n][1], 'k1': idx1[n][2],
                                   'i2': idx2[n][0], 'j2': idx2[n][1], 'k2': idx2[n][2]}
                                  for n in range(len(idx1))]
                return True, pd.DataFrame(match_location)

        # --- Transposed comparison ---
        n_outer2, n_inner2 = _varying_dims(perm_lb, perm_ub)
        # Transpose only makes sense if swapping dims gives the right shape
        if n_outer1 == n_inner2 and n_inner1 == n_outer2:
            grid2 = pts2.reshape(n_outer2, n_inner2, 3)
            pts2_T = grid2.transpose(1, 0, 2).reshape(-1, 3)
            if pts1.shape == pts2_T.shape:
                diffs = np.linalg.norm(pts1 - pts2_T, axis=1)
                if diffs.max() < tol:
                    idx2_arr = np.array(_extract_face_indices(perm_lb, perm_ub))
                    idx2_2d = idx2_arr.reshape(n_outer2, n_inner2, 3)
                    idx2_T = idx2_2d.transpose(1, 0, 2).reshape(-1, 3)
                    match_location = [{'i1': idx1[n][0], 'j1': idx1[n][1], 'k1': idx1[n][2],
                                       'i2': int(idx2_T[n][0]), 'j2': int(idx2_T[n][1]), 'k2': int(idx2_T[n][2])}
                                      for n in range(len(idx1))]
                    return True, pd.DataFrame(match_location)

    return False, pd.DataFrame(columns=['i1','j1','k1','i2','j2','k2'])


def find_matching_blocks(block1:Block,block2:Block,block1_outer:List[Face], block2_outer:List[Face],tol:float=1E-6):  
    """Takes two blocks and finds all matching pairs

    Args:
        block1 (Block): Any plot3d Block that is not the same as block2
        block2 (Block): Any plot3d Block that is not the same as block1
        block1_outer (List[Face]): outer faces for block 1. 
        block2_outer (List[Face]): Outer faces for block 2
        tol (float, Optional): tolerance to use. Defaults to 1E-6
    
    Note:
        This function was changed to be given an input of outer faces for block 1 and block 2. Outer faces can change and we should use the updated value
    Returns:
        (tuple): containing
            - **df** (pandas.DataFrame): corners of matching pair as block1_corners,block2_corners ([imin,jmin,kmin],[imax,jmax,kmax]), ([imin,jmin,kmin],[imax,jmax,kmax])
            - **block1_outer** (List[Face]):
            - **block2_outer** (List[Face]): 
    """
    # Check to see if outer face of block 1 matches any of the outer faces of block 2
    block_match_indices = list()

    block1_split_faces = list()
    block2_split_faces = list() 
    # Create a dataframe for block1 and block 2 inner matches, add to df later
    # df,split_faces1,split_faces2 = get_face_intersection(block1_outer[3],block2_outer[4],block1,block2,tol=1E-6)

    # Checks the nodes of the outer faces to see if any of them match 
    match = True
    while match:
        match = False
        for p in range(len(block1_outer)):
            block1_face = block1_outer[p]
            for q in range(len(block2_outer)):
                block2_face = block2_outer[q]
                df, split_faces1, split_faces2 = get_face_intersection(block1_face,block2_face,block1,block2,tol)
                if len(df)>0:   # the number of intersection points has to be more than 4
                    # if not block1_face in block1MatchingFace and not block2_face in block2MatchingFace:
                    block_match_indices.append(df)
                    block1_split_faces.extend(split_faces1)
                    block2_split_faces.extend(split_faces2)
                    match = True
                    break
            if match:
                break
        if match:
            block1_outer.pop(p) # type: ignore
            block2_outer.pop(q) # type: ignore
            block1_outer.extend(block1_split_faces)
            block2_outer.extend(block2_split_faces)
            block1_split_faces.clear()
            block2_split_faces.clear()

    return block_match_indices, block1_outer, block2_outer # Remove duplicates using set and list 

def select_multi_dimensional(T:np.ndarray,dim1:tuple,dim2:tuple, dim3:tuple):
    """Takes a block (T) and selects X,Y,Z from the block given a face's dimensions
        theres really no good way to do this in python 
        
    Args:
        T (np.ndarray): arbitrary array so say a full matrix containing X
        dim1 (tuple): 20,50 this selects X in the i direction from i=20 to 50
        dim2 (tuple): 40,60 this selects X in the j direction from j=40 to 60
        dim3 (tuple): 10,20 this selects X in the k direction from k=10 to 20

    Returns:
        np.ndarray: returns X or Y or Z given some range of I,J,K

    """
    if dim1[0] == dim1[1]:
        return T[ dim1[0], dim2[0]:dim2[1]+1, dim3[0]:dim3[1]+1 ]
    if dim2[0] == dim2[1]:
        return T[ dim1[0]:dim1[1]+1, dim2[0], dim3[0]:dim3[1]+1 ]
    if dim3[0] == dim3[1]:
        return T[ dim1[0]:dim1[1]+1, dim2[0]:dim2[1]+1, dim3[0] ]
    
    return T[dim1[0]:dim1[1], dim2[0]:dim2[1], dim3[0]:dim3[1]]

def get_face_intersection(face1:Face,face2:Face,block1:Block,block2:Block,tol:float=1E-6):
    """Get the index of the intersection between two faces located on two different blocks.

    Three-step approach:
      Step 1 (fast): Same-size faces — try 8 permutations (4 direction + 4 transposed).
      Step 2 (subregion): Different-size faces — find smaller face's corners on larger
             face to identify subregion, then try 8 permutations on subregion.
      Step 3 (fallback): Per-point geometric matching for any remaining cases.

    Args:
        face1 (Face): An exterior face
        face2 (Face): An exterior face from a different block
        block1 (Block): block containing face1
        block2 (Block): block containing face2
        tol (float): matching tolerance

    Returns:
        (Tuple): containing

            - (pandas.DataFrame): dataframe with matches. Columns = i1, j1, k1, i2, j2, k2
            - (List[Face]): any split faces from block 1
            - (List[Face]): any split faces from block 2
    """
    df = pd.DataFrame(columns=['i1','j1','k1','i2','j2','k2'])
    split_faces1 = list()
    split_faces2 = list()

    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]

    lb1 = [I1[0], J1[0], K1[0]]
    ub1 = [I1[1], J1[1], K1[1]]
    lb2 = [I2[0], J2[0], K2[0]]
    ub2 = [I2[1], J2[1], K2[1]]

    n1 = _face_point_count(lb1, ub1)
    n2 = _face_point_count(lb2, ub2)

    # ── Step 1: Same-size fast path (8 permutations: 4 direct + 4 transposed) ──
    if n1 == n2 and n1 >= 4:
        matched, df_result = _try_permutations_with_transpose(
            block1, lb1, ub1, block2, lb2, ub2, tol)
        if matched:
            return df_result, split_faces1, split_faces2

    # ── Step 2: Subregion path (different-size faces) ──
    # Find smaller face's diagonal corners on the larger face to identify
    # the matching subregion, then try 8 permutations on that subregion.
    if n1 != n2 and n1 >= 4 and n2 >= 4:
        if n1 <= n2:
            # face1 is smaller, find its corners on face2
            small_block, small_lb, small_ub = block1, lb1, ub1
            large_block, large_lb, large_ub = block2, lb2, ub2
            small_is_face1 = True
        else:
            # face2 is smaller, find its corners on face1
            small_block, small_lb, small_ub = block2, lb2, ub2
            large_block, large_lb, large_ub = block1, lb1, ub1
            small_is_face1 = False

        # Find smaller face's lb corner on the larger face
        x_lb = small_block.X[small_lb[0], small_lb[1], small_lb[2]]
        y_lb = small_block.Y[small_lb[0], small_lb[1], small_lb[2]]
        z_lb = small_block.Z[small_lb[0], small_lb[1], small_lb[2]]
        corner1 = _find_corner_on_face(large_block, large_lb, large_ub,
                                       x_lb, y_lb, z_lb, tol)

        # Find smaller face's ub corner on the larger face
        x_ub = small_block.X[small_ub[0], small_ub[1], small_ub[2]]
        y_ub = small_block.Y[small_ub[0], small_ub[1], small_ub[2]]
        z_ub = small_block.Z[small_ub[0], small_ub[1], small_ub[2]]
        corner2 = _find_corner_on_face(large_block, large_lb, large_ub,
                                       x_ub, y_ub, z_ub, tol)

        if corner1 is not None and corner2 is not None:
            sub_lb = corner1
            sub_ub = corner2
            n_sub = _face_point_count(sub_lb, sub_ub)
            n_small = _face_point_count(small_lb, small_ub)

            if n_sub == n_small:
                if small_is_face1:
                    # Hold face1 (small) constant, permute subregion on face2 (large)
                    matched, df_result = _try_permutations_with_transpose(
                        small_block, small_lb, small_ub,
                        large_block, sub_lb, sub_ub, tol)
                else:
                    # Hold face2 (small) constant, permute subregion on face1 (large)
                    # df columns will have i1=face2, i2=face1 — swap after
                    matched, df_result = _try_permutations_with_transpose(
                        small_block, small_lb, small_ub,
                        large_block, sub_lb, sub_ub, tol)

                if matched:
                    if not small_is_face1:
                        # Swap: i1↔i2, j1↔j2, k1↔k2 since small was face2
                        df_result = df_result.rename(columns={
                            'i1': '_i2', 'j1': '_j2', 'k1': '_k2',
                            'i2': 'i1', 'j2': 'j1', 'k2': 'k1',
                        }).rename(columns={
                            '_i2': 'i2', '_j2': 'j2', '_k2': 'k2',
                        })
                    df = df_result
                    # Fall through to split face logic below

    # ── Step 3: Fallback — per-point geometric matching ──
    # Handles edge cases where Steps 1-2 don't find a match.
    if len(df) == 0:
        match_location = list()

        X1 = select_multi_dimensional(block1.X, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
        Y1 = select_multi_dimensional(block1.Y, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
        Z1 = select_multi_dimensional(block1.Z, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))

        X2 = select_multi_dimensional(block2.X, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
        Y2 = select_multi_dimensional(block2.Y, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
        Z2 = select_multi_dimensional(block2.Z, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))

        if I1[0] == I1[1]:
            combo = product(range(X1.shape[0]), range(X1.shape[1]))
            for c in combo:
                p, q = c
                x = X1[p,q]; y = Y1[p,q]; z = Z1[p,q]
                block2_match_location = point_match(x, y, z, X2, Y2, Z2, tol)
                if sum(block2_match_location) != -2:
                    p2 = int(block2_match_location[0])
                    q2 = int(block2_match_location[1])
                    if I2[0]==I2[1]:
                        match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
            df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)

        elif J1[0] == J1[1]:
            combo = product(range(X1.shape[0]), range(X1.shape[1]))
            for c in combo:
                p, q = c
                x = X1[p,q]; y = Y1[p,q]; z = Z1[p,q]
                block2_match_location = point_match(x, y, z, X2, Y2, Z2, tol)
                if sum(block2_match_location) != -2:
                    p2 = int(block2_match_location[0])
                    q2 = int(block2_match_location[1])
                    if I2[0]==I2[1]:
                        match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
            df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)

        elif K1[0] == K1[1]:
            combo = product(range(X1.shape[0]), range(X1.shape[1]))
            for c in combo:
                p, q = c
                x = X1[p,q]; y = Y1[p,q]; z = Z1[p,q]
                block2_match_location = point_match(x, y, z, X2, Y2, Z2, tol)
                if sum(block2_match_location) != -2:
                    p2 = int(block2_match_location[0])
                    q2 = int(block2_match_location[1])
                    if I2[0]==I2[1]:
                        match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
            df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)

    # ── Split face checking (applies to Steps 2 and 3) ──
    if len(df)>=4:
        if (__check_edge(df)):
            df = pd.DataFrame()     # If it's an edge
        else:                       # not edge
            # Filter match increasing - This keeps uniqueness
            if I1[0]==I1[1]:
                df = __filter_block_increasing(df,'j1')
                df = __filter_block_increasing(df,'k1')
            elif J1[0]==J1[1]:
                df = __filter_block_increasing(df,'i1')
                df = __filter_block_increasing(df,'k1')
            elif K1[0]==K1[1]:
                df = __filter_block_increasing(df,'i1')
                df = __filter_block_increasing(df,'j1')

            if I2[0]==I2[1]:
                df = __filter_block_increasing(df,'j2')
                df = __filter_block_increasing(df,'k2')
            elif J2[0]==J2[1]:
                df = __filter_block_increasing(df,'i2')
                df = __filter_block_increasing(df,'k2')
            elif K2[0]==K2[1]:
                df = __filter_block_increasing(df,'i2')
                df = __filter_block_increasing(df,'j2')

            # Do a final check after doing all these checks
            if len(df)>=4:       # Greater than 4 because match can occur with simply 4 corners but the interior doesn't match.
                # Check for Split faces
                ## Block 1
                main_face = create_face_from_diagonals(block1,[I1[0],J1[0],K1[0]],[I1[1],J1[1],K1[1]])
                ilb, jlb, klb = df['i1'].min(), df['j1'].min(), df['k1'].min()
                iub, jub, kub = df['i1'].max(), df['j1'].max(), df['k1'].max()
                if int(ilb==iub) + int(jlb==jub) + int(klb==kub)==1:
                    split_faces1 = split_face(main_face,block1,ilb=ilb,iub=iub,jlb=jlb,jub=jub,klb=klb,kub=kub)
                    [s.set_block_index(face1.blockIndex) for s in split_faces1]
                    [s.set_face_id(face1.id) for s in split_faces1]

                ## Block 2
                main_face = create_face_from_diagonals(block2,[I2[0],J2[0],K2[0]],[I2[1],J2[1],K2[1]])
                ilb, jlb, klb = df['i2'].min(), df['j2'].min(), df['k2'].min()
                iub, jub, kub = df['i2'].max(), df['j2'].max(), df['k2'].max()
                if int(ilb==iub) + int(jlb==jub) + int(klb==kub)==1:
                    split_faces2 = split_face(main_face,block2,ilb=ilb,iub=iub,jlb=jlb,jub=jub,klb=klb,kub=kub)
                    [s.set_block_index(face2.blockIndex) for s in split_faces2]
                    [s.set_face_id(face2.id) for s in split_faces2]

    else:
        df = pd.DataFrame() # set df to empty dataframe
    return df, split_faces1, split_faces2

def __filter_block_increasing(df:pd.DataFrame,key1:str):
    """Filters dataframe results of get_face_intersection to make sure both key1 is increasing. 
        When searching through a plot3D we check based on the planes e.g. const i, j, or k 
        values will be removed if they are not

    Args:
        df (pd.DataFrame): DataFrame containing matching points 
        key1 (str): column that you want to be in increasing order 

    Returns:
        pd.DataFrame: sorted dataframe
    """        

    '''
        Sometimes there's a match on 2 edges and we do not want to keep that 
            | face1 | face2 | face1  | 
        Above shows face 2 touching face 1 at 2 edges. this is not a match. 
    '''
    if len(df)==0:
        return df
    
    key1_vals = list(df[key1].unique()) # get the unique values 
    key1_vals.sort()
    key1_vals_to_use = list()

    if len(key1_vals)<=1:
        return pd.DataFrame() # Returning an empty dataframe. This solves the condition where you have edge matching 

    for i in range(len(key1_vals)-1):
        if (key1_vals[i+1] - key1_vals[i])==1: # Remove
            key1_vals_to_use.append(key1_vals[i])
    # Look backwards 
    if (key1_vals[-1] - key1_vals[-2])==1: # Remove
            key1_vals_to_use.append(key1_vals[-1])
    df = df[df[key1].isin(key1_vals_to_use)]        
    return df

def __check_edge(df:pd.DataFrame):
    """ Check if the results of get_face_intersection is an edge instead of a face.  
        if it's an edge then both intersecting blocks are connected by an edge on both blocks 

    Args:
        df (pd.DataFrame): dataframe containing columns i1, j1, k1, i2, j2, k2

    Returns:
        boolean: True = It is an edge, False = not edge 
    """
    face1_diagonal = [(df['i1'].min(),df['j1'].min(),df['k1'].min()),(df['i1'].max(),df['j1'].max(),df['k1'].max()) ]
    face2_diagonal = [(df['i2'].min(),df['j2'].min(),df['k2'].min()), (df['i2'].max(),df['j2'].max(),df['k2'].max())]
    edge1 = face1_diagonal[0]
    edge2 = face1_diagonal[1]
    edge_matches = 0 
    for i in range(3):
        if edge1[i]==edge2[i]:
            edge_matches+=1
    if edge_matches<2:
        return False
    else:
        return True
    
def combinations_of_nearest_blocks(blocks:List[Block],nearest_nblocks:int=4):
    """Returns the indices of the nearest 6 blocks based on their centroid

    Args:
        block (Block): block you are interested in
        blocks (List[Block]): list of all your blocks

    Returns:
        List[Tuple[int,int]]: combinations of nearest blocks 
    """
    # Pick a block get centroid of all outer faces        
    centroids = np.array([(b.cx,b.cy,b.cz) for b in blocks])
    distance_matrix = np.zeros((centroids.shape[0],centroids.shape[0]))+10000
    # Build a matrix 
    for i in range(centroids.shape[0]):
        for j in range(centroids.shape[0]):
            if i!=j:
                dx = centroids[i,0]-centroids[j,0]
                dy = centroids[i,1]-centroids[j,1]
                dz = centroids[i,2]-centroids[j,2]
                distance_matrix[i,j] = np.sqrt(dx*dx+dy*dy+dz*dz)
    
    # Now that we have this matrix, we sort the distances by rows and pick the closest 8 blocks, can use 4 but 8 might be safer
    new_combos = list()
    for i in range(len(blocks)): # For block i
        indices = np.argsort(distance_matrix[i,:])
        for j in indices[:nearest_nblocks]:
            if distance_matrix[i,j] < 10000:
                new_combos.append((i,j))
    return new_combos
    
def connectivity_fast(blocks:List[Block]):
    """Reduces the size of the blocks by a factor of the minimum gcd. This speeds up finding the connectivity 

    Args:
        blocks (List[Block]): Lists of blocks you want to find the connectivity for

    Returns:
        (List[Dict]): All matching faces formatted as a list of { 'block1': {'block_index', 'lb', 'ub'} }
        (List[Dict]): All exterior surfaces formatted as a list of { 'block_index', 'lb', 'ub', 'id' }
        
    """
    gcd_array = list()
    # Find the gcd of all the blocks 
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))        
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    print(f"gcd to use {gcd_to_use}")
    new_blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)

    # Find Connectivity 
    face_matches, outer_faces_formatted = connectivity(new_blocks)
    # scale it up
    for i in range(len(face_matches)):
        for side in ['block1', 'block2']:
            face_matches[i][side]['lb'] = [x * gcd_to_use for x in face_matches[i][side]['lb']]
            face_matches[i][side]['ub'] = [x * gcd_to_use for x in face_matches[i][side]['ub']]
    for j in range(len(outer_faces_formatted)):
        outer_faces_formatted[j]['lb'] = [x * gcd_to_use for x in outer_faces_formatted[j]['lb']]
        outer_faces_formatted[j]['ub'] = [x * gcd_to_use for x in outer_faces_formatted[j]['ub']]
    return face_matches, outer_faces_formatted

def connectivity(blocks:List[Block]):
    """Returns a dictionary outlining the connectivity of the blocks along with any exterior surfaces 

    Args:
        blocks (List[Block]): List of all blocks in multi-block plot3d mesh

    Returns:
        (List[Dict]): All matching faces formatted as a list of { 'block1': {'block_index', 'lb', 'ub'} }
        (List[Dict]): All exterior surfaces formatted as a list of { 'block_index', 'lb', 'ub', 'id' }

    """

    outer_faces = list()      
    face_matches = list()
    matches_to_remove = list()
    temp = [get_outer_faces(b) for b in blocks]
    block_outer_faces = [t[0] for t in temp]
    combos = combinations_of_nearest_blocks(blocks,6) # Find the 6 nearest Blocks and search through all that. 
    # df_matches, blocki_outerfaces, blockj_outerfaces = find_matching_blocks(blocks[0],blocks[28],
    #                                                                         [block_outer_faces[0][-1]],
    #                                                                         [block_outer_faces[28][1]],1E-12)    # This function finds partial matches between blocks

    t = trange(len(combos))    
    for indx in t:     # block i        
        i,j = combos[indx]
        t.set_description(f"Checking connections block {i} with {j}")
        # Takes 2 blocks, gets the matching faces exterior faces of both blocks 
        df_matches, blocki_outerfaces, blockj_outerfaces = find_matching_blocks(blocks[i],blocks[j],block_outer_faces[i],block_outer_faces[j])    # This function finds partial matches between blocks
        [o.set_block_index(i) for o in blocki_outerfaces]
        [o.set_block_index(j) for o in blockj_outerfaces]
        block_outer_faces[i] = blocki_outerfaces
        block_outer_faces[j] = blockj_outerfaces
        # Update connectivity for blocks with matching faces 
        if (len(df_matches)>0):
            for df in df_matches:
                matches_to_remove.append(create_face_from_diagonals(blocks[i],
                    [df['i1'].min(),df['j1'].min(),df['k1'].min()],
                    [df['i1'].max(),df['j1'].max(),df['k1'].max()]))
                matches_to_remove[-1].set_block_index(i)

                matches_to_remove.append(create_face_from_diagonals(blocks[j],
                    [df['i2'].min(),df['j2'].min(),df['k2'].min()],
                    [df['i2'].max(),df['j2'].max(),df['k2'].max()]))
                matches_to_remove[-1].set_block_index(j)
                
                face1 = matches_to_remove[-2]
                face2 = matches_to_remove[-1]

                # Derive lb/ub directly from the DataFrame traversal order.
                # iloc[0] = first point in face1's traversal, iloc[-1] = last.
                # This correctly handles cross-axis matches where corner-based
                # matching would lose the traversal order.
                lb1_out = [int(df.iloc[0]['i1']), int(df.iloc[0]['j1']), int(df.iloc[0]['k1'])]
                ub1_out = [int(df.iloc[-1]['i1']), int(df.iloc[-1]['j1']), int(df.iloc[-1]['k1'])]
                lb2_out = [int(df.iloc[0]['i2']), int(df.iloc[0]['j2']), int(df.iloc[0]['k2'])]
                ub2_out = [int(df.iloc[-1]['i2']), int(df.iloc[-1]['j2']), int(df.iloc[-1]['k2'])]

                # Compute the orientation vector mapping face1's axes
                # to face2's axes.  orientation[d] is the 1-indexed
                # face2 axis corresponding to face1's axis d.
                # Direction is already encoded in lb/ub.
                orientation = _compute_orientation(df, lb1_out, ub1_out)

                temp = {
                    'block1': {
                        'block_index': i,
                        'lb': lb1_out,
                        'ub': ub1_out,
                        'id': face1.id
                    },
                    'block2': {
                        'block_index': j,
                        'lb': lb2_out,
                        'ub': ub2_out,
                        'id': face2.id
                    },
                    'orientation': orientation,
                    'match': df
                }
                face_matches.append(temp)

    # Update Outer Faces 
    [outer_faces.extend(o) for o in block_outer_faces] # all the outer faces 
    outer_faces = list(set(outer_faces))    # Get most unique
        
    outer_faces = [o for o in outer_faces if o not in matches_to_remove]
    # Remove any outer faces that may have been found by mistake
    # Check I,J,K if J and K are the same with another outer face, select the face with shorter I 
    outer_faces_to_remove = list() 
    for i in range(len(blocks)):
        block_outerfaces = [o for o in outer_faces if o.BlockIndex == i]
        for o in block_outerfaces:
            IJK = np.array([o.IMIN,o.JMIN,o.KMIN,o.IMAX,o.JMAX,o.KMAX])
            for o2 in block_outerfaces:
                IJK2 = np.array([o2.IMIN,o2.JMIN,o2.KMIN,o2.IMAX,o2.JMAX,o2.KMAX])
                if sum((IJK-IJK2)==0) == 5: # [0,0,0,40,100,0] (outer) [0,0,0,56,100,0] (outer) -> remove the longer face
                    if (o2.diagonal_length>o.diagonal_length):
                        outer_faces_to_remove.append(o2)
                    else:
                        outer_faces_to_remove.append(o)

    outer_faces = [o for o in outer_faces if o not in outer_faces_to_remove]


    # Find self-matches: Do any faces of, for example, block1 match another face in block 1
    for i in range(len(blocks)):
        _,self_matches = get_outer_faces(blocks[i]) 
        for match in self_matches: # Append to face matches
            face_matches.append({'block1':{
                                            'block_index':i,
                                            'lb':[int(match[0].I.min()),int(match[0].J.min()),int(match[0].K.min())],
                                            'ub':[int(match[0].I.max()),int(match[0].J.max()),int(match[0].K.max())]
                                        },
                                    'block2':{
                                            'block_index':i,
                                            'lb':[int(match[1].I.min()),int(match[1].J.min()),int(match[1].K.min())],
                                            'ub':[int(match[1].I.max()),int(match[1].J.max()),int(match[1].K.max())]
                                        },
                                    'match':pd.DataFrame([{
                                            'block_index':i,
                                            'lb':[int(match[0].I.min()),int(match[0].J.min()),int(match[0].K.min())],
                                            'ub':[int(match[0].I.max()),int(match[0].J.max()),int(match[0].K.max())]
                                        },{
                                            'block_index':i,
                                            'lb':[int(match[1].I.min()),int(match[1].J.min()),int(match[1].K.min())],
                                            'ub':[int(match[1].I.max()),int(match[1].J.max()),int(match[1].K.max())]
                                        }])
                                    })

    # Update the outer faces
    outer_faces_formatted = list() # This will contain 
    id = 1 
    for face in outer_faces:
        outer_faces_formatted.append({
                            'lb':[int(min(face.I)), int(min(face.J)), int(min(face.K))],
                            'ub':[int(max(face.I)), int(max(face.J)), int(max(face.K))],
                            'id':id, 'block_index':face.BlockIndex })
        id += 1

    return face_matches, outer_faces_formatted  

def face_matches_to_dict(face1:Face, face2:Face,block1:Block,block2:Block):
    """Makes sure the diagonal of face 1 match the diagonal of face 2

    Args:
        face1 (Face): Face 1 with block index 
        face2 (Face): Face 2 with block index
        block1 (Block): Block 1 
        block2 (Block): Block 2 
    
    Returns:
        (dict): dictionary describing the corner matches 

    """
    match = {
            'block1':{
                            'block_index':face1.BlockIndex,
                            'lb':[-1,-1,-1],  # Lower Corner
                            'ub':[-1,-1,-1],   # Upper Corner
                            'id':face1.id
                        },
                'block2':{
                            'block_index':face2.BlockIndex,
                            'lb':[-1,-1,-1],  # Lower Corner
                            'ub':[-1,-1,-1],   # Upper Corner
                            'id':face2.id
                        }
                }
            
    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]

    # Search for corners     
    x1_l = block1.X[I1[0],J1[0],K1[0]]    # lower corner of block 1 
    y1_l = block1.Y[I1[0],J1[0],K1[0]]
    z1_l = block1.Z[I1[0],J1[0],K1[0]]
    # Matches which corner in block 2 
    search_results = list()
    for p in I2: 
        for q in J2:
            for r in K2:
                x2 = block2.X[p,q,r]  
                y2 = block2.Y[p,q,r]
                z2 = block2.Z[p,q,r]
                dx = x2-x1_l; dy = y2-y1_l; dz = z2 -z1_l                
                search_results.append({'I':p,'J':q,'K':r,'d':math.sqrt(dx*dx + dy*dy + dz*dz)})    
    df = pd.DataFrame(search_results)
    df = df.sort_values(by=['d'])
    match['block1']['lb'] = [face1.IMIN, face1.JMIN, face1.KMIN]
    match['block2']['lb'] = [int(df.iloc[0]['I']), int(df.iloc[0]['J']), int(df.iloc[0]['K'])]
    
    # Search for corners     
    x1_u = block1.X[I1[1],J1[1],K1[1]]    # lower corner of block 1 
    y1_u = block1.Y[I1[1],J1[1],K1[1]]
    z1_u = block1.Z[I1[1],J1[1],K1[1]]
    # Matches which corner in block 2 
    search_results = list()
    for p in I2: 
        for q in J2:
            for r in K2:
                x2 = block2.X[p,q,r]  
                y2 = block2.Y[p,q,r]
                z2 = block2.Z[p,q,r]
                dx = x2-x1_u; dy = y2-y1_u; dz = z2 -z1_u
                search_results.append({'I':p,'J':q,'K':r,'d':math.sqrt(dx*dx + dy*dy + dz*dz)})
    df = pd.DataFrame(search_results)
    df = df.sort_values(by=['d'])
    match['block1']['ub'] = [face1.IMAX, face1.JMAX, face1.KMAX]
    match['block2']['ub'] = [int(df.iloc[0]['I']), int(df.iloc[0]['J']), int(df.iloc[0]['K'])]
    return match


