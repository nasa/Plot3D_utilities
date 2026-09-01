"""Connectivity detection for multi-block structured (Plot3D) meshes.

The algorithm runs in three phases:

1. **Phase 1 -- Candidate pairing**: Block pairs whose axis-aligned bounding
   boxes (AABBs) overlap are identified by :func:`candidate_neighbor_pairs`.
   Only these pairs proceed to expensive point matching.

2. **Phase 2 -- Face matching**: For each candidate pair,
   :func:`find_matching_blocks` compares every outer face of block *i* against
   every outer face of block *j*.  Each comparison
   (:func:`get_face_intersection`) tries three strategies in order:

   * *Step 1 (same-size fast path)* -- if both faces have the same number of
     points, try all 8 permutations via
     :func:`_try_permutations_with_transpose`.
   * *Step 2 (sub-region)* -- for different-size faces, locate the smaller
     face's diagonal corners on the larger face to extract a sub-region, then
     try the same 8 permutations.
   * *Step 3 (fallback)* -- per-point geometric matching for any remaining
     cases.

3. **Phase 3 -- Fresh-face validation**: Outer faces left unmatched after
   Phase 2 (because a partner block's face pool was consumed by an earlier
   pair) are re-checked against fresh (un-consumed) outer faces of
   neighbouring blocks.

Orientation system
------------------
When two faces match, the module records *how* their local (u, v) axes
relate.  Because each of the two varying axes can independently be reversed
and the two axes can be swapped, there are 2 x 2 x 2 = **8 distinct
orientations**, encoded as 2x2 signed permutation matrices in
:data:`PERMUTATION_MATRICES`.

The index into that array is computed with the bit formula::

    index = u_reversed | (v_reversed << 1) | (swapped << 2)

where *u_reversed* / *v_reversed* are booleans indicating whether the
traversal direction along that axis is flipped between the two faces, and
*swapped* indicates whether the u and v axes of face 1 map to v and u of
face 2 (a cross-plane match).  Indices 0--3 are the four direct (non-swapped)
permutations; indices 4--7 are the four transposed (swapped) permutations.

JSON export convention
~~~~~~~~~~~~~~~~~~~~~~
The exported ``permutation_index`` follows the **diagonal (lb/ub)** convention:

- **In-plane** matches (perm 0--3): ``permutation_index`` is **-1** because the
  traversal direction is fully encoded in block2's lb/ub ordering.
- **Cross-plane** matches (perm 4--7): ``permutation_index`` is the actual index
  because lb/ub alone cannot represent an axis swap.

The ``permutation_matrix`` field always contains the actual 2x2 matrix.
"""

from .block import Block
from .blockfunctions import reduce_blocks, compute_min_gcd, scale_face_bounds, constant_axis as _constant_axis
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


PERMUTATION_MATRICES = np.array([
    [[ 1,  0], [ 0,  1]],   # 0: identity
    [[-1,  0], [ 0,  1]],   # 1: u reversed
    [[ 1,  0], [ 0, -1]],   # 2: v reversed
    [[-1,  0], [ 0, -1]],   # 3: both reversed
    [[ 0,  1], [ 1,  0]],   # 4: swapped
    [[ 0, -1], [ 1,  0]],   # 5: swap + u reversed
    [[ 0,  1], [-1,  0]],   # 6: swap + v reversed
    [[ 0, -1], [-1,  0]],   # 7: swap + both reversed
], dtype=np.int8)
"""8 canonical 2x2 signed permutation matrices.

Bit encoding: ``index = u_reversed | (v_reversed << 1) | (swapped << 2)``
"""


def _orient_vec_to_permutation(orient_vec: list, lb1: list, ub1: list,
                               lb2: list, ub2: list) -> Tuple[int, str]:
    """Convert an orientation vector to a ``(permutation_index, plane)`` tuple.

    This function translates the high-level 3-axis orientation vector produced
    by :func:`_compute_orientation` into the compact representation used by
    :data:`PERMUTATION_MATRICES`.  It works as follows:

    1. Identify the constant axis on each face to determine the two *varying*
       axes for face1 (``u1, v1``) and the two axes they map to on face2
       (``u2, v2``).
    2. Detect whether the axes are **swapped** (i.e. face1's u maps to a
       different positional axis on face2 than its own).
    3. Compare the traversal direction (sign of ``ub - lb``) along each
       varying axis between the two faces to detect **reversal**.
    4. Combine the three boolean flags into a permutation index using the
       :data:`PERMUTATION_MATRICES` bit-encoding formula::

           index = int(u_reversed) | (int(v_reversed) << 1) | (int(swapped) << 2)

       Indices 0--3 are non-swapped; indices 4--7 are swapped (cross-axis).

    The *plane* string is ``'in-plane'`` when both faces share the same
    constant axis (e.g. both are I-constant) and ``'cross-plane'`` otherwise
    (e.g. one is I-constant and the other is K-constant).

    Parameters
    ----------
    orient_vec : list of int
        1-indexed face2 axis per face1 axis, as returned by
        :func:`_compute_orientation`.
    lb1, ub1 : list of int
        Lower-bound and upper-bound ``[i, j, k]`` corners of face1.
    lb2, ub2 : list of int
        Lower-bound and upper-bound ``[i, j, k]`` corners of face2.

    Returns
    -------
    tuple of (int, str)
        ``(permutation_index, plane)`` where *permutation_index* is in
        ``0..7`` and *plane* is ``'in-plane'`` or ``'cross-plane'``.
    """
    ca1 = _constant_axis(lb1, ub1)
    ca2 = _constant_axis(lb2, ub2)
    plane = 'in-plane' if ca1 == ca2 else 'cross-plane'

    varying1 = [d for d in range(3) if d != ca1]
    if len(varying1) != 2:
        return 0, plane

    u1, v1 = varying1[0], varying1[1]
    u2 = orient_vec[u1] - 1  # 0-indexed
    v2 = orient_vec[v1] - 1

    swapped = (u2 != u1) or (v2 != v1)

    step1 = lambda d: 1 if ub1[d] >= lb1[d] else -1
    step2 = lambda d: 1 if ub2[d] >= lb2[d] else -1

    u_reversed = step1(u1) != step2(u2)
    v_reversed = step1(v1) != step2(v2)

    perm_idx = int(u_reversed) | (int(v_reversed) << 1) | (int(swapped) << 2)
    return perm_idx, plane


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
    """Compute the orientation vector that maps face1's axes to face2's axes.

    Returns a 3-element list ``orientation`` where ``orientation[d]`` is the
    **1-indexed** face2 axis that corresponds to face1 axis *d* (1=I, 2=J,
    3=K).  For example, ``[3, 1, 2]`` means face1-I maps to face2-K, face1-J
    maps to face2-I, and face1-K maps to face2-J.

    Direction (forward / reversed) is **not** encoded in the orientation
    vector -- it is already captured by the lb/ub traversal order of each
    face.

    The mapping is derived from the point-match DataFrame *df* by inspecting
    which face2 index column changes when a single face1 axis is incremented:

    * Row 0 vs Row 1 reveals the face2 axis for face1's *inner* varying axis.
    * Row 0 vs Row ``inner_size`` reveals the axis for face1's *outer* varying
      axis.
    * The remaining (third) face2 axis is assigned to face1's constant axis.

    Parameters
    ----------
    df : pandas.DataFrame
        Point-match table with columns ``i1, j1, k1, i2, j2, k2``.  Rows
        are ordered by face1's I-J-K nested traversal.
    lb1, ub1 : list of int
        Lower-bound and upper-bound [i, j, k] corners of face1 (determines
        which axis is constant and the traversal order).

    Returns
    -------
    list of int
        ``[o_I, o_J, o_K]`` -- 1-indexed face2 axis per face1 axis.
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
    """Test all 8 permutations (4 direct + 4 transposed) to match face1 to face2.

    A structured face has one constant axis and two varying axes, giving four
    possible traversal-direction combinations (forward/reverse for each
    varying axis).  These are the **4 direct permutations**, generated by
    :func:`_generate_face2_permutations`.

    For cross-plane matches (where the two faces have *different* constant
    axes), the outer and inner loop dimensions of face2 must be transposed
    before comparison.  Applying the transpose to each of the 4 direct
    permutations yields the **4 transposed permutations**, for a total of 8.

    The function extracts face1's points once and then iterates over the 4
    ``(lb, ub)`` direction variants of face2.  For each variant it performs:

    * **Direct comparison** -- flatten both grids to (N, 3) and check if
      every point pair is within *tol*.
    * **Transposed comparison** -- reshape face2's points to its 2-D grid
      shape, transpose the two spatial axes, re-flatten, and compare.

    The first permutation that produces a full match (all point distances
    below *tol*) is accepted and the function returns immediately.

    Parameters
    ----------
    block1 : Block
        Block containing face1.
    lb1, ub1 : list of int
        Lower/upper ``[i, j, k]`` corners of face1.
    block2 : Block
        Block containing face2.
    lb2, ub2 : list of int
        Lower/upper ``[i, j, k]`` corners of face2.
    tol : float
        Point-matching tolerance (Euclidean distance).

    Returns
    -------
    tuple of (bool, pandas.DataFrame)
        ``(True, df)`` on success, where *df* has columns
        ``i1, j1, k1, i2, j2, k2`` mapping every matched point.
        ``(False, empty_df)`` if no permutation produces a match.
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
    
    return T[dim1[0]:dim1[1]+1, dim2[0]:dim2[1]+1, dim3[0]:dim3[1]+1]

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

            # Reject matches where matched points don't cover the full
            # matched sub-face area.  Two blocks that share only edges
            # (e.g. O-grid SS and PS sharing LE/TE lines) can pass the
            # edge check above because the matched points span two
            # separate edges, making the diagonal look like a face.
            # Verify that matched point count == expected sub-face area.
            if len(df) >= 4:
                ilb1, jlb1, klb1 = int(df['i1'].min()), int(df['j1'].min()), int(df['k1'].min())
                iub1, jub1, kub1 = int(df['i1'].max()), int(df['j1'].max()), int(df['k1'].max())
                matched_area = _face_point_count([ilb1, jlb1, klb1], [iub1, jub1, kub1])
                if matched_area > 0 and len(df) < matched_area:
                    df = pd.DataFrame()  # Not a face — only partial (edge) coverage

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

    # With only 2 unique values, contiguity is trivially satisfied — keep all.
    # This handles small faces (e.g. 2 nodes wide after GCD reduction) matching
    # a large face where the matching indices may span a wide gap (e.g. [0, 113])
    # but are still a valid match. The __check_edge() call upstream has already
    # verified this isn't a degenerate edge.
    if len(key1_vals) == 2:
        return df

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
    
def candidate_neighbor_pairs(blocks:List[Block], tol:float=1e-6):
    """Returns candidate block pairs whose AABBs overlap or nearly touch.

    This replaces the former centroid-distance approach which only considered
    the N nearest blocks and could miss neighbours for L-shaped or elongated
    geometries.  AABB overlap is both more robust and more correct.

    Args:
        blocks (List[Block]): list of all your blocks
        tol (float): AABB expansion tolerance

    Returns:
        List[Tuple[int,int]]: candidate (i, j) pairs with i < j
    """
    n = len(blocks)
    # Precompute AABBs: [xmin, xmax, ymin, ymax, zmin, zmax]
    aabbs = np.empty((n, 6), dtype=np.float64)
    for i, b in enumerate(blocks):
        aabbs[i, 0] = b.X.min()
        aabbs[i, 1] = b.X.max()
        aabbs[i, 2] = b.Y.min()
        aabbs[i, 3] = b.Y.max()
        aabbs[i, 4] = b.Z.min()
        aabbs[i, 5] = b.Z.max()

    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = aabbs[i], aabbs[j]
            if (a[1] + tol >= b[0] and b[1] + tol >= a[0] and
                a[3] + tol >= b[2] and b[3] + tol >= a[2] and
                a[5] + tol >= b[4] and b[5] + tol >= a[4]):
                pairs.append((i, j))
    return pairs


def combinations_of_nearest_blocks(blocks:List[Block],nearest_nblocks:int=4):
    """Returns the indices of the nearest N blocks based on their centroid.

    .. deprecated::
        Use :func:`candidate_neighbor_pairs` instead for AABB-based pairing.

    Args:
        blocks (List[Block]): list of all your blocks
        nearest_nblocks (int): number of nearest blocks to consider

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

def _phase3_overlaps_existing(bi, lb1, ub1, bj, lb2, ub2, face_matches):
    """Return True if a candidate Phase-3 match's region overlaps any
    existing ``face_matches`` record on the same block pair.

    Overlap guard used inside :func:`connectivity` Phase 3 to prevent
    re-emission of regions that Phase 2 already claimed. On meshes with
    an O-grid seam, ``__filter_block_increasing`` drops seam-adjacent
    points during Phase 2 and leaves a residue that Phase 3 then
    re-matches against the fresh neighbor face pool, producing a
    wrap-around record whose bbox already overlaps two clean Phase 2
    sub-faces on the same block pair. The per-face dedup key in Phase 3
    (``(bi, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)``) does not catch this
    because it only tracks which ``block1`` faces have been processed
    in Phase 3.

    The check is per-dimension on both sides, with ``lb``/``ub``
    normalised to (min, max) per axis so direction does not matter.
    """
    def _ranges_overlap(a_lo, a_hi, b_lo, b_hi):
        return not (a_hi < b_lo or b_hi < a_lo)

    def _normalise(lb, ub):
        return ([min(lb[d], ub[d]) for d in range(len(lb))],
                [max(lb[d], ub[d]) for d in range(len(lb))])

    cand_lb1, cand_ub1 = _normalise(lb1, ub1)
    cand_lb2, cand_ub2 = _normalise(lb2, ub2)

    for m in face_matches:
        m1 = m.get('block1', {})
        m2 = m.get('block2', {})
        mbi = m1.get('block_index')
        mbj = m2.get('block_index')
        if mbi is None or mbj is None:
            continue
        if mbi == bi and mbj == bj:
            mlb1, mub1 = _normalise(m1['lb'], m1['ub'])
            mlb2, mub2 = _normalise(m2['lb'], m2['ub'])
        elif mbi == bj and mbj == bi:
            mlb1, mub1 = _normalise(m2['lb'], m2['ub'])
            mlb2, mub2 = _normalise(m1['lb'], m1['ub'])
        else:
            continue
        side1 = all(
            _ranges_overlap(cand_lb1[d], cand_ub1[d], mlb1[d], mub1[d])
            for d in range(3)
        )
        side2 = all(
            _ranges_overlap(cand_lb2[d], cand_ub2[d], mlb2[d], mub2[d])
            for d in range(3)
        )
        if side1 and side2:
            return True
    return False

def connectivity_fast(blocks:List[Block], use_minmax:bool=False):
    """Find connectivity by GCD-reducing blocks first for speed.

    Computes the minimum GCD across all block dimensions, reduces all blocks
    uniformly, runs :func:`connectivity`, then scales bounds back up.

    Face match dicts follow the diagonal (lb/ub) convention:
    ``permutation_index`` is **-1** for in-plane, actual index for cross-plane.

    Args:
        blocks (List[Block]): List of blocks to find connectivity for.
        use_minmax (bool): If True, normalise lb/ub to strict min/max
            order (IMIN,JMIN,KMIN → IMAX,JMAX,KMAX) and recompute the
            permutation matrix accordingly.  Default is False (traversal
            order).

    Returns:
        (List[Dict]): Face matches with orientation info.
        (List[Dict]): Outer (non-connected) faces.
    """
    gcd_to_use = compute_min_gcd(blocks)
    print(f"gcd to use {gcd_to_use}")
    new_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # Find Connectivity
    face_matches, outer_faces_formatted = connectivity(new_blocks)
    # scale it up
    scale_face_bounds(face_matches, gcd_to_use)
    scale_face_bounds(outer_faces_formatted, gcd_to_use)
    if use_minmax:
        face_matches = normalize_face_matches(face_matches)
    return face_matches, outer_faces_formatted

def connectivity(blocks:List[Block]):
    """Returns a dictionary outlining the connectivity of the blocks along with any exterior surfaces.

    Each face match dict includes an ``orientation`` sub-dict with:

    - ``permutation_index``: **-1** for in-plane matches (direction encoded in
      lb/ub), or the actual index (4-7) for cross-plane matches.
    - ``plane``: ``'in-plane'`` or ``'cross-plane'``.
    - ``permutation_matrix``: the actual 2x2 signed permutation matrix.

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
    combos = candidate_neighbor_pairs(blocks) # AABB overlap pairs (i < j)

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
                perm_idx, plane = _orient_vec_to_permutation(
                    orientation, lb1_out, ub1_out, lb2_out, ub2_out)

                # In-plane: direction fully encoded in lb/ub → export -1.
                # Cross-plane: lb/ub can't encode axis swap → export actual index.
                export_perm = -1 if plane == 'in-plane' else perm_idx

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
                    'orientation': {
                        'permutation_index': export_perm,
                        'plane': plane,
                        'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
                    },
                    'match': df
                }
                face_matches.append(temp)

    # ===== PHASE 3: Fresh-face validation for remaining outer faces =====
    # Some outer faces remain unmatched because the matching block's face pool
    # was consumed by an earlier combo in Phase 2.  Re-check each remaining
    # outer face against *fresh* (un-consumed) outer faces of overlapping blocks.
    print("Phase 3: Fresh-face validation...")
    neighbors = [[] for _ in range(len(blocks))]
    for i_combo, j_combo in combos:
        neighbors[i_combo].append(j_combo)
        neighbors[j_combo].append(i_combo)

    fresh_all = [get_outer_faces(b)[0] for b in blocks]
    phase3_keys = set()
    phase3_partial_matches: dict = {}
    phase3_count = 0

    for bi in range(len(blocks)):
        for face in block_outer_faces[bi]:
            face_key = (bi, face.IMIN, face.JMIN, face.KMIN, face.IMAX, face.JMAX, face.KMAX)
            if face_key in phase3_keys:
                continue
            for bj in neighbors[bi]:
                for fresh_face in fresh_all[bj]:
                    df, _, _ = get_face_intersection(face, fresh_face, blocks[bi], blocks[bj])
                    if len(df) < 4:
                        continue
                    if __check_edge(df):
                        continue

                    # Apply axis filter
                    I1 = [face.IMIN, face.IMAX]
                    J1 = [face.JMIN, face.JMAX]
                    K1 = [face.KMIN, face.KMAX]
                    I2 = [fresh_face.IMIN, fresh_face.IMAX]
                    J2 = [fresh_face.JMIN, fresh_face.JMAX]
                    K2 = [fresh_face.KMIN, fresh_face.KMAX]

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

                    if len(df) < 4:
                        continue

                    # Build match record
                    lb1_out = [int(df.iloc[0]['i1']), int(df.iloc[0]['j1']), int(df.iloc[0]['k1'])]
                    ub1_out = [int(df.iloc[-1]['i1']), int(df.iloc[-1]['j1']), int(df.iloc[-1]['k1'])]
                    lb2_out = [int(df.iloc[0]['i2']), int(df.iloc[0]['j2']), int(df.iloc[0]['k2'])]
                    ub2_out = [int(df.iloc[-1]['i2']), int(df.iloc[-1]['j2']), int(df.iloc[-1]['k2'])]

                    # Skip if this region overlaps an existing Phase 2
                    # match on the same block pair — happens around
                    # O-grid seams where __filter_block_increasing drops
                    # seam-adjacent points and leaves a residue that
                    # Phase 3 re-matches against the fresh neighbor face.
                    if _phase3_overlaps_existing(
                        bi, lb1_out, ub1_out, bj, lb2_out, ub2_out,
                        face_matches,
                    ):
                        continue

                    orientation = _compute_orientation(df, lb1_out, ub1_out)
                    perm_idx, plane = _orient_vec_to_permutation(
                        orientation, lb1_out, ub1_out, lb2_out, ub2_out)
                    export_perm = -1 if plane == 'in-plane' else perm_idx

                    face_matches.append({
                        'block1': {
                            'block_index': bi,
                            'lb': lb1_out,
                            'ub': ub1_out,
                        },
                        'block2': {
                            'block_index': bj,
                            'lb': lb2_out,
                            'ub': ub2_out,
                        },
                        'orientation': {
                            'permutation_index': export_perm,
                            'plane': plane,
                            'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
                        },
                        'match': df,
                    })
                    phase3_keys.add(face_key)
                    # Track the matched sub-region so we can split the
                    # original face and keep the unmatched leftover as
                    # an outer face. Without this, partial Phase 3
                    # matches silently discard the rest of the face,
                    # which then can't participate in periodicity.
                    phase3_partial_matches.setdefault(face_key, []).append(
                        (lb1_out, ub1_out)
                    )
                    phase3_count += 1
                    # Don't break — continue checking other neighbors' faces
                    # for split-face matches

    print(f"  Phase 3 done ({phase3_count} new matches)")

    # Remove Phase 3 matched faces from outer faces, but for partial
    # matches (matched sub-region smaller than the original face)
    # add the unmatched leftover splits back so they can still be
    # detected as outer faces / periodic candidates.
    for bi in range(len(blocks)):
        new_outer = []
        for f in block_outer_faces[bi]:
            fk = (bi, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)
            if fk not in phase3_keys:
                new_outer.append(f)
                continue
            # face was matched in Phase 3 — emit leftovers from each
            # partial match region
            partials = phase3_partial_matches.get(fk, [])
            leftovers = [f]
            for lb_m, ub_m in partials:
                next_leftovers = []
                for left in leftovers:
                    # Skip if match region doesn't overlap this leftover
                    if (ub_m[0] < left.IMIN or lb_m[0] > left.IMAX or
                        ub_m[1] < left.JMIN or lb_m[1] > left.JMAX or
                        ub_m[2] < left.KMIN or lb_m[2] > left.KMAX):
                        next_leftovers.append(left)
                        continue
                    # Clip match region to leftover, then split
                    ilb = max(lb_m[0], left.IMIN); iub = min(ub_m[0], left.IMAX)
                    jlb = max(lb_m[1], left.JMIN); jub = min(ub_m[1], left.JMAX)
                    klb = max(lb_m[2], left.KMIN); kub = min(ub_m[2], left.KMAX)
                    splits = split_face(left, blocks[bi],
                                        ilb=ilb, jlb=jlb, klb=klb,
                                        iub=iub, jub=jub, kub=kub)
                    for s in splits:
                        s.set_block_index(bi)
                    next_leftovers.extend(splits)
                leftovers = next_leftovers
            new_outer.extend(leftovers)
        block_outer_faces[bi] = new_outer

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


def normalize_face_matches(face_matches: list) -> list:
    """Convert face match lb/ub to strict min/max order.

    By default the library encodes traversal direction in lb/ub ordering
    (lb is not necessarily < ub).  This function normalises every face
    match so that ``lb = [IMIN, JMIN, KMIN]`` and ``ub = [IMAX, JMAX,
    KMAX]`` on **both** sides, and recomputes the orientation /
    permutation matrix accordingly.

    The resulting PM satisfies::

        Face_A(IMIN,JMIN,KMIN -> IMAX,JMAX,KMAX) * PM
            = Face_B(IMIN,JMIN,KMIN -> IMAX,JMAX,KMAX)

    Parameters
    ----------
    face_matches : list of dict
        Face match dicts as returned by :func:`connectivity_fast` or
        :func:`connectivity`.  Each dict must have ``block1`` and
        ``block2`` sub-dicts with ``lb`` and ``ub`` keys.  If a ``match``
        key (point-match DataFrame) is present it is used to recompute
        the orientation; otherwise the orientation is recomputed from the
        new min/max bounds.

    Returns
    -------
    list of dict
        A **new** list with the same structure, but lb/ub in min/max
        order and orientation updated.
    """
    out = []
    for fm in face_matches:
        fm = deepcopy(fm)
        lb1 = fm['block1']['lb']
        ub1 = fm['block1']['ub']
        lb2 = fm['block2']['lb']
        ub2 = fm['block2']['ub']

        new_lb1 = [min(lb1[d], ub1[d]) for d in range(3)]
        new_ub1 = [max(lb1[d], ub1[d]) for d in range(3)]
        new_lb2 = [min(lb2[d], ub2[d]) for d in range(3)]
        new_ub2 = [max(lb2[d], ub2[d]) for d in range(3)]

        fm['block1']['lb'] = new_lb1
        fm['block1']['ub'] = new_ub1
        fm['block2']['lb'] = new_lb2
        fm['block2']['ub'] = new_ub2

        # Recompute orientation if we have the point-match DataFrame
        # with the expected columns (i1,j1,k1,i2,j2,k2)
        df = fm.get('match')
        has_point_match = (df is not None and isinstance(df, pd.DataFrame)
                          and len(df) > 1 and 'i2' in df.columns)
        if has_point_match:
            orientation = _compute_orientation(df, new_lb1, new_ub1)
            perm_idx, plane = _orient_vec_to_permutation(
                orientation, new_lb1, new_ub1, new_lb2, new_ub2)
            export_perm = -1 if plane == 'in-plane' else perm_idx
            fm['orientation'] = {
                'permutation_index': export_perm,
                'plane': plane,
                'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
            }
        elif 'orientation' in fm:
            # No DataFrame — recompute from min/max bounds.
            # With min/max ordering, all axes go forward, so reversal
            # flags depend only on the axis mapping (no direction flip).
            ca1 = _constant_axis(new_lb1, new_ub1)
            ca2 = _constant_axis(new_lb2, new_ub2)
            plane = 'in-plane' if ca1 == ca2 else 'cross-plane'

            # Build identity-like orientation: face1 axis d -> face2 axis d
            # unless cross-plane, in which case use the old orientation's
            # permutation matrix to infer the mapping.
            old_perm = fm['orientation'].get('permutation_matrix')
            if old_perm is not None and plane == 'cross-plane':
                # Preserve the cross-plane axis swap from original
                fm['orientation']['plane'] = plane
            else:
                # In-plane with min/max ordering: both sides go forward,
                # so the PM is identity (no reversal, no swap).
                perm_idx = 0
                fm['orientation'] = {
                    'permutation_index': -1,
                    'plane': plane,
                    'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
                }

        out.append(fm)
    return out
