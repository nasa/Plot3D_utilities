from typing import List, Dict, Tuple, Optional
from itertools import combinations_with_replacement, permutations
import numpy as np
import warnings
from .block import Block
from .blockfunctions import rotate_block, reduce_blocks, compute_min_gcd, scale_face_bounds
from .face import Face
from .facefunctions import outer_face_dict_to_list,match_faces_dict_to_list, create_face_from_diagonals, find_bounding_faces, split_face
from .connectivity import (
    get_face_intersection, _compute_orientation, _orient_vec_to_permutation,
    PERMUTATION_MATRICES, _face_point_count,
    revalidate_full_resolution, demote_to_outer,
    _declared_perm_idx, _failure_severity, _describe_failure,
)
from .geometry import coincidence_count
from .permutation import patch_from_bounds
from .correspondence import certify_correspondence, certify_permutation, MappingFailure, CertifiedMapping
import pandas as pd
from math import cos, radians, sin, sqrt, acos
from copy import deepcopy
from tqdm import trange, tqdm
from scipy.spatial import cKDTree

def periodicity_fast(blocks:List[Block],outer_faces:List[Dict[str,int]], matched_faces:List[Dict[str,int]], periodic_direction:str='k', rotation_axis:str='x',nblades:int=55):
    """Finds the connectivity of blocks when they are rotated by an angle defined by the number of blades. Only use this if your mesh is of an annulus. 
        This function reduces the size of the blocks by a factor of the minimum gcd. This speeds up finding the connectivity 

    Args:
        blocks (List[Block]): List of blocks that will be scanned for perodicity
        outer_faces (List[Dict[str,int]]): List of outer faces for each block as a dictionary format. You can get this from connectivity
        matched_faces (List[Dict[str,int]]): List of matched faces from connectivity. Matched faces was added so that it's always removed from outer faces
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
        nblades (int): Number of blades to consider, this affects the rotation angle. 

    Returns:
        (Tuple): containing

        - **periodic_faces_export** (List[Dict[str,int]]):  This is list of all the surfaces/faces that match when rotated by an angle formatted as a dictionary.
        - **outer_faces_export** (List[Dict[str,int]]): These are the list of outer faces that are not periodic formatted as a dictionary.
        - **periodic_faces** (List[Tuple[Face,Face]]): - This is a list of Face objects that are connected to each other organized as a list of tuples: [Face1, Face2] where Face 1 will contain the block number and the diagonals [IMIN,JMIN,KMIN,IMAX,JMAX,KMAX]. Example: blk: 1 [168,0,0,268,100,0].
        - **outer_faces_all** (List[Face]): This is a list of outer faces save as a list of Faces
        
    """
    gcd_to_use = compute_min_gcd(blocks)
    print(f"gcd to use {gcd_to_use}")
    new_blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)
    matched_faces = deepcopy(matched_faces)
    scale_face_bounds(matched_faces, gcd_to_use, divide=True)
    outer_faces = deepcopy(outer_faces)
    scale_face_bounds(outer_faces, gcd_to_use, divide=True)

    # Find Periodicity
    periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all = periodicity(new_blocks, outer_faces, matched_faces, periodic_direction, rotation_axis, nblades)
    # scale it up
    scale_face_bounds(periodic_faces_export, gcd_to_use)

    for i in range(len(periodic_faces)):
        periodic_faces[i][0].I *= gcd_to_use
        periodic_faces[i][0].J *= gcd_to_use
        periodic_faces[i][0].K *= gcd_to_use

        periodic_faces[i][1].I *= gcd_to_use
        periodic_faces[i][1].J *= gcd_to_use
        periodic_faces[i][1].K *= gcd_to_use

    scale_face_bounds(outer_faces_export, gcd_to_use)

    for j in range(len(outer_faces_all)):
        outer_faces_all[j].I *= gcd_to_use
        outer_faces_all[j].J *= gcd_to_use
        outer_faces_all[j].K *= gcd_to_use

    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all

def create_rotation_matrix(rotation_angle:float, rotation_axis:str="x"):
    """Creates a rotation matrix given an angle and axis 

    Args:
        rotation_angle (float): Rotation angle in radians
        rotation_axis (str, optional): Axis of rotation "x", "y", or "z". Defaults to "x".

    Returns:
        np.ndarray: 3x3 rotation matrix 
    """
    
    if rotation_axis=='x':
        rotation_matrix = np.array([[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]])

    elif rotation_axis=='y':
        rotation_matrix = np.array([[cos(rotation_angle),0,sin(rotation_angle)],
                            [0,1,0],
                            [-sin(rotation_angle),0,cos(rotation_angle)]])
    elif rotation_axis=='z':
        rotation_matrix = np.array([[cos(rotation_angle),-sin(rotation_angle), 0],
                            [sin(rotation_angle),cos(rotation_angle), 0],
                            [0, 0, 1]])
    else:
        raise ValueError(f"rotation_axis must be 'x', 'y', or 'z', got '{rotation_axis}'")

    return rotation_matrix

def _compute_periodic_lb_ub_orientation(
    blk1: 'Block', lb1: list, ub1: list,
    blk2: 'Block', lb2_orig: list, ub2_orig: list,
    shift_axis: Optional[int] = None, shift_amount: float = 0.0
) -> Tuple[list, list, List[int]]:
    """Compute corrected lb2, ub2, and orientation for a periodic face pair.

    Shifts face1 by shift_amount along shift_axis, then uses KDTree lookup
    to find which face2 indices correspond to face1's lb and ub corners.
    Returns corrected lb2/ub2 that produce matching traversal order, plus
    the orientation vector.

    Args:
        blk1: Block containing face 1.
        lb1: Lower-bound [i, j, k] of face 1.
        ub1: Upper-bound [i, j, k] of face 1.
        blk2: Block containing face 2.
        lb2_orig: Original lower-bound [i, j, k] of face 2.
        ub2_orig: Original upper-bound [i, j, k] of face 2.
        shift_axis: Axis index (0/1/2) to shift face1 along, or None.
        shift_amount: Distance to shift face1 along shift_axis.

    Returns:
        (corrected_lb2, corrected_ub2, orientation): Corrected diagonal
        corners for face 2 and 3-element orientation vector.
    """
    from scipy.spatial import cKDTree

    def _get_point(blk, i, j, k):
        return np.array([blk.X[i, j, k], blk.Y[i, j, k], blk.Z[i, j, k]])

    # face1 lb and ub corners (shifted)
    p1_lb = _get_point(blk1, lb1[0], lb1[1], lb1[2])
    p1_ub = _get_point(blk1, ub1[0], ub1[1], ub1[2])
    if shift_axis is not None:
        p1_lb[shift_axis] += shift_amount
        p1_ub[shift_axis] += shift_amount

    lo2 = [min(lb2_orig[d], ub2_orig[d]) for d in range(3)]
    hi2 = [max(lb2_orig[d], ub2_orig[d]) for d in range(3)]

    # Build KDTree of all face2 grid points
    indices2 = []
    coords2 = []
    for i in range(lo2[0], hi2[0] + 1):
        for j in range(lo2[1], hi2[1] + 1):
            for k in range(lo2[2], hi2[2] + 1):
                indices2.append([i, j, k])
                coords2.append(_get_point(blk2, i, j, k))
    coords2 = np.array(coords2)
    indices2 = np.array(indices2)
    tree2 = cKDTree(coords2)

    # Find face2 indices matching face1 lb and ub corners
    _, idx_lb = tree2.query(p1_lb)
    _, idx_ub = tree2.query(p1_ub)
    corrected_lb2 = indices2[idx_lb].tolist()
    corrected_ub2 = indices2[idx_ub].tolist()

    # Compute orientation: step along each face1 axis, find which face2 axis changes
    dims1 = [abs(ub1[d] - lb1[d]) + 1 for d in range(3)]
    step1 = [1 if ub1[d] >= lb1[d] else -1 for d in range(3)]
    cdims2 = [abs(corrected_ub2[d] - corrected_lb2[d]) + 1 for d in range(3)]

    orientation = [0, 0, 0]
    for d1 in range(3):
        if dims1[d1] == 1:
            for d2 in range(3):
                if cdims2[d2] == 1:
                    orientation[d1] = d2 + 1
                    break
        else:
            next_idx1 = list(lb1)
            next_idx1[d1] += step1[d1]
            p1_next = _get_point(blk1, next_idx1[0], next_idx1[1], next_idx1[2])
            if shift_axis is not None:
                p1_next[shift_axis] += shift_amount
            _, idx_next = tree2.query(p1_next)
            face2_next = indices2[idx_next]
            for d2 in range(3):
                if face2_next[d2] != corrected_lb2[d2] and cdims2[d2] > 1:
                    orientation[d1] = d2 + 1
                    break

    # Fill any missing entries (sanity)
    if 0 in orientation:
        used = set(orientation) - {0}
        missing_d1 = [d for d in range(3) if orientation[d] == 0]
        missing_d2 = list({1, 2, 3} - used)
        for d1, d2 in zip(missing_d1, missing_d2):
            orientation[d1] = d2

    return corrected_lb2, corrected_ub2, orientation


def _shift_transform(axis_idx: int, shift_amount: float):
    """Build a `certify_correspondence` ``transform`` for a pure per-pair
    axis translation.

    `translational_periodicity` relates two faces by ``point_L + shift ≈
    point_U`` along `axis_idx` (see `_compute_periodic_lb_ub_orientation`'s
    own ``p1_lb[shift_axis] += shift_amount`` convention, which this
    mirrors). `certify_correspondence` compares ``grid_a`` against
    ``transform(grid_b)``, so bringing U's grid into L's frame means
    subtracting the shift back off, not adding it.
    """
    def _apply(pts: np.ndarray) -> np.ndarray:
        out = np.array(pts, dtype=float, copy=True)
        out[..., axis_idx] -= shift_amount
        return out
    return _apply


def _locate_and_certify_periodic_patch(
    fL: Face, blkL: Block, fU: Face, blkU: Block,
    axis_idx: int, shift_amt: float, tol: float,
) -> Optional[Tuple[list, list, list, list, 'CertifiedMapping']]:
    """Locate the sub-range on face U overlapping face L's FULL extent under
    a pure per-pair axis translation, then certify the pair node-for-node.

    `translational_periodicity`'s coverage-fraction tests (`faces_match`'s
    `touches_by_nodes`/orthogonal precheck, and the oblique-fallback's
    footprint-overlap gate) only ever answer "do these faces touch enough to
    count as periodic" -- they never pin down an explicit claimed sub-patch,
    so nothing downstream can be node-for-node certified. This locates that
    sub-patch and certifies it, mirroring `connectivity.get_face_intersection`
    Step 2's "find the smaller face's corners on the larger face" technique,
    specialized to a known pure-axis translation: because the shift is along
    `axis_idx` only, the two in-plane coordinates are untouched by it, so
    locating the sub-region is a plain nearest-node corner lookup (the
    existing `_compute_periodic_lb_ub_orientation` helper already does
    exactly this) -- no general corner search is required. What Step 2 does
    with `_find_corner_on_face`'s exact-match search, this does with a
    KDTree nearest-neighbor query instead, because the match is only
    provisional until certification below re-checks it in full.

    Callers must only invoke this once their own coverage-fraction precheck
    already indicates ``fL`` is (essentially) fully covered by ``fU`` --
    ``fL``'s own full index range becomes the claimed patch on the L side,
    unconditionally. `certify_correspondence` is what actually confirms or
    rejects this: it checks every node of the claimed patch (not just the
    two corners used to locate it), so a proposal accidentally accepted by
    the coverage-fraction test (e.g. two interior nodes with swapped roles,
    still individually coincident somewhere in the other face's point cloud)
    is caught here even though it fooled the coverage test.

    Returns ``(lb1, ub1, lb2, ub2, certified)`` on success. Returns ``None``
    on `correspondence.MappingFailure` (or a degenerate located region that
    is not a valid single-constant-axis `Patch`) -- callers must treat this
    exactly like "no match", never fall back to the uncertified location.
    """
    lb1 = [fL.IMIN, fL.JMIN, fL.KMIN]
    ub1 = [fL.IMAX, fL.JMAX, fL.KMAX]
    lb2_full = [fU.IMIN, fU.JMIN, fU.KMIN]
    ub2_full = [fU.IMAX, fU.JMAX, fU.KMAX]

    lb2, ub2, _orient = _compute_periodic_lb_ub_orientation(
        blkL, lb1, ub1, blkU, lb2_full, ub2_full,
        shift_axis=axis_idx, shift_amount=shift_amt)

    try:
        patch1 = patch_from_bounds(fL.blockIndex, lb1, ub1)
        patch2 = patch_from_bounds(fU.blockIndex, lb2, ub2)
    except ValueError:
        # Corners snapped to a degenerate region (not a single-constant-axis
        # box) -- not a valid patch to certify.
        return None

    try:
        certified = certify_correspondence(
            blkL, patch1, blkU, patch2, tol,
            transform=_shift_transform(axis_idx, shift_amt))
    except MappingFailure:
        return None

    return lb1, ub1, lb2, ub2, certified


def _build_periodic_export(df: pd.DataFrame, periodic_faces_temp: list,
                           rotation_sign: int = 1) -> dict:
    """Build export dict from DataFrame with orientation, consistent with connectivity().

    Derives lb/ub from the DataFrame traversal order (first/last row) and computes
    the orientation vector using _compute_orientation. Uses the -1 convention for
    in-plane permutation_index.

    Args:
        df: DataFrame with columns i1,j1,k1,i2,j2,k2 from face matching.
        periodic_faces_temp: List of [face1, face2] Face pairs.
        rotation_sign: Which pitch direction produced this match. ``+1`` means
            face1 rotated by *+angle* landed on face2; ``-1`` means it took
            *-angle*. Negative matches are stored swapped -- see below.

    Returns:
        Face match dict with block1, block2 sub-dicts and orientation.
    """
    if rotation_sign < 0:
        # This pair matched only when face1 was rotated by *-angle*, i.e.
        # R(-a) @ face1 == face2, equivalently face1 == R(+a) @ face2. A
        # connectivity.json carries a single global transformation_matrix and
        # declares "Face_B = transformation_matrix @ Face_A", so that promise
        # holds for every pair only if all pairs are stored in the same
        # rotational direction. Swap the two sides here so the stored block1
        # is always the one you rotate *forward* to reach block2; otherwise a
        # consumer trusting the declared convention mis-rotates these pairs by
        # twice the pitch.
        df = df.rename(columns={'i1': 'i2', 'j1': 'j2', 'k1': 'k2',
                                'i2': 'i1', 'j2': 'j1', 'k2': 'k1'})
        # Rows arrived in the *old* face1's nested traversal order. lb/ub and
        # _compute_orientation both read that ordering off the row sequence,
        # so re-sort into the new face1's ascending I-J-K traversal -- which
        # is also the "block1 is always ascending" half of the lb/ub convention.
        df = df.sort_values(['i1', 'j1', 'k1']).reset_index(drop=True)
        periodic_faces_temp = [periodic_faces_temp[1], periodic_faces_temp[0]]

    lb1 = [int(df.iloc[0]['i1']), int(df.iloc[0]['j1']), int(df.iloc[0]['k1'])]
    ub1 = [int(df.iloc[-1]['i1']), int(df.iloc[-1]['j1']), int(df.iloc[-1]['k1'])]
    lb2 = [int(df.iloc[0]['i2']), int(df.iloc[0]['j2']), int(df.iloc[0]['k2'])]
    ub2 = [int(df.iloc[-1]['i2']), int(df.iloc[-1]['j2']), int(df.iloc[-1]['k2'])]
    orientation = _compute_orientation(df, lb1, ub1)
    perm_idx, plane = _orient_vec_to_permutation(orientation, lb1, ub1, lb2, ub2)
    export_perm = -1 if plane == 'in-plane' else perm_idx
    return {
        'block1': {'block_index': periodic_faces_temp[0].BlockIndex,
                    'lb': lb1, 'ub': ub1, 'id': periodic_faces_temp[0].id},
        'block2': {'block_index': periodic_faces_temp[1].BlockIndex,
                    'lb': lb2, 'ub': ub2, 'id': periodic_faces_temp[1].id},
        'orientation': {
            'permutation_index': export_perm,
            'plane': plane,
            'permutation_matrix': PERMUTATION_MATRICES[perm_idx].tolist(),
        }
    }


def periodicity(blocks:List[Block],outer_faces:List[Dict[str,int]], matched_faces:List[Dict[str,int]], periodic_direction:str='k', rotation_axis:str='x',nblades:int=55):
    """This function is used to check for periodicity of the other faces rotated about an axis. Use periodicity_fast instead. Periodicity_fast calls this function after reducing the size of the mesh.
        The way it works is to find faces of a constant i,j, or k value

    Args:
        blocks (List[Block]): List of blocks that will be scanned for perodicity
        outer_faces (List[Dict[str,int]]): List of outer faces for each block as a dictionary format. You can get this from connectivity
        matched_faces (List[Dict[str,int]]): List of matched faces from connectivity. Matched faces was added so that it's always removed from outer faces
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
        nblades (int): Number of blades to consider, this affects the rotation angle. 

    Returns:
        (Tuple): containing
            
            - **periodic_faces_export** (List[Dict[str,int]]):  This is list of all the surfaces/faces that match when rotated by an angle formatted as a dictionary.
            - **outer_faces_export** (List[Dict[str,int]]): These are the list of outer faces that are not periodic formatted as a dictionary.
            - **periodic_faces** (List[Tuple[Face,Face]]): - This is a list of Face objects that are connected to each other organized as a list of tuples: [Face1, Face2] where Face 1 will contain the block number and the diagonals [IMIN,JMIN,KMIN,IMAX,JMAX,KMAX]. Example: blk: 1 [168,0,0,268,100,0].
            - **outer_faces_all** (List[Face]): This is a list of outer faces save as a list of Faces

    """
    
    rotation_angle = radians(360.0/nblades)
    rotation_matrix1 = create_rotation_matrix(rotation_angle,rotation_axis)
    rotation_matrix2 = create_rotation_matrix(-rotation_angle,rotation_axis)
    
    # Check periodic within a block 
    periodic_found = True
    
    # Here we make a list of all the outer faces
    periodic_faces = list()      # This is the output of the code 
    periodic_faces_export = list() 
    outer_faces_all = outer_face_dict_to_list(blocks,outer_faces)
    matched_faces_all = match_faces_dict_to_list(blocks,matched_faces)

    split_faces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list
    outer_faces_to_remove = list()  # Integer list of which outer surfaces to remove
    while periodic_found:
        periodic_found = False
        outer_faces_to_remove = list()
        outer_face_combos = list(combinations_with_replacement(range(len(outer_faces_all)),2))
        t = trange(len(outer_face_combos))
        for i in t: 
            # Check if surfaces are periodic with each other
            face1_indx = outer_face_combos[i][0]
            face2_indx = outer_face_combos[i][1]
            face1 = outer_faces_all[face1_indx]
            face2 = outer_faces_all[face2_indx]
            t.set_description(f"Checking connections block {face1.blockIndex} with {face2.blockIndex}")

            dir_map = {'i': 0, 'j': 1, 'k': 2}
            target_const = dir_map[periodic_direction.lower()]

            if face1.const_type == target_const and face2.const_type == target_const:
                block1_rotated = rotate_block(blocks[face1.blockIndex], rotation_matrix1)
                block2 = blocks[face2.blockIndex]
                rotation_sign = 1
                df, periodic_faces_temp, split_faces_temp = __periodicity_check__(
                    face1, face2, block1_rotated, block2)

                if len(periodic_faces_temp) == 0:
                    block1_rotated = rotate_block(blocks[face1.blockIndex], rotation_matrix2)
                    rotation_sign = -1
                    df, periodic_faces_temp, split_faces_temp = __periodicity_check__(
                        face1, face2, block1_rotated, block2)

                if len(periodic_faces_temp) > 0:
                    outer_faces_to_remove.append(face1)
                    outer_faces_to_remove.append(face2)
                    outer_faces_to_remove.append(periodic_faces_temp[0])
                    outer_faces_to_remove.append(periodic_faces_temp[1])
                    periodic_faces.append(periodic_faces_temp)
                    periodic_faces_export.append(_build_periodic_export(df, periodic_faces_temp, rotation_sign))
                    split_faces.extend(split_faces_temp)
                    periodic_found = True
                    break

        if (periodic_found):
            outer_faces_to_remove = list(set(outer_faces_to_remove))
            outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]
            if len(split_faces)>0:
                outer_faces_all.extend(split_faces)
                split_faces.clear()

    # This is an added check to make sure all periodic faces are in the outer_faces_to_remove
    for p in periodic_faces:
        outer_faces_to_remove.append(p[0])
        outer_faces_to_remove.append(p[1])

    for m in matched_faces_all:
        outer_faces_to_remove.append(m)

    outer_faces_to_remove = list(set(outer_faces_to_remove))    # Use only unique values
    outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]    # remove from outer faces 
    # remove any duplicate periodic face pairs 
    indx_to_remove = list()
    for i in range(len(periodic_faces)):
        for j in range(i+1,len(periodic_faces)):
            if periodic_faces[i][0] == periodic_faces[j][0]:
                if periodic_faces[i][1] == periodic_faces[j][1]:
                    indx_to_remove.append(j)
            if periodic_faces[i][1] == periodic_faces[j][0]:
                if periodic_faces[i][0] == periodic_faces[j][1]:
                    indx_to_remove.append(j)
    

    periodic_faces_export = [periodic_faces_export[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    periodic_faces = [periodic_faces[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    # Export periodic faces and outer faces
    outer_faces_export = list() 

    for o in outer_faces_all:
        outer_faces_export.append(o.to_dict())
                        
    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all


def rotated_periodicity(blocks:List[Block], matched_faces:List[Dict[str,int]], outer_faces:List[Dict[str,int]], rotation_angle:float, rotation_axis:str = "x", ReduceMesh:bool=True, use_minmax:bool=False, tol:float=1E-4):
    """Finds the peridocity/connectivity by rotating a block. This is a bit different from "periodicity_fast" where you specify the periodic direction. This method doesn't care about the direction as long as the angle you specify results in a match between the Left Face and the Right Face. I would use this instead.

    Example 1:              
        L      RL           R
        | blk1 || Copy blk1 |
        | blk2 || Copy blk2 |
        | blk3 || Copy blk3 |
        Rotates the set of blocks by an angle and checks the matching surfaces for R and L. 

    Args:
        blocks (List[Block]): List of blocks for a particular geometry. Do not duplicate the geometry and pass it in!
        matched_faces (List[Dict[str,int]]): List of matched faces from connectivity. Matched faces was added so that it's always removed from outer faces
        outer_faces (List[Dict[str,int]]): List of outer faces in dictionary form
        rotation_angle (float): rotation angle in between geometry in degrees.
        rotation_axis (str, Optional): "x", "y", or "z" 
        ReduceMesh (bool, Optional): True, reduces the mesh for faster matching
    
    periodic_faces, outer_faces_export, _, _ = rotated_periodicity(blocks,face_matches, outer_faces, rotation_angle=rotation_angle, rotation_axis = "x")
    
    Replaces:
        Is the same as 

        periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity_fast(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)
        and
        periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)


    Returns:
        (Tuple): containing
            
            - **periodic_faces_export** (List[Dict[str,int]]):  This is list of all the surfaces/faces that match when rotated by an angle formatted as a dictionary.
            - **outer_faces_export** (List[Dict[str,int]]): These are the list of outer faces that are not periodic formatted as a dictionary.
            - **periodic_faces** (List[Tuple[Face,Face]]): - This is a list of Face objects that are connected to each other organized as a list of tuples: [Face1, Face2] where Face 1 will contain the block number and the diagonals [IMIN,JMIN,KMIN,IMAX,JMAX,KMAX]. Example: blk: 1 [168,0,0,268,100,0].
            - **outer_faces_all** (List[Face]): This is a list of outer faces save as a list of Faces
    """
    gcd_to_use = 1
    # Retain a reference to the original full-resolution blocks BEFORE any
    # reduction below reassigns `blocks` -- A5 full-resolution re-validation
    # (at the end of this function) must certify against these, never
    # against the GCD-reduced copy.
    full_res_blocks = blocks
    if ReduceMesh:
        gcd_to_use = compute_min_gcd(blocks)
        blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)

    rotation_matrix_pos = create_rotation_matrix(radians(rotation_angle),rotation_axis)
    rotation_matrix_neg = create_rotation_matrix(-radians(rotation_angle),rotation_axis)

    _rotated_cache_pos: dict = {}
    _rotated_cache_neg: dict = {}
    def get_rotated(block_index, sign=+1):
        cache = _rotated_cache_pos if sign > 0 else _rotated_cache_neg
        rmat = rotation_matrix_pos if sign > 0 else rotation_matrix_neg
        if block_index not in cache:
            cache[block_index] = rotate_block(blocks[block_index], rmat)
        return cache[block_index]

    # Check periodic within a block 
    periodic_found = True
    
    # Here we make a list of all the outer faces
    periodic_faces = list()      # This is the output of the code 
    periodic_faces_export = list() 

    outer_faces_all = outer_face_dict_to_list(blocks,outer_faces,gcd_to_use)
    matched_faces_all = match_faces_dict_to_list(blocks,matched_faces,gcd_to_use)

    split_faces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list
    non_matching = list()
    outer_faces_to_remove = list()  # Integer list of which outer surfaces to remove
    while periodic_found:
        periodic_found = False
        outer_faces_to_remove = list()
        outer_face_combos = list(permutations(range(len(outer_faces_all)),2))
        outer_face_combos = list(set(outer_face_combos) - set(non_matching)) # removes face combinations already checked
        t = trange(len(outer_face_combos))
        for i in t:
            # Check if surfaces are periodic with each other
            face1_indx = outer_face_combos[i][0]
            face2_indx = outer_face_combos[i][1]
            face1 = outer_faces_all[face1_indx]
            face2 = outer_faces_all[face2_indx]

            if (face1.IMIN == face1.IMAX) and (face2.IMIN == face2.IMAX) or \
                (face1.JMIN == face1.JMAX) and (face2.JMIN == face2.JMAX) or \
                (face1.KMIN == face1.KMAX) and (face2.KMIN == face2.KMAX):

                # Rotate Block 1 -> Check periodicity -> if not periodic -> Rotate Block 1 opposite direction -> Check periodicity
                block2 = blocks[face2.blockIndex]
                t.set_description(f"Blk {face1.blockIndex} <-> {face2.blockIndex} | found {len(periodic_faces)}")
                # Try +rotation first
                block1_rotated = get_rotated(face1.blockIndex, sign=+1)
                rotation_sign = 1
                df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2, tol=tol)
                # If no match, retry with -rotation. Without this we miss
                # matches whose pitch direction is opposite to the (face1,
                # face2) ordering tried — the old `periodicity` function
                # tried both directions; `rotated_periodicity` did not, so
                # any pair filtered out of the reverse ordering by the
                # stale-index `non_matching` set was silently dropped.
                if len(periodic_faces_temp) == 0:
                    block1_rotated = get_rotated(face1.blockIndex, sign=-1)
                    rotation_sign = -1
                    df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2, tol=tol)
                
                if len(periodic_faces_temp) > 0:
                    outer_faces_to_remove.append(face1)
                    outer_faces_to_remove.append(face2)
                    outer_faces_to_remove.append(periodic_faces_temp[0])
                    outer_faces_to_remove.append(periodic_faces_temp[1])
                    periodic_faces.append(periodic_faces_temp)
                    periodic_faces_export.append(_build_periodic_export(df, periodic_faces_temp, rotation_sign))
                    split_faces.extend(split_faces_temp)
                    periodic_found = True
                    break
            else: 
                non_matching.append((face1_indx,face2_indx))

        if (periodic_found):
            outer_faces_to_remove = list(set(outer_faces_to_remove))
            outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]
            if len(split_faces)>0:
                outer_faces_all.extend(split_faces)
                split_faces.clear()

    # Free rotated block copies no longer needed after the matching loop
    del _rotated_cache_pos, _rotated_cache_neg

    # This is an added check to make sure all periodic faces are in the outer_faces_to_remove
    for p in periodic_faces:
        outer_faces_to_remove.append(p[0])
        outer_faces_to_remove.append(p[1])

    for m in matched_faces_all:
        outer_faces_to_remove.append(m)

    outer_faces_to_remove = list(set(outer_faces_to_remove))    # Use only unique values
    outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]    # remove from outer faces

    # remove any duplicate periodic face pairs
    indx_to_remove = list()
    for i in range(len(periodic_faces)):
        for j in range(i+1,len(periodic_faces)):
            if periodic_faces[i][0] == periodic_faces[j][0]:
                if periodic_faces[i][1] == periodic_faces[j][1]:
                    indx_to_remove.append(j)
            if periodic_faces[i][1] == periodic_faces[j][0]:
                if periodic_faces[i][0] == periodic_faces[j][1]:
                    indx_to_remove.append(j)


    periodic_faces_export = [periodic_faces_export[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    periodic_faces = [periodic_faces[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    # Export periodic faces and outer faces
    outer_faces_export = list()

    for o in outer_faces_all:
        outer_faces_export.append(o.to_dict())

    # scale it up
    scale_face_bounds(periodic_faces_export, gcd_to_use)
    scale_face_bounds(outer_faces_export, gcd_to_use)

    for i in range(len(periodic_faces)):
        periodic_faces[i][0].I *= gcd_to_use
        periodic_faces[i][0].J *= gcd_to_use
        periodic_faces[i][0].K *= gcd_to_use

        periodic_faces[i][1].I *= gcd_to_use
        periodic_faces[i][1].J *= gcd_to_use
        periodic_faces[i][1].K *= gcd_to_use

    for j in range(len(outer_faces_all)):
        outer_faces_all[j].I *= gcd_to_use
        outer_faces_all[j].J *= gcd_to_use
        outer_faces_all[j].K *= gcd_to_use

    # A5: full-resolution re-validation after GCD reduction. GCD reduction
    # can occasionally produce a proposal that looks valid on the coarse
    # grid but does not actually hold node-for-node at full resolution.
    # Re-certify every scaled-up proposal against the ORIGINAL,
    # full-resolution `full_res_blocks` (never the reduced `blocks`) before
    # returning it. Both the +angle and -angle rotations are offered as
    # candidate transforms: `_build_periodic_export` already normalizes
    # every stored pair so block1 rotated *forward* reaches block2, but
    # trying both here is cheap insurance against that convention rather
    # than a second, independent assumption to keep in sync.
    if ReduceMesh and gcd_to_use > 1:
        t_fwd = lambda pts, R=rotation_matrix_pos: pts @ R.T
        t_bwd = lambda pts, R=rotation_matrix_neg: pts @ R.T
        kept, rejected = revalidate_full_resolution(
            full_res_blocks, periodic_faces_export, transforms=[t_fwd, t_bwd],
            tol=tol, stage="rotated_periodicity",
        )
        if rejected:
            # Keep the parallel Face-object-form return channel
            # (`periodic_faces`/`outer_faces_all`) in sync with the
            # dict-form channel: `periodic_faces_export` and `periodic_faces`
            # have been built and filtered in lockstep everywhere above, so
            # they are still index-aligned here -- match rejected dicts back
            # to their Face-tuple by identity (revalidate_full_resolution
            # never copies the proposal dicts it is given).
            rejected_ids = {id(r) for r in rejected}
            new_periodic_faces = []
            for export_rec, face_pair in zip(periodic_faces_export, periodic_faces):
                if id(export_rec) in rejected_ids:
                    outer_faces_all.append(face_pair[0])
                    outer_faces_all.append(face_pair[1])
                else:
                    new_periodic_faces.append(face_pair)
            periodic_faces = new_periodic_faces
        periodic_faces_export = kept
        demote_to_outer(outer_faces_export, rejected)

    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all

def _median_inplane_spacing(face: Face, block: Block) -> Optional[float]:
    """Median edge length on the face (in-plane), used to derive an
    adaptive per-pair tolerance in `translational_periodicity`.

    Returns None when the face has too few in-plane points to measure any
    spacing from (a degenerate 1xN or 0xN footprint) — callers must treat
    this as "no tolerance derivable", never fabricate a value.
    """
    I0, I1, J0, J1, K0, K1 = face.IMIN, face.IMAX, face.JMIN, face.JMAX, face.KMIN, face.KMAX
    X, Y, Z = block.X, block.Y, block.Z
    if face.const_type == 0:  # I const → vary (J,K)
        i = I0
        x = X[i, J0:J1+1, K0:K1+1]; y = Y[i, J0:J1+1, K0:K1+1]; z = Z[i, J0:J1+1, K0:K1+1]
    elif face.const_type == 1:  # J const → vary (I,K)
        j = J0
        x = X[I0:I1+1, j, K0:K1+1]; y = Y[I0:I1+1, j, K0:K1+1]; z = Z[I0:I1+1, j, K0:K1+1]
    else:  # K const → vary (I,J)
        k = K0
        x = X[I0:I1+1, J0:J1+1, k]; y = Y[I0:I1+1, J0:J1+1, k]; z = Z[I0:I1+1, J0:J1+1, k]
    s = []
    if x.shape[0] > 1:
        dx = np.diff(x, axis=0); dy = np.diff(y, axis=0); dz = np.diff(z, axis=0)
        s.append(np.sqrt(dx*dx + dy*dy + dz*dz))
    if x.shape[1] > 1:
        dx = np.diff(x, axis=1); dy = np.diff(y, axis=1); dz = np.diff(z, axis=1)
        s.append(np.sqrt(dx*dx + dy*dy + dz*dz))
    if not s:
        return None
    return float(np.median(np.concatenate([v.ravel() for v in s])))


def _combine_pair_spacings(sA: Optional[float], sB: Optional[float]) -> Optional[float]:
    """Combine two per-face in-plane spacing estimates into one adaptive
    absolute tolerance for `translational_periodicity`.

    Both None → None (no data on either side, tolerance is not derivable).
    Exactly one None → falls back to the other side's value (never averages
    a None with a number). Otherwise ~3% of the larger spacing, floored at
    1e-4.
    """
    if sA is None and sB is None:
        return None
    s = max(v for v in (sA, sB) if v is not None)
    return max(0.03 * s, 1e-4)


def translational_periodicity(
    blocks: List[Block],
    outer_faces: List[Dict[str,int]],
    delta: Optional[float] = None,
    translational_direction: str = "z",
    node_tol_xyz: Optional[float] = None,        # global override; if None we compute per-pair adaptively
    min_shared_frac: float = 0.02,
    min_shared_abs: int = 4,
    stride_u: int = 1,
    stride_v: int = 1,
    use_minmax: bool = False,
    ) -> Tuple[ 
               List[Dict[str, Dict[str, int]]], 
               List[Tuple[Face, Face, Dict[str, str]]], 
               List[Dict[str,int]]
            ]:
    """
    Detect translational periodicity between block faces along a given axis.

    This function takes a set of outer block faces and attempts to identify 
    periodic counterparts across the domain in the specified translational 
    direction ('x', 'y', or 'z'). It works by:

      1. **Bounding faces:** Uses `find_bounding_faces` to identify candidate 
         lower/upper faces for the given axis.
      2. **Grid reduction:** Reduces blocks to their greatest common divisor 
         (GCD) resolution to make indexing consistent across blocks.
      3. **Shifting:** Creates shifted copies of all blocks in both positive 
         and negative directions along the periodic axis.
      4. **Precheck in orthogonal plane:** Uses a fast projection test 
         (orthogonal to the periodic axis) to determine whether two faces 
         could possibly match. This greatly reduces false negatives when 
         spacing/tolerances differ slightly.
      5. **Node-based match:** Calls `Face.touches_by_nodes` on candidate 
         pairs to check shared node positions, with an adaptive tolerance 
         based on the in-plane spacing of each face.
      6. **Pairing:** Records each valid pair of periodic faces, their 
         IJK index mappings (min→min or min→max), and removes matched faces 
         from the outer-face list.
      7. **Scaling back:** Rescales reduced indices back to the original grid 
         spacing so results are consistent with input block resolution.

    Args:
        blocks (List[Block]): List of blocks.
        outer_faces (List[Dict[str,int]]): Outer faces represented as
            dictionaries (with lb and ub lists).
        delta (float, optional): Periodicity spacing along the chosen axis. 
            If None, it is inferred from the global block min/max extent.
        translational_direction (str, optional): Axis to check ('x','y','z'). 
            Default is 'z'.
        node_tol_xyz (float, optional): Absolute coordinate tolerance for 
            node-matching. If None, tolerance is computed adaptively based 
            on median in-plane spacing of candidate faces.
        min_shared_frac (float, optional): Minimum fraction of nodes that must 
            overlap for two faces to be considered periodic. Default 0.02.
        min_shared_abs (int, optional): Minimum absolute number of shared nodes. 
            Default 4.
        stride_u (int, optional): Subsampling stride along the first face index 
            direction. Default 1 (no skipping).
        stride_v (int, optional): Subsampling stride along the second face index 
            direction. Default 1 (no skipping).

    Returns:
        Tuple[
            List[Dict[str, Dict[str,int]]], 
            List[Tuple[Face, Face, Dict[str,str]]], 
            List[Dict[str,int]]
        ]:
            - **periodic_faces_export**: Export-ready dictionaries describing 
              each periodic pair (block indices, face extents, index mapping, 
              and match mode).
            - **periodic_pairs**: Matched periodic face pairs as `Face` objects 
              with IJK mapping.
            - **outer_faces_remaining**: Updated list of outer faces with 
              periodic ones removed (preserving any existing `id` fields).

    Notes:
        - Works for periodicity in **x**, **y**, or **z** directions.
        - The adaptive tolerance makes the method robust to small spacing 
          differences between blocks.
        - The orthogonal-plane precheck avoids expensive node comparisons 
          when faces clearly do not align.
    """
    # 0) lower/upper via your finder (dicts at original scale)
    lower_connected_faces, upper_connected_faces, _, _ = find_bounding_faces(
        blocks, outer_faces, translational_direction, "both"
    )

    axis = translational_direction.lower().strip()
    assert axis in ("x","y","z")

    # 1) GCD reduce
    gcd_to_use = compute_min_gcd(blocks)

    lower_faces_r = outer_face_dict_to_list(blocks, lower_connected_faces, gcd_to_use)
    upper_faces_r = outer_face_dict_to_list(blocks, upper_connected_faces, gcd_to_use)
    blocks_r = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # 2) Δ along axis (if not provided)
    if axis == "x":
        a_min = min(b.X.min() for b in blocks_r); a_max = max(b.X.max() for b in blocks_r)
    elif axis == "y":
        a_min = min(b.Y.min() for b in blocks_r); a_max = max(b.Y.max() for b in blocks_r)
    else:
        a_min = min(b.Z.min() for b in blocks_r); a_max = max(b.Z.max() for b in blocks_r)
    d_axis = (a_max - a_min) if (delta is None) else float(delta)

    # 3) Shifted copies
    def shift_blocks(bb: List[Block], amount: float) -> List[Block]:
        cp = deepcopy(bb)
        for b in cp:
            b.shift(amount, axis)
        return cp

    blocks_up = shift_blocks(blocks_r, +d_axis)
    blocks_dn = shift_blocks(blocks_r, -d_axis)

    def B(which: str, idx: int) -> Block:
        return {"orig": blocks_r, "up": blocks_up, "dn": blocks_dn}[which][idx]

    # 4) Helpers for adaptive tolerance (module-level: see
    #    `_median_inplane_spacing`/`_combine_pair_spacings` above)
    def _pair_tol(fA: Face, fB: Face) -> Optional[float]:
        """Adaptive absolute tolerance per pair (use global override if provided).

        Returns None when neither face has enough in-plane points to
        estimate a spacing from — callers must skip the pair rather than
        match it against a fabricated tolerance.
        """
        if node_tol_xyz is not None:
            return float(node_tol_xyz)
        sA = _median_inplane_spacing(fA, B("orig", fA.BlockIndex))
        sB = _median_inplane_spacing(fB, B("orig", fB.BlockIndex))
        return _combine_pair_spacings(sA, sB)

    # 5) General orthogonal-plane precheck (works for x/y/z periodicity)
    def _orthogonal_precheck(fA: Face, fB: Face, bA: Block, bB: Block,
                             d_axis_local: float, tol: float, axis_local: str) -> bool:
        """
        Shift face A along 'axis_local' by d_axis_local, then compare projections onto the
        orthogonal plane within tolerance. Requires both absolute and fractional overlap.
        """
        PA = fA.grid_points(bA, stride_u=1, stride_v=1)
        PB = fB.grid_points(bB, stride_u=1, stride_v=1)
        if PA.size == 0 or PB.size == 0:
            return False

        if axis_local == "x":
            PA[:,0] += d_axis_local
            projA, projB = PA[:,1:], PB[:,1:]          # (y,z)
        elif axis_local == "y":
            PA[:,1] += d_axis_local
            projA, projB = PA[:,[0,2]], PB[:,[0,2]]    # (x,z)
        else:  # "z"
            PA[:,2] += d_axis_local
            projA, projB = PA[:,:2], PB[:,:2]          # (x,y)

        shared = coincidence_count(projA, projB, tol)
        return shared >= max(min_shared_abs, int(min_shared_frac * min(len(projA), len(projB))))

    # 6) Node-sharing matcher using per-pair tol + precheck
    def faces_match(fL: Face, fU: Face) -> Tuple[bool, str]:
        bl, bu = fL.BlockIndex, fU.BlockIndex
        tol_pair = _pair_tol(fL, fU)
        if tol_pair is None:
            return False, ""

        # Fast precheck on orthogonal plane (lower up vs upper orig)
        if _orthogonal_precheck(fL, fU, B("orig", bl), B("orig", bu), d_axis, tol_pair, axis):
            return True, f"{axis}_precheck_lower_up"

        # lower moved up vs upper orig
        if fL.touches_by_nodes(fU, B("up", bl), B("orig", bu),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "lower_up_vs_upper_orig"

        # lower orig vs upper moved down
        if fL.touches_by_nodes(fU, B("orig", bl), B("dn", bu),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "lower_orig_vs_upper_dn"

        # Symmetric guards
        if _orthogonal_precheck(fU, fL, B("orig", bu), B("orig", bl), d_axis, tol_pair, axis):
            return True, f"{axis}_precheck_upper_up"
        if fU.touches_by_nodes(fL, B("up", bu), B("orig", bl),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "upper_up_vs_lower_orig"
        if fU.touches_by_nodes(fL, B("orig", bu), B("dn", bl),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "upper_orig_vs_lower_dn"

        return False, ""

    # 7) Index mapping
    def mapping_minmax(fA: Face, fB: Face) -> Dict[str, str]:
        out = {}
        for ax in ("I","J","K"):
            Amin, Amax = getattr(fA, ax+"MIN"), getattr(fA, ax+"MAX")
            Bmin, Bmax = getattr(fB, ax+"MIN"), getattr(fB, ax+"MAX")
            if (Amin == Bmin) and (Amax == Bmax): out[ax] = "min->min"
            elif (Amin == Bmax) and (Amax == Bmin): out[ax] = "min->max"
            else:
                d_mm = abs(Amin-Bmin)+abs(Amax-Bmax)
                d_mM = abs(Amin-Bmax)+abs(Amax-Bmin)
                out[ax] = "min->min" if d_mm <= d_mM else "min->max"
        return out

    # 8) Centroid-sorted pairing (avoids edge-sharing false matches)
    #    For each lower face, compute shifted centroid and try upper faces
    #    in order of centroid proximity (nearest first).
    lower_pool = list(dict.fromkeys(lower_faces_r))
    upper_pool = list(dict.fromkeys(upper_faces_r))
    periodic_pairs_r: List[Tuple[Face, Face, Dict[str,str]]] = []
    periodic_export: List[Dict[str, Dict[str,int]]] = []

    axis_idx = {"x": 0, "y": 1, "z": 2}[axis]

    def _face_key(f: Face) -> Tuple:
        return (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)

    def _record_match(
        fL: Face, fU: Face, mode: str,
        loc: Tuple[list, list, list, list, 'CertifiedMapping'],
    ) -> None:
        """Build the export record + pair entry for one certified matched pair.

        ``loc`` is the ``(lb1, ub1, lb2, ub2, certified)`` tuple returned by
        `_locate_and_certify_periodic_patch`: ``lb2``/``ub2`` are the
        located (not merely fU's raw full-face) sub-range, and the
        certified permutation/plane -- not a corners-only orientation
        vector -- drive the orientation fields, consistent with how
        `correspondence.CertifiedMapping.plane` already matches
        `connectivity._orient_vec_to_permutation`'s 'in-plane'/'cross-plane'
        convention (see `correspondence._plane_for`).
        """
        m = mapping_minmax(fL, fU)
        periodic_pairs_r.append((fL, fU, m))
        lb1, ub1, lb2, ub2, certified = loc
        export_perm = -1 if certified.plane == 'in-plane' else certified.permutation_index
        periodic_export.append({
            "block1": {"block_index": fL.BlockIndex,
                       "lb": lb1, "ub": ub1},
            "block2": {"block_index": fU.BlockIndex,
                       "lb": lb2, "ub": ub2},
            "orientation": {
                "permutation_index": export_perm,
                "plane": certified.plane,
                "permutation_matrix": PERMUTATION_MATRICES[certified.permutation_index].tolist(),
            },
            "mapping": m,
            "mode": mode
            }) # type: ignore

    # Pre-compute centroids for all faces (in-plane XY after shift)
    def _centroid_inplane(f: Face, shifted: bool = False) -> np.ndarray:
        pts = f.grid_points(B("orig", f.BlockIndex), stride_u=1, stride_v=1)
        c = pts.mean(axis=0)
        if shifted:
            c[axis_idx] += d_axis
        # Project out the periodic axis
        return np.delete(c, axis_idx)

    upper_centroids = np.array([_centroid_inplane(fU) for fU in upper_pool])
    consumed_keys: set = set()

    for fL in list(lower_pool):
        cL = _centroid_inplane(fL, shifted=True)
        # Sort upper candidates by centroid distance
        dists = np.linalg.norm(upper_centroids - cL, axis=1)
        order = np.argsort(dists)

        matched_j = -1
        matched_mode = ""
        matched_loc = None
        for rank in order:
            j = int(rank)
            fU = upper_pool[j]
            ok, mode = faces_match(fL, fU)
            if not ok:
                continue
            blk1_orig = B("orig", fL.BlockIndex)
            blk2_orig = B("orig", fU.BlockIndex)
            arr1 = [blk1_orig.X, blk1_orig.Y, blk1_orig.Z][axis_idx]
            arr2 = [blk2_orig.X, blk2_orig.Y, blk2_orig.Z][axis_idx]
            p1_val = arr1[fL.IMIN, fL.JMIN, fL.KMIN]
            p2_val = arr2[fU.IMIN, fU.JMIN, fU.KMIN]
            shift_amt = d_axis if p1_val < p2_val else -d_axis
            tol_pair = _pair_tol(fL, fU)
            if tol_pair is None:
                continue
            # A coverage-fraction pass (`faces_match`) is only provisional --
            # locate the actual overlapping sub-patch and certify it
            # node-for-node before accepting. A candidate that fails here is
            # treated exactly like a `faces_match` rejection: keep trying
            # the next-nearest upper face rather than accepting the
            # uncertified location.
            loc = _locate_and_certify_periodic_patch(
                fL, blk1_orig, fU, blk2_orig, axis_idx, shift_amt, tol_pair)
            if loc is None:
                continue
            matched_j = j
            matched_mode = mode
            matched_loc = loc
            break

        if matched_j >= 0:
            fU = upper_pool[matched_j]
            _record_match(fL, fU, matched_mode, matched_loc)
            consumed_keys.add(_face_key(fL))
            consumed_keys.add(_face_key(fU))
            # Remove matched upper face from pool and centroids
            upper_pool.pop(matched_j)
            upper_centroids = np.delete(upper_centroids, matched_j, axis=0)

    # 8b) Oblique fallback — bladed-cascade pitch boundaries.
    #
    # The bounding-face candidate selection above only sees faces lying at
    # the GLOBAL axis extremes, i.e. flat constant-coordinate pitch planes.
    # A bladed cascade's pitch boundaries are oblique (they follow the
    # metal-angle inlet/outlet extensions and hug the O-grid), so the
    # lower/upper pools come back empty and the legacy path finds nothing
    # even on an exactly-periodic mesh.
    #
    # This pass considers ALL still-unmatched outer faces as candidates,
    # exploiting two properties of a pure translation along `axis`:
    #   1. The pair's centroids COINCIDE in the orthogonal plane.
    #   2. The translation Δ is the centroid difference along `axis`
    #      (no global-extent assumption — that heuristic is wrong for
    #      cascades where the domain is wider than one pitch).
    # Candidate pairs are sorted by orthogonal-centroid distance and
    # verified with a full quantized-node intersection under the per-pair
    # shift, so false positives require two faces that genuinely coincide
    # node-for-node after translation — which IS translational periodicity.
    all_faces_r = list(dict.fromkeys(
        outer_face_dict_to_list(blocks, outer_faces, gcd_to_use)
    ))
    remaining_faces = [f for f in all_faces_r
                       if _face_key(f) not in consumed_keys]

    if len(remaining_faces) >= 2:
        pts_cache: Dict[Tuple, np.ndarray] = {}

        def _pts(f: Face) -> np.ndarray:
            k = _face_key(f)
            if k not in pts_cache:
                pts_cache[k] = f.grid_points(
                    B("orig", f.BlockIndex), stride_u=1, stride_v=1)
            return pts_cache[k]

        def _orth_map(f: Face, tol: float) -> Optional[Dict[Tuple[int, int], float]]:
            """Map quantized orthogonal-plane coords -> axis coordinate.

            A face is a valid translational-periodic candidate for `axis`
            only if it is SINGLE-VALUED along the axis over its orthogonal
            footprint (a height field) — true of pitch boundaries, false
            of faces whose normal is orthogonal to the axis (inlet/outlet,
            hub/shroud: their footprint collapses to a line) and of
            wrap-around walls (a blade wall has SS and PS at the same
            (orth1, orth2) with different axis values). Those returned a
            degenerate footprint that produced false matches; reject them
            by counting multi-valued collisions.
            """
            P = _pts(f)
            orth = np.delete(P, axis_idx, axis=1)
            # Intentionally left as round(x/tol) bucketing (not
            # geometry.coincidence_count): this is a cheap PRUNE for
            # footprint-quality classification, not a coincidence decision.
            # The real accept/reject test is the distance-based
            # coincidence_count check in _match_pair below, which is already
            # fixed for the bin-boundary bug.
            keys = np.round(orth / tol).astype(np.int64)
            out: Dict[Tuple[int, int], float] = {}
            ax_vals = P[:, axis_idx]
            n_multi = 0
            for (k1, k2), av in zip(map(tuple, keys), ax_vals):
                prev = out.get((k1, k2))
                if prev is None:
                    out[(k1, k2)] = av
                elif abs(av - prev) > tol:
                    n_multi += 1
            # Footprint-quality filters (a cheap PRUNE — final acceptance is
            # the 99 % Δ-agreement test in _match_pair). Calibrated on the
            # tgs-py cascade mesh:
            #   ratio = unique footprint keys / nodes. Faces whose normal is
            #     ⊥ to the axis collapse to a line (inlet/outlet ≈ 0.06,
            #     constant-z hub/shroud ≈ 0.18); pitch faces ≥ 0.70.
            #   mfrac = multi-valued collisions / keys. Wrap-around blade
            #     walls see SS+PS at the same footprint (≈ 0.34+); pitch
            #     faces stay ≤ 0.19 (steep LE/TE wrap only).
            ratio = len(out) / max(len(P), 1)
            mfrac = n_multi / max(len(out), 1)
            if ratio < 0.5 or mfrac > 0.30:
                return None  # not a height field over the orthogonal plane
            return out

        def _match_pair(fA: Face, fB: Face, tol: Optional[float]):
            """Try to match A onto B by a pure axis translation.

            Returns ``(n_common, frac_small, d_pair)`` or None. Works for
            PARTIAL containment too (e.g. a split inlet-extension face
            matching the corresponding piece of an unsplit pitch face) —
            the shared orthogonal footprint defines the overlap, and Δ is
            the median axis offset over that footprint. The smaller side
            must be essentially fully covered (≥ 95 % of its footprint —
            quantization collisions in steep LE/TE wrap strips cost a few
            % even on exactly-periodic meshes): discovery has to be
            conservative because the result is consumed as an exact
            index-mapped interface by the solver.

            ``tol=None`` means no tolerance could be derived for this pair
            (e.g. a degenerate face) — declines to match rather than
            comparing against a missing tolerance.
            """
            if tol is None:
                return None
            mA = _orth_map(fA, tol)
            mB = _orth_map(fB, tol)
            if mA is None or mB is None:
                return None
            # Footprint overlap gate + Δ estimation. The per-key map's
            # "first stored value" is ambiguous where quantization
            # collapses steep LE/TE strips, so it is used ONLY to gate
            # overlap and estimate Δ (median — robust to those keys).
            common = mA.keys() & mB.keys()
            n_small = min(len(mA), len(mB))
            if len(common) < max(min_shared_abs, int(0.5 * n_small)):
                return None
            diffs = np.array([mB[k] - mA[k] for k in common])
            d_pair = float(np.median(diffs))
            if abs(d_pair) <= tol:
                return None  # coincident along axis — interface, not periodic
            if delta is not None and abs(abs(d_pair) - abs(delta)) > 10 * tol:
                return None  # caller pinned the pitch; reject other shifts
            # Verification: quantized 3D node intersection under the
            # estimated shift. Identical geometry quantizes identically on
            # both sides, so this has none of the first-stored ambiguity —
            # an exactly-periodic overlap scores ~100 %.
            PA = _pts(fA).copy()
            PA[:, axis_idx] += d_pair
            PB = _pts(fB)
            # Dedup by quantized position (a coarse tol-grid is fine here --
            # this only decides which points are "the same node" within one
            # side, it is not the coincidence decision between sides). The
            # actual accept/reject coincidence test below is distance-based.
            _, idxA = np.unique(np.round(PA / tol), axis=0, return_index=True)
            _, idxB = np.unique(np.round(PB / tol), axis=0, return_index=True)
            QA = PA[idxA]
            QB = PB[idxB]
            shared = coincidence_count(QA, QB, tol)
            n_3d_small = min(len(QA), len(QB))
            need = max(min_shared_abs, int(0.95 * n_3d_small))
            if shared < need:
                return None
            return int(shared), shared / max(n_3d_small, 1), d_pair

        # Test all remaining pairs; faces may match PARTIALLY (an unsplit
        # pitch face can host several smaller counterparts), so only the
        # fully-covered (smaller) side of each accepted pair is consumed.
        pair_hits = []  # (frac_small, n_common, ia, ib, d_pair, tol_pair)
        for ia in range(len(remaining_faces)):
            for ib in range(ia + 1, len(remaining_faces)):
                fA, fB = remaining_faces[ia], remaining_faces[ib]
                tol_pair = _pair_tol(fA, fB)
                if tol_pair is None:
                    continue
                hit = _match_pair(fA, fB, tol_pair)
                if hit is not None:
                    n_common, frac_small, d_pair = hit
                    pair_hits.append(
                        (frac_small, n_common, ia, ib, d_pair, tol_pair))

        # Best-coverage pairs first; consume fully-matched sides. A face
        # is only retired once it is essentially fully covered (an
        # unsplit pitch face can host several smaller counterparts).
        pair_hits.sort(key=lambda t: (-t[0], -t[1]))
        fully_used: set = set()
        for frac_small, n_common, ia, ib, d_pair, tol_pair in pair_hits:
            fA, fB = remaining_faces[ia], remaining_faces[ib]
            kA, kB = _face_key(fA), _face_key(fB)
            if kA in fully_used or kB in fully_used:
                continue
            # block1 = the smaller (contained) face; Δ maps block1 → block2.
            # frac_small (3D-intersection coverage of the smaller side)
            # ≥ 0.95 retires that side from further pairing.
            if len(_pts(fA)) <= len(_pts(fB)):
                small, large, shift_amt = fA, fB, d_pair
            else:
                small, large, shift_amt = fB, fA, -d_pair
            # `_match_pair`'s footprint-overlap gate is a PRUNE (its own
            # docstring: verification is "the distance-based coincidence_count
            # check"), not full node-for-node certification -- it never
            # located an explicit sub-patch on the (possibly larger) other
            # side. Do that now and certify before accepting.
            blk_small = B("orig", small.BlockIndex)
            blk_large = B("orig", large.BlockIndex)
            loc = _locate_and_certify_periodic_patch(
                small, blk_small, large, blk_large, axis_idx, shift_amt, tol_pair)
            if loc is None:
                continue
            _record_match(small, large, f"{axis}_oblique_pair", loc)
            if frac_small >= 0.95:
                fully_used.add(_face_key(small))

    # 9) scale back up
    scale_face_bounds(periodic_export, gcd_to_use)

    periodic_pairs: List[Tuple[Face, Face, Dict[str,str]]] = []
    for (fL, fU, m) in periodic_pairs_r:
        gL = deepcopy(fL); gU = deepcopy(fU)
        gL.I *= gcd_to_use; gL.J *= gcd_to_use; gL.K *= gcd_to_use
        gU.I *= gcd_to_use; gU.J *= gcd_to_use; gU.K *= gcd_to_use
        periodic_pairs.append((gL, gU, m))

    # 9b) A5: full-resolution re-validation after GCD reduction. Unlike
    # `rotated_periodicity`'s single global forward/backward rotation, each
    # pair here carries its own per-pair axis shift, so a shared
    # `transforms` list doesn't fit -- re-derive each proposal's own
    # `shift_amt`/tolerance from the ORIGINAL full-resolution `blocks`
    # (never `blocks_r`, the internally-reduced copy) and certify inline
    # rather than over-generalizing `revalidate_full_resolution` for this
    # one caller.
    if gcd_to_use > 1:
        def _pair_tol_full(fA: Face, blkA: Block, fB: Face, blkB: Block) -> float:
            """Per-pair tolerance for full-resolution re-validation, mirroring
            `_pair_tol` above but sourced from the FULL-RESOLUTION blocks.
            Falls back to the same 1e-4 floor `_combine_pair_spacings` uses
            when neither face has enough in-plane points to measure -- a
            proposal that reached this point already matched at reduced
            resolution, so re-validation must not silently skip it for lack
            of a derivable tolerance.
            """
            if node_tol_xyz is not None:
                return float(node_tol_xyz)
            sA = _median_inplane_spacing(fA, blkA)
            sB = _median_inplane_spacing(fB, blkB)
            combined = _combine_pair_spacings(sA, sB)
            return combined if combined is not None else 1e-4

        revalidated_export: List[Dict[str, Dict[str, int]]] = []
        revalidated_pairs: List[Tuple[Face, Face, Dict[str,str]]] = []
        rejected: List[dict] = []
        rejected_failures: List[MappingFailure] = []

        for exp, pair in zip(periodic_export, periodic_pairs):
            b1, b2 = exp['block1'], exp['block2']
            bi1, bi2 = b1['block_index'], b2['block_index']
            block1_full = blocks[bi1]
            block2_full = blocks[bi2]
            try:
                patch1 = patch_from_bounds(bi1, b1['lb'], b1['ub'])
                patch2 = patch_from_bounds(bi2, b2['lb'], b2['ub'])
            except ValueError:
                rejected.append(exp)
                continue

            fL_full, fU_full = pair[0], pair[1]
            tol_pair = _pair_tol_full(fL_full, block1_full, fU_full, block2_full)

            # Recover this proposal's own +/- d_axis shift sign the same
            # way it was originally derived in step 8: compare the axis
            # coordinate of each patch's lb corner on the (unshifted)
            # full-resolution blocks. `d_axis` is a physical distance, so
            # it is unaffected by GCD reduction and reusable as-is here.
            arr1 = [block1_full.X, block1_full.Y, block1_full.Z][axis_idx]
            arr2 = [block2_full.X, block2_full.Y, block2_full.Z][axis_idx]
            p1_val = arr1[tuple(b1['lb'])]
            p2_val = arr2[tuple(b2['lb'])]
            shift_amt = d_axis if p1_val < p2_val else -d_axis
            transform = _shift_transform(axis_idx, shift_amt)

            declared_idx = _declared_perm_idx(
                exp.get('orientation'), b1['lb'], b1['ub'], b2['lb'], b2['ub'])
            try:
                if declared_idx is not None:
                    certify_permutation(
                        block1_full, patch1, block2_full, patch2,
                        declared_idx, tol_pair, transform=transform)
                else:
                    certify_correspondence(
                        block1_full, patch1, block2_full, patch2,
                        tol_pair, transform=transform)
            except MappingFailure as exc:
                rejected.append(exp)
                rejected_failures.append(exc)
                continue

            revalidated_export.append(exp)
            revalidated_pairs.append(pair)

        periodic_export = revalidated_export
        periodic_pairs = revalidated_pairs

        if rejected:
            # Rebind the local `outer_faces` name (never mutate the
            # caller's own list) so the un-touched step 10 filter below
            # naturally keeps demoted proposals as outer faces.
            outer_faces = list(outer_faces)
            demote_to_outer(outer_faces, rejected)
            if rejected_failures:
                worst = max(rejected_failures, key=_failure_severity)
                worst_desc = _describe_failure(worst)
            else:
                worst_desc = "no diagnostic available"
            warnings.warn(
                f"translational_periodicity: demoted {len(rejected)} "
                f"proposal(s) that failed full-resolution re-certification "
                f"(worst: {worst_desc})",
                RuntimeWarning, stacklevel=2)

    # 10) remove periodic from outer_faces (keep 'id' on remaining)
    periodic_keys = set()
    for rec in periodic_export:
        for side in ("block1","block2"):
            bi = rec[side]["block_index"]
            key = (bi, tuple(rec[side]["lb"]), tuple(rec[side]["ub"]))
            periodic_keys.add(key)

    outer_faces_remaining = []
    for o in outer_faces:
        key = (o["block_index"], tuple(o["lb"]), tuple(o["ub"]))
        if key not in periodic_keys:
            outer_faces_remaining.append(o)

    if use_minmax:
        from .connectivity import normalize_face_matches
        periodic_export = normalize_face_matches(periodic_export)

    return periodic_export, periodic_pairs, outer_faces_remaining

def linear_real_transform(face1:Face,face2:Face) -> Tuple:
    """Computes the rotation angle from Face1 to Face2. This can be used to check if the faces are periodic 
        This function assumes the rotation axis is in the "x" direction. This is good for faces within the same block 

    Reference:
        - Linear Real Transforms (M_ccMBMesh.F, computeLRT)
        
    Args:
        Face1 (Face): Face to rotate
        Face2 (Face): Face to rotate to

    Returns:
        (tuple): tuple containing:

            - **ang** (float): rotation angle
            - **rotation_matrix** (numpy.ndarray): Rotation matrix 3x3
        
    """

    cTo3,cTo1 = face1.get_corners() 
    cFrom3,cFrom1 = face2.get_corners()

    dTo  = np.array(cTo3).transpose() - np.array(cTo1).transpose()                      # difference in corner points = diagonal vector for Face 1
    ldTo=np.sqrt(np.sum(dTo*dTo))
    if ldTo > 0:
        dTo=dTo/ldTo
    
    dFrom = np.array(cFrom3).transpose() - np.array(cFrom1).transpose()                 # difference in corner points = diagonal vector for Face 2
    ldFrom = np.sqrt(np.sum(dFrom*dFrom))
    if( ldFrom > 0 ):
        dFrom=dFrom/ldFrom
    
    dotprod = np.sum(dTo * dFrom)
    
    if( abs(dotprod-1) < 1E-10 ): # Case of no rotation
        ang = 0

        rotation_matrix  = np.zeros(shape=(3,3))
    else:
        #Compute the angle of rotation  
        yz_mag_to = sqrt(dTo[1]*dTo[1]+dTo[2]*dTo[2])
        yz_mag_from = sqrt(dFrom[1]*dFrom[1]+dFrom[2]*dFrom[2])
        if yz_mag_to < 1e-15 or yz_mag_from < 1e-15:
            cosAng = 1.0
            sinAng = 0.0
        else:
            cosAng=(dTo[1]*dFrom[1]+dTo[2]*dFrom[2])/(yz_mag_to*yz_mag_from)
            sinAng=(dTo[2]*dFrom[1]-dTo[1]*dFrom[2])/(yz_mag_to*yz_mag_from)
        ang=acos(cosAng)

        rotation_matrix = [ [1, 0, 0],
                            [0, cosAng, -sinAng],
                            [0, sinAng, cosAng] ]
        if( sinAng < 0 ):
            ang*=-1
    return ang, rotation_matrix

def _extract_face_points(face: Face, block: Block):
    """Extract all (x,y,z) points and their (i,j,k) indices from a face.

    Returns:
        pts (np.ndarray): (N, 3) array of xyz coordinates
        ijk (List[Tuple[int,int,int]]): parallel list of (i,j,k) indices
    """
    rows = []
    ijk = []
    if face.IMIN == face.IMAX:
        ic = face.IMIN
        for j in range(face.JMIN, face.JMAX + 1):
            for k in range(face.KMIN, face.KMAX + 1):
                rows.append([block.X[ic, j, k], block.Y[ic, j, k], block.Z[ic, j, k]])
                ijk.append((ic, j, k))
    elif face.JMIN == face.JMAX:
        jc = face.JMIN
        for i in range(face.IMIN, face.IMAX + 1):
            for k in range(face.KMIN, face.KMAX + 1):
                rows.append([block.X[i, jc, k], block.Y[i, jc, k], block.Z[i, jc, k]])
                ijk.append((i, jc, k))
    elif face.KMIN == face.KMAX:
        kc = face.KMIN
        for i in range(face.IMIN, face.IMAX + 1):
            for j in range(face.JMIN, face.JMAX + 1):
                rows.append([block.X[i, j, kc], block.Y[i, j, kc], block.Z[i, j, kc]])
                ijk.append((i, j, kc))
    return np.array(rows), ijk


def __periodicity_check__(face1:Face, face2:Face,block1:Block,block2:Block,tol:float=1E-4):
    """Check if two faces are periodic using cKDTree geometric matching.

    Uses a spatial tree to robustly match rotated coordinates, handling
    different-sized faces and floating-point error from rotation matrices.

    Args:
        face1 (Face): An arbitrary face
        face2 (Face): An arbitrary face
        block1 (Block): block 1 corresponding to face 1 (already rotated)
        block2 (Block): block 2 corresponding to face 2
        tol (float): Matching tolerance for point distances

    Returns:
        (tuple): containing

            - **df** (pandas.DataFrame): Point matches with columns i1,j1,k1,i2,j2,k2
            - **periodic_surface** (List[Face]): Faces that are periodic
            - **split_surfaces** (List[Face]): Split faces to be treated as outer faces
    """
    periodic_faces = list()
    split_faces_out = list()
    swapped = False
    if (face2.diagonal_length < face1.diagonal_length):
        temp = deepcopy(face1)
        face1 = deepcopy(face2)
        face2 = temp

        temp_block = deepcopy(block1)
        block1 = deepcopy(block2)
        block2 = temp_block
        swapped = True

    # Extract xyz points with index tracking
    pts1, ijk1 = _extract_face_points(face1, block1)
    pts2, ijk2 = _extract_face_points(face2, block2)

    if len(pts1) == 0 or len(pts2) == 0:
        return pd.DataFrame(), periodic_faces, split_faces_out

    # cKDTree matching: find face1 points on face2
    tree = cKDTree(pts2)
    dists, indices = tree.query(pts1, k=1)
    mask = dists < tol

    # Build match DataFrame
    match_records = []
    for idx in range(len(pts1)):
        if mask[idx]:
            i1, j1, k1 = ijk1[idx]
            i2, j2, k2 = ijk2[indices[idx]]
            match_records.append({'i1': i1, 'j1': j1, 'k1': k1,
                                  'i2': i2, 'j2': j2, 'k2': k2})

    df = pd.DataFrame(match_records, columns=['i1','j1','k1','i2','j2','k2'])

    if len(df) >= 4:
        # Check it's not just an edge
        n_const = (int(df['i1'].min() == df['i1'].max()) +
                   int(df['j1'].min() == df['j1'].max()) +
                   int(df['k1'].min() == df['k1'].max()))
        if n_const >= 2:
            # It's an edge or degenerate, not a face
            return pd.DataFrame(), periodic_faces, split_faces_out

        ilb1, jlb1, klb1 = int(df['i1'].min()), int(df['j1'].min()), int(df['k1'].min())
        iub1, jub1, kub1 = int(df['i1'].max()), int(df['j1'].max()), int(df['k1'].max())
        ilb2, jlb2, klb2 = int(df['i2'].min()), int(df['j2'].min()), int(df['k2'].min())
        iub2, jub2, kub2 = int(df['i2'].max()), int(df['j2'].max()), int(df['k2'].max())

        # Completeness check: the KDTree nearest-neighbor match must cover
        # every point of the claimed sub-patch, not just a partial subset
        # (e.g. two faces that only share an edge/corner strip can still
        # produce >=4 non-collinear point pairs and pass the edge check
        # above). Mirrors connectivity.py's get_face_intersection, which
        # applies the exact same standard: reject iff the matched-point
        # count is less than the sub-patch's index-space area.
        matched_area = _face_point_count([ilb1, jlb1, klb1], [iub1, jub1, kub1])
        if matched_area > 0 and len(df) < matched_area:
            return pd.DataFrame(), periodic_faces, split_faces_out

        # Full node-for-node certification: confirm every node of the
        # claimed sub-patch corresponds under exactly one of the 8
        # structured permutations, not just the KDTree's nearest-neighbor
        # proposal. `block1` is already the pre-rotated copy that produced
        # the KDTree match above (the caller rotates one side before
        # calling this function; the swap earlier in this function keeps
        # `block1`/`block2` and `pts1`/`pts2` in lockstep), so both patches
        # are already expressed in the same frame -- no further transform.
        patch1 = patch_from_bounds(face1.blockIndex, [ilb1, jlb1, klb1], [iub1, jub1, kub1])
        patch2 = patch_from_bounds(face2.blockIndex, [ilb2, jlb2, klb2], [iub2, jub2, kub2])
        try:
            certify_correspondence(block1, patch1, block2, patch2, tol, transform=None)
        except MappingFailure:
            return pd.DataFrame(), periodic_faces, split_faces_out

        # Create Face objects from matched region
        f1 = create_face_from_diagonals(block1,
            [ilb1, jlb1, klb1],
            [iub1, jub1, kub1])
        f1.set_block_index(face1.blockIndex)
        f1.set_face_id(face1.id)

        f2 = create_face_from_diagonals(block2,
            [ilb2, jlb2, klb2],
            [iub2, jub2, kub2])
        f2.set_block_index(face2.blockIndex)
        f2.set_face_id(face2.id)

        # Split faces for partial matches (matched region smaller than original face)
        # Face 1 splits
        ilb1, jlb1, klb1 = int(df['i1'].min()), int(df['j1'].min()), int(df['k1'].min())
        iub1, jub1, kub1 = int(df['i1'].max()), int(df['j1'].max()), int(df['k1'].max())
        if (ilb1 != face1.IMIN or iub1 != face1.IMAX or
            jlb1 != face1.JMIN or jub1 != face1.JMAX or
            klb1 != face1.KMIN or kub1 != face1.KMAX):
            if int(ilb1==iub1) + int(jlb1==jub1) + int(klb1==kub1) == 1:
                main_face1 = create_face_from_diagonals(block1,
                    [face1.IMIN, face1.JMIN, face1.KMIN],
                    [face1.IMAX, face1.JMAX, face1.KMAX])
                sf1 = split_face(main_face1, block1,
                    ilb=ilb1, jlb=jlb1, klb=klb1, iub=iub1, jub=jub1, kub=kub1)
                [s.set_block_index(face1.blockIndex) for s in sf1]
                [s.set_face_id(face1.id) for s in sf1]
                split_faces_out.extend(sf1)

        # Face 2 splits
        ilb2, jlb2, klb2 = int(df['i2'].min()), int(df['j2'].min()), int(df['k2'].min())
        iub2, jub2, kub2 = int(df['i2'].max()), int(df['j2'].max()), int(df['k2'].max())
        if (ilb2 != face2.IMIN or iub2 != face2.IMAX or
            jlb2 != face2.JMIN or jub2 != face2.JMAX or
            klb2 != face2.KMIN or kub2 != face2.KMAX):
            if int(ilb2==iub2) + int(jlb2==jub2) + int(klb2==kub2) == 1:
                main_face2 = create_face_from_diagonals(block2,
                    [face2.IMIN, face2.JMIN, face2.KMIN],
                    [face2.IMAX, face2.JMAX, face2.KMAX])
                sf2 = split_face(main_face2, block2,
                    ilb=ilb2, jlb=jlb2, klb=klb2, iub=iub2, jub=jub2, kub=kub2)
                [s.set_block_index(face2.blockIndex) for s in sf2]
                [s.set_face_id(face2.id) for s in sf2]
                split_faces_out.extend(sf2)

        if swapped:
            df = df.rename(columns={
                'i1': '_i2', 'j1': '_j2', 'k1': '_k2',
                'i2': 'i1',  'j2': 'j1',  'k2': 'k1',
            }).rename(columns={
                '_i2': 'i2', '_j2': 'j2', '_k2': 'k2',
            })
            periodic_faces.append(f2)
            periodic_faces.append(f1)
        else:
            periodic_faces.append(f1)
            periodic_faces.append(f2)

    return df, periodic_faces, split_faces_out


