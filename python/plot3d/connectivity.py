from .block import Block, compute_gcd, reduce_blocks
from .face import Face
from .facefunctions import create_face_from_diagonals, split_face, get_outer_faces
from .utils import (euclidean_distance, enumerate_unique_corners, scale_face_dict_indices,
    face_key, face_grid_dims, divide_face_dict_indices)
import gc
import math
from itertools import product
from tqdm import trange
import numpy as np
from typing import List, Dict, Tuple
from .point_match import point_match


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

def _try_full_face_corner_match(face1, face2, block1, block2,
                                I1, J1, K1, I2, J2, K2, tol):
    """O(1) full face match using corner coordinates only.

    Compares the 4 corner XYZ coordinates of face1 against the 4 corners of
    face2. If all corners match bijectively within *tol*, the faces are a full
    match and we can skip the expensive point-by-point search entirely.

    Handles all 8 possible orientations (flips/transposes) because matching is
    done by spatial proximity, not by index ordering.

    Args:
        face1, face2: Face objects
        block1, block2: Block objects containing the faces
        I1, J1, K1: [min, max] index ranges for face1
        I2, J2, K2: [min, max] index ranges for face2
        tol (float): matching tolerance

    Returns:
        List of dicts with keys i1,j1,k1,i2,j2,k2 (4 corner rows) if full match, or None.
    """
    # 1. Dimension check — face grid sizes must match (allowing transpose)
    d1 = face_grid_dims(I1[0], I1[1], J1[0], J1[1], K1[0], K1[1])
    d2 = face_grid_dims(I2[0], I2[1], J2[0], J2[1], K2[0], K2[1])
    if d1 != d2:
        return None

    # 2. Extract 4 corner (i,j,k,x,y,z) tuples for each face
    def _face_corners(block, I, J, K):
        """Get the 4 corner vertices of a face, based on which axis is constant."""
        corners = []
        if I[0] == I[1]:       # I-constant face
            i = I[0]
            for j in [J[0], J[1]]:
                for k in [K[0], K[1]]:
                    corners.append((i, j, k, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k]))
        elif J[0] == J[1]:     # J-constant face
            j = J[0]
            for i in [I[0], I[1]]:
                for k in [K[0], K[1]]:
                    corners.append((i, j, k, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k]))
        elif K[0] == K[1]:     # K-constant face
            k = K[0]
            for i in [I[0], I[1]]:
                for j in [J[0], J[1]]:
                    corners.append((i, j, k, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k]))
        else:
            return None        # Not a proper face
        return corners

    c1 = _face_corners(block1, I1, J1, K1)
    c2 = _face_corners(block2, I2, J2, K2)
    if c1 is None or c2 is None or len(c1) != 4 or len(c2) != 4:
        return None

    # 3. Bijective corner matching — each corner of face1 must match a
    #    unique corner of face2 within tolerance
    used = [False, False, False, False]
    match_pairs = []       # list of (c1_idx, c2_idx)
    for ci, (i1, j1, k1, x1, y1, z1) in enumerate(c1):
        best_idx = -1
        best_d = float('inf')
        for cj, (i2, j2, k2, x2, y2, z2) in enumerate(c2):
            if used[cj]:
                continue
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            if d < best_d:
                best_d = d
                best_idx = cj
        if best_idx < 0 or best_d >= tol:
            return None        # No match for this corner → not a full face match
        used[best_idx] = True
        match_pairs.append((ci, best_idx))

    # 4. Build match rows
    rows = []
    for ci, cj in match_pairs:
        rows.append({
            'i1': c1[ci][0], 'j1': c1[ci][1], 'k1': c1[ci][2],
            'i2': c2[cj][0], 'j2': c2[cj][1], 'k2': c2[cj][2],
        })
    return rows


def _verify_match_interior(block1, block2, I1, J1, K1, I2, J2, K2, tol):
    """Verify a corner-based match by checking that interior nodes are spatially coincident.

    After corner matching succeeds, samples the center point of face1 and
    checks whether any node on face2 is within tolerance. This rejects false
    positives where corners happen to coincide but interior nodes diverge
    (e.g. curved surfaces that share corners but bow in different directions).

    Mirrors Rust's ``verify_match_interior``.

    Args:
        block1, block2: Block objects containing the faces
        I1, J1, K1: [min, max] index ranges for face1
        I2, J2, K2: [min, max] index ranges for face2
        tol: matching tolerance

    Returns:
        True if the interior check passes (or face is too small to have interior points)
    """
    # Find center index for face1
    ic1 = (I1[0] + I1[1]) // 2
    jc1 = (J1[0] + J1[1]) // 2
    kc1 = (K1[0] + K1[1]) // 2

    # Only verify if center differs from at least one corner
    is_corner = (ic1 in I1) and (jc1 in J1) and (kc1 in K1)
    if is_corner:
        return True  # Center is a corner, nothing extra to verify

    x_c = block1.X[ic1, jc1, kc1]
    y_c = block1.Y[ic1, jc1, kc1]
    z_c = block1.Z[ic1, jc1, kc1]

    # Check if any node on face2 is within tolerance
    X2 = select_multi_dimensional(block2.X, (I2[0], I2[1]), (J2[0], J2[1]), (K2[0], K2[1]))
    Y2 = select_multi_dimensional(block2.Y, (I2[0], I2[1]), (J2[0], J2[1]), (K2[0], K2[1]))
    Z2 = select_multi_dimensional(block2.Z, (I2[0], I2[1]), (J2[0], J2[1]), (K2[0], K2[1]))

    dists_sq = (X2.ravel() - x_c)**2 + (Y2.ravel() - y_c)**2 + (Z2.ravel() - z_c)**2
    return float(np.min(dists_sq)) < tol * tol


def get_face_intersection(face1:Face,face2:Face,block1:Block,block2:Block,tol:float=1E-6):
    """Get the index of the intersection between two faces located on two different blocks
        Face1 needs to be the smaller face.

    Args:
        face1 (Face): An exterior face
        face2 (Face): An exterior face from a different block
        block1 (Block): block containing face1
        block2 (Block): block containing face2
        tol (float): matching tolerance

    Returns:
        (Tuple): containing

            - (List[Dict]): list of match dicts with keys i1,j1,k1,i2,j2,k2
            - (List[Face]): any split faces from block 1
            - (List[Face]): any split faces from block 2
    """

    match_location = list()
    split_faces1 = list()
    split_faces2 = list()

    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]

    # FAST PATH: Try full face corner match (O(1) — no point_match needed)
    full_match = _try_full_face_corner_match(
        face1, face2, block1, block2,
        I1, J1, K1, I2, J2, K2, tol
    )
    if full_match is not None:
        # Verify interior points to reject false positives (matches Rust verify_match_interior)
        if _verify_match_interior(block1, block2, I1, J1, K1, I2, J2, K2, tol):
            return full_match, split_faces1, split_faces2
        # Fall through to slow path if interior check fails

    # SLOW PATH: Point-by-point matching for partial/split faces
    # Grab the points of Face 1 and Face 2
    X1 = select_multi_dimensional(block1.X, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
    Y1 = select_multi_dimensional(block1.Y, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
    Z1 = select_multi_dimensional(block1.Z, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))

    X2 = select_multi_dimensional(block2.X, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
    Y2 = select_multi_dimensional(block2.Y, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
    Z2 = select_multi_dimensional(block2.Z, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))

    # Determine which axis is constant for face1 and build the index mapper
    # Each branch maps (p, q) -> (i1, j1, k1) differently based on the constant axis
    if I1[0] == I1[1]:
        def _map_face1(p, q): return (I1[0], p + J1[0], q + K1[0])
    elif J1[0] == J1[1]:
        def _map_face1(p, q): return (p + I1[0], J1[0], q + K1[0])
    elif K1[0] == K1[1]:
        def _map_face1(p, q): return (p + I1[0], q + J1[0], K1[0])
    else:
        # Not a proper face (no constant axis)
        return [], split_faces1, split_faces2

    # Build mappers for face2's constant axis -> (i2, j2, k2)
    face2_mappers = []
    if I2[0] == I2[1]:
        face2_mappers.append(lambda p2, q2: (I2[0], p2 + J2[0], q2 + K2[0]))
    if J2[0] == J2[1]:
        face2_mappers.append(lambda p2, q2: (p2 + I2[0], J2[0], q2 + K2[0]))
    if K2[0] == K2[1]:
        face2_mappers.append(lambda p2, q2: (p2 + I2[0], q2 + J2[0], K2[0]))

    # Single loop for all constant-axis cases
    for p in range(X1.shape[0]):
        for q in range(X1.shape[1]):
            block2_match_location = point_match(X1[p,q], Y1[p,q], Z1[p,q], X2, Y2, Z2, tol)
            if sum(block2_match_location) != -2:
                p2 = int(block2_match_location[0])
                q2 = int(block2_match_location[1])
                i1, j1, k1 = _map_face1(p, q)
                for mapper in face2_mappers:
                    i2, j2, k2 = mapper(p2, q2)
                    match_location.append({"i1": i1, "j1": j1, "k1": k1, "i2": i2, "j2": j2, "k2": k2})

    # Checking for split faces
    if len(match_location) >= 4:
        if _check_edge(match_location):
            match_location = []     # If it's an edge
        else:                       # not edge
            # Filter match increasing - This keeps uniqueness
            if I1[0]==I1[1]:
                match_location = _filter_block_increasing(match_location, 'j1')
                match_location = _filter_block_increasing(match_location, 'k1')
            elif J1[0]==J1[1]:
                match_location = _filter_block_increasing(match_location, 'i1')
                match_location = _filter_block_increasing(match_location, 'k1')
            elif K1[0]==K1[1]:
                match_location = _filter_block_increasing(match_location, 'i1')
                match_location = _filter_block_increasing(match_location, 'j1')

            if I2[0]==I2[1]:
                match_location = _filter_block_increasing(match_location, 'j2')
                match_location = _filter_block_increasing(match_location, 'k2')
            elif J2[0]==J2[1]:
                match_location = _filter_block_increasing(match_location, 'i2')
                match_location = _filter_block_increasing(match_location, 'k2')
            elif K2[0]==K2[1]:
                match_location = _filter_block_increasing(match_location, 'i2')
                match_location = _filter_block_increasing(match_location, 'j2')

            # Do a final check after doing all these checks
            if len(match_location) >= 4:
                # Check for Split faces
                ## Block 1
                main_face = create_face_from_diagonals(block1,imin=I1[0],imax=I1[1], jmin=J1[0],jmax=J1[1],kmin=K1[0],kmax=K1[1])
                imin = min(m['i1'] for m in match_location)
                jmin = min(m['j1'] for m in match_location)
                kmin = min(m['k1'] for m in match_location)
                imax = max(m['i1'] for m in match_location)
                jmax = max(m['j1'] for m in match_location)
                kmax = max(m['k1'] for m in match_location)
                if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                    split_faces1 = split_face(main_face,block1,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
                    [s.set_block_index(face1.blockIndex) for s in split_faces1]
                    [s.set_face_id(face1.id) for s in split_faces1]

                ## Block 2
                main_face = create_face_from_diagonals(block2,imin=I2[0],imax=I2[1], jmin=J2[0],jmax=J2[1],kmin=K2[0],kmax=K2[1])
                imin = min(m['i2'] for m in match_location)
                jmin = min(m['j2'] for m in match_location)
                kmin = min(m['k2'] for m in match_location)
                imax = max(m['i2'] for m in match_location)
                jmax = max(m['j2'] for m in match_location)
                kmax = max(m['k2'] for m in match_location)
                if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                    split_faces2 = split_face(main_face,block2,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
                    [s.set_block_index(face2.blockIndex) for s in split_faces2]
                    [s.set_face_id(face2.id) for s in split_faces2]
            else:
                match_location = []

    else:
        match_location = []
    return match_location, split_faces1, split_faces2

def _filter_block_increasing(matches: List[Dict], key1: str) -> List[Dict]:
    """Filters match results to ensure key1 values are monotonically increasing.

    Removes rows where the specified key doesn't form a contiguous increasing
    sequence. This handles the edge-matching case where face2 touches face1
    at two disconnected edges.

    When there are exactly 2 unique values, they are always kept regardless of
    gap size. This handles the case where a small face (e.g. 2 nodes wide after
    GCD reduction) matches a large face — the matching indices on the large face
    may span a wide gap (e.g. [0, 113]) but are still a valid match. The
    _check_edge() call upstream has already verified this isn't a degenerate edge.

    Args:
        matches: list of match dicts with keys i1,j1,k1,i2,j2,k2
        key1: key that should be in increasing order

    Returns:
        Filtered list of match dicts
    """
    if len(matches) == 0:
        return matches

    key1_vals = sorted(set(m[key1] for m in matches))

    if len(key1_vals) <= 1:
        return []  # Edge matching condition

    # With only 2 unique values, contiguity is trivially satisfied — keep all.
    if len(key1_vals) == 2:
        return matches

    key1_vals_to_use = set()
    for i in range(len(key1_vals) - 1):
        if (key1_vals[i + 1] - key1_vals[i]) == 1:
            key1_vals_to_use.add(key1_vals[i])
    # Look backwards
    if (key1_vals[-1] - key1_vals[-2]) == 1:
        key1_vals_to_use.add(key1_vals[-1])

    return [m for m in matches if m[key1] in key1_vals_to_use]


def _check_edge(matches: List[Dict]) -> bool:
    """Check if the match results represent an edge instead of a face.

    Args:
        matches: list of match dicts with keys i1,j1,k1,i2,j2,k2

    Returns:
        True if the intersection is an edge, False if it's a face
    """
    i1_min = min(m['i1'] for m in matches)
    j1_min = min(m['j1'] for m in matches)
    k1_min = min(m['k1'] for m in matches)
    i1_max = max(m['i1'] for m in matches)
    j1_max = max(m['j1'] for m in matches)
    k1_max = max(m['k1'] for m in matches)

    edge_matches = int(i1_min == i1_max) + int(j1_min == j1_max) + int(k1_min == k1_max)
    return edge_matches >= 2
    
def candidate_neighbor_pairs(blocks:List[Block], tol:float=1e-6):
    """Return (i, j) pairs of block indices whose bounding boxes overlap or touch.

    Two blocks can only share a face if their axis-aligned bounding boxes (AABBs)
    overlap or touch. This replaces the old centroid-distance approach which could
    miss neighbors for irregularly shaped blocks (L-shaped, elongated, etc.).

    Args:
        blocks (List[Block]): list of all your blocks
        tol (float): AABB expansion tolerance. Blocks whose bounding boxes are
            within this distance of touching are still considered candidates.

    Returns:
        List[Tuple[int,int]]: candidate block pairs (i, j) with i < j
    """
    n = len(blocks)
    # Precompute axis-aligned bounding boxes: [xmin, xmax, ymin, ymax, zmin, zmax]
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
            # Check AABB overlap/touch with tolerance on each axis
            if (aabbs[i, 1] + tol >= aabbs[j, 0] and
                aabbs[j, 1] + tol >= aabbs[i, 0] and
                aabbs[i, 3] + tol >= aabbs[j, 2] and
                aabbs[j, 3] + tol >= aabbs[i, 2] and
                aabbs[i, 5] + tol >= aabbs[j, 4] and
                aabbs[j, 5] + tol >= aabbs[i, 4]):
                pairs.append((i, j))
    return pairs
    
def connectivity_fast(blocks:List[Block]):
    """GCD-accelerated connectivity detection.

    Down-samples all blocks by their minimum GCD, runs the three-phase
    connectivity algorithm (see :func:`connectivity`), then scales indices
    back to the original resolution.

    Args:
        blocks (List[Block]): Lists of blocks you want to find the connectivity for

    Returns:
        (List[Dict]): All matching faces formatted as a list of
            ``{ 'block1': {'block_index', 'IMIN', ...}, 'block2': ... }``
        (List[Dict]): All exterior surfaces formatted as a list of
            ``{ 'block_index', 'IMIN', ..., 'id' }``

    """
    gcd_to_use = compute_gcd(blocks)
    print(f"gcd to use {gcd_to_use}")
    new_blocks = reduce_blocks([b.copy() for b in blocks],gcd_to_use)

    # Find Connectivity
    face_matches, outer_faces_formatted = connectivity(new_blocks)
    # scale it up
    scale_face_dict_indices(face_matches, gcd_to_use, nested_sides=['block1', 'block2'])
    scale_face_dict_indices(outer_faces_formatted, gcd_to_use)
    return face_matches, outer_faces_formatted

def connectivity(blocks:List[Block]):
    """Detect face connectivity between blocks using a three-phase algorithm.

    **Phase 1-2** (in ``find_matching_blocks``): For each candidate block pair
    (identified by AABB overlap), faces are matched iteratively. Phase 1 uses
    O(1) corner matching for full 1:1 faces. Phase 2 handles partial/split
    faces via point-by-point intersection, creating remnant sub-faces that
    re-enter the pool until convergence.

    **Phase 3** (fresh-face validation): After Phase 1-2 converge, remaining
    unmatched faces are re-checked against *fresh* outer faces of neighbor
    blocks. Axis-aligned bounding boxes (AABBs) computed from **all face
    nodes** (not just corners) are used for fast pre-filtering before calling
    ``get_face_intersection``. This catches edge cases where curved or skewed
    faces have interior extents beyond their corner coordinates.

    Args:
        blocks (List[Block]): List of all blocks in multi-block plot3d mesh

    Returns:
        (List[Dict]): All matching faces formatted as a list of
            ``{ 'block1': {'block_index', 'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX'}, 'block2': ... }``
        (List[Dict]): All exterior surfaces formatted as a list of
            ``{ 'block_index', 'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX', 'id' }``

    """

    face_matches = list()
    # Use a set of (block_index, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX) for fast removal lookups
    matches_to_remove_keys = set()
    temp = [get_outer_faces(b) for b in blocks]
    block_outer_faces = [t[0] for t in temp]
    del temp
    combos = candidate_neighbor_pairs(blocks) # Find all block pairs whose bounding boxes touch/overlap

    gc_interval = max(1, len(combos) // 20)  # gc every ~5% of progress
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
            for match_rows in df_matches:
                face1 = create_face_from_diagonals(block=blocks[i],
                    imin=min(m['i1'] for m in match_rows), jmin=min(m['j1'] for m in match_rows), kmin=min(m['k1'] for m in match_rows),
                    imax=max(m['i1'] for m in match_rows), jmax=max(m['j1'] for m in match_rows), kmax=max(m['k1'] for m in match_rows))
                face1.set_block_index(i)
                matches_to_remove_keys.add((i, face1.IMIN, face1.JMIN, face1.KMIN, face1.IMAX, face1.JMAX, face1.KMAX))

                face2 = create_face_from_diagonals(block=blocks[j],
                    imin=min(m['i2'] for m in match_rows), jmin=min(m['j2'] for m in match_rows), kmin=min(m['k2'] for m in match_rows),
                    imax=max(m['i2'] for m in match_rows), jmax=max(m['j2'] for m in match_rows), kmax=max(m['k2'] for m in match_rows))
                face2.set_block_index(j)
                matches_to_remove_keys.add((j, face2.IMIN, face2.JMIN, face2.KMIN, face2.IMAX, face2.JMAX, face2.KMAX))

                face_matches.append(face_matches_to_dict(face1,face2,blocks[i],blocks[j]))
        # Periodic garbage collection to prevent memory buildup
        if indx % gc_interval == 0:
            gc.collect()

    # ===== PHASE 3: Fresh-face validation for remaining outer faces =====
    # Some outer faces remain unmatched because the matching block's face pool
    # was consumed by an earlier combo in Phases 1-2. Re-check each remaining
    # outer face against *fresh* (un-consumed) outer faces of overlapping blocks.
    # Build neighbor adjacency from combos
    neighbors = [[] for _ in range(len(blocks))]
    for (ci, cj) in combos:
        neighbors[ci].append(cj)
        neighbors[cj].append(ci)

    # Collect remaining outer faces (not yet matched)
    remaining = []
    for face_list in block_outer_faces:
        for o in face_list:
            key = (o.BlockIndex, o.IMIN, o.JMIN, o.KMIN, o.IMAX, o.JMAX, o.KMAX)
            if key not in matches_to_remove_keys:
                remaining.append(o)

    # Precompute fresh outer faces and their all-node AABBs for every block
    fresh_all = []
    fresh_aabbs = []
    for bi, blk in enumerate(blocks):
        fresh_faces, _ = get_outer_faces(blk)
        fresh_all.append(fresh_faces)
        aabbs_for_block = []
        for f in fresh_faces:
            pts = f.grid_points(blk)
            if len(pts) > 0:
                aabbs_for_block.append([
                    pts[:, 0].min(), pts[:, 0].max(),
                    pts[:, 1].min(), pts[:, 1].max(),
                    pts[:, 2].min(), pts[:, 2].max(),
                ])
            else:
                aabbs_for_block.append([0, 0, 0, 0, 0, 0])
        fresh_aabbs.append(aabbs_for_block)

    phase3_matched_keys = set()
    tol_pre = 0.01
    for face in remaining:
        face_key = (face.BlockIndex, face.IMIN, face.JMIN, face.KMIN, face.IMAX, face.JMAX, face.KMAX)
        if face_key in phase3_matched_keys:
            continue
        bi = face.BlockIndex

        # Compute AABB from all nodes of this face
        pts = face.grid_points(blocks[bi])
        if len(pts) == 0:
            continue
        fxn, fxx = pts[:, 0].min(), pts[:, 0].max()
        fyn, fyx = pts[:, 1].min(), pts[:, 1].max()
        fzn, fzx = pts[:, 2].min(), pts[:, 2].max()

        for bj in neighbors[bi]:
            for fi, ff in enumerate(fresh_all[bj]):
                # AABB pre-check using all-node AABBs
                gaabb = fresh_aabbs[bj][fi]
                if (fxx + tol_pre < gaabb[0] or gaabb[1] + tol_pre < fxn or
                    fyx + tol_pre < gaabb[2] or gaabb[3] + tol_pre < fyn or
                    fzx + tol_pre < gaabb[4] or gaabb[5] + tol_pre < fzn):
                    continue

                match_rows, _, _ = get_face_intersection(face, ff, blocks[bi], blocks[bj])
                if len(match_rows) > 0:
                    face1 = create_face_from_diagonals(
                        block=blocks[bi],
                        imin=min(m['i1'] for m in match_rows), jmin=min(m['j1'] for m in match_rows), kmin=min(m['k1'] for m in match_rows),
                        imax=max(m['i1'] for m in match_rows), jmax=max(m['j1'] for m in match_rows), kmax=max(m['k1'] for m in match_rows)
                    )
                    face1.set_block_index(bi)
                    face2 = create_face_from_diagonals(
                        block=blocks[bj],
                        imin=min(m['i2'] for m in match_rows), jmin=min(m['j2'] for m in match_rows), kmin=min(m['k2'] for m in match_rows),
                        imax=max(m['i2'] for m in match_rows), jmax=max(m['j2'] for m in match_rows), kmax=max(m['k2'] for m in match_rows)
                    )
                    face2.set_block_index(bj)
                    face_matches.append(face_matches_to_dict(face1, face2, blocks[bi], blocks[bj]))
                    matches_to_remove_keys.add((bi, face1.IMIN, face1.JMIN, face1.KMIN, face1.IMAX, face1.JMAX, face1.KMAX))
                    matches_to_remove_keys.add((bj, face2.IMIN, face2.JMIN, face2.KMIN, face2.IMAX, face2.JMAX, face2.KMAX))
                    phase3_matched_keys.add(face_key)
                    # Don't break — continue checking other neighbors for split-face matches

    del fresh_all, fresh_aabbs, neighbors, remaining

    # Update Outer Faces — use dict-based dedup instead of list(set(...))
    outer_faces_dict = dict()
    for face_list in block_outer_faces:
        for o in face_list:
            key = (o.BlockIndex, o.IMIN, o.JMIN, o.KMIN, o.IMAX, o.JMAX, o.KMAX)
            if key not in matches_to_remove_keys:
                outer_faces_dict[key] = o
    del block_outer_faces  # Free memory
    del matches_to_remove_keys

    # Remove any outer faces that may have been found by mistake
    # Check I,J,K if J and K are the same with another outer face, select the face with shorter I
    outer_faces_to_remove_keys = set()
    # Group faces by block index for efficient comparison
    block_face_groups = dict()
    for key, o in outer_faces_dict.items():
        bi = o.BlockIndex
        if bi not in block_face_groups:
            block_face_groups[bi] = []
        block_face_groups[bi].append((key, o))

    for bi, face_pairs in block_face_groups.items():
        for idx_a in range(len(face_pairs)):
            key_a, o = face_pairs[idx_a]
            if key_a in outer_faces_to_remove_keys:
                continue
            ijk = (o.IMIN, o.JMIN, o.KMIN, o.IMAX, o.JMAX, o.KMAX)
            for idx_b in range(idx_a + 1, len(face_pairs)):
                key_b, o2 = face_pairs[idx_b]
                if key_b in outer_faces_to_remove_keys:
                    continue
                ijk2 = (o2.IMIN, o2.JMIN, o2.KMIN, o2.IMAX, o2.JMAX, o2.KMAX)
                # Count matching indices
                matching = sum(1 for a, b in zip(ijk, ijk2) if a == b)
                if matching == 5:  # 5 of 6 indices match -> remove the longer face
                    if o2.diagonal_length > o.diagonal_length:
                        outer_faces_to_remove_keys.add(key_b)
                    else:
                        outer_faces_to_remove_keys.add(key_a)

    outer_faces = [o for key, o in outer_faces_dict.items() if key not in outer_faces_to_remove_keys]
    del outer_faces_dict, block_face_groups, outer_faces_to_remove_keys


    # Find self-matches: Do any faces of, for example, block1 match another face in block 1
    for i in range(len(blocks)):
        _,self_matches = get_outer_faces(blocks[i])
        for match in self_matches: # Append to face matches
            face_matches.append({'block1':{
                                            'block_index':i,'IMIN':int(match[0].I.min()),'JMIN':int(match[0].J.min()),'KMIN':int(match[0].K.min()),
                                            'IMAX':int(match[0].I.max()),'JMAX':int(match[0].J.max()),'KMAX':int(match[0].K.max())
                                        },
                                    'block2':{
                                            'block_index':i,'IMIN':int(match[1].I.min()),'JMIN':int(match[1].J.min()),'KMIN':int(match[1].K.min()),
                                            'IMAX':int(match[1].I.max()),'JMAX':int(match[1].J.max()),'KMAX':int(match[1].K.max())
                                        }
                                    })

    # Update the outer faces
    outer_faces_formatted = list() # This will contain 
    id = 1 
    for face in outer_faces:        
        outer_faces_formatted.append({ 'IMIN':min(face.I), 'JMIN':min(face.J), 'KMIN':min(face.K),
                            'IMAX':max(face.I), 'JMAX':max(face.J), 'KMAX':max(face.K),
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
                            'IMIN':-1,'JMIN':-1,'KMIN':-1,  # Lower Corner
                            'IMAX':-1,'JMAX':-1,'KMAX':-1,   # Upper Corner
                            'id':face1.id
                        },
                'block2':{
                            'block_index':face2.BlockIndex,
                            'IMIN':-1,'JMIN':-1,'KMIN':-1,  # Lower Corner
                            'IMAX':-1,'JMAX':-1,'KMAX':-1,   # Upper Corner
                            'id':face2.id
                        }
                }
            
    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]

    # Search for lower corner match
    x1_l = block1.X[I1[0],J1[0],K1[0]]
    y1_l = block1.Y[I1[0],J1[0],K1[0]]
    z1_l = block1.Z[I1[0],J1[0],K1[0]]
    best_d = float('inf')
    best_ijk = (I2[0], J2[0], K2[0])
    for p in I2:
        for q in J2:
            for r in K2:
                dx = block2.X[p,q,r] - x1_l
                dy = block2.Y[p,q,r] - y1_l
                dz = block2.Z[p,q,r] - z1_l
                d = dx*dx + dy*dy + dz*dz
                if d < best_d:
                    best_d = d
                    best_ijk = (p, q, r)
    match['block1']['IMIN'] = face1.IMIN
    match['block1']['JMIN'] = face1.JMIN
    match['block1']['KMIN'] = face1.KMIN
    match['block2']['IMIN'] = best_ijk[0]
    match['block2']['JMIN'] = best_ijk[1]
    match['block2']['KMIN'] = best_ijk[2]

    # Search for upper corner match
    x1_u = block1.X[I1[1],J1[1],K1[1]]
    y1_u = block1.Y[I1[1],J1[1],K1[1]]
    z1_u = block1.Z[I1[1],J1[1],K1[1]]
    best_d = float('inf')
    best_ijk = (I2[0], J2[0], K2[0])
    for p in I2:
        for q in J2:
            for r in K2:
                dx = block2.X[p,q,r] - x1_u
                dy = block2.Y[p,q,r] - y1_u
                dz = block2.Z[p,q,r] - z1_u
                d = dx*dx + dy*dy + dz*dz
                if d < best_d:
                    best_d = d
                    best_ijk = (p, q, r)
    match['block1']['IMAX'] = face1.IMAX
    match['block1']['JMAX'] = face1.JMAX
    match['block1']['KMAX'] = face1.KMAX
    match['block2']['IMAX'] = best_ijk[0]
    match['block2']['JMAX'] = best_ijk[1]
    match['block2']['KMAX'] = best_ijk[2]
    return match


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
    # Compute GCD and reduce blocks
    gcd_to_use = compute_gcd(blocks)
    reduced_blocks = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # Scale down face_matches indices by GCD
    scaled_matches = [{k: (dict(v) if isinstance(v, dict) else v) for k, v in fm.items()} for fm in face_matches]
    divide_face_dict_indices(scaled_matches, gcd_to_use, nested_sides=['block1', 'block2'])

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
        x1_l = block1.X[b1['IMIN'], b1['JMIN'], b1['KMIN']]
        y1_l = block1.Y[b1['IMIN'], b1['JMIN'], b1['KMIN']]
        z1_l = block1.Z[b1['IMIN'], b1['JMIN'], b1['KMIN']]

        x1_u = block1.X[b1['IMAX'], b1['JMAX'], b1['KMAX']]
        y1_u = block1.Y[b1['IMAX'], b1['JMAX'], b1['KMAX']]
        z1_u = block1.Z[b1['IMAX'], b1['JMAX'], b1['KMAX']]

        # Enumerate unique corners of block2's face
        unique_corners = enumerate_unique_corners(
            [b2['IMIN'], b2['IMAX']], [b2['JMIN'], b2['JMAX']], [b2['KMIN'], b2['KMAX']])

        # Check stored diagonal first
        x2_l = block2.X[b2['IMIN'], b2['JMIN'], b2['KMIN']]
        y2_l = block2.Y[b2['IMIN'], b2['JMIN'], b2['KMIN']]
        z2_l = block2.Z[b2['IMIN'], b2['JMIN'], b2['KMIN']]

        x2_u = block2.X[b2['IMAX'], b2['JMAX'], b2['KMAX']]
        y2_u = block2.Y[b2['IMAX'], b2['JMAX'], b2['KMAX']]
        z2_u = block2.Z[b2['IMAX'], b2['JMAX'], b2['KMAX']]

        d_lower = euclidean_distance(x1_l, y1_l, z1_l, x2_l, y2_l, z2_l)
        d_upper = euclidean_distance(x1_u, y1_u, z1_u, x2_u, y2_u, z2_u)

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

                dl = euclidean_distance(x1_l, y1_l, z1_l,
                                       block2.X[il, jl, kl], block2.Y[il, jl, kl], block2.Z[il, jl, kl])
                du = euclidean_distance(x1_u, y1_u, z1_u,
                                       block2.X[iu, ju, ku], block2.Y[iu, ju, ku], block2.Z[iu, ju, ku])

                if dl < best_d_lower:
                    best_d_lower = dl
                if du < best_d_upper:
                    best_d_upper = du

                if dl < tol and du < tol:
                    corrected = {k: (dict(v) if isinstance(v, dict) else v) for k, v in face_matches[idx].items()}
                    corrected['block2']['IMIN'] = il * gcd_to_use
                    corrected['block2']['JMIN'] = jl * gcd_to_use
                    corrected['block2']['KMIN'] = kl * gcd_to_use
                    corrected['block2']['IMAX'] = iu * gcd_to_use
                    corrected['block2']['JMAX'] = ju * gcd_to_use
                    corrected['block2']['KMAX'] = ku * gcd_to_use
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
                  f"lower=({b1_orig['IMIN']},{b1_orig['JMIN']},{b1_orig['KMIN']}) "
                  f"upper=({b1_orig['IMAX']},{b1_orig['JMAX']},{b1_orig['KMAX']})")
            print(f"  block2 (block_index={b2_orig['block_index']}): "
                  f"lower=({b2_orig['IMIN']},{b2_orig['JMIN']},{b2_orig['KMIN']}) "
                  f"upper=({b2_orig['IMAX']},{b2_orig['JMAX']},{b2_orig['KMAX']})")
            print(f"  block1 lower xyz = ({x1_l:.6e}, {y1_l:.6e}, {z1_l:.6e})")
            print(f"  block1 upper xyz = ({x1_u:.6e}, {y1_u:.6e}, {z1_u:.6e})")
            print(f"  Closest block2 corner dist to block1 lower: {best_d_lower:.6e}")
            print(f"  Closest block2 corner dist to block1 upper: {best_d_upper:.6e}")
            mismatched.append(face_matches[idx])

    return verified, mismatched
