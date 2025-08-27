from typing import Dict, List, Optional, Tuple
from .listfunctions import unique_pairs
from .block import Block, reduce_blocks
from .face import Face 
from copy import deepcopy
import numpy.typing as npt
import numpy as np
import math

def faces_match(face1: Tuple[npt.NDArray, npt.NDArray, npt.NDArray],face2: Tuple[npt.NDArray, npt.NDArray, npt.NDArray],tol: float = 1e-12) -> Tuple[bool, Optional[Tuple[bool, bool]]]:
    """
    Compare two block faces and return whether they match and the flip required on face2 to match face1.
    Returns (True, (flip_ud, flip_lr)) if matching, otherwise (False, None).
    """
    def get_corners(X, Y, Z):
        return np.array([
            [X[0, 0], Y[0, 0], Z[0, 0]],
            [X[0, -1], Y[0, -1], Z[0, -1]],
            [X[-1, 0], Y[-1, 0], Z[-1, 0]],
            [X[-1, -1], Y[-1, -1], Z[-1, -1]],
        ])

    X1, Y1, Z1 = face1
    X2, Y2, Z2 = face2

    if X1.shape != X2.shape:
        return False, None

    corners1 = get_corners(X1, Y1, Z1)
    for flip_ud in [False, True]:
        for flip_lr in [False, True]:
            X2f, Y2f, Z2f = X2.copy(), Y2.copy(), Z2.copy()
            if flip_ud:
                X2f, Y2f, Z2f = np.flip(X2f, axis=0), np.flip(Y2f, axis=0), np.flip(Z2f, axis=0)
            if flip_lr:
                X2f, Y2f, Z2f = np.flip(X2f, axis=1), np.flip(Y2f, axis=1), np.flip(Z2f, axis=1)

            corners2 = get_corners(X2f, Y2f, Z2f)
            diffs = np.linalg.norm(corners1 - corners2, axis=1)
            if np.all(diffs <= tol):
                return True, (flip_ud, flip_lr)

    return False, None
def find_matching_faces(block1, block2, tol=1e-8):
    """
    Returns a tuple (face1_name, face2_name, flip_flags) if a matching face is found.
    Otherwise returns (None, None, None).
    """
    faces1 = block1.get_faces()
    faces2 = block2.get_faces()
    for face1_name, face1_data in faces1.items():
        for face2_name, face2_data in faces2.items():
            match, flip_flags = faces_match(face1_data, face2_data, tol=tol)
            if match:
                return face1_name, face2_name, flip_flags
    return None, None, None


def get_outer_faces(block1:Block):
    """Get the outer faces of a block

    Args:
        block1 (Block): A plot3D block

    Returns:
        List[Face]: Non matching faces of the block 
        List[(Face,Face)]: Matching faces inside the block 
    """
    I = [0,block1.IMAX-1]               # Python index starts at 0, need to subtract 1 for it to get the i,j,k
    J = [0,block1.JMAX-1]
    K = [0,block1.KMAX-1]
    # Create the outer faces        
    faces = list()
    face = Face(4)
    i=I[0]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    i=I[1]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    j=J[0]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    j=J[1]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    k=K[0]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    k=K[1]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    faces.append(face)

    # Check if faces match each other
    matching = list()
    non_matching = list()
    for i in range(len(faces)):
        matchFound = False
        for j in range(len(faces)):
            if (i!=j and faces[i].vertices_equals(faces[j])):
                matching.append((i,j))
                matchFound = True
        if not matchFound:
            non_matching.append(faces[i]) # these are guaranteed to be exterior 
    matching = list(unique_pairs(matching))
    matching = [(faces[i],faces[j]) for i,j in matching]
    
    # Make sure normals do not intersect 
    # block_center_to_face_center =  block1.cx
    return non_matching, matching # these should be the outer faces
    
def create_face_from_diagonals(block:Block,imin:int,jmin:int,kmin:int,imax:int,jmax:int,kmax:int) -> Face:
    """Creates a face on a block given a the diagonals defined as (IMIN,JMIN,KMIN), (IMAX, JMAX, KMAX)

    Args:
        block (Block): Block to create a face on 
        imin (int): Lower Corner IMIN
        jmin (int): Lower Corner JMIN
        kmin (int): Lower Corner KMIN
        imax (int): Upper Corner IMAX
        jmax (int): Upper Corner JMAX
        kmax (int): Upper Corner

    Returns:
        (Face): Face created from diagonals 
    """
    newFace = Face(4)           # This is because two of the corners either imin or imax can be equal
    if imin==imax:
        i = imin
        for j in [jmin,jmax]:
            for k in [kmin,kmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    elif jmin==jmax:
        j = jmin
        for i in [imin,imax]:
            for k in [kmin,kmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    elif kmin==kmax:
        k = kmin
        for i in [imin,imax]:
            for j in [jmin,jmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    return newFace


AxisMap = {"x": 0, "y": 1, "z": 2}

def _face_axis_extreme(face: Face, axis: str) -> Tuple[float, float]:
    """Return (vmin, vmax) of the face along axis ∈ {'x','y','z'} using stored vertices."""
    if axis == "x":
        arr = face.x[:face.nvertex]
    elif axis == "y":
        arr = face.y[:face.nvertex]
    else:
        arr = face.z[:face.nvertex]
    return float(np.min(arr)), float(np.max(arr))

def _global_axis_extreme(blocks: List[Block], axis: str) -> Tuple[float, float]:
    """Global (min, max) along axis across all blocks."""
    vals = []
    idx = AxisMap[axis]
    for b in blocks:
        if idx == 0: vals.append(b.X.reshape(-1))
        elif idx == 1: vals.append(b.Y.reshape(-1))
        else: vals.append(b.Z.reshape(-1))
    cat = np.concatenate(vals, axis=0)
    return float(cat.min()), float(cat.max())

def _select_seed_faces(outer_faces: List[Face], blocks: List[Block],
                       axis: str, side: str, tol_abs: float) -> List[Face]:
    """Pick all outer faces whose face extreme equals the global extreme within tol."""
    gmin, gmax = _global_axis_extreme(blocks, axis)
    target = gmin if side == "min" else gmax
    seeds: List[Face] = []
    for f in outer_faces:
        fmin, fmax = _face_axis_extreme(f, axis)
        face_ext = fmin if side == "min" else fmax
        if abs(face_ext - target) <= tol_abs:
            seeds.append(f)
    return seeds

def _bfs_collect_boundary(seed_faces: List[Face],
                          all_outer_faces: List[Face],
                          blocks: List[Block],
                          axis: str,
                          side: str,
                          tol_abs: float,
                          node_tol_xyz: float,
                          min_shared_abs: int = 2,
                          min_shared_frac: float = 0.005) -> List[Face]:
    """
    Flood-fill over outer faces: neighbors share nodes (not just centroids)
    and lie on the same extreme plane (within tol_abs).
    """
    # index faces by (block, IMIN..KMAX) for visited bookkeeping
    def key(f: Face) -> Tuple[int,int,int,int,int,int,int]:
        return (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)

    # all faces that lie on the extreme plane
    on_plane = []
    for f in all_outer_faces:
        fmin, fmax = _face_axis_extreme(f, axis)
        v = fmin if side == "min" else fmax
        # Require the *entire* face to be on/above(below) the plane within tol to avoid picking X side
        gmin, gmax = _global_axis_extreme(blocks, axis)
        plane = gmin if side == "min" else gmax
        # Face must touch the plane and not protrude past it more than tol_abs
        touch_plane = abs(v - plane) <= tol_abs
        # The opposite extreme should not "cross" the plane the wrong way
        opp = fmax if side == "min" else fmin
        not_past = (opp - plane <= tol_abs) if side == "min" else (plane - opp <= tol_abs)
        if touch_plane and not_past:
            on_plane.append(f)

    pool = {key(f): f for f in on_plane}
    q = [f for f in seed_faces if key(f) in pool]
    visited = set()
    result: List[Face] = []

    while q:
        cur = q.pop()
        kcur = key(cur)
        if kcur in visited:
            continue
        visited.add(kcur)
        result.append(cur)

        # find neighbors by node-sharing (fast and reliable for structured grids)
        for k2, cand in list(pool.items()):
            if k2 in visited or k2 == kcur:
                continue
            if cur.touches_by_nodes(cand,
                                    blocks[cur.BlockIndex], blocks[cand.BlockIndex],
                                    tol_xyz=node_tol_xyz,
                                    min_shared_abs=min_shared_abs,
                                    min_shared_frac=min_shared_frac,
                                    stride_u=1, stride_v=1):
                q.append(cand)

    return result

def find_bounding_faces(blocks: List[Block],
                        outer_faces: List[Dict[str,int]],
                        direction: str = "z",
                        side: str = "both",
                        tol_rel: float = 1e-8,
                        node_tol_xyz: float = 1e-6) -> Tuple[List[Dict[str,int]], List[Dict[str,int]],List[Face], List[Face]]:
    """
    Find *outer* bounding faces at the global min/max of a given direction ('x','y','z').
    Uses node-sharing BFS to gather continuous boundary patches across blocks.

    Args:
        blocks: list of Block
        outer_faces: optional precomputed outer faces (dict form)
        direction: 'x','y','z'
        side: 'both' | 'min' | 'max'
        tol_rel: relative tolerance against global extreme value
        node_tol_xyz: absolute node matching tolerance (geometry units)

    Returns:
        lower_connected_faces_export, upper_connected_faces_export, lower_connected_faces, upper_connected_faces
    """
    # 1) Reduce by GCD so grids line up
    gcd_array = [math.gcd(b.IMAX-1, math.gcd(b.JMAX-1, b.KMAX-1)) for b in blocks]
    gcd_to_use = min(gcd_array)
    blocks_r = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # 2) Build outer face list at reduced resolution
    if len(outer_faces) == 0:
        outer_faces_all: List[Face] = []
        for bi, b in enumerate(blocks_r):
            outs, _ = get_outer_faces(b)
            for o in outs:
                o.set_block_index(bi)
            outer_faces_all.extend(outs)
    else:
        outer_faces_all = outer_face_dict_to_list(blocks_r, outer_faces, gcd_to_use)

    # 3) Absolute tolerance for plane selection
    gmin, gmax = _global_axis_extreme(blocks_r, direction)
    # scale tolerance by magnitude so meshes near origin work too
    tol_abs = max(1.0, abs(gmin) + abs(gmax)) * tol_rel

    # 4) Seeds on min/max planes
    want_min = (side in ("min", "both"))
    want_max = (side in ("max", "both"))

    lower_connected_faces: List[Face] = []
    upper_connected_faces: List[Face] = []

    if want_min:
        seeds_min = _select_seed_faces(outer_faces_all, blocks_r, direction, "min", tol_abs)
        lower_connected_faces = _bfs_collect_boundary(
            seed_faces=seeds_min,
            all_outer_faces=outer_faces_all,
            blocks=blocks_r,
            axis=direction,
            side="min",
            tol_abs=tol_abs,
            node_tol_xyz=node_tol_xyz,
            min_shared_abs=2,          # edge-wise continuity OK
            min_shared_frac=0.005      # small overlap fraction enough to chain
        )

    if want_max:
        seeds_max = _select_seed_faces(outer_faces_all, blocks_r, direction, "max", tol_abs)
        upper_connected_faces = _bfs_collect_boundary(
            seed_faces=seeds_max,
            all_outer_faces=outer_faces_all,
            blocks=blocks_r,
            axis=direction,
            side="max",
            tol_abs=tol_abs,
            node_tol_xyz=node_tol_xyz,
            min_shared_abs=2,
            min_shared_frac=0.005
        )

    # 5) Scale indices back up to the original grid
    #    (Make copies so we don't mutate faces held elsewhere.)
    def _rescale_faces(faces: List[Face]) -> List[Face]:
        out: List[Face] = []
        for f in faces:
            g = Face()
            g.x = f.x.copy(); g.y = f.y.copy(); g.z = f.z.copy()
            g.I = (f.I * gcd_to_use).astype(f.I.dtype)
            g.J = (f.J * gcd_to_use).astype(f.J.dtype)
            g.K = (f.K * gcd_to_use).astype(f.K.dtype)
            g.nvertex = f.nvertex
            g.cx, g.cy, g.cz = f.cx, f.cy, f.cz
            g.blockIndex = f.blockIndex
            g.id = f.id
            out.append(g)
        # de-duplicate by index ranges
        uniq = {}
        for f in out:
            k = (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)
            uniq[k] = f
        return list(uniq.values())

    lower_connected_faces = _rescale_faces(lower_connected_faces)
    upper_connected_faces = _rescale_faces(upper_connected_faces)

    lower_connected_faces_export = [f.to_dict() for f in lower_connected_faces]
    upper_connected_faces_export = [f.to_dict() for f in upper_connected_faces]

    return (lower_connected_faces_export,
            upper_connected_faces_export,
            lower_connected_faces,
            upper_connected_faces)

def find_closest_block(blocks:List[Block],x:np.ndarray,y:np.ndarray,z:np.ndarray,centroid:np.ndarray,translational_direction:str="x",minvalue:bool=True):
    """Find the closest block to an extreme in the x,y, or z direction and returns the targeting point. 
    Target point is the reference point where we want the closest block and the closest face 

    Args:
        x (np.ndarray): x coordinate of all the blocks' centroid
        y (np.ndarray): y coordinate of all the blocks' centroid
        z (np.ndarray): z coordinate of all the blocks' centroid
        centroid (np.ndarray): centroid (cx,cy,cz)
        translational_direction (str, optional): _description_. Defaults to "x".
        minvalue (bool, optional): _description_. Defaults to True.

    Returns:
        (tuple): containing

            *selected block index* (int): index of closest block
            *target_x* (float): this is the x value where selected block is closest to
            *target_y* (float): this is the y value where selected block is closest to
            *target_z* (float): this is the z value where selected block is closest to
    """
    cx = centroid[0]; cy = centroid[1]; cz = centroid[2]
    target_x = cx; target_y = cy; target_z = cz
    if translational_direction=="x":
        xmins = [b.X.min() for b in blocks]
        xmaxes = [b.X.max() for b in blocks]
        xmin = min(xmins)
        xmax = max(xmaxes)
        dx = xmax - xmin
        if minvalue:
            target_x = xmin - dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))            
        else:
            target_x = xmax + dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))            
        target_y= cy
        target_z = cz
    elif translational_direction=="y":
        ymins = [b.Y.min() for b in blocks]
        ymaxes = [b.Y.max() for b in blocks]
        ymin = min(ymins)
        ymax = max(ymaxes)
        dy = ymax-ymin
        if minvalue:
            target_y = ymin - dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        else:
            target_y = ymax + dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        target_x = cx
        target_z = cz
    else: #  translational_direction=="z":
        zmins = [b.Z.min() for b in blocks]
        zmaxes = [b.Z.max() for b in blocks]
        zmin = min(zmins)
        zmax = max(zmaxes)
        dz = zmax - zmin
        if minvalue:
            target_z = zmin - dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))            
        else:
            target_z = zmax + dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))
        target_x = cx
        target_y = cy

    return selected_block_indx,target_x,target_y,target_z 


def split_face(face_to_split:Face, block:Block,imin:int,jmin:int,kmin:int,imax:int,jmax:int,kmax:int):
    """Splits a face with another face within the same block 
        picture the split as a two rectangles inside each other

    Args:
        face_to_split (Face): Face on the block to be split 
        block (Block): Block the split is occuring on 
        imin (int): IMIN index of the split (diagonals)
        jmin (int): JMIN index of the split (diagonals)
        kmin (int): KMIN index of the split (diagonals) 
        imax (int): IMAX index of the split
        jmax (int): JMAX index of the split
        kmax (int): KMAX index of the split 

    :: 

                    left face    top face     right face
         ________       __          __            __
        |   __   |     |  |        |__|          |  |     __
        |  |__|  |     |  |         __           |  |    |__|  face_to_split/center face 
        |________|     |__|        |__|          |__|
                                bottom face

    Returns:
        [List[Faces]]: List of unique faces from the split 
    """
    center_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmin,kmax=kmax)

    if kmin == kmax:
        # In the picture above Horizontal = i, vertical = j
        left_face = create_face_from_diagonals(block,
                imin=face_to_split.IMIN,imax=imin,
                jmin=face_to_split.JMIN,jmax=face_to_split.JMAX,
                kmin=kmin, kmax=kmax)

        
        right_face = create_face_from_diagonals(block,
                imin=imax, imax=face_to_split.IMAX,
                jmin=face_to_split.JMIN, jmax=face_to_split.JMAX,
                kmin=kmin, kmax=kmax)

        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmax,jmax=face_to_split.JMAX,
            kmin=kmin,kmax=kmax)
        
        bottom_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=face_to_split.JMIN,jmax=jmin,
            kmin=kmin,kmax=kmax)  

    elif (imin==imax):
        # In the picture above Horizontal = j, vertical = k
        left_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=face_to_split.JMIN, jmax=jmin,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        right_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmax, jmax=face_to_split.JMAX,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmax,kmax=face_to_split.KMAX)

        bottom_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=kmin)

    elif (jmin==jmax):
        # In the picture above Horizontal = i, vertical = k 
        left_face = create_face_from_diagonals(block,
            imin=face_to_split.IMIN,imax=imin,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        right_face = create_face_from_diagonals(block,
            imin=imax,imax=face_to_split.IMAX,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)
        
        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmax,kmax=face_to_split.KMAX)
        
        bottom_face = create_face_from_diagonals(block,
            imin=imin, imax=imax,
            jmin=jmin, jmax=jmax,
            kmin=face_to_split.KMIN, kmax=kmin)
    
    faces = [top_face,bottom_face,left_face,right_face]
    faces = [f for f in faces if not f.isEdge and not f.index_equals(center_face)] # Remove edges
    [f.set_block_index(face_to_split.blockIndex) for f in faces] 
    return faces 

def find_face_nearest_point(faces:List[Face], x:float,y:float,z:float):
    """Find a face nearest to a given point

    Args:
        blocks (List[Block]): List of blocks
        faces (List[Face]): List of faces
        x (float): x coordinate of a reference point
        y (float): y coordinate of a reference point
        z (float): z coordinate of a reference point
    """
    n = list(range(len(faces)))
    dv = list()
    for i in n:
        dx = x-faces[i].cx
        dy = y-faces[i].cy
        dz = z-faces[i].cz
        dv.append(math.sqrt(dx*dx + dy*dy + dz*dz))
    face_index = np.argmin(np.array(dv))
    return face_index

def outer_face_dict_to_list(blocks:List[Block],outer_faces:List[Dict[str,int]],gcd:int=1) -> List[Face]:
    """Converts a list of dictionary face representations to a list of faces. Use this only for outer faces 

    Args:
        blocks (List[Block]): List of blocks
        outer_faces (List[Dict[str,int]]): List of outer faces represented as a dictionary 
        gcd (int, optional): Greatst common divisor. Defaults to 1.

    Returns:
        List[Face]: List of Face objects 
    """
    outer_faces_all = list()
    for o in outer_faces:
        face = create_face_from_diagonals(blocks[o['block_index']], int(o['IMIN']/gcd), int(o['JMIN']/gcd), 
            int(o['KMIN']/gcd), int(o['IMAX']/gcd), int(o['JMAX']/gcd), int(o['KMAX']/gcd))
        if 'id' in o.keys():
            face.id = o['id']
        face.set_block_index(o['block_index'])
        outer_faces_all.append(face)

    return outer_faces_all

def match_faces_dict_to_list(blocks:List[Block],matched_faces:List[Dict[str,int]],gcd:int=1):
    """Converts a list of dictionaries representing matched faces to a list of Faces

    Args:
        blocks (List[Block]): List of blocks 
        matched_faces (List[Dict[str,int]]): List of matched faces represented as a dictionary 
        gcd (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    matched_faces_all = list() 
    for _,m in enumerate(matched_faces):
        face1 = create_face_from_diagonals(blocks[m['block1']['block_index']], 
                            int(m['block1']['IMIN']/gcd), int(m['block1']['JMIN']/gcd), int(m['block1']['KMIN']/gcd), 
                            int(m['block1']['IMAX']/gcd), int(m['block1']['JMAX']/gcd), int(m['block1']['KMAX']/gcd))
        face2 = create_face_from_diagonals(blocks[m['block2']['block_index']], 
                            int(m['block2']['IMIN']/gcd), int(m['block2']['JMIN']/gcd), int(m['block2']['KMIN']/gcd), 
                            int(m['block2']['IMAX']/gcd), int(m['block2']['JMAX']/gcd), int(m['block2']['KMAX']/gcd))
        face1.set_block_index(m['block1']['block_index'])
        if 'id' in m['block1'].keys():
            face1.id = m['block1']['id']
        face2.set_block_index(m['block2']['block_index'])
        if 'id' in m['block2'].keys():
            face2.id = m['block2']['id']
        matched_faces_all.append(face1)
        matched_faces_all.append(face2)
    return matched_faces_all