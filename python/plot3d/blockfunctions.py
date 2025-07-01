from copy import deepcopy
from itertools import combinations, product
import math
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from tqdm import trange
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_matching_faces
from .block import Block, reduce_blocks
from .face import Face
import tqdm
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
import numpy.typing as npt 

def rotate_block(block,rotation_matrix:np.ndarray) -> Block:
    """Rotates a block by a rotation matrix 

    Args:
        rotation_matrix (np.ndarray): 3x3 rotation matrix 

    Returns:
        Block: returns a new rotated block 
    """
    X = block.X.copy()
    Y = block.Y.copy()
    Z = block.Z.copy()
    points = np.zeros(shape=(3,block.IMAX*block.JMAX*block.KMAX))
    indx = 0
    for i in range(block.IMAX):
        for j in range(block.JMAX):
            for k in range(block.KMAX):
                points[0,indx] = block.X[i,j,k]
                points[1,indx] = block.Y[i,j,k]
                points[2,indx] = block.Z[i,j,k]
                indx+=1
    points_rotated = np.matmul(rotation_matrix,points)
    indx=0
    for i in range(block.IMAX):
        for j in range(block.JMAX):
            for k in range(block.KMAX):
                X[i,j,k] = points_rotated[0,indx]
                Y[i,j,k] = points_rotated[1,indx]
                Z[i,j,k] = points_rotated[2,indx]
                indx+=1
                
    return Block(X,Y,Z)
def get_outer_bounds(blocks:List[Block]):
    """Get outer bounds for a set of blocks

    Args:
        blocks (List[Block]): Blocks defining your shape

    Returns:
        (Tuple) containing: 

            **xbounds** (Tuple[float,float]): xmin,xmax
            **ybounds** (Tuple[float,float]): ymin,ymax
            **zbounds** (Tuple[float,float]): zmin,zmax
    """
    xbounds = [blocks[0].X.min(),blocks[0].X.max()]
    ybounds = [blocks[0].Y.min(),blocks[0].Y.max()]
    zbounds = [blocks[0].Z.min(),blocks[0].Z.max()]
    
    for i in range(1,len(blocks)):
        xmin = blocks[i].X.min()
        xmax = blocks[i].X.max()

        ymin = blocks[i].Y.min()
        ymax = blocks[i].Y.max()

        zmin = blocks[i].Z.min()
        zmax = blocks[i].Z.max()

        if xmin<xbounds[0]:
            xbounds[0] = xmin
        elif xmax>xbounds[1]:
            xbounds[1] = xmax
        
        if ymin<ybounds[0]:
            ybounds[0] = ymin
        elif ymax>ybounds[1]:
            ybounds[1] = ymax

        if zmin<zbounds[0]:
            zbounds[0] = zmin
        elif zmax>zbounds[1]:
            zbounds[1] = zmax
    
    return tuple(xbounds),tuple(ybounds),tuple(zbounds)

def block_connection_matrix(blocks:List[Block],outer_faces:List[Dict[str,int]]=[],tol:float=1E-8):
    """Creates a matrix representing how block edges are connected to each other 

    Args:
        blocks (List[Block]): List of blocks that describe the Plot3D mesh
        outer_faces (List[Dict[str,int]], optional): List of outer faces remaining from connectivity. Useful if you are interested in finding faces that are exterior to the block. Also useful if you combine outerfaces with match faces, this will help identify connections by looking at split faces. Defaults to [].
        tol (float, optional): Matching tolerance to look for when comparing face centroids.
    
    Returns:
        (Tuple): containing

            *connectivity* (np.ndarray): integer matrix defining how the blocks are connected to each other
            *connectivity_i* (np.ndarray): integer matrix defining connectivity of all blocks where IMAX=IMIN
            *connectivity_j* (np.ndarray): integer matrix defining connectivity of all blocks where JMAX=JMIN
            *connectivity_k* (np.ndarray): integer matrix defining connectivity of all blocks where KMAX=KMIN
            
    """
    # Reduce the size of the blocks by the GCD 
    gcd_array = list()    
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)

    # Face to List 
    outer_faces_all = list()
    for o in outer_faces:
        face = create_face_from_diagonals(blocks[o['block_index']], int(o['IMIN']/gcd_to_use), int(o['JMIN']/gcd_to_use), 
            int(o['KMIN']/gcd_to_use), int(o['IMAX']/gcd_to_use), int(o['JMAX']/gcd_to_use), int(o['KMAX']/gcd_to_use))
        face.set_block_index(o['block_index'])
        if "id" in o:
            face.id = o['id']
        outer_faces_all.append(face)

    outer_faces = outer_faces_all

    n = len(blocks)
    connectivity = np.eye(n,dtype=np.int8)
    combos = list(combinations(range(n),2))    
    for indx in (pbar:=trange(len(combos))):
        i,j = combos[indx]
        pbar.set_description(f"Building block to block connectivity matrix: checking {i}")
        b1 = blocks[i]

        if len(outer_faces)==0:                     # Get the outerfaces to search
            b1_outer_faces,_ = get_outer_faces(b1)
        else:
            b1_outer_faces = [o for o in outer_faces if o.BlockIndex == i] # type: ignore
        
        if i != j and connectivity[i,j]!=-1:
            b2 = blocks[j]

            if len(outer_faces)==0:                 # Get the outerfaces to search
                b2_outer_faces,_ = get_outer_faces(b2)
            else:
                b2_outer_faces = [o for o in outer_faces if o.BlockIndex == j]                 # type: ignore

            # Check to see if any of the outer faces of the blocks match   
            connection_found=False             
            for f1 in b1_outer_faces:
                for f2 in b2_outer_faces:
                    if (f1.is_connected(f2,tol)):   # type: ignore # Check if face centroid is the same
                        connectivity[i,j] = 1       # Default block to block connection matrix 
                        connectivity[j,i] = 1
                        connection_found=True
                        # c = np.sum(connectivity[i,:]==1)
                        # print(f"block {i} connections {c}")
                        break
                if connection_found:
                    break
            if not connection_found:
                connectivity[i,j] = -1
                connectivity[j,i] = -1      
    return connectivity
def plot_blocks(blocks):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    markers = ['o', 's']  # circle, square

    for i,b in enumerate(blocks):
        x_flat = b.X.ravel()
        y_flat = b.Y.ravel()
        z_flat = b.Z.ravel()
        ax.scatter(
                x_flat, y_flat, z_flat,
                s=1, alpha=0.4, # type: ignore
                marker=markers[i % len(markers)],
                label=f'Block {i}'
            )    
        ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z') # type: ignore
    ax.set_title('3D Block Grid')
    plt.show()

def transform_block_to_match_face(block, face_name_src, face_name_dst):
    """
    Given a block and two face names, transform the block so that face_name_src
    aligns with face_name_dst in orientation.
    
    Returns a new transformed Block instance.
    """
    X, Y, Z = block.X.copy(), block.Y.copy(), block.Z.copy()

    # Face alignment logic
    flip_map = {
        ('imin', 'imax'): 0,  # flip i
        ('imax', 'imin'): 0,
        ('jmin', 'jmax'): 1,
        ('jmax', 'jmin'): 1,
        ('kmin', 'kmax'): 2,
        ('kmax', 'kmin'): 2,
    }

    flip_axis = flip_map.get((face_name_src, face_name_dst))
    if flip_axis is not None:
        X = np.flip(X, axis=flip_axis)
        Y = np.flip(Y, axis=flip_axis)
        Z = np.flip(Z, axis=flip_axis)

    # Note: for true rotation (e.g., rotated faces), we’d need to handle transposes or axis permutations.
    # You can expand this logic if needed.

    return Block(X, Y, Z)

def combine_8_blocks(blocks: List[Block], tol: float = 1e-8, max_tries: int = 4) -> Tuple[List[Block], List[int]]:
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
                face1, face2 = find_matching_faces(blk_a, blk_b, tol=tol)
                if face1 is not None:
                    try:
                        merged = combine_blocks(blk_a, blk_b, tol=tol)
                        found = True
                        break
                    except Exception as e:
                        print(f"⚠️ Failed to merge blocks {i} and {j}: {e}")

            if found:
                new_merged.append(merged)
                skip.update([i, j])
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


def combine_blocks(block1, block2, tol=1e-8):
    """
    Combine block1 and block2 by matching and aligning one of their faces.
    Supports all 6 face pairs, including same-direction ones (e.g., jmax-jmax).
    Returns a new merged Block with geometry preserved.
    """
    face1, face2 = find_matching_faces(block1, block2, tol=tol)
    if face1 is None:
        print("No matching faces found between the blocks.")
        return block1

    print(f"Blocks connected: block1.{face1} matches block2.{face2}")

    axis_map = {
        ('imax', 'imin'): 0, ('imin', 'imax'): 0, ('imin', 'imin'): 0, ('imax', 'imax'): 0,
        ('jmax', 'jmin'): 1, ('jmin', 'jmax'): 1, ('jmin', 'jmin'): 1, ('jmax', 'jmax'): 1,
        ('kmax', 'kmin'): 2, ('kmin', 'kmax'): 2, ('kmin', 'kmin'): 2, ('kmax', 'kmax'): 2,
    }

    key = (face1, face2)
    if key not in axis_map:
        raise NotImplementedError(f"Merge not supported for face pair: {key}")
    axis = axis_map[key]

    # Transform block2 to align with block1
    aligned_block2 = transform_block_to_match_face(block2, face2, face1)

    # Slice off overlapping face to avoid duplication
    slicer = [slice(None), slice(None), slice(None)]
    if face2 in ['imin', 'jmin', 'kmin']:
        slicer[axis] = slice(1, None)  # Keep tail of block2
    else:
        slicer[axis] = slice(0, -1)  # Keep head of block2

    X2s = aligned_block2.X[tuple(slicer)]
    Y2s = aligned_block2.Y[tuple(slicer)]
    Z2s = aligned_block2.Z[tuple(slicer)]

    X = np.concatenate([block1.X, X2s], axis=axis)
    Y = np.concatenate([block1.Y, Y2s], axis=axis)
    Z = np.concatenate([block1.Z, Z2s], axis=axis)

    return Block(X, Y, Z)

def checkCollinearity(v1:npt.NDArray, v2:npt.NDArray):
    # Calculate their cross product
    cross_P = np.cross(v1,v2) 
 
    # Check if their cross product
    # is a NULL Vector or not
    if (cross_P[0] == 0 and
        cross_P[1] == 0 and
        cross_P[2] == 0):
        return True
    else:
        return False
def calculate_outward_normals(block:Block):
    # Calculate Normals
    X = block.X
    Y = block.Y
    Z = block.Z
    imax = block.IMAX
    jmax = block.JMAX
    kmax = block.KMAX 
    # IMAX - Normal should be out of the page        
    # Normals I direction: IMIN https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
    x = [X[0,0,0],X[0,jmax,0],X[0,0,kmax]] 
    y = [Y[0,0,0],Y[0,jmax,0],Y[0,0,kmax]]
    z = [Z[0,0,0],Z[0,jmax,0],Z[0,0,kmax]]
    u = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]]) 
    v = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_imin = np.cross(v1,v2) # type: ignore
    
    # Normals I direction: IMAX
    x = [X[imax,0,0],X[imax,jmax,0],X[imax,0,kmax]] 
    y = [Y[imax,0,0],Y[imax,jmax,0],Y[imax,0,kmax]]
    z = [Z[imax,0,0],Z[imax,jmax,0],Z[imax,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_imax = np.cross(v1,v2)

    # Normals J direction: JMIN
    x = [X[0,0,0],X[imax,0,0],X[0,0,kmax]] 
    y = [Y[0,0,0],Y[imax,0,0],Y[0,0,kmax]]
    z = [Z[0,0,0],Z[imax,0,0],Z[0,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_jmin = np.cross(v1,v2)

    # Normals J direction: JMAX
    x = [X[0,jmax,0],X[imax,jmax,0],X[0,jmax,kmax]] 
    y = [Y[0,jmax,0],Y[imax,jmax,0],Y[0,jmax,kmax]]
    z = [Z[0,jmax,0],Z[imax,jmax,0],Z[0,jmax,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_jmax = np.cross(v1,v2)

    # Normals K direction: KMIN
    x = [X[imax,0,0],X[0,jmax,0],X[0,0,0]] 
    y = [Y[imax,0,0],Y[0,jmax,0],Y[0,0,0]]
    z = [Z[imax,0,0],Z[0,jmax,0],Z[0,0,0]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_kmin = np.cross(v1,v2)

    # Normals K direction: KMAX
    x = [X[imax,0,kmax],X[0,jmax,kmax],X[0,0,kmax]] 
    y = [Y[imax,0,kmax],Y[0,jmax,kmax],Y[0,0,kmax]]
    z = [Z[imax,0,kmax],Z[0,jmax,kmax],Z[0,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_kmax = np.cross(v1,v2)

    return n_imin,n_jmin,n_kmin,n_imax,n_jmax,n_kmax
def split_blocks(blocks:List[Block],gcd:int=4):
    """Split blocks but also keep greatest common divisor

    Args:
        blocks (List[]): _description_
        gcd (int, optional): _description_. Defaults to 4.
    """
    pass

def common_neighbor(G: nx.Graph, a: int, b: int, exclude: Set[int]) -> int:
    """
    Return a node that is connected to both `a` and `b` and not in `exclude`.
    """
    for n in G.neighbors(a):
        if n in exclude:
            continue
        if G.has_edge(n, b):
            return n
    return None

def build_connectivity_graph(connectivities: List[List[Dict]]) -> nx.Graph:
    """
    Build an undirected graph from a list of face-to-face block connectivities.
    Each edge connects two block indices.
    """
    G = nx.Graph()
    for pair in connectivities:
        block1 = pair['block1']['block_index'] # type: ignore
        block2 = pair['block2']['block_index'] # type: ignore
        G.add_edge(block1, block2)
    return G

def combine_2x2x2_cubes(
    blocks: List[Block],
    connectivities: List[List[Dict]],
    tol: float = 1e-8
) -> List[Tuple[Block, Set[int]]]:
    """
    Find and combine all non-overlapping 2x2x2 cube groups of blocks
    using face connectivity data and return the merged components.

    Parameters
    ----------
    blocks : List[Block]
        All input block objects.
    connectivities : List of face match metadata pairs
        Face connectivity between blocks.
    tol : float
        Face match tolerance.

    Returns
    -------
    List[Tuple[Block, Set[int]]]
        A list of merged Block objects and their source block indices.
    """
    from itertools import combinations

    used = set()
    merged_groups = []
    G = build_connectivity_graph(connectivities)

    remaining_indices = list(range(len(blocks)))

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

            neighbors = list(G.neighbors(seed_index))
            found_group = None

            for a, b, c in combinations(neighbors, 3):
                if any(G.has_edge(u, v) for u, v in combinations([a, b, c], 2)):
                    continue

                base_corner = {seed_index, a, b, c}

                ab = common_neighbor(G, a, b, exclude=base_corner)
                if ab is None: continue
                bc = common_neighbor(G, b, c, exclude=base_corner | {ab})
                if bc is None: continue
                ac = common_neighbor(G, a, c, exclude=base_corner | {ab, bc})
                if ac is None: continue

                abc = common_neighbor(G, ab, bc, exclude=base_corner | {ac})
                if abc is None or not G.has_edge(abc, ac):
                    continue

                candidate_group = base_corner | {ab, bc, ac, abc}
                if candidate_group & used or candidate_group & new_used:
                    continue

                found_group = candidate_group
                break

            if not found_group:
                i += 1
                continue

            group_block_list = [blocks[k] for k in sorted(found_group)]
            index_mapping = {i: orig_idx for i, orig_idx in enumerate(sorted(found_group))}

            try:
                partial_merges, local_indices = combine_8_blocks(group_block_list, tol=tol)
                for merged_block in partial_merges:
                    merged_group = {index_mapping[i] for i in local_indices}
                    merged_groups.append((merged_block, merged_group))
                    new_used.update(merged_group)
                merged_this_round = True
            except Exception as e:
                print(f"⚠️ Skipping group {found_group} due to error: {e}")
                i += 1
                continue

            # Update remaining_indices and restart inner loop
            remaining_indices = [idx for idx in remaining_indices if idx not in new_used]
            i = 0

        # Apply used indices after full round
        used.update(new_used)

        if not merged_this_round or len(remaining_indices) == before_len:
            print("✅ No further merges possible. Appending unmerged blocks.")
            for idx in remaining_indices:
                merged_groups.append((blocks[idx], {idx}))
            break

    return merged_groups

def combine_nxnxn_cubes(
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
                partial_merges, local_indices = combine_8_blocks(group_block_list, tol=tol)
                for merged_block in partial_merges:
                    merged_group = {index_mapping[i] for i in local_indices}
                    merged_groups.append((merged_block, merged_group))
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