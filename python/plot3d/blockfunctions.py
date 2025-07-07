from copy import deepcopy
from itertools import combinations, product
import math
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from tqdm import trange
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_matching_faces,faces_match
from .block import Block, reduce_blocks
from .write import write_plot3D
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
    gcd_array = list()    
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    blocks = reduce_blocks(deepcopy(blocks),4)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    markers = ['o', 's']  # alternate between circle and square

    for i, b in enumerate(blocks):
        color = f"C{i % 10}"
        X, Y, Z = b.X, b.Y, b.Z
        ax.scatter(X.ravel(), Y.ravel(), Z.ravel(), s=20, alpha=0.4,
                   marker=markers[i % len(markers)], label=f'Block {i}', color=color)

        # Draw lines along i-direction (axis 0)
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                ax.plot(X[:, j, k], Y[:, j, k], Z[:, j, k], color=color, linewidth=0.8, alpha=0.6)

        # Draw lines along j-direction (axis 1)
        for i_ in range(X.shape[0]):
            for k in range(X.shape[2]):
                ax.plot(X[i_, :, k], Y[i_, :, k], Z[i_, :, k], color=color, linewidth=0.8, alpha=0.6)

        # Draw lines along k-direction (axis 2)
        for i_ in range(X.shape[0]):
            for j_ in range(X.shape[1]):
                ax.plot(X[i_, j_, :], Y[i_, j_, :], Z[i_, j_, :], color=color, linewidth=0.8, alpha=0.6)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Block Grid with Connected Lines')
    ax.legend()
    plt.tight_layout()
    plt.show()
    

def transform_block_to_match_face(block, face_src: str, face_dst: str) -> Block:
    """
    Transform `block` so that `face_src` aligns with `face_dst`.

    Handles both axis permutation and flipping.
    """
    X, Y, Z = block.X.copy(), block.Y.copy(), block.Z.copy()

    def face_to_axis(face: str) -> int:
        if face.startswith('i'): return 0
        if face.startswith('j'): return 1
        if face.startswith('k'): return 2
        raise ValueError(f"Invalid face name: {face}")

    axis_src = face_to_axis(face_src)
    axis_dst = face_to_axis(face_dst)

    axes = [0, 1, 2]
    # Step 1: transpose if axis_src != axis_dst
    if axis_src != axis_dst:
        axes[axis_src], axes[axis_dst] = axes[axis_dst], axes[axis_src]
        X = np.transpose(X, axes)
        Y = np.transpose(Y, axes)
        Z = np.transpose(Z, axes)
        # After transpose, axis_dst is now at position axis_src
        flip_axis = axis_src
    else:
        flip_axis = axis_dst

    # Step 2: flip if min ↔ max
    flip_needed = (
        (face_src.endswith('min') and face_dst.endswith('max')) or
        (face_src.endswith('max') and face_dst.endswith('min'))
    )
    if flip_needed:
        X = np.flip(X, axis=flip_axis)
        Y = np.flip(Y, axis=flip_axis)
        Z = np.flip(Z, axis=flip_axis)

    return Block(X, Y, Z)


def combine_blocks(blocks: List[Block], tol: float = 1e-8, max_tries: int = 4) -> Tuple[List[Block], List[int]]:
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
                        merged = combine_2_blocks(blk_a, blk_b, tol=tol)
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


def standardize_block_orientation(block:Block):
    """Standardizes the orientation of a block so that its physical coordinates increase
    consistently along each of the indexing axes:

        - X increases along the i-axis
        - Y increases along the j-axis
        - Z increases along the k-axis

    This ensures consistent face orientation and alignment across multiple blocks,
    especially useful when merging, visualizing, or exporting grids. The function
    checks the dominant physical component (X, Y, or Z) along each axis, and flips
    the block along that axis if the component decreases.

    Parameters:
        block (Block): The input block to be standardized.

    Returns:
        Block: A new block instance with all three axes oriented consistently.
    """

    X, Y, Z = block.X.copy(), block.Y.copy(), block.Z.copy()
    i_center = X.shape[0] // 2
    j_center = X.shape[1] // 2
    k_center = X.shape[2] // 2

    # Check i-direction
    dx_i = X[-1, j_center, k_center] - X[0, j_center, k_center]
    if dx_i < 0:
        X = np.flip(X, axis=0)
        Y = np.flip(Y, axis=0)
        Z = np.flip(Z, axis=0)

    # Check j-direction
    dy_j = Y[i_center, -1, k_center] - Y[i_center, 0, k_center]
    if dy_j < 0:
        X = np.flip(X, axis=1)
        Y = np.flip(Y, axis=1)
        Z = np.flip(Z, axis=1)

    # Check k-direction
    dz_k = Z[i_center, j_center, -1] - Z[i_center, j_center, 0]
    if dz_k < 0:
        X = np.flip(X, axis=2)
        Y = np.flip(Y, axis=2)
        Z = np.flip(Z, axis=2)

    return Block(X, Y, Z)

def fix_physical_direction(block1: Block, axis: int) -> Block:
    """Ensures that the dominant physical direction (X, Y, or Z) along the given index axis (i, j, or k)
    increases from the start to the end of the block. If it decreases, the block is flipped along
    the specified axis to maintain consistent stacking and orientation.

    This is useful when merging blocks where one block may be geometrically flipped (e.g., Z increases
    as j decreases), causing discontinuities or 'jumps' at the interface.

    Parameters:
        block (Block): The input block whose orientation should be checked and potentially flipped.
        axis (int): The index axis (0 for i, 1 for j, 2 for k) used for stacking and physical direction check.

    Returns:
        Block: A new Block instance with corrected orientation along the specified axis if needed.
    """
    X, Y, Z = block1.X.copy(), block1.Y.copy(), block1.Z.copy()
    shape = X.shape
    imid, jmid, kmid = shape[0] // 2, shape[1] // 2, shape[2] // 2

    if axis == 0:
        line = np.array([X[:, jmid, kmid], Y[:, jmid, kmid], Z[:, jmid, kmid]])
    elif axis == 1:
        line = np.array([X[imid, :, kmid], Y[imid, :, kmid], Z[imid, :, kmid]])
    else:  # axis == 2
        line = np.array([X[imid, jmid, :], Y[imid, jmid, :], Z[imid, jmid, :]])

    # Compute average delta of each spatial component
    dz = line[2][-1] - line[2][0]
    dy = line[1][-1] - line[1][0]
    dx = line[0][-1] - line[0][0]

    # Find dominant spatial direction
    deltas = [dx, dy, dz]
    dominant_axis = np.argmax(np.abs(deltas))
    dominant_delta = deltas[dominant_axis]

    if dominant_delta < 0:
        # print(f"Flipping block1 along axis {axis} to align dominant direction (axis {dominant_axis})")
        X = np.flip(X, axis=axis)
        Y = np.flip(Y, axis=axis)
        Z = np.flip(Z, axis=axis)

    return Block(X, Y, Z)


def combine_2_blocks(block1, block2, tol=1e-8):
    """
    Combine block1 and block2 by matching and aligning one of their faces.
    Corrects both index and physical direction mismatches.
    """
    face1, face2, flip_flags = find_matching_faces(block1, block2, tol=tol)

    if face1 is None or flip_flags is None:
        print("No matching faces or face mismatch.")
        return block1

    # print(f"Blocks connected: block1.{face1} matches block2.{face2}")
    flip_ud, flip_lr = flip_flags
    X2, Y2, Z2 = block2.X.copy(), block2.Y.copy(), block2.Z.copy()

    # Determine axis of stacking
    axis_map = {
        ('imax', 'imin'): 0, ('imin', 'imax'): 0, ('imin', 'imin'): 0, ('imax', 'imax'): 0,
        ('jmax', 'jmin'): 1, ('jmin', 'jmax'): 1, ('jmin', 'jmin'): 1, ('jmax', 'jmax'): 1,
        ('kmax', 'kmin'): 2, ('kmin', 'kmax'): 2, ('kmin', 'kmin'): 2, ('kmax', 'kmax'): 2,
    }
    key = (face1, face2)
    if key not in axis_map:
        raise NotImplementedError(f"Merge not supported for face pair: {key}")
    axis = axis_map[key]

    # Flip block1 for correct stacking direction
    block1 = fix_physical_direction(block1, axis)
    block2 = fix_physical_direction(block2, axis)

    # Step 1: Flip block2 for face alignment
    if face2 in ['imin', 'imax']:
        if flip_ud:
            X2 = np.flip(X2, axis=1)
            Y2 = np.flip(Y2, axis=1)
            Z2 = np.flip(Z2, axis=1)
        if flip_lr:
            X2 = np.flip(X2, axis=2)
            Y2 = np.flip(Y2, axis=2)
            Z2 = np.flip(Z2, axis=2)
    elif face2 in ['jmin', 'jmax']:
        if flip_ud:
            X2 = np.flip(X2, axis=0)
            Y2 = np.flip(Y2, axis=0)
            Z2 = np.flip(Z2, axis=0)
        if flip_lr:
            X2 = np.flip(X2, axis=2)
            Y2 = np.flip(Y2, axis=2)
            Z2 = np.flip(Z2, axis=2)
    elif face2 in ['kmin', 'kmax']:
        if flip_ud:
            X2 = np.flip(X2, axis=0)
            Y2 = np.flip(Y2, axis=0)
            Z2 = np.flip(Z2, axis=0)
        if flip_lr:
            X2 = np.flip(X2, axis=1)
            Y2 = np.flip(Y2, axis=1)
            Z2 = np.flip(Z2, axis=1)

    # Step 2: Slice off overlapping face from block2
    if face2 == 'imin':
        X2s, Y2s, Z2s = X2[1:,:,:], Y2[1:,:,:], Z2[1:,:,:]
    elif face2 == 'imax':
        X2s, Y2s, Z2s = X2[:-1,:,:], Y2[:-1,:,:], Z2[:-1,:,:]
    elif face2 == 'jmin':
        X2s, Y2s, Z2s = X2[:,1:,:], Y2[:,1:,:], Z2[:,1:,:]
    elif face2 == 'jmax':
        X2s, Y2s, Z2s = X2[:,:-1,:], Y2[:,:-1,:], Z2[:,:-1,:]
    elif face2 == 'kmin':
        X2s, Y2s, Z2s = X2[:,:,1:], Y2[:,:,1:], Z2[:,:,1:]
    elif face2 == 'kmax':
        X2s, Y2s, Z2s = X2[:,:,:-1], Y2[:,:,:-1], Z2[:,:,:-1]
    else:
        raise ValueError(f"Unexpected face2: {face2}")

    # Step 3: Concatenate in correct order
    if face2 in ['imin', 'jmin', 'kmin']:
        X = np.concatenate([block1.X, X2s], axis=axis)
        Y = np.concatenate([block1.Y, Y2s], axis=axis)
        Z = np.concatenate([block1.Z, Z2s], axis=axis)
    else:
        X = np.concatenate([X2s, block1.X], axis=axis)
        Y = np.concatenate([Y2s, block1.Y], axis=axis)
        Z = np.concatenate([Z2s, block1.Z], axis=axis)

    # Step 4: Standardize and return
    merged = Block(X, Y, Z)
    write_plot3D('block1.xyz',[block1,block2,merged],False)
    return standardize_block_orientation(merged)




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
                partial_merges, local_indices = combine_blocks(group_block_list, tol=tol)
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
                partial_merges, local_indices = combine_blocks(group_block_list, tol=tol)
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