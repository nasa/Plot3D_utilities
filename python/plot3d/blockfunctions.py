from copy import deepcopy
from itertools import combinations, product
import math
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from tqdm import trange
from .facefunctions import create_face_from_diagonals, get_outer_faces
from .block import Block, combine_8_blocks,find_matching_faces
from .face import Face
import tqdm
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  


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

def reduce_blocks(blocks:List[Block],factor:int):
    """reduce the blocks by a factor of (factor)

    Args:
        blocks (List[Block]): list of blocks to reduce in size
        factor (int, optional): Number of indicies to skip . Defaults to 2.

    Returns:
        [type]: [description]
    """
    for i in range(len(blocks)):
        blocks[i].X = blocks[i].X[::factor,::factor,::factor]
        blocks[i].Y = blocks[i].Y[::factor,::factor,::factor]
        blocks[i].Z = blocks[i].Z[::factor,::factor,::factor]
        blocks[i].IMAX,blocks[i].JMAX,blocks[i].KMAX = blocks[i].X.shape
    return blocks

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

def combine_spatial_group_from_connectivity(blocks:List[Block],
    seed_index: int,
    connectivities: List[List[Dict]],
    already_used: Optional[Set[int]] = None
) -> Tuple[Optional[Block], Set[int]]:
    """
    Combine a group of `count` connected blocks starting from `seed_index` using connectivity data.

    Args:
        seed_index: Index of the starting block.
        blocks_by_index: Dict from block_index to Block object.
        connectivities: List of connectivity pairs (each pair = two dicts with block_index and face extents).
        count: Number of blocks to group (4, 8, or 26).
        already_used: Optional set of block indices already grouped.

    Returns:
        A merged Block and the set of block indices used.
    """
    if already_used is None:
        already_used = set()

    if seed_index in already_used:
        return None, set() # type: ignore

    # Build graph from connectivity info
    G = build_connectivity_graph(connectivities)
    neighbors = list(G.neighbors(seed_index))
    full_group = []
    for i in range(len(neighbors)):
        for j in range(i+1, len(neighbors)):
            for k in range(j+1, len(neighbors)):
                n1, n2, n3 = neighbors[i], neighbors[j], neighbors[k]

                # Ensure they are not directly connected to each other (orthogonal assumption)
                if G.has_edge(n1, n2) or G.has_edge(n1, n3) or G.has_edge(n2, n3):
                    continue

                base_corner = {seed_index, n1, n2, n3}

                # Try to find the 4 diagonal opposites
                ab = common_neighbor(G, n1, n2, exclude=base_corner)
                if ab is None:
                    continue
                bc = common_neighbor(G, n2, n3, exclude=base_corner.union({ab}))
                if bc is None:
                    continue
                ac = common_neighbor(G, n1, n3, exclude=base_corner.union({ab, bc}))
                if ac is None:
                    continue

                abc = common_neighbor(G, ab, bc, exclude=base_corner.union({ac}))
                if abc is None or not G.has_edge(abc, ac):
                    continue

                full_group = base_corner.union({ab, bc, ac, abc})


    if len(full_group) == 8:
        visited_blocks = [blocks[b] for b in full_group]
        merged = combine_8_blocks(visited_blocks)
        visited_blocks.append(merged)
        plot_blocks(visited_blocks)
    print('check')