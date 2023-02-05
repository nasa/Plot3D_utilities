from copy import deepcopy
from itertools import combinations
import math
import numpy as np
from typing import Dict, List
from tqdm import trange
from .facefunctions import create_face_from_diagonals
from .block import Block
from .face import Face

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

def block_connection_matrix(blocks:List[Block],outer_faces:List[Dict[str,int]]=[]):
    """Creates a matrix representing how block edges are connected to each other 

    Args:
        blocks (List[Block]): _description_
        outer_faces (List[Dict[str,int]], optional): List of outer faces remaining from connectivity. Useful if you are interested in finding faces that are exterior to the block. Also useful if you combine outerfaces with match faces, this will help identify connections by looking at split faces. Defaults to [].

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
            b1_outer_faces = [o for o in outer_faces if o.BlockIndex == i]
        
        if i != j and connectivity[i,j]!=-1:
            b2 = blocks[j]

            if len(outer_faces)==0:                 # Get the outerfaces to search
                b2_outer_faces,_ = get_outer_faces(b2)
            else:
                b2_outer_faces = [o for o in outer_faces if o.BlockIndex == j]                

            # Check to see if any of the outer faces of the blocks match   
            connection_found=False             
            for f1 in b1_outer_faces:
                for f2 in b2_outer_faces:
                    if (f1.is_connected(f2)): # Check if face centroid is the same
                        connectivity[i,j] = 1       # Default block to block connection matrix 
                        connectivity[j,i] = 1
                        connection_found=True
                        break
                if connection_found:
                    break
            if not connection_found:
                connectivity[i,j] = -1
                connectivity[j,i] = -1
        # c = np.sum(connectivity[i,:]==1)
        # print(f"block {i} connections {c}")
    return connectivity
