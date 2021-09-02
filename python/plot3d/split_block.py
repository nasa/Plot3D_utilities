'''
Split blocks is combination of help from Dave Rigby and Tim Beach 
'''

from face import Face
import face
from .block import Block
from typing import List
from enum import Enum
from math import gcd, sqrt 
import numpy as np 

class Direction(Enum):
    i = 0
    j = 1
    k = 2

def max_aspect_ratio(X:np.ndarray,Y:np.ndarray,Z:np.ndarray,ix:int,jx:int,kx:int):
    """Finds the maximum cell aspect ratio

    Args:
        X (np.ndarray): 3 dimensional Array containing X values. X[i,j,k]
        Y (np.ndarray): 3 dimensional Array containing X values. X[i,j,k]
        Z (np.ndarray): 3 dimensional Array containing X values. X[i,j,k]

    Returns:
        float: Maximum value of aspect ratio 
    """
    [ix,jx,kx]  = X.shape # ix, jx, kx are the max values like IMAX, JMAX, KMAX
    ix-=1; jx-=1; kx-=1 # Python is 0 index so need to subtract 1 from ix, iy, kx to reference the last value 

    # Each column is a corner 
    # [corner 1, corner 2, corner 3, corner 4]
    i1 = [0, ix, 0, 0] 
    j1 = [0, 0, jx, 0]
    k1 = [0, 0, 0, kx]

    # Appending more corners []
    i1 = i1 + [0,ix,ix,ix]
    j1 = j1 + [jx,0,jx,jx]
    k1 = k1 + [kx,kx,0,kx]

    # i2, j2, k2 are the first nodes diagonally from the block corners 
    i2 = [1, ix-1, 1, 1]
    j2 = [1, 1, jx-1, 1]
    k2 = [1, 1, 1, kx-1]

    i2 = i2 + [1, ix-1, ix-1, ix-1]
    j2 = j2 + [jx-1, 1, jx-1, jx-1]
    k2 = k2 + [kx-1, kx-1, 1, kx-1]
    
    ds = list()
    for n in range(len(i2)):
        ds.append(sqrt((X(i2[n],j1[n],k1[n])-X(i1[n],j1[n],k1[n]))**2 +
                    (Y(i2[n],j1[n],k1[n])-Y(i1[n],j1[n],k1[n]))**2 +  
                    (Z(i2[n],j1[n],k1[n])-Z(i1[n],j1[n],k1[n]))**2)
                )
        ds.append(sqrt((X(i1[n],j2[n],k1[n])-X(i1[n],j1[n],k1[n]))**2 +  
                    (Y(i1[n],j2[n],k1[n])-Y(i1[n],j1[n],k1[n]))**2 +  
                    (Z(i1[n],j2[n],k1[n])-Z(i1[n],j1[n],k1[n]))**2)
                )
        ds.append(sqrt((X(i1[n],j1[n],k2[n])-X(i1[n],j1[n],k1[n]))**2 + 
                    (Y(i1[n],j1[n],k2[n])-Y(i1[n],j1[n],k1[n]))**2 + 
                    (Z(i1[n],j1[n],k2[n])-Z(i1[n],j1[n],k1[n]))**2)
                )
    aspect = [0,0,0]
    if ds[0]>0:
        aspect[0] = max(ds[1],ds[2])/ds[0]
    elif(ds[1]>0):
        aspect[1] = max(ds[2],ds[0])/ds[1]
    elif(ds[2]>0):
        aspect[2] = max(ds[0],ds[1])/ds[2]

    return max(aspect) 

"""Divide the matched block according to a direction 

    Args:
        block (Block): [description]
        n_cells (int): [description]
        direction (Direction): [description]
        face_to_match (Face): This is the new face on a different block that we will match.
    
    Returns:
        
    """
def __divide_matched_block(block:Block,istep:int, jstep:int, kstep:int,face_to_split:Face,face_to_match:Face):
    """Divides the matched block with consideration to it's connectivity. face_to_split is the face on the block. face_to_match is the new face from a different block. 
    
    This function is only called if:   
        Connectivity block 1 and 2 starting (I,J,K) <= block 1 split (I, J, or K <= Connectivity block 1 and 2 ending (I,J,K)

    Args:
        block (Block): [description]
        istep (int): [description]
        jstep (int): [description]
        kstep (int): [description]
        face_to_split (Face): [description]
        face_to_match (Face): [description]
    """
    # Check if to advance in i forward or backward by looking at the verticies 
    
    face_to_match.x

def split_blocks(blocks:List[Block], ncells_per_block:int, face_matches:dict, outer_faces_formatted:dict,direction:Direction,split_matching_faces:bool=False):
    """Split an array of blocks based on number of cells per block. This takes into account the outer faces and matching faces in the block
        So say you had a block with 2M cells and you wanted to split into blocks containing around 400K cells. This function will allow you to do that. 

        Wisdom from Dave Rigby:
            For example, for radial equilibrium we must integrate across the span.  Some codes (GlennHT used to) would want a single block across the entire span.  In that case you would want some additional control.
 
            Another example might be if you would like a block to include the entire boundary layer.  In that case you might introduce an aspect ratio control.
            
    Args:
        blocks (List[Block]): List of blocks
        ncells_per_block (int): number of cells desired per block 

        face_matches (dict): This dictionary describes how the face of a block matches the face of another block or within itself. The format is {'block1': {'block_index','IMIN','IMAX','JMIN','JMAX','KMIN','KMAX'}}

        outer_faces_formatted (dict): [description]
        direction (Direction): direction to split the blocks in. Direction.(i,j,k)
        split_matching_faces (bool): choose whether or not to split where the faces are matching. Defaults to False.

    Returns:

    """  

    new_blocks = list()
    new_face_matches = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        total_cells = block.IMAX*block.JMAX*block.KMAX
        # Use greatest common divsor to maintain multi-grid so say the entire block is divisible by 4 then we want to maintain than for all the splits! 
        greatest_common_divisor = gcd(block.IMAX, block.JMAX, block.KMAX) # Gets the maximum number of partitions that we can make for this given block 
        if direction == Direction.i: 
            iprev = 0
            # In order to get close to the number of cells per block, we need to control how many steps of the greatest_common_divisor to advance so for example if you have a multigrid mesh that has gcd of 16 (fine) => 8 (coarse) => 4 (coarser) => 2 (coarsest) and you want 400K cells per block then JMAX*KMAX*gcd*some_factor has to be close to 400K cells
            max_divisions_in_direction = total_cells/(block.JMAX*block.KMAX*greatest_common_divisor)
            step_size = round(ncells_per_block/max_divisions_in_direction)
            cells_per_step = step_size * max_divisions_in_direction 
            for i in range(0,block.IMAX,step=step_size):
                X = block.X[iprev:i-1,:,:]      # New X, Y, Z splits 
                Y = block.Y[iprev:i-1,:,:]
                Z = block.Z[iprev:i-1,:,:]

                # Check face matches 
                for m in face_matches:
                    for block_name in ['block1', 'block2']:
                        if m[block_name]['block_index'] == block_indx: 
                            if m[block_name]['block_index']['IMIN'] < iprev: 
                                # if the split is larger than one of the matching blocks 
                                if split_matching_faces: # if we are splitting matching faces
                                    # Logic to divide up the connectivity file 
                                    
                                    iprev = i
                            else:
                                iprev = i
                
                # Check outer faces. Keep same surface ID but do the split
                        
                
                # There needs to be an adjustment when we overlap a matched face, we don't want to break the match or if we do we rematch it perfectly. So before appending the block, need to check face_matches and outer faces to see if we are splitting it right

                new_blocks.append(Block(X,Y,Z))
                
    # Check if any of the splits is in the face_matches with a different block or outer_faces of each block 
