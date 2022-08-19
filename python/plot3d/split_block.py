'''
Split blocks is combination of help from Dave Rigby and Tim Beach 
'''

from .face import Face
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


def __step_search(total_cells:int,greatest_common_divisor:int,ncells_per_block:int,denominator:float,direction:str='forward'):
    """Searches for the right step size that leads to block splits with the same greatest common divisor. Same greatest common denominator means you will always have the same multi-grid when you solve. 

    Args:
        total_cells (int): [description]
        greatest_common_divisor (int): [description]
        ncells_per_block (int): [description]
        denominator (float): If the direction is "i" then denominator is JMAX*KMAX
        direction (str, optional): 'forward' or 'backward' This chooses to increment step direction by +1 or -1. Defaults to 'forward'.

    Returns:
        [type]: [description]
    """
    step_size = round(ncells_per_block/denominator)     # initial Guess
    initial_guess = step_size
    Number_of_Cells_Remaining = total_cells % (step_size*denominator)
    IJK_MAX_Remainder = Number_of_Cells_Remaining/denominator
    increment = -1
    if direction.lower() == 'forward':
        increment=1 
    while ( (step_size %greatest_common_divisor != 0) or 
            ((IJK_MAX_Remainder-1)%greatest_common_divisor !=0) and 
            step_size>initial_guess/2 and 
            step_size<initial_guess*1.5):
        if (step_size %greatest_common_divisor == 0) and (IJK_MAX_Remainder-1)%greatest_common_divisor ==0:
            break
        step_size+=increment 
        Number_of_Cells_Remaining = total_cells % ((step_size)*denominator)
        IJK_MAX_Remainder = Number_of_Cells_Remaining/denominator
    
    if step_size %greatest_common_divisor != 0:             # This checks if each of the split blocks is divisible by gcd
        if (IJK_MAX_Remainder-1)%greatest_common_divisor !=0:   #  This checks the remaining/final block to see if it is divisible by gcd
            return -1
    return step_size

def split_blocks(blocks:List[Block], ncells_per_block:int,direction:Direction=None):
    """Split blocks is used to divide an array of blocks based on number of cells per block. This code maintains the greatest common denominator of the parent block. Number of cells per block is simply an estimate of how many you want. The actual number will change to meet the greatest common denominator (GCD). GCD of 4 means multigrid of 3 e.g. grid/4 (coarse), 2 (fine), and 1 (finest). If a direction is not specified then for each block the longest index either i,j, or k is used. 
    

        Wisdom from Dave Rigby:
            For example, for radial equilibrium we must integrate across the span.  Some codes (GlennHT used to) would want a single block across the entire span.  In that case you would want some additional control.
 
            Another example might be if you would like a block to include the entire boundary layer.  In that case you might introduce an aspect ratio control.
            
    Args:
        blocks (List[Block]): List of blocks
        ncells_per_block (int): number of cells desired per block 
        direction (Direction): direction to split the blocks in. Direction.(i,j,k). Defaults to None. None means it will pick the direction for you based on which is greater IMAX, JMAX, or KMAX

    Returns:
        Blocks (List[Block]): list of blocks split in the specified direction 
    """ 
    
    direction_to_use = direction # store the user input variable 
    
    new_blocks = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        total_cells = block.IMAX*block.JMAX*block.KMAX

        if direction==None: 
            indx = np.argmin(np.array([block.IMAX,block.JMAX,block.KMAX]))
            if indx == 0:
                direction_to_use=Direction.i
            elif indx == 1:
                direction_to_use=Direction.j
            elif indx == 2:
                direction_to_use=Direction.k

        if total_cells>ncells_per_block:
            # Use greatest common divsor to maintain multi-grid so say the entire block is divisible by 4 then we want to maintain than for all the splits! 
            greatest_common_divisor =gcd(block.IMAX-1, gcd(block.JMAX-1, block.KMAX-1)) # Gets the maximum number of partitions that we can make for this given block
            if direction_to_use == Direction.i: 
                
                # In order to get close to the number of cells per block, we need to control how many steps of the greatest_common_divisor to advance so for example if you have a multigrid mesh that has gcd of 16 (fine) => 8 (coarse) => 4 (coarser) => 2 (coarsest) and you want 400K cells per block then JMAX*KMAX*gcd*some_factor has to be close to 400K cells
                denominator = block.JMAX*block.KMAX
                step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='backward')
                if step_size==-1:
                    step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='forward')
                if step_size==-1:
                    assert('no valid step size found, do you have multi-block? gcd > 1')
                # step_size-1 is the IMAX of the sub_blocks e.g. 0 to 92 this shows IMAX=93, (93-1) % 4 = 0 (good)
                
                iprev = 0
                for i in range(step_size,block.IMAX,step_size):
                    if (i+1) > block.IMAX:
                        break

                    X = block.X[iprev:i+1,:,:]      # New X, Y, Z splits 
                    Y = block.Y[iprev:i+1,:,:]      # This indexes to iprev:i so if iprev=2 and i = 10 it will go from 2 to 9
                    Z = block.Z[iprev:i+1,:,:]
                    iprev=i                     # Blocks have to share the same face, Pick the previous face
                    new_blocks.append(Block(X,Y,Z))

                # Check for remainder
                if i+1 < block.IMAX:
                    # Add remainder to last block
                    X = block.X[i:,:,:] # New X, Y, Z splits 
                    Y = block.Y[i:,:,:]
                    Z = block.Z[i:,:,:]
                    new_blocks.append(Block(X,Y,Z))
                
            elif direction_to_use == Direction.j:                
                denominator = block.IMAX*block.KMAX
                step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='backward')
                if step_size==-1:
                    step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='forward')
                if step_size==-1:
                    assert('no valid step size found, do you have multi-block? gcd > 1')
                jprev = 0
                for j in range(step_size,block.JMAX,step_size):
                    if (j+1) > block.IMAX:
                        break
                    X = block.X[:,jprev:j,:]      # New X, Y, Z splits 
                    Y = block.Y[:,jprev:j,:]
                    Z = block.Z[:,jprev:j,:]
                    jprev=j
                    new_blocks.append(Block(X,Y,Z))

                # Check for remainder
                if j+1 < block.JMAX:
                    # Add remainder to last block
                    X = block.X[:,j:,:] # New X, Y, Z splits 
                    Y = block.Y[:,j:,:]
                    Z = block.Z[:,j:,:]
                    new_blocks.append(Block(X,Y,Z))
            else:
                denominator = block.IMAX*block.JMAX
                step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='backward')
                if step_size==-1:
                    step_size = __step_search(total_cells,greatest_common_divisor,ncells_per_block,denominator,direction='forward')
                if step_size==-1:
                    assert('no valid step size found, do you have multi-block? gcd > 1')
                kprev = 0
                for k in range(step_size,block.KMAX,step_size):
                    if (k+1) > block.KMAX:
                        break
                    X = block.X[:,:,kprev:k+1]     # New X, Y, Z splits 
                    Y = block.Y[:,:,kprev:k+1]
                    Z = block.Z[:,:,kprev:k+1]
                    kprev=k
                    new_blocks.append(Block(X,Y,Z))

                # Check for remainder
                if k+1 < block.KMAX:
                   # Add remainder to last block
                    X = block.X[:,:,k:] # New X, Y, Z splits 
                    Y = block.Y[:,:,k:]
                    Z = block.Z[:,:,k:]                    
                    new_blocks.append(Block(X,Y,Z)) # replace it 
        else:
            new_blocks.append(block)
    return new_blocks

        