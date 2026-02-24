"""Utilities for splitting structured multi-block grids into smaller blocks.

This module provides functions and enumerations used to divide Plot3D structured
blocks along a chosen index direction (i, j, or k) while preserving the
greatest common divisor (GCD) of the parent block's cell counts.  Maintaining
the GCD ensures that all child blocks remain compatible with the same multi-grid
hierarchy used by the solver.

Credits: Dave Rigby and Tim Beach.
"""

from .face import Face
from .block import Block
from typing import List
from enum import Enum
from math import gcd, sqrt 
import numpy as np 

class Direction(Enum):
    """Enumeration of the three structured-grid index directions.

    Attributes:
        i: First index direction (IMAX axis).
        j: Second index direction (JMAX axis).
        k: Third index direction (KMAX axis).
    """

    i = 0
    j = 1
    k = 2

def max_aspect_ratio(X:np.ndarray,Y:np.ndarray,Z:np.ndarray,ix:int,jx:int,kx:int):
    """Compute the maximum cell aspect ratio sampled at the block corners.

    The aspect ratio at each corner is estimated by comparing the edge lengths
    along the i, j, and k directions to their neighbours one cell inward.  The
    function samples all eight block corners and returns the largest ratio found.

    Args:
        X (np.ndarray): 3-D array of x-coordinates indexed as ``X[i, j, k]``.
        Y (np.ndarray): 3-D array of y-coordinates indexed as ``Y[i, j, k]``.
        Z (np.ndarray): 3-D array of z-coordinates indexed as ``Z[i, j, k]``.
        ix (int): Size of the i-dimension (unused; inferred from ``X.shape``).
        jx (int): Size of the j-dimension (unused; inferred from ``X.shape``).
        kx (int): Size of the k-dimension (unused; inferred from ``X.shape``).

    Returns:
        float: Maximum aspect ratio across all sampled corner cells.
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
        ds.append(sqrt((X[i2[n],j1[n],k1[n]]-X[i1[n],j1[n],k1[n]])**2 +
                    (Y[i2[n],j1[n],k1[n]]-Y[i1[n],j1[n],k1[n]])**2 +
                    (Z[i2[n],j1[n],k1[n]]-Z[i1[n],j1[n],k1[n]])**2)
                )
        ds.append(sqrt((X[i1[n],j2[n],k1[n]]-X[i1[n],j1[n],k1[n]])**2 +
                    (Y[i1[n],j2[n],k1[n]]-Y[i1[n],j1[n],k1[n]])**2 +
                    (Z[i1[n],j2[n],k1[n]]-Z[i1[n],j1[n],k1[n]])**2)
                )
        ds.append(sqrt((X[i1[n],j1[n],k2[n]]-X[i1[n],j1[n],k1[n]])**2 +
                    (Y[i1[n],j1[n],k2[n]]-Y[i1[n],j1[n],k1[n]])**2 +
                    (Z[i1[n],j1[n],k2[n]]-Z[i1[n],j1[n],k1[n]])**2)
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
    """Find a step size along one index direction that preserves the block GCD.

    Starting from an initial guess of ``round(ncells_per_block / denominator)``,
    the function walks forward or backward by one until it finds a step size
    where both the step itself and the final remainder block are divisible by
    ``greatest_common_divisor``.  This guarantees every child block shares the
    same multi-grid hierarchy as the parent.

    Args:
        total_cells (int): Total number of cells in the parent block
            (``IMAX * JMAX * KMAX``).
        greatest_common_divisor (int): GCD of ``(IMAX-1)``, ``(JMAX-1)``, and
            ``(KMAX-1)`` for the parent block.
        ncells_per_block (int): Target number of cells per child block.
        denominator (float): Product of the two index dimensions *not* being
            split (e.g. ``JMAX * KMAX`` when splitting along i).
        direction (str, optional): Search direction; ``'forward'`` increments
            the step size by +1 and ``'backward'`` decrements by -1.
            Defaults to ``'forward'``.

    Returns:
        int: A valid step size along the split direction, or ``-1`` if no
        valid step size is found within the search bounds.
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
    """Divide a list of structured blocks so that each child has roughly the target cell count.

    The split preserves the greatest common divisor (GCD) of the parent block's
    cell counts so that all child blocks remain compatible with the same
    multi-grid hierarchy.  For example, a GCD of 4 supports three multigrid
    levels: grid/4 (coarse), grid/2 (medium), and grid/1 (finest).

    ``ncells_per_block`` is a target, not a guarantee.  The actual count is
    adjusted to the nearest value that satisfies the GCD constraint.  Blocks
    whose total cell count is already at or below ``ncells_per_block`` are
    passed through unchanged.

    Note:
        Wisdom from Dave Rigby: additional controls such as integrating across
        the full span (radial equilibrium) or keeping the entire boundary layer
        inside a single block may require overriding the automatic direction
        selection via the ``direction`` argument.

    Args:
        blocks (List[Block]): Input list of structured Plot3D blocks to split.
        ncells_per_block (int): Approximate target number of cells per child
            block after splitting.
        direction (Direction, optional): Index direction along which to split
            every block.  When ``None`` (the default), the direction with the
            largest index size (IMAX, JMAX, or KMAX) is chosen automatically
            for each block.

    Returns:
        List[Block]: New list of blocks after splitting.  Blocks that were
        already smaller than ``ncells_per_block`` are included unmodified.

    Raises:
        ValueError: If no valid GCD-preserving step size can be found for a
            block in either the forward or backward search direction.
    """
    
    direction_to_use = direction # store the user input variable 
    
    new_blocks = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        total_cells = block.IMAX*block.JMAX*block.KMAX

        if direction==None: 
            indx = np.argmax(np.array([block.IMAX,block.JMAX,block.KMAX]))
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
                    raise ValueError('no valid step size found, do you have multi-block? gcd > 1')
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
                    raise ValueError('no valid step size found, do you have multi-block? gcd > 1')
                jprev = 0
                for j in range(step_size,block.JMAX,step_size):
                    if (j+1) > block.JMAX:
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
                    raise ValueError('no valid step size found, do you have multi-block? gcd > 1')
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

        