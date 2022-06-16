import numpy as np 
import math 
from tqdm import trange
from typing import List

class Block:
    """Plot3D Block definition
    """
    def __init__(self, X:np.ndarray,Y:np.ndarray,Z:np.ndarray):
        """Initializes the block using all the X,Y,Z coordinates of the block

        Args:
            X (np.ndarray): All the X coordinates (i,j,k)
            Y (np.ndarray): All the Y coordinates (i,j,k)
            Z (np.ndarray): All the Z coordinates (i,j,k)

        """
        self.IMAX, self.JMAX, self.KMAX = X.shape
        self.X = X
        self.Y = Y
        self.Z = Z
        # Centroid 
        self.cx = np.mean(X) 
        self.cy = np.mean(Y)
        self.cz = np.mean(Z)
    
    def cylindrical(self):
        """Converts the block to cylindrical coordinate system. The rotation axis is assumed to be "x" direction
        """
        self.r = np.sqrt(self.Z*self.Z + self.Y*self.Y)
        self.theta = np.arctan2(self.Y,self.Z)

    def cell_volumes(self):
        """Compute volume of all cells

        Returns:
            numpy.ndarray: volume of all cells

        Reference:
            Davies, D.E. and Salmond, D.J., "Calculation of the Volume of a General Hexahedron for Flow Predicitons", AIAA Journal, vol. 23, No. 6, pp. 954-956, June 1985.  It is (supposedly) exact for a hexahedral whose faces are bi-linear surfaces (i.e., the simplest surface that can be fit through the four nodes defining the face).
        """
        X = self.X
        Y = self.Y
        Z = self.Z
        a = [np.zeros(shape=(self.IMAX,self.JMAX,self.KMAX))  for _ in range(9)]
        # face csi=const 
        for k in range(1,self.KMAX):
            for j in range(1,self.JMAX):
                for i in range(self.IMAX):          # csi-const
                    dx1 = X[i,j,k-1] - X[i,j-1,k]
                    dy1 = Y[i,j,k-1] - Y[i,j-1,k]
                    dz1 = Z[i,j,k-1] - Z[i,j-1,k]

                    dx2 = X[i,j,k] - X[i,j-1,k-1]
                    dy2 = Y[i,j,k] - Y[i,j-1,k-1]
                    dz2 = Z[i,j,k] - Z[i,j-1,k-1]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[0][i,j,k] =  ax*0.5
                    a[1][i,j,k]  = ay*0.5
                    a[2][i,j,k]  = az*0.5

        
        for k in range(1,self.KMAX):
            for j in range(self.JMAX):              # face eta=const 
                for i in range(1,self.IMAX):
                    dx1 = X[i,j,k] - X[i-1,j,k-1]
                    dy1 = Y[i,j,k] - Y[i-1,j,k-1]
                    dz1 = Z[i,j,k] - Z[i-1,j,k-1]

                    dx2 = X[i,j,k-1] - X[i-1,j,k]
                    dy2 = Y[i,j,k-1] - Y[i-1,j,k]
                    dz2 = Z[i,j,k-1] - Z[i-1,j,k]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[3][i,j,k] = ax*0.5
                    a[4][i,j,k] = ay*0.5
                    a[5][i,j,k] = az*0.5

        # face zit=const 
        for k in range(self.KMAX):                  # zit=const
            for j in range(1,self.JMAX):            
                for i in range(1,self.IMAX):
                    dx1 = X[i,j,k] - X[i-1,j-1,k]
                    dy1 = Y[i,j,k] - Y[i-1,j-1,k]
                    dz1 = Z[i,j,k] - Z[i-1,j-1,k]

                    dx2 = X[i-1,j,k] - X[i,j-1,k]
                    dy2 = Y[i-1,j,k] - Y[i,j-1,k]
                    dz2 = Z[i-1,j,k] - Z[i,j-1,k]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[6][i,j,k] = ax*0.5
                    a[7][i,j,k] = ay*0.5
                    a[8][i,j,k] = az*0.5

        cf = np.zeros(shape=(6,3))
        v = np.zeros(shape=(self.IMAX,self.JMAX,self.KMAX))
        
        for k in trange(1,self.KMAX,desc='Calculating the volumes'):
            for j in range(1,self.JMAX):            
                for i in range(1,self.IMAX):
                    cf[0,0] = X[i-1,j-1,k-1] + X[i-1,j-1,k] + X[i-1,j,k-1] + X[i-1,j,k]
                    cf[0,1] = Y[i-1,j-1,k-1] + Y[i-1,j-1,k] + Y[i-1,j,k-1] + Y[i-1,j,k]
                    cf[0,2] = Z[i-1,j-1,k-1] + Z[i-1,j-1,k] + Z[i-1,j,k-1] + Z[i-1,j,k]
                    cf[1,0] = X[i,j-1,k-1] + X[i,j-1,k] + X[i,j,k-1] + X[i,j,k]
                    cf[1,1] = Y[i,j-1,k-1] + Y[i,j-1,k] + Y[i,j,k-1] + Y[i,j,k]
                    cf[1,2] = Z[i,j-1,k-1] + Z[i,j-1,k] + Z[i,j,k-1] + Z[i,j,k]
                    cf[2,0] = X[i-1,j-1,k-1] + X[i-1,j-1,k] + X[i,j-1,k-1] + X[i,j-1,k]
                    cf[2,1] = Y[i-1,j-1,k-1] + Y[i-1,j-1,k] + Y[i,j-1,k-1] + Y[i,j-1,k]
                    cf[2,2] = Z[i-1,j-1,k-1] + Z[i-1,j-1,k] + Z[i,j-1,k-1] + Z[i,j-1,k]
                    cf[3,0] = X[i-1,j,k-1] + X[i-1,j,k] + X[i,j,k-1] + X[i,j,k]
                    cf[3,1] = Y[i-1,j,k-1] + Y[i-1,j,k] + Y[i,j,k-1] + Y[i,j,k]
                    cf[3,2] = Z[i-1,j,k-1] + Z[i-1,j,k] + Z[i,j,k-1] + Z[i,j,k]
                    cf[4,0] = X[i-1,j-1,k-1] + X[i-1,j,k-1] + X[i,j-1,k-1] + X[i,j,k-1]
                    cf[4,1] = Y[i-1,j-1,k-1] + Y[i-1,j,k-1] + Y[i,j-1,k-1] + Y[i,j,k-1]
                    cf[4,2] = Z[i-1,j-1,k-1] + Z[i-1,j,k-1] + Z[i,j-1,k-1] + Z[i,j,k-1]
                    cf[5,0] = X[i-1,j-1,k] + X[i-1,j,k] + X[i,j-1,k] + X[i,j,k]
                    cf[5,1] = Y[i-1,j-1,k] + Y[i-1,j,k] + Y[i,j-1,k] + Y[i,j,k]
                    cf[5,2] = Z[i-1,j-1,k] + Z[i-1,j,k] + Z[i,j-1,k] + Z[i,j,k]

                    vol12=0
                    for n in range(0,2): # n = 0,1
                        for l in range(0,3): # l = 0,1,2
                            vol12 += math.pow(-1,n+1) *(               
                                        +cf[n,l]*a[l][i-1+n,j,k]
                                        +cf[2+n,l]*a[3+l][i,j-1+n,k]
                                        +cf[4+n,l]*a[6+l][i,j,k-1+n])
                    v[i,j,k]= vol12/12
        return v

    
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

