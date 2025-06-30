import numpy as np 
import math 
from tqdm import trange
from typing import List
from copy import deepcopy
import numpy.typing as npt 

class Block:
    """Plot3D Block definition
    """
    X: npt.NDArray
    Y: npt.NDArray
    Z: npt.NDArray
    IMAX:int
    JMAX:int
    KMAX:int
    def __init__(self, X:npt.NDArray,Y:npt.NDArray,Z:npt.NDArray):
        """Initializes the block using all the X,Y,Z coordinates of the block

        Args:
            X (npt.NDArray): All the X coordinates (i,j,k)
            Y (npt.NDArray): All the Y coordinates (i,j,k)
            Z (npt.NDArray): All the Z coordinates (i,j,k)

        """
        self.IMAX,self.JMAX,self.KMAX = X.shape; 
        self.X = X
        self.Y = Y
        self.Z = Z
        # Centroid 
        self.cx = np.mean(X) 
        self.cy = np.mean(Y)
        self.cz = np.mean(Z)
        
    def __repr__(self):
        return f"({self.IMAX},{self.JMAX},{self.KMAX}"
    
    def scale(self,factor:float):
        """Scales a mesh by a certain factor 

        Args:
            factor (float): _description_
        """
        
        self.X *= factor
        self.Y *= factor
        self.Z *= factor

    def shift(self,shift_amount:float,direction:str="z"):
        """shifts the blocks by a certain amount

        Args:
            shift_amount (float): _description_
            direction (str, optional): _description_. Defaults to "z".
        """
        if direction.lower() == 'z':
            self.Z +=shift_amount
        elif direction.lower() == 'y':
            self.Y +=shift_amount
        elif direction.lower() == 'x':
            self.X +=shift_amount

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
    
    def get_faces(self):
        """
        Returns a dictionary of the six faces of the block.
        Each face is a tuple of (X_face, Y_face, Z_face).
        """
        return {
            'imin': (self.X[0,:,:], self.Y[0,:,:], self.Z[0,:,:]),
            'imax': (self.X[-1,:,:], self.Y[-1,:,:], self.Z[-1,:,:]),
            'jmin': (self.X[:,0,:], self.Y[:,0,:], self.Z[:,0,:]),
            'jmax': (self.X[:,-1,:], self.Y[:,-1,:], self.Z[:,-1,:]),
            'kmin': (self.X[:,:,0], self.Y[:,:,0], self.Z[:,:,0]),
            'kmax': (self.X[:,:,-1], self.Y[:,:,-1], self.Z[:,:,-1]),
        }
        
    @property
    def size(self)->int:
        """returns the total number of nodes 

        Returns:
            int: number of nodes
        """
        return self.IMAX*self.JMAX*self.KMAX


def faces_match(face1, face2, tol=1e-12):
    """
    Quickly compare two block faces by checking if their corners and diagonals match geometrically.
    Faces are (X, Y, Z) arrays of shape (m, n).
    """
    def extract_key_points(X, Y, Z):
        # Flattened index corners: top-left, top-right, bottom-left, bottom-right
        corners = [
            (X[0,0], Y[0,0], Z[0,0]),
            (X[0,-1], Y[0,-1], Z[0,-1]),
            (X[-1,0], Y[-1,0], Z[-1,0]),
            (X[-1,-1], Y[-1,-1], Z[-1,-1])
        ]
        # Add diagonals: [0,0] to [-1,-1] and [0,-1] to [-1,0]
        diag = [
            ((X[0,0] + X[-1,-1]) / 2, (Y[0,0] + Y[-1,-1]) / 2, (Z[0,0] + Z[-1,-1]) / 2),
            ((X[0,-1] + X[-1,0]) / 2, (Y[0,-1] + Y[-1,0]) / 2, (Z[0,-1] + Z[-1,0]) / 2),
        ]
        return np.array(corners + diag)

    X1, Y1, Z1 = face1
    X2, Y2, Z2 = face2

    if X1.shape != X2.shape:
        return False

    pts1 = extract_key_points(X1, Y1, Z1)
    pts2 = extract_key_points(X2, Y2, Z2)

    # Sort points by coordinates to align order-insensitive comparison
    pts1_sorted = np.array(sorted(pts1.tolist()))
    pts2_sorted = np.array(sorted(pts2.tolist()))

    d = np.linalg.norm(pts1_sorted - pts2_sorted, axis=1)
    return np.all(d <= tol)


def find_matching_faces(block1:Block, block2:Block,tol:float=1E-12):
    """
    Return the matching face pair (name1, name2) if any faces match
    between block1 and block2.
    """
    faces1 = block1.get_faces()
    faces2 = block2.get_faces()
    
    for name1, f1 in faces1.items():
        for name2, f2 in faces2.items():
            if faces_match(f1, f2, tol):
                return name1, name2
    return None, None

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


def combine_8_blocks(blocks: List[Block], tol=1e-8) -> Block:
    """
    Combine exactly 8 structured Block objects into a 2×2×2 cube.
    Groups blocks by shared face match direction before merging.

    Parameters
    ----------
    blocks : List[Block]
        Exactly 8 blocks that form a 2×2×2 cube via face adjacency.
    tol : float
        Tolerance used for face matching.

    Returns
    -------
    Block
        The merged Block.

    Raises
    ------
    RuntimeError
        If not all blocks can be merged.
    """
    if len(blocks) != 8:
        raise ValueError("combine_8_blocks requires exactly 8 Block objects.")

    from itertools import combinations

    remaining = blocks.copy()
    merged = None
    used_ids = set()

    # Start with the first two blocks that match on any face
    for a, b in combinations(remaining, 2):
        face1, face2 = find_matching_faces(a, b, tol=tol)
        if face1 is not None:
            print(f"Starting with face pair: {face1}-{face2}")
            merged = combine_blocks(a, b, tol=tol)
            used_ids.update({id(a), id(b)})
            remaining.remove(a)
            remaining.remove(b)
            break

    if merged is None:
        raise RuntimeError("❌ No initial face match found to begin merge.")

    # Try to keep merging the rest, prioritizing blocks matching the last merged face
    while remaining:
        progress = False
        for blk in remaining:
            face1, face2 = find_matching_faces(merged, blk, tol=tol)
            if face1 is not None:
                merged = combine_blocks(merged, blk, tol=tol)
                used_ids.add(id(blk))
                remaining.remove(blk)
                progress = True
                break

        if not progress:
            # Try pairwise merges within remaining blocks
            for a, b in combinations(remaining, 2):
                face1, face2 = find_matching_faces(a, b, tol=tol)
                if face1 is not None:
                    temp_merge = combine_blocks(a, b, tol=tol)
                    remaining.remove(a)
                    remaining.remove(b)
                    merged = combine_blocks(merged, temp_merge, tol=tol)
                    used_ids.update({id(a), id(b)})
                    progress = True
                    break

        if not progress:
            raise RuntimeError(f"❌ Unable to merge all 8 blocks. {len(remaining)} remain.")

    return merged


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

