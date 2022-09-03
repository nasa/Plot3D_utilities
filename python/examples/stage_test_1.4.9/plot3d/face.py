import itertools
from typing import Dict, List, Tuple
import numpy as np
from numpy.lib import math
from .block import Block

class Face:
    """Defines a Face of a block for example IMIN,JMIN,JMIN to IMAX,JMIN,JMIN
    """
    def __init__(self,nvertex:int=4):
        """ Defines a face using nvertex 4 = quad 3 = triangles 

        Args:
            nvertex (int, optional): Number of vertices. Defaults to 4.
            id (int, optional): A unique index indentifying a face
        """
        self.x = np.zeros(4)
        self.y = np.zeros(4)
        self.z = np.zeros(4)        
        self.I = np.zeros(4,dtype=np.int64)
        self.J = np.zeros(4,dtype=np.int64)
        self.K = np.zeros(4,dtype=np.int64)
        self.cx = 0         # centroid 
        self.cy = 0
        self.cz = 0
        self.nvertex=0
        self.blockIndex = 0 # not really needed except in periodicity 
        
    def to_dict(self):
        """Returns a dictionary representaon of a face
        """
        return {'IMIN':min(self.I), 'JMIN':min(self.J), 'KMIN':min(self.K),
                    'IMAX':max(self.I), 'JMAX':max(self.J), 'KMAX':max(self.K),
                    'id':0, 'block_index':self.blockIndex}

    @property
    def IMIN(self):
        return self.I.min()
    
    @property
    def JMIN(self):
        return self.J.min()
    
    @property
    def KMIN(self):
        return self.K.min()
    
    @property
    def IMAX(self):
        return self.I.max()
    
    @property
    def JMAX(self):
        return self.J.max()
    
    @property
    def KMAX(self):
        return self.K.max()
    
    @property 
    def BlockIndex(self):
        return self.blockIndex
    
    @property
    def isEdge(self):
        """check if the face is actually an edge. This is an edge if two indicies IMIN == IMAX or JMIN=JMAX or KMIN=KMAX

        Returns:
            [bool]: True if face is really an edge 
        """
        return (int(self.IMIN == self.IMAX) + int(self.JMIN == self.JMAX) + int(self.KMIN == self.KMAX)) > 1

    @property
    def isPoint(self):
        """check if the face is actually an edge. This is an edge if two indicies IMIN == IMAX or JMIN=JMAX or KMIN=KMAX

        Returns:
            [type]: True if face is really a point 
        """
        return (int(self.IMIN == self.IMAX) + int(self.JMIN == self.JMAX) + int(self.KMIN == self.KMAX)) > 2

    @property
    def get_val(self,i_val:int,j_val:int,k_val:int):
        """Get the value where key (I,J,K) is equal to val

        Args:            
            i_val (int): value of I
            j_val (int): value of J
            k_val (int): value of K

        Returns:
            [float]: x value
            [float]: y value
            [float]: z value
        """
        
        indx_i = np.where(self.I == i_val).tolist()
        indx_j = np.where(self.J == j_val).tolist()
        indx_k = np.where(self.K == k_val).tolist()
        
        indx_i.extend(indx_j)
        indx_i.extend(indx_k)

        indx = list(set([indx_i]))[0] # Get the common one through a union
        return self.x[indx], self.y[indx], self.z[indx]


    def add_vertex(self, x:float,y:float,z:float, i:int, j:int, k:int):
        """Add vertex to define a face

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
            z (float): z-coordinate
            i (int): i-th index of the coordinates (x,y,z)
            j (int): j-th index of the coordinates (x,y,z)
            k (int): k-th index of the coordinates (x,y,z)
        """
        
        self.x[self.nvertex] = x
        self.y[self.nvertex] = y 
        self.z[self.nvertex] = z 
        self.I[self.nvertex] = i
        self.J[self.nvertex] = j
        self.K[self.nvertex] = k        
        self.nvertex+=1    
        if self.nvertex==4:
            self.cx = self.x.mean()
            self.cy = self.y.mean()
            self.cz = self.z.mean()
    @property
    def size(self):
        if self.IMIN==self.IMAX:
            return (self.JMAX- self.JMIN)*(self.KMAX-self.KMIN)
        elif (self.JMIN==self.JMAX):
            return (self.IMAX-self.IMIN)*(self.KMAX-self.KMIN)
        elif (self.KMIN==self.KMAX):
            return (self.IMAX-self.IMIN)*(self.JMAX- self.JMIN)
        else:
            return (self.IMAX-self.IMIN)*(self.JMAX- self.JMIN)*(self.KMAX-self.KMIN)

    def set_block_index(self,val):
        self.blockIndex = val
        
    def __normal__(self):
        """Computes the normal vector of the face 
            not really used but if anyone wants it. 
        """
        if (self.I[0]!=self.I[1]) and (self.I[0]!=self.I[2]):
            indx = np.argsort(self.I)
        elif (self.J[0]!=self.J[1]) and (self.J[0]!=self.J[2]):
            indx = np.argsort(self.J)
        elif (self.K[0]!=self.K[1]) and (self.K[0]!=self.K[2]):
            indx = np.argsort(self.K)

        self.x = self.x[indx]
        self.y = self.y[indx]
        self.z = self.z[indx]
        self.I = self.I[indx]
        self.J = self.J[indx]
        self.K = self.K[indx]
        x1 = self.x[1]-self.x[0]; y1 = self.y[1]-self.y[0]; z1 = self.z[1]-self.z[0]
        x2 = self.x[2]-self.x[0]; y2 = self.y[2]-self.y[0]; z2 = self.z[2]-self.z[0]
        nx = y1*z2-y2*z1; ny = -1*(x1*z2-x2*z1); nz = x1*y2-x2*y1
        self.nx = nx
        self.ny = ny
        self.nz = nz 

    def match_indices(self,f):
        """Check to see if two faces are the same. Checks to see if any of vertices x,y,z match
            Normally this is used by Face1==Face2

        Args:
            f (Face): another face

        Returns:
            List[(int,int)]: list of indicies where there's a match. 
        """
        matched_vertices = list()
        tol = 1E-6
        matchedIndices = list()
        for i in range(self.nvertex):
            for j in range(f.nvertex):
                dx = abs(self.x[i] - f.x[j])
                dy = abs(self.y[i] - f.y[j])
                dz = abs(self.z[i] - f.z[j])
                if (dx<tol and dy<tol and dz<tol and (j not in matched_vertices)):
                    matchedIndices.append([i,j])
                    matched_vertices.append(j) # This vertex has been matched, remove from list
                    break # each index can only have a single match
        return matchedIndices

    def __eq__(self, f):
        """Check to see if two faces are the same by looking at the I,J,K
        Checks to see if any of vertices x,y,z match

        Args:
            f (Face): another face 

        Returns:
            Boolean: True if faces match, False if no match is found 
        """
        # matchedIndices = self.match_indices(f)
        # (len(matchedIndices)==self.nvertex) and
        return ((self.BlockIndex == f.BlockIndex) 
            and (self.IMIN == f.IMIN) and (self.IMAX == f.IMAX) 
            and (self.JMIN == f.JMIN) and (self.JMAX == f.JMAX) 
            and (self.KMIN == f.KMIN) and (self.KMAX == f.KMAX) )
    
    def vertices_equals(self,f):
        """Checks to see if two faces are the same by looking at the vertices

        Args:
            f (Face): Another face

        Returns:
            bool: True = face vertices are equal
        """
        matchedIndices = self.match_indices(f)
        return (len(matchedIndices)==self.nvertex)

    def __ne__(self,f):
        """Checks if two faces are not equal 

        Args:
            f (Face): another face 

        Returns:
            Boolean: True if faces match, False if no match is found 
        """
        match = self.__eq__(f)
        return not match
    
    def index_equals(self,f2):
        """Check to see of the face indices are equal

        Args:
            f2 ([type]): [description]
        """
        if (self.IMIN == f2.IMIN and 
            self.JMIN == f2.JMIN and 
            self.KMIN == f2.KMIN and 
            self.IMAX == f2.IMAX and 
            self.JMAX == f2.JMAX and 
            self.KMAX == f2.KMAX):
            return True
    def __hash__(self):
        if (len(self.I)>0):
            return hash((self.I[0], self.J[0], self.K[0], self.I[-1], self.J[-1], self.K[-1]))
        else:
            return hash((0, 0, 0, 0, 0, 0))
            
    def __str__(self):
        if (len(self.I)>0):
            return 'blk: {:d} [{:d},{:d},{:d},{:d},{:d},{:d}]'.format(self.blockIndex,self.IMIN, self.JMIN, self.KMIN, self.IMAX, self.JMAX, self.KMAX)
        else:
            return 'blk: {:d} [{:d},{:d},{:d},{:d},{:d},{:d}]'.format(self.blockIndex,0,0,0,0,0,0)
    
    def __repr__(self):
        return str(self)
    
    @property
    def diagonal_length(self) -> float:
        """Returns the diagonal length of the face 

        Returns:
            float: diagonal length computed using IMIN, IMAX, JMIN, JMAX, KMIN, KMAX 
        """
        minIndx = 0; maxIndx = 0 
        for indx in range(len(self.I)):
            if self.I[indx] == self.IMIN and self.J[indx] == self.JMIN and self.K[indx] == self.KMIN:
                minIndx = indx
            if self.I[indx] == self.IMAX and self.J[indx] == self.JMAX and self.K[indx] == self.KMAX:
                maxIndx = indx
        dx = self.x[minIndx] - self.x[maxIndx]
        dy = self.y[minIndx] - self.y[maxIndx]
        dz = self.z[minIndx] - self.z[maxIndx]
        return math.sqrt(dx*dx + dy*dy + dz*dz)

    def get_corners(self) -> Tuple:
        """Get the corners defined by (IMIN,JMIN,KMIN), (IMAX,JMAX,KMAX),
        
        Returns:
            Tuple: containing

                - **(x,y,z)** (float,float,float): at IMIN,JMIN,KMIN
                - **(x,y,z)** (float,float,float): at IMAX,JMAX,KMAX
        
        Reference: 
            - GlennHT source code https://gitlab.grc.nasa.gov/lte-turbo/GlennHT/-/blob/master/src/M_ccMBMesh.F function computeLRT

        """
        minIndx = 0; maxIndx = 0 
        for indx in range(len(self.I)):
            if self.I[indx] == self.IMIN and self.J[indx] == self.JMIN and self.K[indx] == self.KMIN:
                minIndx = indx
            if self.I[indx] == self.IMAX and self.J[indx] == self.JMAX and self.K[indx] == self.KMAX:
                maxIndx = indx
        return (self.x[minIndx],self.y[minIndx], self.z[minIndx]),(self.x[maxIndx],self.y[maxIndx], self.z[maxIndx])
        


def create_face_from_diagonals(block:Block,imin:int,jmin:int,kmin:int,imax:int,jmax:int,kmax:int):
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