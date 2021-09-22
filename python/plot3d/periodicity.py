from typing import List, Dict, Tuple
from itertools import combinations, product, permutations
import numpy as np
from numpy.core.shape_base import block
from .block import Block, rotate_block
from .face import Face, split_face
from .connectivity import get_face_intersection
from .write import write_plot3D
from math import cos, radians, sin, sqrt, acos, radians
from copy import deepcopy
from tqdm import trange


def create_face(block:Block,imin:int,imax:int,jmin:int,jmax:int,kmin:int,kmax:int) -> Face:
    """Creates a face/surface from IJK bounds. Face = either i is constant or j constant or k constant. 

    Args:
        block (Block): Block object containing X,Y,Z as 3 dimensional arrays 
        imin (int): Minimum I-index
        imax (int): Maximum I-index
        jmin (int): Minimum J-index
        jmax (int): Maximum J-index
        kmin (int): Minimum K-index
        kmax (int): Maximum K-index

    Returns:
        Face: Face object with vertices of I,J,K
    """
    f = Face(4)
    if imin==imax:
        for j in [jmin,jmax]:
            for k in [kmin,kmax]:
                f.add_vertex(block.X[imin,j,k], block.Y[imin,j,k], block.Z[imin,j,k],imin,j,k)
    if jmin==jmax:
        for i in [imin,imax]:
            for k in [kmin,kmax]:
                f.add_vertex(block.X[i,jmin,k], block.Y[i,jmin,k], block.Z[i,jmin,k],i,jmin,k)
    if kmin==kmax:
        for i in [imin,imax]:
            for j in [jmin,jmax]:
                f.add_vertex(block.X[i,j,kmin], block.Y[i,j,kmin], block.Z[i,j,kmin],i,j,kmin)
    return f

# def save_face(face_list:List[Face],block_indices:List[int]):
#     """Converts a list of faces and block indices to a dictionary 

#     Args:
#         face_list (List[Face]): [description]
#         block_indices (List[int]): [description]

#     Returns:
#         (List[dict]): Dictionary specifying the properties of the face 
#     """
#     temp = list()
#     for i,f in enumerate(face_list):
#         temp.append({'block_indx':block_indices[i], 'IMIN':f.IMIN, 'IMAX':f.IMAX, 'JMIN':f.JMIN, 'JMAX':f.JMAX,'KMIN':f.KMIN, 'KMAX':f.KMAX}) 
#     return temp 
def find_periodicity(blocks:List[Block],outer_faces:List, periodic_direction:str='k', rotation_axis:str='x',nblades:int=55):
    """This function is used to check for periodicity of the other faces rotated about an axis 
        The way it works is to find faces of a constant i,j, or k value

    Args:
        blocks (List[Block]): [description]
        outer_faces (List): [description]
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
    """
    

    rotation_angle = radians(360/nblades)
    if rotation_axis=='x':
        rotation_matrix1 = np.array([[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]])
        rotation_matrix2 = np.array([[1,0,0],
                            [0,cos(-rotation_angle),-sin(-rotation_angle)],
                            [0,sin(-rotation_angle),cos(-rotation_angle)]])
    elif rotation_axis=='y':
        rotation_matrix1 = np.array([[cos(rotation_angle),0,sin(rotation_angle)],
                            [0,1,0],
                            [-sin(rotation_angle),0,cos(rotation_angle)]])
        rotation_matrix2 = np.array([[cos(-rotation_angle),0,sin(-rotation_angle)],
                            [0,1,0],
                            [-sin(-rotation_angle),0,cos(-rotation_angle)]])
    elif rotation_axis=='z':
        rotation_matrix1 = np.array([[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]])
        rotation_matrix2 = np.array([[cos(-rotation_angle),-sin(-rotation_angle),0],
                            [sin(-rotation_angle),cos(-rotation_angle),0],
                            [0,0,1]])
    outer_faces_to_keep = list()
    # Check periodic within a block 
    periodic_found = True
    
    # Here we make a list of all the outer faces
    outer_faces_all = list() 
    periodic_faces = list()      # This is the output of the code 

    for o in outer_faces:
        for s in o['surfaces']:
            s['block_indx'] = o['index']
            face = create_face(blocks[o['index']], s['IMIN'], s['IMAX'], s['JMIN'], s['JMAX'], s['KMIN'], s['KMAX'])
            face.set_block_index(o['index'])
            outer_faces_all.append(face)

    split_faces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list
    while periodic_found:
        periodic_found = False
        outer_faces_to_remove = list()  # Integer list of which outher surfaces to remove
        outer_face_combos = list(combinations(range(len(outer_faces_all)),2))
        t = trange(len(outer_face_combos))
        for i in t: 
            # Check if surfaces are periodic with each other
            face1_indx = outer_face_combos[i][0]
            face2_indx = outer_face_combos[i][1]
            face1 = outer_faces_all[face1_indx]
            face2 = outer_faces_all[face2_indx]
            t.set_description(f"Checking connections block {face1.blockIndex} with {face2.blockIndex}")

            if periodic_direction.lower() == "i":
                
                if (face1.IMIN == face1.IMAX) and (face2.IMIN == face2.IMAX):
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                
                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            periodic_faces.extend(periodic_faces_temp)
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        periodic_faces.extend(periodic_faces_temp)
                        split_faces.extend(split_faces_temp)
                        periodic_found = True

            elif periodic_direction.lower() == "j":
                
                if (face1.JMIN == face1.JMAX) and (face2.JMIN == face2.JMAX): 
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)

                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            periodic_faces.extend(periodic_faces_temp)
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        periodic_faces.extend(periodic_faces_temp)
                        split_faces.extend(split_faces_temp)
                        periodic_found = True

            elif periodic_direction.lower() == 'k':       # constant k between surfaces 
                
                if (face1.KMIN == face1.KMAX) and (face2.KMIN == face2.KMAX): 
                    # Rotate Block 1 -> Check periodicity -> if not periodic -> Rotate Block 1 opposite direction -> Check periodicity
                    #   Rotate Block 1
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    #   Check periodicity
                    periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)

                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            periodic_faces.append(periodic_faces_temp)
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        periodic_faces.append(periodic_faces_temp)
                        split_faces.extend(split_faces_temp)
                        periodic_found = True

        if (periodic_found):
            [outer_faces_all.remove(p) for p in outer_faces_to_remove if p in outer_faces_all]
            if len(split_faces)>0:
                outer_faces_all.extend(split_faces)
                split_faces.clear()

    # Export periodic faces and outer faces
    periodic_faces_export = list() 
    outer_faces_export = list() 

    for f in periodic_faces:
        periodic_faces_export.append({
                                'block1':{
                                            'index':f[0].blockIndex,
                                            'IMIN':f[0].IMIN,'JMIN':f[0].JMIN,'KMIN':f[0].KMIN,
                                            'IMAX':f[0].IMAX,'JMAX':f[0].JMAX,'KMAX':f[0].KMAX
                                        },
                                'block2':{
                                            'index':f[1].blockIndex,
                                            'IMIN':f[1].IMIN,'JMIN':f[1].JMIN,'KMIN':f[1].KMIN,
                                            'IMAX':f[1].IMAX,'JMAX':f[1].JMAX,'KMAX':f[1].KMAX
                                        },
                                })

    for o in outer_faces_all:
        outer_faces_export.append(o.to_dict())
                        

    return periodic_faces_export, outer_faces_export
                        


def linear_real_transform(face1:Face,face2:Face) -> Tuple:
    """Computes the rotation angle from Face1 to Face2. This can be used to check if the faces are periodic 
        This function assumes the rotation axis is in the "x" direction. This is good for faces within the same block 

    Reference:
        - Linear Real Transforms from GlennHT https://gitlab.grc.nasa.gov/lte-turbo/GlennHT/-/blob/master/src/M_ccMBMesh.F See computeLRT
        
    Args:
        Face1 (Face): Face to rotate
        Face2 (Face): Face to rotate to

    Returns:
        (tuple): tuple containing:

            - **ang** (float): rotation angle
            - **rotation_matrix** (numpy.ndarray): Rotation matrix 3x3
        
    """

    cTo3,cTo1 = face1.get_corners() 
    cFrom3,cFrom1 = face2.get_corners()

    dTo  = np.array(cTo3).transpose() - np.array(cTo1).transpose()                      # difference in corner points = diagonal vector for Face 1
    ldTo=np.sqrt(np.sum(dTo*dTo))
    if ldTo > 0:
        dTo=dTo/ldTo
    
    dFrom = np.array(cFrom3).transpose() - np.array(cFrom1).transpose()                 # difference in corner points = diagonal vector for Face 2
    ldFrom = np.sqrt(np.sum(dFrom*dFrom))
    if( ldFrom > 0 ):
        dFrom=dFrom/ldFrom
    
    dotprod = np.sum(dTo * dFrom)
    
    if( abs(dotprod-1) < 1E-10 ): # Case of no rotation
        ang = 0

        rotation_matrix  = np.zeros(shape=(3,3))
    else:
        #Compute the angle of rotation  
        cosAng=(dTo[1]*dFrom[1]+dTo[2]*dFrom[2])/sqrt(dTo[1]*dTo[1]+dTo[2]*dTo[2])/sqrt(dFrom[1]*dFrom[1]+dFrom[2]*dFrom[2])
        sinAng=(dTo[2]*dFrom[1]-dTo[1]*dFrom[2])/sqrt(dTo[1]*dTo[1]+dTo[2]*dTo[2])/sqrt(dFrom[1]*dFrom[1]+dFrom[2]*dFrom[2])
        ang=acos(cosAng)

        rotation_matrix = [ [1, 0, 0],
                            [0, cosAng, -sinAng],
                            [0, sinAng, cosAng] ]
        if( sinAng < 0 ):
            ang*=-1
    return ang, rotation_matrix

def __periodicity_check__(face1:Face, face2:Face,block1:Block,block2:Block):
    """General function to find periodicity within a given block. 
    
    Steps:
        - 1: Take the face with the shorter diagonal. 
        - 2: Rotate the shorter face by angle 360/nblades.  
        - 3: Check to see if faces intersect

    Args:
        face1 (Face): An arbitrary face 
        face2 (Face): An arbitrary face 
        blocks (List[Block]): List of all blocks
        rotation_matrix (np.ndarray): rotation matrix
        face1_indx (int): Index of face 1 inside the big array of faces. This is added to the list of faces to remove if periodicity is found
        face2_indx (int): Index of face 1 inside the big array of faces. This is added to the list of faces to remove if periodicity is found
        face1_block_indx (int): what block index face 1 is located in 
        face2_block_indx (int): what block index face 2 is located in 

    Returns:
        (tuple): containing

            - **periodic_surface** (List[Face]):  These are faces that are periodic 
            - **split_surfaces** (List[Face]): Some blocks may have periodic faces with other blocks. But the faces may need to be split so say you pair a small face with a larger face. The split surfaces should be treated as an outer face  
            - **outer_faces_to_remove** (List[int]): Indicies containing [face1_indx, face2_indx]. These are not needed anymore and split_blocks should be added to the list of outer faces 

    """
    outer_faces_to_remove = list() 
    if (face2.diagonal_length < face1.diagonal_length): # switch so that face 2 is always longer
        temp_face = deepcopy(face1)     # Swap face 1 with face 2
        face1 = deepcopy(face1)
        face2 = temp_face
        
        temp = deepcopy(block1)
        block1 = block2
        block2 = temp 
        
    df,split_face1,split_face2 = get_face_intersection(face1,face2,block1,block2)
    
    periodic_faces = list()
    split_faces = list()
    if len(df)>4:
        f1 = create_face(block1,imin=df['i1'].min(),jmin=df['j1'].min(),kmin=df['k1'].min(), imax=df['i1'].max(),jmax=df['j1'].max(),kmax=df['k1'].max())
        f1.set_block_index(face1.blockIndex)
        periodic_faces.append(f1)

        f2 = create_face(block2,imin=df['i2'].min(),jmin=df['j2'].min(),kmin=df['k2'].min(), imax=df['i2'].max(),jmax=df['j2'].max(),kmax=df['k2'].max())
        f2.set_block_index(face2.blockIndex)
        periodic_faces.append(f2)

        # periodic_faces_dict = {
        #         'block1':
        #             {
        #                 'index':face1.blockIndex,'IMIN':df['i1'].min(),'JMIN':df['j1'].min(),'KMIN':df['k1'].min(),
        #                 'IMAX':df['i1'].max(),'JMAX':df['j1'].max(),'KMAX':df['k1'].max()
        #             },
        #         'block2':
        #             {
        #                 'index':face2.blockIndex,'IMIN':df['i2'].min(),'JMIN':df['j2'].min(),'KMIN':df['k2'].min(),
        #                 'IMAX':df['i2'].max(),'JMAX':df['j2'].max(),'KMAX':df['k2'].max()
        #             },
        #         'match':df
        #     }
        split_faces.extend(split_face1)
        split_faces.extend(split_face2)


    return periodic_faces,split_faces