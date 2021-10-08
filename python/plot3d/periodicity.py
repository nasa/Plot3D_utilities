from typing import List, Dict, Tuple
from itertools import combinations_with_replacement
import numpy as np
from .block import Block, rotate_block
from .face import Face, create_face_from_diagonals, split_face
from .connectivity import get_face_intersection, face_matches_to_dict
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
def periodicity(blocks:List[Block],outer_faces:List[Face], matched_faces:List[Face], periodic_direction:str='k', rotation_axis:str='x',nblades:int=55):
    """This function is used to check for periodicity of the other faces rotated about an axis 
        The way it works is to find faces of a constant i,j, or k value

    Args:
        blocks (List[Block]): [description]
        outer_faces (List): [description]
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
    """
    

    rotation_angle = radians(360.0/nblades)
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
        rotation_matrix2 = np.array([[1,0,0],
                            [0,cos(-rotation_angle),-sin(-rotation_angle)],
                            [0,sin(-rotation_angle),cos(-rotation_angle)]])
    # Check periodic within a block 
    periodic_found = True
    
    # Here we make a list of all the outer faces
    outer_faces_all = list() 
    matched_faces_all = list()
    periodic_faces = list()      # This is the output of the code 
    periodic_faces_export = list() 
    for o in outer_faces:
        face = create_face(blocks[o['block_index']], o['IMIN'], o['IMAX'], o['JMIN'], o['JMAX'], o['KMIN'], o['KMAX'])
        face.set_block_index(o['block_index'])
        outer_faces_all.append(face)
    for m in matched_faces:
        face1 = create_face(blocks[m['block1']['block_index']], m['block1']['IMIN'], m['block1']['IMAX'], m['block1']['JMIN'], m['block1']['JMAX'], m['block1']['KMIN'], m['block1']['KMAX'])
        face2 = create_face(blocks[m['block2']['block_index']], m['block2']['IMIN'], m['block2']['IMAX'], m['block2']['JMIN'], m['block2']['JMAX'], m['block2']['KMIN'], m['block2']['KMAX'])
        face1.set_block_index(m['block1']['block_index'])
        face2.set_block_index(m['block2']['block_index'])
        matched_faces_all.append(face1)
        matched_faces_all.append(face2)
    # Debug
    # debug_per1 = create_face(blocks[2], 0, 40, 0, 100,0,0) 
    # debug_per1.set_block_index(2)
    # debug_per2 = create_face(blocks[2], 0, 40, 0, 100,52,52) 
    # debug_per2.set_block_index(2)
    # End debug 

    split_faces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list
    while periodic_found:
        periodic_found = False        
        outer_faces_to_remove = list()  # Integer list of which outher surfaces to remove
        outer_face_combos = list(combinations_with_replacement(range(len(outer_faces_all)),2))
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
                    # Rotate Block 1 -> Check periodicity -> if not periodic -> Rotate Block 1 opposite direction -> Check periodicity
                    #   Rotate Block 1
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    #   Check periodicity
                    df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)

                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            outer_faces_to_remove.append(periodic_faces_temp[0])    # Make sure periodic faces are also removed from outer faces during the loop
                            outer_faces_to_remove.append(periodic_faces_temp[1])
                            periodic_faces.append(periodic_faces_temp)
                            periodic_faces_export.append(face_matches_to_dict(face1,face2,block1_rotated,block2))
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                            break
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        outer_faces_to_remove.append(periodic_faces_temp[0])
                        outer_faces_to_remove.append(periodic_faces_temp[1])
                        periodic_faces.append(periodic_faces_temp)
                        periodic_faces_export.append(face_matches_to_dict(face1,face2,block1_rotated,block2))
                        split_faces.extend(split_faces_temp)
                        periodic_found = True
                        break

            elif periodic_direction.lower() == "j":
                
                if (face1.JMIN == face1.JMAX) and (face2.JMIN == face2.JMAX): 
                    # Rotate Block 1 -> Check periodicity -> if not periodic -> Rotate Block 1 opposite direction -> Check periodicity
                    #   Rotate Block 1
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    #   Check periodicity
                    df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)

                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            outer_faces_to_remove.append(periodic_faces_temp[0])
                            outer_faces_to_remove.append(periodic_faces_temp[1])
                            periodic_faces.append(periodic_faces_temp)
                            periodic_faces_export.append(face_matches_to_dict(face1,face2,block1_rotated,block2))
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                            break
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        outer_faces_to_remove.append(periodic_faces_temp[0])
                        outer_faces_to_remove.append(periodic_faces_temp[1])
                        periodic_faces.append(periodic_faces_temp)
                        periodic_faces_export.append(face_matches_to_dict(face1,face2,block1_rotated,block2))
                        split_faces.extend(split_faces_temp)
                        periodic_found = True
                        break

            elif periodic_direction.lower() == 'k':       # constant k between surfaces 
                 if (face1.KMIN == face1.KMAX) and (face2.KMIN == face2.KMAX): 
                    # Rotate Block 1 -> Check periodicity -> if not periodic -> Rotate Block 1 opposite direction -> Check periodicity
                    #   Rotate Block 1
                    block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix1)
                    block2 = blocks[face2.blockIndex]
                    #   Check periodicity
                    df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)

                    if len(periodic_faces_temp) == 0:
                        block1_rotated = rotate_block(blocks[face1.blockIndex],rotation_matrix2)
                        block2 = blocks[face2.blockIndex]
                        df, periodic_faces_temp, split_faces_temp = __periodicity_check__(face1,face2,block1_rotated, block2)
                        if len(periodic_faces_temp)>0:  # Save the data
                            outer_faces_to_remove.append(face1)
                            outer_faces_to_remove.append(face2)
                            outer_faces_to_remove.append(periodic_faces_temp[0])
                            outer_faces_to_remove.append(periodic_faces_temp[1])
                            periodic_faces.append(periodic_faces_temp)
                            periodic_faces_export.append(face_matches_to_dict(periodic_faces_temp[0],periodic_faces_temp[1],block1_rotated,block2))
                            split_faces.extend(split_faces_temp)
                            periodic_found = True
                            break
                    else:   # Save the data
                        outer_faces_to_remove.append(face1)
                        outer_faces_to_remove.append(face2)
                        outer_faces_to_remove.append(periodic_faces_temp[0])
                        outer_faces_to_remove.append(periodic_faces_temp[1])
                        periodic_faces.append(periodic_faces_temp)
                        periodic_faces_export.append(face_matches_to_dict(periodic_faces_temp[0],periodic_faces_temp[1],block1_rotated,block2))
                        split_faces.extend(split_faces_temp)
                        periodic_found = True
                        break

        if (periodic_found):
            outer_faces_to_remove = list(set(outer_faces_to_remove))
            outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]
            if len(split_faces)>0:
                outer_faces_all.extend(split_faces)
                split_faces.clear()

    # This is an added check to make sure all periodic faces are in the outer_faces_to_remove
    for p in periodic_faces:
        outer_faces_to_remove.append(p[0])
        outer_faces_to_remove.append(p[1])

    for m in matched_faces_all:
        outer_faces_to_remove.append(m)

    outer_faces_to_remove = list(set(outer_faces_to_remove))    # Use only unique values
    outer_faces_all = [p for p in outer_faces_all if p not in outer_faces_to_remove]    # remove from outer faces 
    # remove any duplicate periodic face pairs 
    indx_to_remove = list()
    for i in range(len(periodic_faces)):
        for j in range(i+1,len(periodic_faces)):
            if periodic_faces[i][0] == periodic_faces[j][0]:
                if periodic_faces[i][1] == periodic_faces[j][1]:
                    indx_to_remove.append(j)
            if periodic_faces[i][1] == periodic_faces[j][0]:
                if periodic_faces[i][0] == periodic_faces[j][1]:
                    indx_to_remove.append(j)
    

    periodic_faces_export = [periodic_faces_export[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    periodic_faces = [periodic_faces[i] for i in range(len(periodic_faces)) if i not in indx_to_remove]
    # Export periodic faces and outer faces
    outer_faces_export = list() 

    for o in outer_faces_all:
        outer_faces_export.append(o.to_dict())
                        
    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all
                        


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

            - **df** (pandas.Dataframe): List of point matches for periodic surfaces 
            - **periodic_surface** (List[Face]):  These are faces that are periodic 
            - **split_surfaces** (List[Face]): Some blocks may have periodic faces with other blocks. But the faces may need to be split so say you pair a small face with a larger face. The split surfaces should be treated as an outer face  

    """
    
    periodic_faces = list()
    split_faces = list()
    swapped = False
    if (face2.diagonal_length < face1.diagonal_length): # switch so that face 2 is always longer
        temp = deepcopy(face1)
        face1 = deepcopy(face2)
        face2 = temp

        temp_block = deepcopy(block1)
        block1 = deepcopy(block2)
        block2 = temp_block        
        swapped = True
    
    df,split_face1,split_face2 = get_face_intersection(face1,face2,block1,block2)

    if len(df)>4:
        f1 = create_face(block1,imin=df['i1'].min(),jmin=df['j1'].min(),kmin=df['k1'].min(), imax=df['i1'].max(),jmax=df['j1'].max(),kmax=df['k1'].max())
        f1.set_block_index(face1.blockIndex)

        f2 = create_face(block2,imin=df['i2'].min(),jmin=df['j2'].min(),kmin=df['k2'].min(), imax=df['i2'].max(),jmax=df['j2'].max(),kmax=df['k2'].max())
        f2.set_block_index(face2.blockIndex)

        split_faces.extend(split_face1)
        split_faces.extend(split_face2)
        if swapped:
            periodic_faces.append(f2)
            periodic_faces.append(f1)
        else:
            periodic_faces.append(f1)
            periodic_faces.append(f2)

    return df,periodic_faces,split_faces