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
i

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


def find_periodicity(blocks:List[Block],outer_faces:List, periodic_direction:str='k', rotation_axis:str='x',nblades:int=55):
    """This function is used to check for periodicity of the other faces rotated about an axis 
        The way it works is to find faces of a constant i,j, or k value

    Args:
        blocks (List[Block]): [description]
        outer_faces (List): [description]
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
    """
    
    periodic_surfaces = list()      # This is the output of the code 
    rotation_angle = radians(360/nblades)
    
    if rotation_axis=='x':
        rotation_matrix1 = [[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]]
        rotation_matrix2 = [[1,0,0],
                            [0,cos(-rotation_angle),-sin(-rotation_angle)],
                            [0,sin(-rotation_angle),cos(-rotation_angle)]]
    elif rotation_axis=='y':
        rotation_matrix1 = [[cos(rotation_angle),0,sin(rotation_angle)],
                            [0,1,0],
                            [-sin(rotation_angle),0,cos(rotation_angle)]]
        rotation_matrix2 = [[cos(-rotation_angle),0,sin(-rotation_angle)],
                            [0,1,0],
                            [-sin(-rotation_angle),0,cos(-rotation_angle)]]
    elif rotation_axis=='z':
        rotation_matrix1 = [[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]]
        rotation_matrix2 = [[cos(-rotation_angle),-sin(-rotation_angle),0],
                            [sin(-rotation_angle),cos(-rotation_angle),0],
                            [0,0,1]]
    outer_faces_to_keep = list()
    # Check periodic within a block 
    block_match = True

    block0_periodic = create_face(blocks[0],0,120,0,100,32,32)
    block1_periodic = create_face(blocks[3],0,72,0,100,52,52)
    
    # Here we make a list of all the outer faces
    outer_faces_all = list() 
    for o in outer_faces:
        for s in o['surfaces']:
            s['block_indx'] = o['index']
            outer_faces_all.append(s)

    while block_match:
        block_match = False
        outer_faces_to_remove = list()  # Integer list of which outher surfaces to remove
        split_faces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list
        outer_face_combos = list(combinations(range(len(outer_faces_all)),2))
        for face1_indx, face2_indx in outer_face_combos:        
            # Check if surfaces are periodic with each other
            face1 = create_face(blocks[outer_faces_all[face1_indx]['block_indx']], outer_faces_all[face1_indx]['IMIN'], outer_faces_all[face1_indx]['IMAX'], 
                                                    outer_faces_all[face1_indx]['JMIN'], outer_faces_all[face1_indx]['JMAX'], 
                                                    outer_faces_all[face1_indx]['KMIN'], outer_faces_all[face1_indx]['KMAX'])

            face2 = create_face(blocks[outer_faces_all[face2_indx]['block_indx']], outer_faces_all[face2_indx]['IMIN'], outer_faces_all[face2_indx]['IMAX'], 
                                            outer_faces_all[face2_indx]['JMIN'], outer_faces_all[face2_indx]['JMAX'], 
                                            outer_faces_all[face2_indx]['KMIN'], outer_faces_all[face2_indx]['KMAX'])

            is_periodic = False
            if periodic_direction.lower() == "i":
                if (face1['IMIN'] == face1['IMAX']) and (face2['IMIN'] == face2['IMAX']): 
                    is_periodic, split_surf, faces_to_remove, angle, rot_mat = __periodicity_check__(face1,face2,blocks,rotation_matrix,face1_indx, face2_indx, outer_faces_all[face2_indx]['block_indx'],outer_faces_all[face2_indx]['block_indx'])
            elif periodic_direction.lower() == "j":
                if (face1['JMIN'] == face1['JMAX']) and (face2['JMIN'] == face2['JMAX']): 
                    is_periodic, split_surf,faces_to_remove, angle, rot_mat = __periodicity_check__(face1,face2,blocks,rotation_matrix,face1_indx, face2_indx, outer_faces_all[face2_indx]['block_indx'],outer_faces_all[face2_indx]['block_indx']))
            else:  # periodic_direction.lower() == 'k':       # constant k between surfaces 
                if (face1['KMIN'] == face1['KMAX']) and (face2['KMIN'] == face2['KMAX']): 
                    rotation_matrix = [[1,0,0],
                                       [0,cos(rotation_angle),-sin(rotation_angle)],
                                       [0,sin(rotation_angle),cos(rotation_angle)]]
                    is_periodic, split_surf, faces_to_remove, angle, rot_mat = __periodicity_check__(face1,face2,blocks,rotation_matrix,face1_indx, face2_indx, outer_faces_all[face2_indx]['block_indx'],outer_faces_all[face2_indx]['block_indx']))
            
            if is_periodic:
                periodic_surfaces.append(is_periodic)
                outer_faces_to_remove.extend(faces_to_remove)
                split_surfaces.extend(split_surf)
                periodic_angle = angle
                rotation_matrix = rot_mat
                block_match = True
                # break
                
        # if block_match:
        #     break
    # House keeping - Remove outer surfaces and add in the split surfaces 
    block['surfaces'] = [block['surfaces'][o] for o in range(len(block['surfaces'])) if o not in outer_faces_to_remove]
    outer_faces[indx]['surfaces'] = block['surfaces']
    for s in split_surfaces:
        outer_faces[indx]['surfaces'].append({'IMIN':s.IMIN,'JMIN':s.JMIN,'KMIN':s.KMIN,'IMAX':s.IMAX,'JMAX':s.JMAX,'KMAX':s.KMAX,'id':0})
    
    block_outer_faces_to_remove = list()
    if rotation_matrix is not None:
        # Check periodic outer_faces from block to block 
        block_combos = list(combinations(range(len(outer_faces)),2))
        for b1,b2 in block_combos: 
            # Get outer surfaces for block 1 and block 2
            # https://www.geeksforgeeks.org/python-program-to-get-all-unique-combinations-of-two-lists/ 
            for b1_surf_indx in range(len(outer_faces[b1]['surfaces'])):
                for b2_surf_indx in range(len(outer_faces[b2]['surfaces'])):
                    # Creates two surfaces from b1_surf and b2_surf
                    b1_surf = outer_faces[b1]['surfaces'][b1_surf_indx]
                    b2_surf = outer_faces[b2]['surfaces'][b2_surf_indx]

                    surface1 = create_face(blocks[b1], b1_surf['IMIN'], b1_surf['IMAX'], 
                                                            b1_surf['JMIN'], b1_surf['JMAX'], 
                                                            b1_surf['KMIN'], b1_surf['KMAX'])

                    surface2 = create_face(blocks[b2], b2_surf['IMIN'], b2_surf['IMAX'], 
                                                            b2_surf['JMIN'], b2_surf['JMAX'], 
                                                            b2_surf['KMIN'], b2_surf['KMAX'])
                    
                    
                    if periodic_direction == "i" and (b1_surf['IMIN'] == b1_surf['IMAX']) and (b2_surf['IMIN'] == b2_surf['IMAX']):
                        block1_rotated = rotate_block(blocks[b1], rotation_matrix) 
                        df, _, _ = get_face_intersection(surface1,surface2,block1_rotated,blocks[b2],tol=1E-5)
                            
                        if len(df)==0:                        # Try with second block rotated
                            # Step 2: if Step 1 does not yield a match then rotate surface 2 and check for intersections 
                            block2_rotated = rotate_block(blocks[b2], rotation_matrix)
                            df, _, _ = get_face_intersection(surface1,surface2,blocks[b1],block2_rotated,tol=1E-5)
                           
                        if len(df)>0:
                            df['j1'] += surface1.JMIN           
                            df['k1'] += surface1.KMIN
                            df['j2'] += surface2.JMIN
                            df['k2'] += surface2.KMIN
                            periodic_surfaces.append({
                                                        'block1':{
                                                                'index':b1,
                                                                'IMIN':df.iloc[0]['i1'],'JMIN':df.iloc[0]['j1'],'KMIN':df.iloc[0]['k1'],
                                                                'IMAX':df.iloc[-1]['i1'],'JMAX':df.iloc[-1]['j1'],'KMAX':df.iloc[-1]['k1']
                                                                },
                                                        'block2':{
                                                                'index':b2,
                                                                'IMIN':df.iloc[0]['i2'],'JMIN':df.iloc[0]['j2'],'KMIN':df.iloc[0]['k2'],
                                                                'IMAX':df.iloc[-1]['i2'],'JMAX':df.iloc[-1]['j2'],'KMAX':df.iloc[-1]['k2']
                                                                },
                                                        })
                            # query faces to remove
                            block_outer_faces_to_remove.append({'block_indx':b1,'surface_indx':b1_surf_indx})
                            block_outer_faces_to_remove.append({'block_indx':b2,'surface_indx':b2_surf_indx})

                    elif periodic_direction == "j" and (b1_surf['JMIN'] == b1_surf['JMAX']) and (b2_surf['JMIN'] == b2_surf['JMAX']):
                        block1_rotated = rotate_block(blocks[b1], rotation_matrix) # Rotate the block 
                        df, _, _ = get_face_intersection(surface1,surface2,block1_rotated,blocks[b2],tol=1E-5)

                        if len(df)==0:                        # Try with second block rotated
                            # Step 2: if Step 1 does not yield a match then rotate surface 2 and check for intersections 
                            block2_rotated = rotate_block(blocks[b2], rotation_matrix)
                            df, _, _ = get_face_intersection(surface1,surface2,blocks[b1],block2_rotated,tol=1E-5)

                        if len(df)>0:
                            df['i1'] += surface1.IMIN
                            df['k1'] += surface1.KMIN
                            df['i2'] += surface2.IMIN
                            df['k2'] += surface2.KMIN
                            periodic_surfaces.append({
                                                        'block1':{
                                                            'index':b1,
                                                            'IMIN':df.iloc[0]['i1'],'JMIN':df.iloc[0]['j1'],'KMIN':df.iloc[0]['k1'],
                                                            'IMAX':df.iloc[-1]['i1'],'JMAX':df.iloc[-1]['j1'],'KMAX':df.iloc[-1]['k1']
                                                            },
                                                        'block2':{
                                                                'index':b2,
                                                                'IMIN':df.iloc[0]['i2'],'JMIN':df.iloc[0]['j2'],'KMIN':df.iloc[0]['k2'],
                                                                'IMAX':df.iloc[-1]['i2'],'JMAX':df.iloc[-1]['j2'],'KMAX':df.iloc[-1]['k2']
                                                            },
                                                        })
                            # query faces to remove
                            block_outer_faces_to_remove.append({'block_indx':b1,'surface_indx':b1_surf_indx})
                            block_outer_faces_to_remove.append({'block_indx':b2,'surface_indx':b2_surf_indx})

                    # If the k's are constant then rotate the block by a rotation matrix
                    elif periodic_direction == "k" and (b1_surf['KMIN'] == b1_surf['KMAX']) and (b2_surf['KMIN'] == b2_surf['KMAX']):
                        block1_rotated = rotate_block(blocks[b1], rotation_matrix) # Rotate block 1
                        
                        # write_plot3D('test_rotated.xyz',[block1_rotated],binary=True) # debug purposes
                        # Step 1: Try rotating surface 1 and check for intersections 
                        df, _, _ = get_face_intersection(surface1,surface2,block1_rotated,blocks[b2],tol=1E-5)
                        if len(df)==0:
                            block2_rotated = rotate_block(blocks[b2], rotation_matrix) # Rotate block 2
                            df, _, _ = get_face_intersection(surface1,surface2,blocks[b1],block2_rotated,tol=1E-5)
                        
                        if len(df)>0:
                            df['i1'] += surface1.IMIN
                            df['j1'] += surface1.JMIN
                            df['i2'] += surface2.IMIN
                            df['j2'] += surface2.JMIN
                            # Here we save the periodicity 
                            periodic_surfaces.append({
                                                        'block1':{
                                                            'index':b1,'IMIN':df.iloc[0]['i1'],'JMIN':df.iloc[0]['j1'],'KMIN':df.iloc[0]['k1'],
                                                            'IMAX':df.iloc[-1]['i1'],'JMAX':df.iloc[-1]['j1'],'KMAX':df.iloc[-1]['k1']
                                                            },
                                                        'block2':{
                                                                'index':b2,'IMIN':df.iloc[0]['i2'],'JMIN':df.iloc[0]['j2'],'KMIN':df.iloc[0]['k2'],
                                                                'IMAX':df.iloc[-1]['i2'],'JMAX':df.iloc[-1]['j2'],'KMAX':df.iloc[-1]['k2']
                                                            },
                                                        })
                            # query faces to remove
                            block_outer_faces_to_remove.append({'block_indx':b1,'surface_indx':b1_surf_indx})
                            block_outer_faces_to_remove.append({'block_indx':b2,'surface_indx':b2_surf_indx})
            
        # Lets remove the outer faces
        for b in range(len(outer_faces)):
            outer_faces_to_keep.append({'index':b, 'surfaces':[]})
            for s in range(len(outer_faces[b]['surfaces'])):                
                lets_remove = [True for to_remove in block_outer_faces_to_remove if to_remove['block_indx'] == b and to_remove['surface_indx'] == s] # Return true if we need to remove
                if not lets_remove: # if lets_remove is false then we want to run this statement to keep 
                    outer_faces_to_keep[b]['surfaces'].append(outer_faces[b]['surfaces'][s])
    return periodic_surfaces, outer_faces_to_keep
                        


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

def __periodicity_check__(face1:Face, face2:Face,blocks:List[Block],rotation_matrix:np.ndarray,face1_index:int, face2_index:int,face1_block_indx:int, face2_block_indx:int):
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
        face1_index (int): Index of face 1 inside the big array of faces. This is added to the list of faces to remove if periodicity is found
        face2_index (int): Index of face 1 inside the big array of faces. This is added to the list of faces to remove if periodicity is found
        face1_block_indx (int): what block index face 1 is located in 
        face2_block_indx (int): what block index face 2 is located in 

    Returns:
        (tuple): containing

            - **periodic_surface** (List[dict]):  
            - **split_surfaces** (List[dict]):  
            - **outer_faces_to_remove** (List[dict]): 
            - **angle** (float):
            - **rotation_matrix** (np.ndarray):

    """
    periodic_surface = None
    split_surfaces = list()
    outer_faces_to_remove = list() 
    if (face1.diagonal_length> face1.diagonal_length):        # switch so that surface 1 is always longer
        temp = deepcopy(face1)
        face1 = deepcopy(face1)
        face2 = temp
        
        temp = face1_indx
        face1_index = face2_indx
        face2_index = deepcopy(temp)

        temp = face2_block_indx
        face1_block_indx = face2_block_indx
        face2_block_indx = deepcopy(temp)
        
    if (face1.diagonal_length> face2.diagonal_length): # Surface 1 is longer than surface 2, shrink surface 1
        # Shrink surface 1 but keep the same k value since it's periodic direction is "k"
        surface1_full = deepcopy(face1)
        surface1 = create_face(blocks[block_indx], surface2.IMIN, surface2.IMAX, 
                                    surface2.JMIN, surface2.JMAX, 
                                    surface1.KMIN, surface1.KMAX)


    angle, rotation_matrix = linear_real_transform(surface1,surface2)   # Do linear real transform to get angle and rotation matrix
    # Check the rotation - rotate the corners to see if this is really equal
    corner1,corner2 = surface2.get_corners()
    corner1_rotated = np.matmul(rotation_matrix,np.array(corner1).transpose()) # y[0], y[1], y[2] are the new x1,y1,z1 rotated about x-axis and they should equal x2,y2,z2 if it's periodic
    corner2_rotated = np.matmul(rotation_matrix,np.array(corner2).transpose())
    corner1_to_match,corner2_to_match = surface1.get_corners()

    corner1_match_error = np.sum(np.square(corner1_rotated-corner1_to_match))
    corner2_match_error = np.sum(np.square(corner2_rotated-corner2_to_match))

    if corner1_match_error<1E-10 and corner2_match_error<1E-10: 
        faces = split_face(surface1,blocks[block_indx],
                            surface1_full.IMIN,surface1_full.JMIN,surface1_full.KMIN,
                            surface1_full.IMAX,surface1_full.JMAX,surface1_full.KMAX)
        
        split_surfaces.extend(faces)              # This list will be appended to the outer_faces
        outer_faces_to_remove.append(surface1_indx) # Schedule these two surfaces to be removed from the outer_faces since they match 
        outer_faces_to_remove.append(surface2_indx)
        
        periodic_surface = {'block1':{'index':block_indx,
            'IMIN':surface1.IMIN,'JMIN':surface1.JMIN,'KMIN':surface1.KMIN, 
            'IMAX':surface1.IMAX,'JMAX':surface1.JMAX,'KMAX':surface1.KMAX, },
            'block2':{'index':block_indx,
            'IMIN':surface2.IMIN,'JMIN':surface2.JMIN,'KMIN':surface2.KMIN, 
            'IMAX':surface2.IMAX,'JMAX':surface2.JMAX,'KMAX':surface2.KMAX, } }

    return periodic_surface,split_surfaces, outer_faces_to_remove, angle, rotation_matrix