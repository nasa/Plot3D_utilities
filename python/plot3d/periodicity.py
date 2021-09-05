from typing import List, Dict, Tuple
from itertools import combinations, product, permutations
import numpy as np
from numpy.core.shape_base import block
from .block import Block, rotate_block
from .face import Face, split_face
from .connectivity import get_face_intersection
from .write import write_plot3D
from math import cos, radians, sin, sqrt, acos
from copy import deepcopy

def create_face(block:Block,imin:int,imax:int,jmin:int,jmax:int,kmin:int,kmax:int) -> Face:
    """Creates a face/surface

    Args:
        block (Block): [description]
        imin (int): [description]
        imax (int): [description]
        jmin (int): [description]
        jmax (int): [description]
        kmin (int): [description]
        kmax (int): [description]

    Returns:
        Face: [description]
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


def find_periodicity(blocks:List[Block],outer_faces:List, periodic_direction:str='k', rotation_axis='x'):
    """This function is used to check for periodicity of the other faces rotated about an axis 
        The way it works is to find faces of a constant i,j, or k value

    Args:
        blocks (List[Block]): [description]
        outer_faces (List): [description]
        periodic_direction (str): either i,j,k to look for
        rotation_axis (str): either x,y,z
    """
    
    periodic_surfaces = list()      # This is the output of the code 
    periodic_angle = 0
    rotation_matrix = None
    outer_faces_to_keep = list()
    # Check periodic within a block 
    block_match = True
    while block_match:
        block_match = False
        for indx in range(len(outer_faces)):
            block = outer_faces[indx]
            outer_faces_to_remove = list()  # Integer list of which outher surfaces to remove
            split_surfaces = list()         # List of split but free surfaces, this will be appended to outer_faces_to_remove list

            blk_indx = block['index']
            surfaces = block['surfaces']   # outer surfaces for a given block 
            surface_combos = list(combinations(range(len(surfaces)),2))
            for s1,s2 in surface_combos:
                surf1_index = s1
                surf2_index = s2            
                # Check if surfaces are periodic with each other
                surface1 = create_face(blocks[blk_indx], surfaces[s1]['IMIN'], surfaces[s1]['IMAX'], 
                                                        surfaces[s1]['JMIN'], surfaces[s1]['JMAX'], 
                                                        surfaces[s1]['KMIN'], surfaces[s1]['KMAX'])

                surface2 = create_face(blocks[blk_indx], surfaces[s2]['IMIN'], surfaces[s2]['IMAX'], 
                                                surfaces[s2]['JMIN'], surfaces[s2]['JMAX'], 
                                                surfaces[s2]['KMIN'], surfaces[s2]['KMAX'])
                if periodic_direction.lower() == "i":
                    if (surfaces[s1]['IMIN'] == surfaces[s1]['IMAX']) and (surfaces[s2]['IMIN'] == surfaces[s2]['IMAX']): 
                        periodic, split_surf, faces_to_remove, angle, rot_mat = __periodicity_within_block__(surface1,surface2,blocks,surf1_index, surf2_index, blk_indx)
                        if periodic:
                            periodic_surfaces.append(periodic)
                            outer_faces_to_remove.extend(faces_to_remove)
                            split_surfaces.extend(split_surf)
                            periodic_angle = angle
                            rotation_matrix = rot_mat
                            block_match = True
                        break

                if periodic_direction.lower() == "i":
                    if (surfaces[s1]['IMIN'] == surfaces[s1]['IMAX']) and (surfaces[s2]['IMIN'] == surfaces[s2]['IMAX']): 
                        periodic, split_surf,faces_to_remove, angle, rot_mat = __periodicity_within_block__(surface1,surface2,blocks,surf1_index, surf2_index, blk_indx)
                        if periodic:
                            periodic_surfaces.append(periodic)
                            outer_faces_to_remove.extend(faces_to_remove)
                            split_surfaces.extend(split_surf)
                            periodic_angle = angle
                            rotation_matrix = rot_mat
                            block_match = True
                        break

                if periodic_direction.lower() == 'k':       # constant k between surfaces 
                    if (surfaces[s1]['KMIN'] == surfaces[s1]['KMAX']) and (surfaces[s2]['KMIN'] == surfaces[s2]['KMAX']): 
                        periodic, split_surf, faces_to_remove, angle, rot_mat = __periodicity_within_block__(surface1,surface2,blocks,surf1_index, surf2_index, blk_indx)
                        if periodic:
                            periodic_surfaces.append(periodic)
                            outer_faces_to_remove.extend(faces_to_remove)
                            split_surfaces.extend(split_surf)
                            periodic_angle = angle
                            rotation_matrix = rot_mat
                            block_match = True
                        break
                    
            if block_match:
                break
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
                    b1_surf = outer_faces[b1]['surfaces'][b1_surf_indx]
                    b2_surf = outer_faces[b2]['surfaces'][b2_surf_indx]

                    if (b1_surf['KMIN'] == b1_surf['KMAX']) and (b2_surf['KMIN'] == b2_surf['KMAX']):
                        surface1 = create_face(blocks[b1], b1_surf['IMIN'], b1_surf['IMAX'], 
                                                            b1_surf['JMIN'], b1_surf['JMAX'], 
                                                            b1_surf['KMIN'], b1_surf['KMAX'])

                        surface2 = create_face(blocks[b2], b2_surf['IMIN'], b2_surf['IMAX'], 
                                                            b2_surf['JMIN'], b2_surf['JMAX'], 
                                                            b2_surf['KMIN'], b2_surf['KMAX'])

                        # Step 1: Try rotating surface 1 and check for intersections 
                        block1_rotated = rotate_block(blocks[b1], rotation_matrix)
                        # write_plot3D('test_rotated.xyz',[block1_rotated],binary=True) # debug purposes
                        df, split_face1, split_face2 = get_face_intersection(surface1,surface2,block1_rotated,blocks[b2],tol=1E-5)
                        if len(df)>0:                          # Try with first block rotated
                            if periodic_direction == "i":
                                df['j1'] += surface1.JMIN
                                df['k1'] += surface1.KMIN
                                df['j2'] += surface2.JMIN
                                df['k2'] += surface2.KMIN

                            if periodic_direction == "j":
                                df['i1'] += surface1.IMIN
                                df['k1'] += surface1.KMIN
                                df['i2'] += surface2.IMIN
                                df['k2'] += surface2.KMIN

                            if periodic_direction == "k":
                                df['i1'] += surface1.IMIN
                                df['j1'] += surface1.JMIN
                                df['i2'] += surface2.IMIN
                                df['j2'] += surface2.JMIN
                        elif len(df)==0:                        # Try with second block rotated
                            # Step 2: if Step 1 does not yield a match then rotate surface 2 and check for intersections 
                            block2_rotated = rotate_block(blocks[b2], rotation_matrix)
                            df, split_face1, split_face2 = get_face_intersection(surface1,surface2,blocks[b1],block2_rotated,tol=1E-5)
                            if periodic_direction == "i":
                                df['j1'] += surface1.JMIN
                                df['k1'] += surface1.KMIN
                                df['j2'] += surface2.JMIN
                                df['k2'] += surface2.KMIN

                            if periodic_direction == "j":
                                df['i1'] += surface1.IMIN
                                df['k1'] += surface1.KMIN
                                df['i2'] += surface2.IMIN
                                df['k2'] += surface2.KMIN

                            if periodic_direction == "k":
                                df['i1'] += surface1.IMIN
                                df['j1'] += surface1.JMIN
                                df['i2'] += surface2.IMIN
                                df['j2'] += surface2.JMIN
                            
                        if len(df)>0:
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
        This function assumes the rotation axis is in the "x" direction 

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

def __periodicity_within_block__(surface1:Face, surface2:Face,blocks:List[Block],surface1_indx:int, surface2_indx:int, block_indx:int):
    """General function to find periodicity within a given block. 

    Steps:
        - 1: Take the surface with the shorter diagonal. Rotate the face about the x-axis. We call this surface-1-rotated. 
        - 2: Shorten surface 2 to match the I and J bounds of surface 1.
        - 3: Do a Linear Real Transform (LRT) Check 

    Args:
        surface1 (Face): [description]
        surface2 (Face): [description]
        blocks (List[Block]): List of blocks
        surface1_indx (int): [description]
        surface2_indx (int): [description]

    Returns:
        [type]: [description]
    """
    periodic_surface = None
    split_surfaces = list()
    outer_faces_to_remove = list() 
    if (surface2.diagonal_length> surface1.diagonal_length):        # switch so that surface 1 is always longer
        temp = deepcopy(surface1)
        surface1 = deepcopy(surface2)
        surface2 = temp
        
        temp = surface1_indx
        surf1_index = surface2_indx
        surf2_index = temp

    if (surface1.diagonal_length> surface2.diagonal_length): # Surface 1 is longer than surface 2, shrink surface 1
        # Shrink surface 1 but keep the same k value since it's periodic direction is "k"
        surface1_full = deepcopy(surface1)
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