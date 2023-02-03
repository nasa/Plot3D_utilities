from typing import Dict, List
from .block import Block
from .face import Face 
from tqdm import trange
import math
import numpy as np

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


#! Possible deletion
def find_connected_face(blocks:List[Block], face:Face, faces:List[Dict[str,int]], look_for_linked:bool=True):
    """Takes a face and a list of faces. Searches for any connected faces. 
        Connections will be checked based on shared verticies. 
        If a face shares at least 2 vertices then it's connected. 

    Args:
        blocks (List[Block]): 
        face (Face): This is the face that you want to find something that matches with it. 
        faces (List[Dict[str,int]]): List of faces not including the face you want to match
        look_for_linked (bool, optional): This takes Face 2 which is connected to face 1 and finds any shared vertices for that face too. Defaults to True.

    Returns:
        Tuple containing: 

            *connected_faces* (List[Dict[str,int]]): List of connected faces in dictionary format
            *faces* (List[Dict[str,int]]): List of faces minus any connected face
    """
    if isinstance(faces[0],dict):
        faces = outer_face_dict_to_list(blocks,faces)
    
    faces = [f for f in faces if f!=face]
    connected_faces = list() # F
    faces_to_search = [face]
    match_found = True
    while (match_found):
        match_found = False
        non_match = list()  
        for face in faces_to_search:
            t = trange(len(faces))
            for i in t:
                t.set_description(f"Matches found {len(connected_faces)}")
                # Look for verticies 2 that match
                ind_x = np.isin(face.x, faces[i].x)
                ind_y = np.isin(face.y, faces[i].y)
                ind_z = np.isin(face.z, faces[i].z)
                if np.sum(ind_x)==2 and np.sum(ind_y)==2 and np.sum(ind_z) == 2 and faces_to_search[0].const_type == faces[i].const_type:
                    connected_faces.append(faces[i])
                    match_found = True
                elif np.sum(ind_x)>2 and np.sum(ind_y)==2 and np.sum(ind_z) == 2 and faces_to_search[0].const_type == faces[i].const_type:
                    connected_faces.append(faces[i])
                    match_found = True
                elif np.sum(ind_x)==2 and np.sum(ind_y)>2 and np.sum(ind_z) == 2 and faces_to_search[0].const_type == faces[i].const_type:
                    connected_faces.append(faces[i])
                    match_found = True
                elif np.sum(ind_x)==2 and np.sum(ind_y)==2 and np.sum(ind_z) > 2 and faces_to_search[0].const_type == faces[i].const_type:
                    connected_faces.append(faces[i])
                    match_found = True
                else:
                    if (faces[i] not in non_match) and (faces[i] not in connected_faces):
                        non_match.append(faces[i])
        if look_for_linked:
            faces_to_search = list(set([c for c in connected_faces if c not in faces_to_search]))
            connected_faces = list(set(connected_faces))
        faces = non_match
        
    if len(connected_faces)>0:
        return [c.to_dict() for c in connected_faces], [f.to_dict() for f in faces]
    else:
        return connected_faces, [f.to_dict() for f in faces] # Returns an empty list and original set of faces to search


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

#! Possible future deletion
def find_face(blocks:List[Block],block_index:int, indices:np.ndarray, outer_faces:List[Face]):
    """Finds a particular face inside a list of outer_faces. This assumes you know the indicies and the block ID. 
    This is important because if you accidentally type in the wrong indicies, it doesn't create a random face. 
    You'll know that this face does not exist. 
    Use this function after you've looked at the geometry in paraview to verify the block id and the indicies of the face. 

    Args:
        blocks (List[Block]): List of blocks
        block_index (int): block index the face belongs to
        indices (np.ndarray): List of integers for [IMIN,JMIN,KMIN,IMAX,JMAX,KMAX]
        outer_faces (List[Face]): _description_

    Returns:
        (Tuple): containing
            *face* (Face | None): Either retuns the matching face or None object.
            *outer_face_index* (int): index of outer face to remove
    """
    outer_face_to_match = None
    for index,o in enumerate(outer_faces):
        if o['block_index'] == block_index:
            a = np.array([o['IMIN'], o['JMIN'], o['KMIN'], o['IMAX'], o['JMAX'], o['KMAX']], dtype=int)
            if np.array_equal(a,indices):
                outer_face_to_match = create_face_from_diagonals(blocks[o['block_index']], imin=o['IMIN'], jmin=o['JMIN'], kmin=o['KMIN'], imax=o['IMAX'], jmax=o['JMAX'], kmax=o['KMAX'])
                outer_face_to_match.set_block_index(block_index)
                break
    return outer_face_to_match, index

def find_face_nearest_point(blocks:List[Block], faces:List[Face], x:float,y:float,z:float):
    """Find a face nearest to a given point

    Args:
        blocks (List[Block]): List of blocks
        faces (List[Face]): List of faces
        x (float): x coordinate of a reference point
        y (float): y coordinate of a reference point
        z (float): z coordinate of a reference point
    """
    n = list(range(len(faces)))
    dv = list()
    for i in n:
        dx = x-faces[i].cx
        dy = y-faces[i].cy
        dz = z-faces[i].cz
        dv.append(math.sqrt(dx*dx + dy*dy + dz*dz))
    face_index = np.argmin(np.array(dv))
    return face_index

#! Possible future deletion
def find_face_near_plane(blocks:List[Block],faces:List[Face],axis:str='y',plane_value:float=0):
    """Returns a list of faces near a given plane. So if your plane is in constant y direction, you specify a value of y that the plane is located at
        
    Args:
        blocks (List[Block]): List of blocks
        faces (List[Face]): List of faces
        axis (str, optional): _description_. Defaults to 'y'.
        plane_value (float, optional): value x,y,z where plane is located at. Defaults to 0.

    Returns:
        (List[Face]): All faces connected to at located nearest to the plane value.
    """
    n = list(range(len(faces)))
    dv = list()
    for i in n:
        if axis=="x":
            dv.append(plane_value-faces[i].cx)
        elif axis=="y":
            dv.append(plane_value-faces[i].cy)
        else: 
            dv.append(plane_value-faces[i].cz)
    face_index = np.argmin(np.array(dv))
    faces_near_plane = find_connected_face(blocks,faces[face_index],faces)
    return faces_near_plane



def outer_face_dict_to_list(blocks:List[Block],outer_faces:List[Dict[str,int]],gcd:int=1) -> List[Face]:
    """Converts a list of dictionary face representations to a list of faces. Use this only for outer faces 

    Args:
        blocks (List[Block]): List of blocks
        outer_faces (List[Dict[str,int]]): List of outer faces represented as a dictionary 
        gcd (int, optional): Greatst common divisor. Defaults to 1.

    Returns:
        List[Face]: List of Face objects 
    """
    outer_faces_all = list()
    for o in outer_faces:
        face = create_face_from_diagonals(blocks[o['block_index']], int(o['IMIN']/gcd), int(o['JMIN']/gcd), 
            int(o['KMIN']/gcd), int(o['IMAX']/gcd), int(o['JMAX']/gcd), int(o['KMAX']/gcd))
        face.set_block_index(o['block_index'])
        outer_faces_all.append(face)

    return outer_faces_all

def match_faces_dict_to_list(blocks:List[Block],matched_faces:List[Dict[str,int]],gcd:int=1):
    """Converts a list of dictionaries representing matched faces to a list of Faces

    Args:
        blocks (List[Block]): List of blocks 
        matched_faces (List[Dict[str,int]]): List of matched faces represented as a dictionary 
        gcd (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    matched_faces_all = list() 
    for _,m in enumerate(matched_faces):
        face1 = create_face_from_diagonals(blocks[m['block1']['block_index']], 
                            int(m['block1']['IMIN']/gcd), int(m['block1']['JMIN']/gcd), int(m['block1']['KMIN']/gcd), 
                            int(m['block1']['IMAX']/gcd), int(m['block1']['JMAX']/gcd), int(m['block1']['KMAX']/gcd))
        face2 = create_face_from_diagonals(blocks[m['block2']['block_index']], 
                            int(m['block2']['IMIN']/gcd), int(m['block2']['JMIN']/gcd), int(m['block2']['KMIN']/gcd), 
                            int(m['block2']['IMAX']/gcd), int(m['block2']['JMAX']/gcd), int(m['block2']['KMAX']/gcd))
        face1.set_block_index(m['block1']['block_index'])
        face2.set_block_index(m['block2']['block_index'])
        matched_faces_all.append(face1)
        matched_faces_all.append(face2)
    return matched_faces_all