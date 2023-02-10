from typing import Dict, List
from .listfunctions import unique_pairs
from .block import Block, reduce_blocks
from .face import Face 
from copy import deepcopy
import math
import numpy as np


def get_outer_faces(block1:Block):
    """Get the outer faces of a block

    Args:
        block1 (Block): A plot3D block

    Returns:
        List[Face]: Non matching faces of the block 
        List[(Face,Face)]: Matching faces inside the block 
    """
    I = [0,block1.IMAX-1]               # Python index starts at 0, need to subtract 1 for it to get the i,j,k
    J = [0,block1.JMAX-1]
    K = [0,block1.KMAX-1]
    # Create the outer faces        
    faces = list()
    face = Face(4)
    i=I[0]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    i=I[1]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    j=J[0]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    j=J[1]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    k=K[0]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    
    faces.append(face)
    face = Face(4)
    k=K[1]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    faces.append(face)

    # Check if faces match each other
    matching = list()
    non_matching = list()
    for i in range(len(faces)):
        matchFound = False
        for j in range(len(faces)):
            if (i!=j and faces[i].vertices_equals(faces[j])):
                matching.append((i,j))
                matchFound = True
        if not matchFound:
            non_matching.append(faces[i]) # these are guaranteed to be exterior 
    matching = list(unique_pairs(matching))
    matching = [(faces[i],faces[j]) for i,j in matching]
    
    # Make sure normals do not intersect 
    # block_center_to_face_center =  block1.cx
    return non_matching, matching # these should be the outer faces
    
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

def find_connected_faces(face_to_search:Face,outer_faces:List[Face],connectivity_matrix:np.ndarray,blocks:List[Block]):
    """Recursive program to return all the matching faces. Note faces must have the same I,J,K definition so faces will be matching if for example: Face1 IMIN=IMAX and Face2 IMIN=IMAX and they share a common edge (2 vertices)

    Args:
        face_to_search (Face): This is the face to search for 
        outer_faces (List[Face]): List of outer faces 
        connectivity_matrix (np.ndarray): block connectivity matrix
        searched_faces (List[Face]): Leave this as blank

    Returns:
        List[Face]: list of all faces that connect with face_to_search and it's neighbors. Beware of duplicates.
    """
    all_matching_faces = [face_to_search]
    faces_to_search = [face_to_search]
    faces_searched = []
    while len(faces_to_search)>0:
        matching_faces = list()
        for face_to_search in faces_to_search:
            n1 = face_to_search.normal(blocks[face_to_search.BlockIndex])
            selected_block_indx = face_to_search.BlockIndex    
            connected_block_indices = np.where(connectivity_matrix[selected_block_indx,:]==1)[0]
            faces_to_check = [o for o in outer_faces if o.BlockIndex in connected_block_indices.tolist()]
            
            for f in faces_to_check:
                if (len(face_to_search.match_indices(f))==2):
                    n2 = f.normal(blocks[f.BlockIndex])
                    angle = abs(math.degrees(math.acos(np.dot(n1,n2)/(np.linalg.norm(n1)*np.linalg.norm(n2)))))
                    if angle>90:
                        angle = 180-angle
                    if angle<50:
                        connectivity_matrix[selected_block_indx, f.BlockIndex] = 0
                        connectivity_matrix[f.BlockIndex, selected_block_indx] = 0
                        matching_faces.append(f)        
              
        faces_searched.extend(faces_to_search)
        matching_faces = [m for m in matching_faces if m not in all_matching_faces]
        faces_to_search.clear()
        faces_to_search.extend(matching_faces)
        all_matching_faces.extend(matching_faces)        
        matching_faces.clear()
    all_matching_faces = list(set(all_matching_faces))  
    return all_matching_faces


def find_closest_block(blocks:List[Block],x:np.ndarray,y:np.ndarray,z:np.ndarray,centroid:np.ndarray,translational_direction:str="x",minvalue:bool=True):
    """Find the closest block to an extreme in the x,y, or z direction and returns the targetting point. 
    Target point is the reference point where we want the closest block and the closest face 

    Args:
        x (np.ndarray): x coordinate of all the blocks' centroid
        y (np.ndarray): y coordinate of all the blocks' centroid
        z (np.ndarray): z coordinate of all the blocks' centroid
        centroid (np.ndarray): centroid (cx,cy,cz)
        translational_direction (str, optional): _description_. Defaults to "x".
        minvalue (bool, optional): _description_. Defaults to True.

    Returns:
        (tuple): containing

            *selected block index* (int): index of closest block
            *target_x* (float): this is the x value where selected block is closest to
            *target_y* (float): this is the y value where selected block is closest to
            *target_z* (float): this is the z value where selected block is closest to
    """
    cx = centroid[0]; cy = centroid[1]; cz = centroid[2]
    target_x = cx; target_y = cy; target_z = cz
    if translational_direction=="x":
        xmins = [b.X.min() for b in blocks]
        xmaxes = [b.X.max() for b in blocks]
        xmin = min(xmins)
        xmax = max(xmaxes)
        dx = xmax - xmin
        if minvalue:
            target_x = xmin - dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))            
        else:
            target_x = xmax + dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))            
        target_y= cy
        target_z = cz
    elif translational_direction=="y":
        ymins = [b.Y.min() for b in blocks]
        ymaxes = [b.Y.max() for b in blocks]
        ymin = min(ymins)
        ymax = max(ymaxes)
        dy = ymax-ymin
        if minvalue:
            target_y = ymin - dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        else:
            target_y = ymax + dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        target_x = cx
        target_z = cz
    else: #  translational_direction=="z":
        zmins = [b.Z.min() for b in blocks]
        zmaxes = [b.Z.max() for b in blocks]
        zmin = min(zmins)
        zmax = max(zmaxes)
        dz = zmax - zmin
        if minvalue:
            target_z = zmin - dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))            
        else:
            target_z = zmax + dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))
        target_x = cx
        target_y = cy

    return selected_block_indx,target_x,target_y,target_z 

def find_bounding_faces(blocks:List[Block],connectivity_matrix:np.ndarray,outer_faces:List[Dict[str,int]], direction:str="z"):
    """Finds bounding faces in x,y,z direction so think of this as the xmin and xmax faces of a block.
    This is useful for translational periodicity where you want the bounds to check

    Args:
        blocks (List[Block]): List of blocks
        outer_faces (List[Dict[str,int]]): outer faces as a list of dictionaries 
        direction (str, optional): Direction to search for. Defaults to "z".

    Returns:
        (tuple) containing:
        
            *lower_connected_faces_export* (List[Dict[str,int]]): Export ready version of lower/left connected faces
            *upper_connected_faces_export* (List[Dict[str,int]]): Export ready version of upper/right connected faces
            *lower_connected_faces* (List[face]): List of lower/left connected faces
            *upper_connected_faces* (List[face]): List of upper/right connected faces
    """
    gcd_array = list()
    # Find the gcd of all the blocks
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)    
    
    xyz_array = np.array([(b.cx, b.cy, b.cz) for b in blocks]); 
    
    cx = np.mean(xyz_array[:,0]) # Centroid 
    cy = np.mean(xyz_array[:,1])
    cz = np.mean(xyz_array[:,2])
    x = xyz_array[:,0]; y = xyz_array[:,1]; z = xyz_array[:,2]

    if len(outer_faces) == 0:
        for i,b in enumerate(blocks):
            outer,_ = get_outer_faces(b)
            [o.set_block_index(i) for o in outer]
            outer_faces.extend(outer)
        outer_faces_all = outer_faces
    else:
        outer_faces_all = outer_face_dict_to_list(blocks,outer_faces,gcd_to_use)

    # Find closest face to the min value 
    selected_block_index,tx,ty,tz = find_closest_block(blocks,x,y,z,np.array([cx,cy,cz]),direction,minvalue=True)
    faces = [f for f in outer_faces_all if f.BlockIndex == selected_block_index]    # Pick all the faces for the selected block 
    min_face_indx = np.argmin(np.array([np.sqrt((f.cx-tx)**2+(f.cy-ty)**2+(f.cz-tz)**2) for f in faces]))
    min_face = faces[min_face_indx]     # Just need this
    
    selected_block_index,tx,ty,tz = find_closest_block(blocks,x,y,z,np.array([cx,cy,cz]),direction,minvalue=False)
    faces = [f for f in outer_faces_all if f.BlockIndex == selected_block_index]    # Pick all the faces for the selected block 
    max_face_indx = np.argmin(np.array([np.sqrt((f.cx-tx)**2+(f.cy-ty)**2+(f.cz-tz)**2) for f in faces]))
    max_face = faces[max_face_indx]     # Also need this

    # Search face connectivity in block connectivity matrix 
    connectivity_matrix = connectivity_matrix - np.eye(connectivity_matrix.shape[0],dtype=np.int8) # turn off connections with itself 

    outer_faces_all = [o for o in outer_faces_all if o.BlockIndex != min_face.BlockIndex]
    outer_faces_all = [o for o in outer_faces_all if o.BlockIndex != max_face.BlockIndex]
    
    # !print("Recursively searching for connected faces")
    lower_connected_faces = find_connected_faces(min_face,outer_faces_all,connectivity_matrix,blocks)
    upper_connected_faces = find_connected_faces(max_face,outer_faces_all,connectivity_matrix,blocks)
    
    lower_connected_faces.append(min_face)
    upper_connected_faces.append(max_face)
    # Unscale faces
    lower_connected_faces = list(set(lower_connected_faces))
    upper_connected_faces = list(set(upper_connected_faces))

    for l in lower_connected_faces:
        l.I*=gcd_to_use
        l.J*=gcd_to_use
        l.K*=gcd_to_use
    for u in upper_connected_faces:
        u.I*=gcd_to_use
        u.J*=gcd_to_use
        u.K*=gcd_to_use

    lower_connected_faces_export = [l.to_dict() for l in lower_connected_faces]
    upper_connected_faces_export = [u.to_dict() for u in upper_connected_faces]
    return lower_connected_faces_export,upper_connected_faces_export, lower_connected_faces, upper_connected_faces

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