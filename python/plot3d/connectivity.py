from .block import Block
from .blockfunctions import reduce_blocks
from .face import Face
from .facefunctions import create_face_from_diagonals, split_face, get_outer_faces
import math 
from itertools import product, combinations
from tqdm import trange
import numpy as np 
import pandas as pd
from typing import List
import math
from .point_match import point_match
from copy import deepcopy


def find_matching_blocks(block1:Block,block2:Block,block1_outer:List[Face], block2_outer:List[Face],tol:float=1E-6):  
    """Takes two blocks and finds all matching pairs

    Args:
        block1 (Block): Any plot3d Block that is not the same as block2
        block2 (Block): Any plot3d Block that is not the same as block1
        block1_outer (List[Face]): outer faces for block 1. 
        block2_outer (List[Face]): Outer faces for block 2
        tol (float, Optional): tolerance to use. Defaults to 1E-6
    
    Note:
        This function was changed to be given an input of outer faces for block 1 and block 2. Outer faces can change and we should use the updated value
    Returns:
        (tuple): containing
            - **df** (pandas.DataFrame): corners of matching pair as block1_corners,block2_corners ([imin,jmin,kmin],[imax,jmax,kmax]), ([imin,jmin,kmin],[imax,jmax,kmax])
            - **block1_outer** (List[Face]):
            - **block2_outer** (List[Face]): 
    """
    # Check to see if outer face of block 1 matches any of the outer faces of block 2
    block_match_indices = list()

    block1_split_faces = list()
    block2_split_faces = list() 
    # Create a dataframe for block1 and block 2 inner matches, add to df later
    # df,split_faces1,split_faces2 = get_face_intersection(block1_outer[3],block2_outer[4],block1,block2,tol=1E-6)

    # Checks the nodes of the outer faces to see if any of them match 
    match = True
    while match:
        match = False
        for p in range(len(block1_outer)):
            block1_face = block1_outer[p]
            for q in range(len(block2_outer)):
                block2_face = block2_outer[q]
                df, split_faces1, split_faces2 = get_face_intersection(block1_face,block2_face,block1,block2,tol)
                if len(df)>0:   # the number of intersection points has to be more than 4
                    # if not block1_face in block1MatchingFace and not block2_face in block2MatchingFace:
                    block_match_indices.append(df)
                    block1_split_faces.extend(split_faces1)
                    block2_split_faces.extend(split_faces2)
                    match = True
                    break
            if match:
                break
        if match:
            block1_outer.pop(p) # type: ignore
            block2_outer.pop(q) # type: ignore
            block1_outer.extend(block1_split_faces)
            block2_outer.extend(block2_split_faces)
            block1_split_faces.clear()
            block2_split_faces.clear()

    return block_match_indices, block1_outer, block2_outer # Remove duplicates using set and list 

def select_multi_dimensional(T:np.ndarray,dim1:tuple,dim2:tuple, dim3:tuple):
    """Takes a block (T) and selects X,Y,Z from the block given a face's dimensions
        theres really no good way to do this in python 
        
    Args:
        T (np.ndarray): arbitrary array so say a full matrix containing X
        dim1 (tuple): 20,50 this selects X in the i direction from i=20 to 50
        dim2 (tuple): 40,60 this selects X in the j direction from j=40 to 60
        dim3 (tuple): 10,20 this selects X in the k direction from k=10 to 20

    Returns:
        np.ndarray: returns X or Y or Z given some range of I,J,K

    """
    if dim1[0] == dim1[1]:
        return T[ dim1[0], dim2[0]:dim2[1]+1, dim3[0]:dim3[1]+1 ]
    if dim2[0] == dim2[1]:
        return T[ dim1[0]:dim1[1]+1, dim2[0], dim3[0]:dim3[1]+1 ]
    if dim3[0] == dim3[1]:
        return T[ dim1[0]:dim1[1]+1, dim2[0]:dim2[1]+1, dim3[0] ]
    
    return T[dim1[0]:dim1[1], dim2[0]:dim2[1], dim3[0]:dim3[1]]

def get_face_intersection(face1:Face,face2:Face,block1:Block,block2:Block,tol:float=1E-6):
    """Get the index of the intersection between two faces located on two different blocks 
        Face1 needs to be the smaller face. 

    Args:
        face1 (Face): An exterior face
        face2 (Face): An exterior face from a different block
        block1 (Block): block containing face1
        block2 (Block): block containing face2
        tol (float): matching tolerance

    Returns:
        (Tuple): containing

            - (pandas.DataFrame): dataframe with matches. Columns = I1, J1, K1, I2, J2, K2
            - (List[Face]): any split faces from block 1
            - (List[Face]): any split faces from block 2 
    """
    
    match_location = list()
    df =pd.DataFrame(columns=['i1','j1','k1','i2','j2','k2'])
    split_faces1 = list()
    split_faces2 = list()
    
    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]
    
    # Grab the points of Face 1 and Face 2
    X1 = select_multi_dimensional(block1.X, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
    Y1 = select_multi_dimensional(block1.Y, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
    Z1 = select_multi_dimensional(block1.Z, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))

    X2 = select_multi_dimensional(block2.X, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
    Y2 = select_multi_dimensional(block2.Y, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
    Z2 = select_multi_dimensional(block2.Z, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))

    # General Search
    if I1[0] == I1[1]: # I is constant in Face 1
        combo = product(range(X1.shape[0]), range(X1.shape[1]))
        for c in combo:
            p, q = c
            x = X1[p,q]
            y = Y1[p,q]
            z = Z1[p,q]
            block2_match_location = point_match(x, y, z, X2, Y2, Z2,tol)
            if sum(block2_match_location)!=-2:
                p2 = int(block2_match_location[0])
                q2 = int(block2_match_location[1])
                # if __edge_match2(df1_edges,df2_edges, p, q, p2, q2):
                if I2[0]==I2[1]:
                    match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})
                if J2[0]==J2[1]:
                    match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                if K2[0]==K2[1]:
                    match_location.append({"i1":I1[0],"j1":p+J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
        df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)

    elif J1[0] == J1[1]: # J is constant in face 1 
        combo = product(range(0,X1.shape[0]), range(0,X1.shape[1]))
        for c in combo:
            p, q = c
            x = X1[p,q]
            y = Y1[p,q]
            z = Z1[p,q]
            block2_match_location = point_match(x, y, z, X2, Y2, Z2,tol)
            if sum(block2_match_location)!=-2:
                p2 = int(block2_match_location[0])
                q2 = int(block2_match_location[1])
                # if __edge_match2(df1_edges,df2_edges, p, q, p2, q2):
                if I2[0]==I2[1]:
                    match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})    # Added an offset because some faces don't start at I=0 or J=0 or K=0 
                if J2[0]==J2[1]:
                    match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                if K2[0]==K2[1]:
                    match_location.append({"i1":p+I1[0],"j1":J1[0],"k1":q+K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
            df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)

    elif K1[0] == K1[1]: # K is constant in face 1 
        combo = product(range(X1.shape[0]), range(X1.shape[1]))
        for c in combo:
            p, q = c
            x = X1[p,q]
            y = Y1[p,q]
            z = Z1[p,q]
            block2_match_location = point_match(x, y, z, X2, Y2, Z2,tol) # pm,qm are the p and q indicies where match occurs 
            if sum(block2_match_location)!=-2:
                p2 = int(block2_match_location[0])
                q2 = int(block2_match_location[1])
                # if __edge_match2(df1_edges,df2_edges, p, q, p2, q2):
                if I2[0]==I2[1]:
                    match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':I2[0],'j2':p2+J2[0],'k2':q2+K2[0]})
                if J2[0]==J2[1]:
                    match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':p2+I2[0],'j2':J2[0],'k2':q2+K2[0]})
                if K2[0]==K2[1]:
                    match_location.append({"i1":p+I1[0],"j1":q+J1[0],"k1":K1[0],'i2':p2+I2[0],'j2':q2+J2[0],'k2':K2[0]})
        df = pd.concat([df, pd.DataFrame(match_location)], ignore_index=True)
    
    # Checking for split faces 
    if len(df)>=4:
        if (__check_edge(df)):
            df = pd.DataFrame()     # If it's an edge
        else:                       # not edge 
            # Filter match increasing - This keeps uniqueness
            if I1[0]==I1[1]:
                df = __filter_block_increasing(df,'j1')
                df = __filter_block_increasing(df,'k1')
            elif J1[0]==J1[1]:
                df = __filter_block_increasing(df,'i1')
                df = __filter_block_increasing(df,'k1')
            elif K1[0]==K1[1]:
                df = __filter_block_increasing(df,'i1')
                df = __filter_block_increasing(df,'j1')
            
            if I2[0]==I2[1]:
                df = __filter_block_increasing(df,'j2')
                df = __filter_block_increasing(df,'k2')
            elif J2[0]==J2[1]:
                df = __filter_block_increasing(df,'i2')
                df = __filter_block_increasing(df,'k2')
            elif K2[0]==K2[1]:
                df = __filter_block_increasing(df,'i2')
                df = __filter_block_increasing(df,'j2')

            # Do a final check after doing all these checks
            if len(df)>=4:       # Greater than 4 because match can occur with simply 4 corners but the interior doesn't match. 
                # Check for Split faces
                ## Block 1
                main_face = create_face_from_diagonals(block1,imin=I1[0],imax=I1[1], jmin=J1[0],jmax=J1[1],kmin=K1[0],kmax=K1[1])
                imin, jmin, kmin = df['i1'].min(), df['j1'].min(), df['k1'].min()
                imax, jmax, kmax = df['i1'].max(), df['j1'].max(), df['k1'].max()
                if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                    split_faces1 = split_face(main_face,block1,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
                    [s.set_block_index(face1.blockIndex) for s in split_faces1]
                    [s.set_face_id(face1.id) for s in split_faces1]

                ## Block 2
                main_face = create_face_from_diagonals(block2,imin=I2[0],imax=I2[1], jmin=J2[0],jmax=J2[1],kmin=K2[0],kmax=K2[1])
                imin, jmin, kmin = df['i2'].min(), df['j2'].min(), df['k2'].min()
                imax, jmax, kmax = df['i2'].max(), df['j2'].max(), df['k2'].max()
                if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                    split_faces2 = split_face(main_face,block2,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
                    [s.set_block_index(face2.blockIndex) for s in split_faces2]
                    [s.set_face_id(face2.id) for s in split_faces2]

    else:
        df = pd.DataFrame() # set df to empty dataframe
    return df, split_faces1, split_faces2

def __filter_block_increasing(df:pd.DataFrame,key1:str):
    """Filters dataframe results of get_face_intersection to make sure both key1 is increasing. 
        When searching through a plot3D we check based on the planes e.g. const i, j, or k 
        values will be removed if they are not

    Args:
        df (pd.DataFrame): DataFrame containing matching points 
        key1 (str): column that you want to be in increasing order 

    Returns:
        pd.DataFrame: sorted dataframe
    """        

    '''
        Sometimes there's a match on 2 edges and we do not want to keep that 
            | face1 | face2 | face1  | 
        Above shows face 2 touching face 1 at 2 edges. this is not a match. 
    '''
    if len(df)==0:
        return df
    
    key1_vals = list(df[key1].unique()) # get the unique values 
    key1_vals.sort()
    key1_vals_to_use = list()

    if len(key1_vals)<=1:
        return pd.DataFrame() # Returning an empty dataframe. This solves the condition where you have edge matching 

    for i in range(len(key1_vals)-1):
        if (key1_vals[i+1] - key1_vals[i])==1: # Remove
            key1_vals_to_use.append(key1_vals[i])
    # Look backwards 
    if (key1_vals[-1] - key1_vals[-2])==1: # Remove
            key1_vals_to_use.append(key1_vals[-1])
    df = df[df[key1].isin(key1_vals_to_use)]        
    return df

def __check_edge(df:pd.DataFrame):
    """ Check if the results of get_face_intersection is an edge instead of a face.  
        if it's an edge then both intersecting blocks are connected by an edge on both blocks 

    Args:
        df (pd.DataFrame): dataframe containing columns i1, j1, k1, i2, j2, k2

    Returns:
        boolean: True = It is an edge, False = not edge 
    """
    face1_diagonal = [(df['i1'].min(),df['j1'].min(),df['k1'].min()),(df['i1'].max(),df['j1'].max(),df['k1'].max()) ]
    face2_diagonal = [(df['i2'].min(),df['j2'].min(),df['k2'].min()), (df['i2'].max(),df['j2'].max(),df['k2'].max())]
    edge1 = face1_diagonal[0]
    edge2 = face1_diagonal[1]
    edge_matches = 0 
    for i in range(3):
        if edge1[i]==edge2[i]:
            edge_matches+=1
    if edge_matches<2:
        return False
    else:
        return True
    
def combinations_of_nearest_blocks(blocks:List[Block],nearest_nblocks:int=4):
    """Returns the indices of the nearest 6 blocks based on their centroid

    Args:
        block (Block): block you are interested in
        blocks (List[Block]): list of all your blocks

    Returns:
        List[Tuple[int,int]]: combinations of nearest blocks 
    """
    # Pick a block get centroid of all outer faces        
    centroids = np.array([(b.cx,b.cy,b.cz) for b in blocks])
    distance_matrix = np.zeros((centroids.shape[0],centroids.shape[0]))+10000
    # Build a matrix 
    for i in range(centroids.shape[0]):
        for j in range(centroids.shape[0]):
            if i!=j:
                dx = centroids[i,0]-centroids[j,0]
                dy = centroids[i,1]-centroids[j,1]
                dz = centroids[i,2]-centroids[j,2]
                distance_matrix[i,j] = np.sqrt(dx*dx+dy*dy+dz*dz)
    
    # Now that we have this matrix, we sort the distances by rows and pick the closest 8 blocks, can use 4 but 8 might be safer
    new_combos = list()
    for i in range(len(blocks)): # For block i
        indices = np.argsort(distance_matrix[i,:])
        for j in indices[:nearest_nblocks]:
            if distance_matrix[i,j] < 10000:
                new_combos.append((i,j))
    return new_combos
    
def connectivity_fast(blocks:List[Block]):
    """Reduces the size of the blocks by a factor of the minimum gcd. This speeds up finding the connectivity 

    Args:
        blocks (List[Block]): Lists of blocks you want to find the connectivity for

    Returns:
        (List[Dict]): All matching faces formatted as a list of { 'block1': {'block_index', 'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX'} }
        (List[Dict]): All exterior surfaces formatted as a list of { 'block_index', 'surfaces': [{'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX', 'ID'}] }
        
    """
    gcd_array = list()
    # Find the gcd of all the blocks 
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))        
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    print(f"gcd to use {gcd_to_use}")
    new_blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)

    # Find Connectivity 
    face_matches, outer_faces_formatted = connectivity(new_blocks)
    # scale it up
    for i in range(len(face_matches)):
        face_matches[i]['block1']['IMIN'] *= gcd_to_use
        face_matches[i]['block1']['JMIN'] *= gcd_to_use
        face_matches[i]['block1']['KMIN'] *= gcd_to_use
        face_matches[i]['block1']['IMAX'] *= gcd_to_use
        face_matches[i]['block1']['JMAX'] *= gcd_to_use
        face_matches[i]['block1']['KMAX'] *= gcd_to_use

        face_matches[i]['block2']['IMIN'] *= gcd_to_use
        face_matches[i]['block2']['JMIN'] *= gcd_to_use
        face_matches[i]['block2']['KMIN'] *= gcd_to_use
        face_matches[i]['block2']['IMAX'] *= gcd_to_use
        face_matches[i]['block2']['JMAX'] *= gcd_to_use
        face_matches[i]['block2']['KMAX'] *= gcd_to_use
    for j in range(len(outer_faces_formatted)):
        outer_faces_formatted[j]['IMIN'] *= gcd_to_use
        outer_faces_formatted[j]['JMIN'] *= gcd_to_use
        outer_faces_formatted[j]['KMIN'] *= gcd_to_use
        outer_faces_formatted[j]['IMAX'] *= gcd_to_use
        outer_faces_formatted[j]['JMAX'] *= gcd_to_use
        outer_faces_formatted[j]['KMAX'] *= gcd_to_use
    return face_matches, outer_faces_formatted

def connectivity(blocks:List[Block]):
    """Returns a dictionary outlining the connectivity of the blocks along with any exterior surfaces 

    Args:
        blocks (List[Block]): List of all blocks in multi-block plot3d mesh

    Returns:
        (List[Dict]): All matching faces formatted as a list of { 'block1': {'block_index', 'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX'} }
        (List[Dict]): All exterior surfaces formatted as a list of { 'block_index', 'surfaces': [{'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX', 'ID'}] }

    """

    outer_faces = list()      
    face_matches = list()
    matches_to_remove = list()
    temp = [get_outer_faces(b) for b in blocks]
    block_outer_faces = [t[0] for t in temp]
    combos = combinations_of_nearest_blocks(blocks,6) # Find the 6 nearest Blocks and search through all that. 
    # df_matches, blocki_outerfaces, blockj_outerfaces = find_matching_blocks(blocks[0],blocks[1],[block_outer_faces[0][0],block_outer_faces[0][1]],[block_outer_faces[1][0],block_outer_faces[1][1]],1E-12)    # This function finds partial matches between blocks

    t = trange(len(combos))    
    for indx in t:     # block i        
        i,j = combos[indx]
        t.set_description(f"Checking connections block {i} with {j}")
        # Takes 2 blocks, gets the matching faces exterior faces of both blocks 
        df_matches, blocki_outerfaces, blockj_outerfaces = find_matching_blocks(blocks[i],blocks[j],block_outer_faces[i],block_outer_faces[j])    # This function finds partial matches between blocks
        [o.set_block_index(i) for o in blocki_outerfaces]
        [o.set_block_index(j) for o in blockj_outerfaces]
        block_outer_faces[i] = blocki_outerfaces
        block_outer_faces[j] = blockj_outerfaces
        # Update connectivity for blocks with matching faces 
        if (len(df_matches)>0):
            for df in df_matches:
                matches_to_remove.append(create_face_from_diagonals(block=blocks[i],imin=df['i1'].min(),jmin=df['j1'].min(),kmin=df['k1'].min(),
                                                                                    imax=df['i1'].max(),jmax=df['j1'].max(),kmax=df['k1'].max()))
                matches_to_remove[-1].set_block_index(i)

                matches_to_remove.append(create_face_from_diagonals(block=blocks[j],imin=df['i2'].min(),jmin=df['j2'].min(),kmin=df['k2'].min(),
                                                                                    imax=df['i2'].max(),jmax=df['j2'].max(),kmax=df['k2'].max()))
                matches_to_remove[-1].set_block_index(j)
                
                face1 = matches_to_remove[-2]
                face2 = matches_to_remove[-1]

                temp = face_matches_to_dict(face1,face2,blocks[i],blocks[j])
                temp['match'] = df
                face_matches.append(temp)

    # Update Outer Faces 
    [outer_faces.extend(o) for o in block_outer_faces] # all the outer faces 
    outer_faces = list(set(outer_faces))    # Get most unique
        
    outer_faces = [o for o in outer_faces if o not in matches_to_remove]
    # Remove any outer faces that may have been found by mistake
    # Check I,J,K if J and K are the same with another outer face, select the face with shorter I 
    outer_faces_to_remove = list() 
    for i in range(len(blocks)):
        block_outerfaces = [o for o in outer_faces if o.BlockIndex == i]
        for o in block_outerfaces:
            IJK = np.array([o.IMIN,o.JMIN,o.KMIN,o.IMAX,o.JMAX,o.KMAX])
            for o2 in block_outerfaces:
                IJK2 = np.array([o2.IMIN,o2.JMIN,o2.KMIN,o2.IMAX,o2.JMAX,o2.KMAX])
                if sum((IJK-IJK2)==0) == 5: # [0,0,0,40,100,0] (outer) [0,0,0,56,100,0] (outer) -> remove the longer face
                    if (o2.diagonal_length>o.diagonal_length):
                        outer_faces_to_remove.append(o2)
                    else:
                        outer_faces_to_remove.append(o)

    outer_faces = [o for o in outer_faces if o not in outer_faces_to_remove]


    # Find self-matches: Do any faces of, for example, block1 match another face in block 1
    for i in range(len(blocks)):
        _,self_matches = get_outer_faces(blocks[i]) 
        for match in self_matches: # Append to face matches 
            face_matches.append({'block1':{
                                            'block_index':i,'IMIN':match[0].I.min(),'JMIN':match[0].J.min(),'KMIN':match[0].K.min(),
                                            'IMAX':match[0].I.max(),'JMAX':match[0].J.max(),'KMAX':match[0].K.max()
                                        },
                                    'block2':{
                                            'block_index':i,'IMIN':match[1].I.min(),'JMIN':match[1].J.min(),'KMIN':match[1].K.min(),
                                            'IMAX':match[1].I.max(),'JMAX':match[1].J.max(),'KMAX':match[1].K.max()
                                        },
                                    'match':pd.DataFrame([{
                                            'block_index':i,'IMIN':match[0].I.min(),'JMIN':match[0].J.min(),'KMIN':match[0].K.min(),
                                            'IMAX':match[0].I.max(),'JMAX':match[0].J.max(),'KMAX':match[0].K.max()
                                        },{
                                            'block_index':i,'IMIN':match[1].I.min(),'JMIN':match[1].J.min(),'KMIN':match[1].K.min(),
                                            'IMAX':match[1].I.max(),'JMAX':match[1].J.max(),'KMAX':match[1].K.max()
                                        }])
                                    })

    # Update the outer faces
    outer_faces_formatted = list() # This will contain 
    id = 1 
    for face in outer_faces:        
        outer_faces_formatted.append({ 'IMIN':min(face.I), 'JMIN':min(face.J), 'KMIN':min(face.K),
                            'IMAX':max(face.I), 'JMAX':max(face.J), 'KMAX':max(face.K),
                            'id':id, 'block_index':face.BlockIndex })
        id += 1

    return face_matches, outer_faces_formatted  

def face_matches_to_dict(face1:Face, face2:Face,block1:Block,block2:Block):
    """Makes sure the diagonal of face 1 match the diagonal of face 2

    Args:
        face1 (Face): Face 1 with block index 
        face2 (Face): Face 2 with block index
        block1 (Block): Block 1 
        block2 (Block): Block 2 
    
    Returns:
        (dict): dictionary describing the corner matches 

    """
    match = {
            'block1':{
                            'block_index':face1.BlockIndex,
                            'IMIN':-1,'JMIN':-1,'KMIN':-1,  # Lower Corner
                            'IMAX':-1,'JMAX':-1,'KMAX':-1,   # Upper Corner
                            'id':face1.id
                        },
                'block2':{
                            'block_index':face2.BlockIndex,
                            'IMIN':-1,'JMIN':-1,'KMIN':-1,  # Lower Corner
                            'IMAX':-1,'JMAX':-1,'KMAX':-1,   # Upper Corner
                            'id':face2.id
                        }
                }
            
    I1 = [face1.IMIN,face1.IMAX]
    J1 = [face1.JMIN,face1.JMAX]
    K1 = [face1.KMIN,face1.KMAX]

    I2 = [face2.IMIN,face2.IMAX]
    J2 = [face2.JMIN,face2.JMAX]
    K2 = [face2.KMIN,face2.KMAX]

    # Search for corners     
    x1_l = block1.X[I1[0],J1[0],K1[0]]    # lower corner of block 1 
    y1_l = block1.Y[I1[0],J1[0],K1[0]]
    z1_l = block1.Z[I1[0],J1[0],K1[0]]
    # Matches which corner in block 2 
    search_results = list()
    for p in I2: 
        for q in J2:
            for r in K2:
                x2 = block2.X[p,q,r]  
                y2 = block2.Y[p,q,r]
                z2 = block2.Z[p,q,r]
                dx = x2-x1_l; dy = y2-y1_l; dz = z2 -z1_l                
                search_results.append({'I':p,'J':q,'K':r,'d':math.sqrt(dx*dx + dy*dy + dz*dz)})    
    df = pd.DataFrame(search_results)
    df = df.sort_values(by=['d'])
    match['block1']['IMIN'] = face1.IMIN
    match['block1']['JMIN'] = face1.JMIN
    match['block1']['KMIN'] = face1.KMIN
    match['block2']['IMIN'] = int(df.iloc[0]['I'])
    match['block2']['JMIN'] = int(df.iloc[0]['J'])
    match['block2']['KMIN'] = int(df.iloc[0]['K'])
    
    # Search for corners     
    x1_u = block1.X[I1[1],J1[1],K1[1]]    # lower corner of block 1 
    y1_u = block1.Y[I1[1],J1[1],K1[1]]
    z1_u = block1.Z[I1[1],J1[1],K1[1]]
    # Matches which corner in block 2 
    search_results = list()
    for p in I2: 
        for q in J2:
            for r in K2:
                x2 = block2.X[p,q,r]  
                y2 = block2.Y[p,q,r]
                z2 = block2.Z[p,q,r]
                dx = x2-x1_u; dy = y2-y1_u; dz = z2 -z1_u
                search_results.append({'I':p,'J':q,'K':r,'d':math.sqrt(dx*dx + dy*dy + dz*dz)})
    df = pd.DataFrame(search_results)
    df = df.sort_values(by=['d'])
    match['block1']['IMAX'] = face1.IMAX
    match['block1']['JMAX'] = face1.JMAX
    match['block1']['KMAX'] = face1.KMAX
    match['block2']['IMAX'] = int(df.iloc[0]['I'])
    match['block2']['JMAX'] = int(df.iloc[0]['J'])
    match['block2']['KMAX'] = int(df.iloc[0]['K'])
    return match
