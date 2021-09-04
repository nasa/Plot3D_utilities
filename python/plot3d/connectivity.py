from pandas.core.algorithms import unique
from .block import Block
from .face import Face, create_face_from_diagonals, split_face
import math 
from itertools import product, combinations
from tqdm import trange
import numpy as np 
from .differencing import find_face_edges
import pandas as pd
from operator import eq 
import tqdm
from typing import List
import math
from .point_match import point_match

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
            if (i!=j and faces[i] == faces[j]):
                matching.append((i,j))
                matchFound = True
        if not matchFound:
            non_matching.append(faces[i]) # these are guaranteed to be exterior 
    matching = list(unique_pairs(matching))
    matching = [(faces[i],faces[j]) for i,j in matching]
    
    # Make sure normals do not intersect 
    # block_center_to_face_center =  block1.cx
    return non_matching, matching # these should be the outer faces

def unique_pairs(listOfItems:list):
    """Checks if an item is not already in the list 
    
    Args:
        listOfItems (list): list of combinations e.g. (1,2),(3,4),(2,1)

    Yields:
        unique pair: [description]
    """
    seen = set()  #use set to keep track of already seen items, sets provide O(1) lookup  
    for x,y in listOfItems:
        if x!=y and (y,x) not in seen:
            seen.add((x,y)) 
            yield x,y


def find_matching_blocks(block1:Block,block2:Block, full_face_match=False):  
    """Takes two blocks and finds all matching pairs

    Args:
        block1 (Block): Any plot3d Block that is not the same as block2
        block2 (Block): Any plot3d Block that is not the same as block1
        full_face_match (bool): use full face matching (Much faster)

    Returns:
        pandas.DataFrame: corners of matching pair as block1_corners,block2_corners ([imin,jmin,kmin],[imax,jmax,kmax]), ([imin,jmin,kmin],[imax,jmax,kmax])
    """
    # Check to see if outer face of block 1 matches any of the outer faces of block 2
    block_match_indices = list()

    block1_outer,_ = get_outer_faces(block1)
    block2_outer,_ = get_outer_faces(block2)
    block1_split_faces = list()
    block2_split_faces = list() 
    # Create a dataframe for block1 and block 2 inner matches, add to df later
    # df,split_faces1,split_faces2 = get_face_intersection(block1_outer[3],block2_outer[4],block1,block2,tol=1E-6)

    # Checks the nodes of the outer faces to see if any of them match 
    block1MatchingFace = list()
    block2MatchingFace = list()
    match = True
    while match:
        match = False
        for p in range(len(block1_outer)):
            block1_face = block1_outer[p]
            for q in range(len(block2_outer)):
                block2_face = block2_outer[q]
                df,split_faces1,split_faces2 = get_face_intersection(block1_face,block2_face,block1,block2,full_face_match)                
                if len(df)>2:
                    block_match_indices.append(df)
                    block1MatchingFace.append(block1_face)
                    block2MatchingFace.append(block2_face)
                    block1_split_faces.extend(split_faces1)
                    block2_split_faces.extend(split_faces2)
                    match = True

        [block1_outer.remove(b) for b in block1MatchingFace if b in block1_outer]
        [block2_outer.remove(b) for b in block2MatchingFace if b in block2_outer]
        if len(block1_split_faces)>0:
            block1_outer.extend(block1_split_faces)
            block1_split_faces.clear()
        if len(block2_split_faces)>0:
            block2_outer.extend(block2_split_faces)
            block2_split_faces.clear()
    
    return block_match_indices, block1_outer, block2_outer # Remove duplicates using set and list 



def get_face_intersection(face1:Face,face2:Face,block1:Block,block2:Block,full_face_match=False,tol:float=1E-6):
    """Get the index of the intersection between two faces located on two different blocks 

    Args:
        face1 (Face): An exterior face
        face2 (Face): An exterior face from a different block
        block1 (Block): block containing face1
        block2 (Block): block containing face2
        full_face_match (bool): use full face matching (Much faster)

    Returns:
        (Tuple): containing

            - (pandas.DataFrame): dataframe with matches. Columns = I1, J1, K1, I2, J2, K2
            - (List[Face]): any split faces from block 1
            - (List[Face]): any split faces from block 2 
    """

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
    
    match_location = list()
    df =pd.DataFrame(columns=['i1','j1','k1','i2','j2','k2'])
    split_faces1 = list()
    split_faces2 = list()
    
    matchedIndicies = face1.match_indices(face2) # face1 == face2, Full face match is automatically checked first. If this is not 4 then code proceeds to partial match
    if len(matchedIndicies)==4: # Checks to see if the two faces corners are actually equal
        for i,j in matchedIndicies:
            i1,j1,k1 = face1.I[i],face1.J[i],face1.K[i]
            i2,j2,k2 = face2.I[j],face2.J[j],face2.K[j]
            match_location.append({"i1":i1,"j1":j1,"k1":k1,'i2':i2,'j2':j2,'k2':k2})
        df = df.append(match_location,ignore_index=True)
         
    elif(not full_face_match): # Check all points interior of the block
        I1 = [face1.IMIN,face1.IMAX]
        J1 = [face1.JMIN,face1.JMAX]
        K1 = [face1.KMIN,face1.KMAX]

        I2 = [face2.IMIN,face2.IMAX]
        J2 = [face2.JMIN,face2.JMAX]
        K2 = [face2.KMIN,face2.KMAX]

        X1 = select_multi_dimensional(block1.X, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
        Y1 = select_multi_dimensional(block1.Y, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))
        Z1 = select_multi_dimensional(block1.Z, (I1[0],I1[1]),(J1[0],J1[1]),(K1[0],K1[1]))

        X2 = select_multi_dimensional(block2.X, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
        Y2 = select_multi_dimensional(block2.Y, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))
        Z2 = select_multi_dimensional(block2.Z, (I2[0],I2[1]),(J2[0],J2[1]),(K2[0],K2[1]))

        # General Search

        if I1[0] == I1[1]: # I is constant in Face 1
            i = I1[0]
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
                        match_location.append({"i1":i,"j1":p,"k1":q,'i2':I2[0],'j2':p2,'k2':q2})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":i,"j1":p,"k1":q,'i2':p2,'j2':J2[0],'k2':q2})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":i,"j1":p,"k1":q,'i2':p2,'j2':q2,'k2':K2[0]})
            df = df.append(match_location,ignore_index=True)

        elif J1[0] == J1[1]: # J is constant in face 1 
            j = J1[0]
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
                        match_location.append({"i1":p,"j1":J1[0],"k1":q,'i2':I2[0],'j2':p2,'k2':q2})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":p,"j1":J1[0],"k1":q,'i2':p2,'j2':J2[0],'k2':q2})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":p,"j1":J1[0],"k1":q,'i2':p2,'j2':q2,'k2':K2[0]})
                df = df.append(match_location,ignore_index=True)

        elif K1[0] == K1[1]: # K is constant in face 1 
            k = K1[0]
            combo = product(range(0,X1.shape[0]), range(0,X1.shape[1]))
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
                        match_location.append({"i1":p,"j1":q,"k1":K1[0],'i2':I2[0],'j2':p2,'k2':q2})
                    if J2[0]==J2[1]:
                        match_location.append({"i1":p,"j1":q,"k1":K1[0],'i2':p2,'j2':J2[0],'k2':q2})
                    if K2[0]==K2[1]:
                        match_location.append({"i1":p,"j1":q,"k1":K1[0],'i2':p2,'j2':q2,'k2':K2[0]})
            df = df.append(match_location,ignore_index=True)
        
        # Checking for split faces 
        if len(df)>0:
            if (len(df)==2 or __check_edge(df)):
                df = pd.DataFrame()
            else:
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

                if len(df)>0:
                    # Check for Split faces
                    ## Block 1
                    main_face = create_face_from_diagonals(block1,imin=I1[0],imax=I1[1], jmin=J1[0],jmax=J1[1],kmin=K1[0],kmax=K1[1])
                    imin, jmin, kmin = df['i1'].min(), df['j1'].min(), df['k1'].min()
                    imax, jmax, kmax = df['i1'].max(), df['j1'].max(), df['k1'].max()
                    if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                        split_faces1 = split_face(main_face,block1,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)

                    ## Block 2
                    main_face = create_face_from_diagonals(block2,imin=I2[0],imax=I2[1], jmin=J2[0],jmax=J2[1],kmin=K2[0],kmax=K2[1])
                    imin, jmin, kmin = df['i2'].min(), df['j2'].min(), df['k2'].min()
                    imax, jmax, kmax = df['i2'].max(), df['j2'].max(), df['k2'].max()
                    if int(imin==imax) + int(jmin==jmax) + int(kmin==kmax)==1:
                        split_faces2 = split_face(main_face,block2,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)

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

    for i in range(len(key1_vals)-1):
        if (key1_vals[i+1] - key1_vals[i])==1: # Remove
            key1_vals_to_use.append(key1_vals[i])
    # Look backwards 
    if (key1_vals[-1] - key1_vals[-2])==1: # Remove
            key1_vals_to_use.append(key1_vals[-1])
    df = df[df[key1].isin(key1_vals_to_use)]        
    return df

def __check_edge(df:pd.DataFrame):
    """Check if the results of get_face_intersection is an edge instead of a face
        if it's an edge then both intersecting blocks are connected by an edge on both blocks 

    Args:
        df (pd.DataFrame): dataframe containing columns i1, j1, k1, i2, j2, k2

    Returns:
        boolean: True = It is an edge, False = not edge 
    """
    
    block = [['i1','j1','k1'],['i2','j2','k2']]
    c=0
    for b in block:        
        for k in b:
            c += int(len(df[k].unique())>1) # example i1 = 1,2,3,4,5
    return c<4 #  If "c" is less than 4 then it's an edge 
    

def connectivity(blocks:List[Block],full_face_match=False):
    """Returns a dictionary outlining the connectivity of the blocks along with any exterior surfaces 

    Args:
        blocks (List[Block]): [description]
        full_face_match (bool): assume the entire face of a block is a perfect match on the opposite block (much faster)

    Returns:
        [list]: All matching faces formatted as a list of { 'block1': {'index', 'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX'} }
        [list]: All exterior surfaces formatted as a list of { 'block_index', 'surfaces': [{'IMIN', 'JMIN','KMIN', 'IMAX','JMAX','KMAX', 'ID'}] }
    """

    def combinations_of_nearest_blocks(blocks:List[Block],nearest_nblocks:int=12):
        """Returns the indices of the nearest 6 blocks based on their centroid

        Args:
            block (Block): block you are interested in
            blocks (List[Block]): list of all your blocks

        Returns:
            [type]: [description]
        """
        combos = list(combinations(range(len(blocks)),2))
        if (len(blocks))<=48:
            return combos
        
        distances = [None] * len(combos)
        for indx,(i,j) in enumerate(combos):
            dx = blocks[i].cx - blocks[j].cx
            dy = blocks[i].cy - blocks[j].cy
            dz = blocks[i].cz - blocks[j].cz
            distances[indx] = math.sqrt(dx*dx + dy*dy + dz*dz)

        df = pd.DataFrame(combos, columns =['block_i', 'block_j'])
        df['distance'] = distances 

        block_i = list(df.block_i.unique())        
        main_df = df.copy(deep=True)
        block_df = list()
        for i in block_i:
            df2 = main_df[main_df['block_i']==i].sort_values('distance')
            if len(df2)>nearest_nblocks:
                d = df2.iloc[nearest_nblocks-1].distance
                block_df.append(main_df[(main_df['block_i']==i) & (main_df['distance'] <= d)])
            else:
                block_df.append(main_df[(main_df['block_i']==i)])
        df = pd.concat(block_df)
        subset = df[['block_i', 'block_j']]
        return [tuple(x) for x in subset.to_numpy()] 
        


    def find_block_index_in_outer_faces(outer_faces:List[Face], block_indx:int):
        """Looks through an array of faces (outer_faces). Outer faces can include multiple blocks
            This function finds the face of corresponding block index

        Args:
            outer_faces ([type]): List of faces
            block_indx (int): block index

        Returns:
            [type]: [description]
        """
        for i in range(len(outer_faces)):
            if outer_faces[i]['block'] == block_indx:
                return i
        return -1 
        
    def update_outer_faces(prev_outer_faces:list, new_outer_faces:list):
        """Takes a previous set of outer faces and a new set of outer faces. 
            If a face id is present in the "previous" but not in the "new" remove that face 

        Args:
            prev_outer_faces (list): previous set of non matching outer faces of a given block 
            new_outer_faces (list): new set of non-matching outer faces of the same block
        """                
        face_intersection = list(set(prev_outer_faces) & set(new_outer_faces))             # Find the intersection 
        return face_intersection                                     
            

    outer_faces = list()      
    face_matches = list()
    # Find the 6 nearest Blocks and search through all that.     
    # combos = list(combinations(range(len(blocks)),2))

    combos = combinations_of_nearest_blocks(blocks)
    t = trange(len(combos))
    for indx in t:     # block i        
        i,j = combos[indx]
        t.set_description(f"Checking connections block {i} with {j}")
        # Takes 2 blocks, gets the matching faces exterior faces of both blocks 
        df_matches, blocki_outerfaces, blockj_outerfaces = find_matching_blocks(blocks[i],blocks[j],full_face_match)
        # Update connectivity for blocks with matching faces 
        if (len(df_matches)>0):
            for df in df_matches:
                face_matches.append({'block1':{
                                            'index':i,'IMIN':df['i1'].min(),'JMIN':df['j1'].min(),'KMIN':df['k1'].min(),
                                            'IMAX':df['i1'].max(),'JMAX':df['j1'].max(),'KMAX':df['k1'].max()
                                        },
                                    'block2':{
                                            'index':j,'IMIN':df['i2'].min(),'JMIN':df['j2'].min(),'KMIN':df['k2'].min(),
                                            'IMAX':df['i2'].max(),'JMAX':df['j2'].max(),'KMAX':df['k2'].max()
                                        },
                                    'match':df})

        # Update Outer Faces 
        ## 01 Find block i and block j
        block_iloc = find_block_index_in_outer_faces(outer_faces,i)
        block_jloc = find_block_index_in_outer_faces(outer_faces,j)
        
        # Update or Append to list of outer faces 
        if block_iloc==-1:      # outer face does not exist 
            outer_faces.append({"block":i, "outer_faces":blocki_outerfaces})                            
        else:                   # Block i is in outer faces 
            outer_faces[block_iloc]['outer_faces'] = update_outer_faces(outer_faces[block_iloc]['outer_faces'], blocki_outerfaces)
        
        if block_jloc==-1:
            outer_faces.append({"block":j, "outer_faces":blockj_outerfaces})
        else:                   # Block j is in outer faces 
            outer_faces[block_jloc]['outer_faces'] = update_outer_faces(outer_faces[block_jloc]['outer_faces'], blockj_outerfaces)
        
    # Find self-matches, do any faces of for example block1 match another face in block 1
    for i in range(len(blocks)):
        _,self_matches = get_outer_faces(blocks[i]) 
        for match in self_matches: # Append to face matches 
            face_matches.append({'block1':{
                                            'index':i,'IMIN':match[0].I.min(),'JMIN':match[0].J.min(),'KMIN':match[0].K.min(),
                                            'IMAX':match[0].I.max(),'JMAX':match[0].J.max(),'KMAX':match[0].K.max()
                                        },
                                    'block2':{
                                            'index':i,'IMIN':match[1].I.min(),'JMIN':match[1].J.min(),'KMIN':match[1].K.min(),
                                            'IMAX':match[1].I.max(),'JMAX':match[1].J.max(),'KMAX':match[1].K.max()
                                        },
                                    'match':pd.DataFrame([{
                                            'index':i,'IMIN':match[0].I.min(),'JMIN':match[0].J.min(),'KMIN':match[0].K.min(),
                                            'IMAX':match[0].I.max(),'JMAX':match[0].J.max(),'KMAX':match[0].K.max()
                                        },{
                                            'index':i,'IMIN':match[1].I.min(),'JMIN':match[1].J.min(),'KMIN':match[1].K.min(),
                                            'IMAX':match[1].I.max(),'JMAX':match[1].J.max(),'KMAX':match[1].K.max()
                                        }])
                                    })



    # Update the outer faces
    outer_faces_formatted = list() # This will contain 
    id = 1 
    for temp in outer_faces:        
        block = {'index':temp['block']}      # Grabs the index 
        surfaces = list()
        for face in temp['outer_faces']:            
            surfaces.append({ 'IMIN':min(face.I), 'JMIN':min(face.J), 'KMIN':min(face.K),
                                'IMAX':max(face.I), 'JMAX':max(face.J), 'KMAX':max(face.K),
                                'id':id })
            id += 1
        block['surfaces'] = surfaces
        outer_faces_formatted.append(block)

    return face_matches, outer_faces_formatted 
