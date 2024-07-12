import numpy as np
import numpy.typing as npt
import networkx as nx 
import itertools as it
from typing import Dict, Tuple, List
import tqdm 

def block_to_graph(IMAX:int,JMAX:int,KMAX:int,offset:int = 0) -> nx.graph.Graph:
    """Converts a block to a graph

    Args:
        IMAX (int): block.IMAX
        JMAX (int): block.JMAX
        KMAX (int): block.KMAX
        offset (int): IMAX*JMAX*KMAX of previous block 

    Returns:
        nx.graph.Graph: networkx graph object 
    """
    G = nx.Graph()
    irange = np.arange(IMAX)
    jrange = np.arange(JMAX)
    krange = np.arange(KMAX)
    kshift = 0 
    for k in range(KMAX): # K slices 
        kshift = IMAX*JMAX*k

        for j in range(JMAX):
            nx.add_star(G, offset+ kshift + IMAX*j + irange,weight=1)

        for i in range(IMAX):
            nx.add_star(G, offset+kshift + i + IMAX*jrange,weight=1)

    for p in range(IMAX*JMAX):
        nx.add_star(G, offset + p + IMAX*JMAX*krange)
    return G

def get_face_vertex_indices(IMIN:int,JMIN:int,KMIN:int,IMAX:int,JMAX:int,KMAX:int,block_size:Tuple[int,int,int]) -> npt.NDArray:
    """Returns an array containing the vertex number of a given face 

    Args:
        IMIN (int): starting I index
        JMIN (int): starting J index
        KMIN (int): starting K index
        IMAX (int): ending I index
        JMAX (int): ending J index
        KMAX (int): ending K index
        block_size (Tuple[int,int,int]): This is the actual IMAX,JMAX,KMAX of the block
    
    Returns:
        npt.NDArray: an array containing all the vertices 
    """

    def create_range(indx1,indx2):
        if indx1<indx2:
            return np.arange(indx1,indx2)
        else:
            return np.arange(indx1-1,indx2-1,-1)
        
    indices = list()
    if IMIN==IMAX:
        jrange = create_range(JMIN,JMAX)
        krange = create_range(KMIN,KMAX)
        if IMIN==block_size[0]: # IMIN is really IMAX
            for j in jrange:
                j_offset = j*block_size[0]    # IMAX * JMAX
                indices.append(j_offset + block_size[0]*block_size[1]*krange + IMAX-1)
        else:
            for j in jrange:
                j_offset = j*block_size[0]
                indices.append(j_offset + block_size[0]*block_size[1]*krange)
        
    elif JMIN == JMAX:
        irange = create_range(IMIN,IMAX)
        krange = create_range(KMIN,KMAX)
        if JMIN==block_size[1]: # JMIN is really JMAX
            for k in krange:
                k_offset = k*block_size[0]*block_size[1] 
                indices.append(k_offset + block_size[0]*(block_size[1]-1)+irange)
        else:                   # JMIN
            for k in krange:
                k_offset = k*block_size[0]*block_size[1]
                indices.append(k_offset + irange)
    else:
        irange = create_range(IMIN,IMAX)
        jrange = create_range(JMIN,JMAX)
        if KMIN == block_size[2]: # KMIN is really KMAX
            offset = (KMIN-1)*block_size[0]*block_size[1] 
        else:
            offset = 0 
        for j in jrange:
            indices.append(offset+block_size[0]*j + irange)
    return np.array(indices).flatten()

def get_starting_vertex(blockIndex:int,block_sizes:List[Tuple[int,int,int]]) -> int:
    """Gets the starting vertex index of the block

    Args:
        blockIndex (int): index of block
        block_sizes (List[Tuple[int,int,int]]): List of all the [[IMAX,JMAX,KMAX]] 

    Returns:
        int: offset
    """
    offset = 0 
    for i in range(blockIndex):
        offset+=block_sizes[i][0]*block_sizes[i][1]*block_sizes[i][2]
    return offset

def add_connectivity_to_graph(G:nx.classes.graph.Graph,block_sizes:List[Tuple[int,int,int]],connectivities:List[Dict[str,int]]) -> nx.graph.Graph:
    """Convert plot3d defined connectivity into additional graph edges 

    Args:
        G (nx.classes.graph.Graph): Giant graph 
        block_sizes (List[Tuple[int,int,int]]): List of all the [[IMAX,JMAX,KMAX]] 
        connectivity (List[Dict[str,int]]): _description_
    
    Returns:
        nx.graph.Graph: networkx graph object with added edges 
    """
    
    for con in tqdm.tqdm(connectivities,"Adding connectivity to Graph"):
        block1_index = con['block1']['block_index']
        block2_index = con['block2']['block_index']
        IMIN1,IMAX1 = con['block1']['IMIN'], con['block1']['IMAX']
        JMIN1,JMAX1 = con['block1']['JMIN'], con['block1']['JMAX']
        KMIN1,KMAX1 = con['block1']['KMIN'], con['block1']['KMAX']
        
        IMIN2,IMAX2 = con['block2']['IMIN'], con['block2']['IMAX']
        JMIN2,JMAX2 = con['block2']['JMIN'], con['block2']['JMAX']
        KMIN2,KMAX2 = con['block2']['KMIN'], con['block2']['KMAX']
        
        # Number of connectivities should match
        face1 = get_face_vertex_indices(IMIN1,JMIN1,KMIN1,IMAX1,JMAX1,KMAX1,block_sizes[block1_index]) + get_starting_vertex(block1_index, block_sizes)    
        face2 = get_face_vertex_indices(IMIN2,JMIN2,KMIN2,IMAX2,JMAX2,KMAX2,block_sizes[block2_index]) + get_starting_vertex(block2_index, block_sizes)
        
        if block1_index!= block2_index:
            nodes_to_add = face1
            nodes_to_replace = face2
            for node_to_add,node_to_replace in tqdm.tqdm(zip(nodes_to_add,nodes_to_replace)):
                G.add_edges_from(
                    it.product(
                        G.neighbors(node_to_add),
                        G.neighbors(node_to_replace)
                        )
                )
                G.remove_node(node_to_replace)
                
        assert len(face1) == len(face2), f"Number of connections from {block1_index} I[{IMIN1},{IMAX1}], J[{JMIN1},{JMAX1}], K[{KMIN1},{KMAX1}] to {block2_index} I[{IMIN2},{IMAX2}], J[{JMIN2},{JMAX2}], K[{KMIN2},{KMAX2}] should match."
        
        for i in range(len(face1)):
            G.add_edge(face1[i],face2[i])
            
    return G

def block_connectivity_to_graph(connectivities:List[Dict[str,int]],block_sizes:List[int],connectivity_multiplier:float=1,block_size_multiplier:float=1) -> nx.graph.Graph:
    """Models the blocks at vertices connected to each other 

    Args:
        connectivities (List[Dict[str,int]]): List of connectivities 
        block_sizes (List[Tuple[int,int,int]]): List of all the [[IMAX,JMAX,KMAX]] 
        connectivity_multiplier (float): amount to weight the block connections by. Defaults to 1.
        block_size_multiplier (float): amount to weight the size of the blocks by. Defaults to 1 
        
    Returns:
        nx.graph.Graph: Graph object 
    """
    block_to_block = list()
    G = nx.Graph()
    for con in tqdm.tqdm(connectivities,"Adding connectivity to Graph"):
        block1_index = con['block1']['block_index']
        block2_index = con['block2']['block_index']
        dI = max(con['block1']['IMAX']-con['block1']['IMIN'],1)
        dJ = max(con['block1']['JMAX']-con['block1']['JMIN'],1)
        dK = max(con['block1']['KMAX']-con['block1']['KMIN'],1)
        edge_weight = dI * dJ * dK *connectivity_multiplier
        # Weight the edges based on number of connections
        block_to_block.append((block1_index,block2_index,{"weight": int(edge_weight)}))
    
    # Adds the nodes and node weights
    for i in range(len(block_sizes)):
        G.add_node(i, weight=block_sizes[i]*block_size_multiplier)
    
    # Adds the connectivity information 
    G.add_edges_from(block_to_block)
    
    return G
    
    