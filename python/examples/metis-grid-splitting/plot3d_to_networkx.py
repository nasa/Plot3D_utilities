import sys, os, pickle
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast,translational_periodicity, write_plot3D, Direction, split_blocks, block_connection_matrix, find_bounding_faces
from plot3d import outer_face_dict_to_list, match_faces_dict_to_list
import numpy as np
import numpy.typing as npt
import networkx as nx 
from typing import Dict, Tuple, List

def block_to_graph(IMAX:int,JMAX:int,KMAX:int,offset:int = 0):
    """Converts a block to a graph

    Args:
        IMAX (int): block.IMAX
        JMAX (int): block.JMAX
        KMAX (int): block.KMAX
        offset (int): IMAX*JMAX*KMAX of previous block 

    Returns:
        G: networkx graph object 
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


def getFaceVertexIndicies(IMIN:int,IMAX:int,JMIN:int,JMAX:int,KMIN:int,KMAX:int,block_size:Tuple[int,int,int]) -> npt.NDArray:
    """Returns an array containing the vertex number of a given face 

    Args:
        IMIN (int): starting I index
        IMAX (int): ending I index
        JMIN (int): starting J index
        JMAX (int): ending J index
        KMIN (int): starting K index
        KMAX (int): ending K index
        block_size (Tuple[int,int,int]): This is the actual IMAX,JMAX,KMAX of the block 
    Returns:
        npt.NDArray: _description_
    """

    def create_range(indx1,indx2):
        if indx1<indx2:
            return np.arange(indx1,indx2)
        else:
            return np.arange(indx2-1,indx1-1,-1)
        
    indices = list()
    if IMIN==IMAX:
        jrange = create_range(JMIN,JMAX)
        krange = create_range(KMIN,KMAX)
        if IMIN==block_size[0]: # IMIN is really IMAX
            for k in krange:
                k_offset = k*block_size[0]*block_size[1]    # IMAX * JMAX
                indices.append(k_offset + block_size[0]*(jrange+1)-1)
        else:
            for k in krange:
                k_offset = k*block_size[0]*block_size[1]
                indices.append(k_offset + block_size[0]*jrange)
        
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
    return indices
    
    
def add_connectivity_to_graph(G:nx.classes.graph.Graph,block_sizes:List[Tuple[int,int,int]],connectivities:List[Dict[str,int]]):
    """_summary_

    Args:
        G (nx.classes.graph.Graph): _description_
        block_sizes (List[Tuple[int,int,int]]): _description_
        connectivity (List[Dict[str,int]]): _description_
    """
    
    
    
            
    for con in connectivities: 
        block1_index = con['block1']['index']
        block2_index = con['block2']['index']
        IMIN1,IMAX1 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        JMIN1,JMAX1 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        KMIN1,KMAX1 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        
        IMIN2,IMAX2 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        JMIN2,JMAX2 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        KMIN2,KMAX2 = con['block1']['index']['IMIN'], con['block1']['index']['IMAX']
        
        # Number of connectivities should match
        nodes1 = (IMAX1-IMIN1)*(JMAX1-JMIN1)*(KMAX1-KMIN1)
        nodes2 = (IMAX2-IMIN2)*(JMAX2-JMIN2)*(KMAX2-KMIN2)
        assert nodes1 == nodes2, f"Number of connections from {block1_index} I[{IMIN1},{IMAX1}], J[{JMIN1},{JMAX1}], K[{KMIN1},{KMAX1}] to {block2_index} I[{IMIN2},{IMAX2}], J[{JMIN2},{JMAX2}], K[{KMIN2},{KMAX2}] should match."
        

        # shift1 = getStartingVIndex(block1_index)
        # shift2 = getStartingVIndex(block2_index)
        
        
        G.add_edge(F,T)

    block_sizes[connectivity]

#%% Example of a 6x6 block 
IMAX = 4
JMAX = 6
KMAX = 3
G1 = block_to_graph(IMAX,JMAX,KMAX)
G2 = block_to_graph(IMAX,JMAX,KMAX,IMAX*JMAX*KMAX)

#%% Test get face vertex indices 
indices_imin_face = getFaceVertexIndicies(0,0,0,JMAX,0,KMAX,(IMAX,JMAX,KMAX))       # Constant IMIN Face
indices_imax_face = getFaceVertexIndicies(IMAX,IMAX,0,JMAX,0,KMAX,(IMAX,JMAX,KMAX)) # Constant IMAX Face

indices_jmin_face = getFaceVertexIndicies(0,IMAX,0,0,0,KMAX,(IMAX,JMAX,KMAX))       # Constant JMIN Face
indices_jmax_face = getFaceVertexIndicies(0,IMAX,JMAX,JMAX,0,KMAX,(IMAX,JMAX,KMAX)) # Constant JMAX Face

indices_kmin_face = getFaceVertexIndicies(0,IMAX,0,JMAX,0,0,(IMAX,JMAX,KMAX))       # Constant KMIN Face
indices_kmax_face = getFaceVertexIndicies(0,IMAX,0,JMAX,KMAX,KMAX,(IMAX,JMAX,KMAX)) # Constant KMAX Face

#%% Test Connectivity 


G = nx.compose_all([G1,G2])
block_sizes = [(IMAX,JMAX,KMAX),(IMAX,JMAX,KMAX)]
# Block 0 and Block 0 share a top face
connectivity = [{
    'block1': 
            {
                'index':0,
                'IMIN':0,'IMAX':6,
                'JMIN':0,'JMAX':6,
                'KMIN':0,'KMAX':0
            },
    'block2': 
            {
                'index':0,
                'IMIN':0,'IMAX':6,
                'JMIN':0,'JMAX':6,
                'KMIN':3,'KMAX':3
            }
    }]

connectivity.append({
    'block1': 
            {
                'index':0,
                'IMIN':6,'IMAX':6,
                'JMIN':0,'JMAX':6,
                'KMIN':0,'KMAX':3
            },
    'block2': 
            {
                'index':1,
                'IMIN':0,'IMAX':0,
                'JMIN':0,'JMAX':6,
                'KMIN':0,'KMAX':3
            }
    })


# IMAX, JMAX, KMAX = block[0].IMAX, block[0].JMAX, block[0].KMAX
graphs = list() 
graphs.append(G)

print('done')
    