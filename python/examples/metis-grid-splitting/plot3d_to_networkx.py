import sys, os, pickle
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast,translational_periodicity, write_plot3D, Direction, split_blocks, block_connection_matrix, find_bounding_faces
from plot3d import outer_face_dict_to_list, match_faces_dict_to_list
import numpy as np
import networkx as nx 
from typing import Dict, Tuple, List

def block_to_graph(IMAX:int,JMAX:int,KMAX:int,offset:int):
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




def add_connectivity_to_graph(G:nx.classes.graph.Graph,block_sizes:List[Tuple[int,int,int]],connectivity:List[Dict[str,int]]):
    """_summary_

    Args:
        G (nx.classes.graph.Graph): _description_
        block_sizes (List[Tuple[int,int,int]]): _description_
        connectivity (List[Dict[str,int]]): _description_
    """
    pass 

# Example of a 6x6 block 
IMAX = 6
JMAX = 6
KMAX = 3
G1 = block_to_graph(IMAX,JMAX,KMAX)
G2 = block_to_graph(IMAX,JMAX,KMAX,IMAX*JMAX*KMAX)

# Block 1 and Block 1 share a face: O-Grid 
connectivity = [
    {'block1': {'index':0,
                'IMIN':0,'IMAX':6,
                'JMIN':0,'JMAX':6,
                'KMIN':0,'KMAX':0}}
]

# IMAX, JMAX, KMAX = block[0].IMAX, block[0].JMAX, block[0].KMAX
graphs = list() 
graphs.append(G)
nx.compose_all(graphs)
print('done')
    