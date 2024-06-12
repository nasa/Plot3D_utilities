import sys
sys.path.insert(0,'../../')
import numpy as np
import networkx as nx 
from plot3d import read_plot3D, write_plot3D
from plot3d.graph import block_to_graph, get_face_vertex_indices, add_connectivity_to_graph




#%% Example of a 6x6 block 
IMAX = 4
JMAX = 6
KMAX = 3
G1 = block_to_graph(IMAX,JMAX,KMAX)
G2 = block_to_graph(IMAX,JMAX,KMAX,IMAX*JMAX*KMAX)

#%% Flatten Test: Create a simple Plot3D Block structure
A = np.arange(IMAX*JMAX*KMAX).reshape((KMAX,JMAX,IMAX))
A = np.transpose(A,[2,1,0]) 
A_flat = A.flatten(order='F')

#%% Test get face vertex indices 
indices_imin_face = get_face_vertex_indices(0,0,0,0,JMAX,KMAX,(IMAX,JMAX,KMAX))       # Constant IMIN Face
indices_imax_face = get_face_vertex_indices(IMAX,0,0,IMAX,JMAX,KMAX,(IMAX,JMAX,KMAX)) # Constant IMAX Face

indices_jmin_face = get_face_vertex_indices(0,0,0,IMAX,0,KMAX,(IMAX,JMAX,KMAX))       # Constant JMIN Face
indices_jmax_face = get_face_vertex_indices(0,JMAX,0,IMAX,JMAX,KMAX,(IMAX,JMAX,KMAX)) # Constant JMAX Face

indices_kmin_face = get_face_vertex_indices(0,0,0,IMAX,JMAX,0,(IMAX,JMAX,KMAX))       # Constant KMIN Face
indices_kmax_face = get_face_vertex_indices(0,0,KMAX,IMAX,JMAX,KMAX,(IMAX,JMAX,KMAX)) # Constant KMAX Face

indices_jmin_face_reverseI = get_face_vertex_indices(IMAX,0,0,0,0,KMAX,(IMAX,JMAX,KMAX))       # Constant JMIN Face, reversing direction of I
#%% Test Connectivity 
G1 = block_to_graph(IMAX,JMAX,KMAX)
G2 = block_to_graph(IMAX,JMAX,KMAX,IMAX*JMAX*KMAX)
G = nx.compose_all([G1,G2])

block_sizes = [(IMAX,JMAX,KMAX),(IMAX,JMAX,KMAX)]

# Block 0 and Block 0 share a top face
interblock_k_connectivity = [{
    'block1': 
            {
                'index':0,
                'IMIN':0,'IMAX':IMAX,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':0,'KMAX':0
            },
    'block2': 
            {
                'index':0,
                'IMIN':0,'IMAX':IMAX,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':KMAX,'KMAX':KMAX
            }
    }]

interblock_i_connectivity = [{
    'block1': 
            {
                'index':0,
                'IMIN':0,'IMAX':0,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':0,'KMAX':KMAX
            },
    'block2': 
            {
                'index':0,
                'IMIN':IMAX,'IMAX':IMAX,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':0,'KMAX':KMAX
            }
}]

block_to_block_connectivity = [{
    'block1': 
            {
                'index':0,
                'IMIN':0,'IMAX':0,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':0,'KMAX':KMAX
            },
    'block2': 
            {
                'index':1,
                'IMIN':IMAX,'IMAX':IMAX,
                'JMIN':0,'JMAX':JMAX,
                'KMIN':0,'KMAX':KMAX
            }
}]
block_sizes=[(IMAX,JMAX,KMAX), (IMAX,JMAX,KMAX)]
G = add_connectivity_to_graph(G,block_sizes,interblock_i_connectivity)
G = add_connectivity_to_graph(G,block_sizes,interblock_k_connectivity)

# Block to Block connectivity 
G1 = block_to_graph(IMAX,JMAX,KMAX)
G2 = block_to_graph(IMAX,JMAX,KMAX,IMAX*JMAX*KMAX)
G = nx.compose_all([G1,G2])
G = add_connectivity_to_graph(G,block_sizes,block_to_block_connectivity)
# IMAX, JMAX, KMAX = block[0].IMAX, block[0].JMAX, block[0].KMAX



