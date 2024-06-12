import sys, os
sys.path.insert(0,'../../')
from plot3d import read_plot3D, write_plot3D, connectivity_fast, periodicity_fast
from plot3d.graph import block_to_graph, get_face_vertex_indices, add_connectivity_to_graph
import numpy as np
import numpy.typing as npt 
import networkx as nx
import pickle

def reindex_vertices(v:npt.NDArray,G:nx.graph.Graph):
    mapping = dict()
    nodes = np.array(list(G.nodes()))
    
    nx.relabel_nodes(G,mapping)
    

if not os.path.exists('./VSPT_Binary.xyz'):
    blocks = read_plot3D('VSPT_ASCII.xyz',binary=False,read_double=False)
    
    face_matches, outer_faces = connectivity_fast(blocks)
    
    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity_fast(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)
    
    face_matches.extend(periodic_surfaces)
    
    block_sizes = [(b.IMAX, b.JMAX, b.KMAX) for b in blocks]
    graphs = list()
    starting_index = 0
    for b in blocks:
        G1 = block_to_graph(b.IMAX,b.JMAX,b.KMAX,starting_index)
        starting_index += b.IMAX*b.JMAX*b.KMAX 
        graphs.append(G1)
    
    G = nx.compose_all(graphs)    
    G = add_connectivity_to_graph(G,block_sizes,face_matches)

    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
                        "face_matches":face_matches, 
                        "outer_faces":outer_faces,
                        "graph":G
                    },f)    
    write_plot3D('VSPT_Binary.xyz',blocks,binary=True)    # Writing plot3D to binary file

else:
    blocks = read_plot3D('VSPT_Binary.xyz',binary=True)
    with open('connectivity.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

