import sys, os
sys.path.insert(0,'../../')
from plot3d import read_plot3D, write_plot3D, connectivity_fast, periodicity_fast
from plot3d.graph import block_to_graph, get_face_vertex_indices, add_connectivity_to_graph
import numpy as np
import numpy.typing as npt 
import networkx as nx
import pickle
from glennht_con import export_to_glennht_conn
import json


def reindex_vertices(v:npt.NDArray,G:nx.graph.Graph):
    mapping = dict()
    nodes = np.array(list(G.nodes()))
    
    nx.relabel_nodes(G,mapping)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

if not os.path.exists('./VSPT_Binary.xyz'):
    blocks = read_plot3D('VSPT_ASCII.xyz',binary=False,read_double=False)
    
    face_matches, outer_faces = connectivity_fast(blocks)
    
    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity_fast(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)
    
    face_matches.extend(periodic_surfaces)
    export_to_glennht_conn(face_matches,outer_faces_to_keep,'vspt')

    
    
    block_sizes = [(b.IMAX, b.JMAX, b.KMAX) for b in blocks]
    graphs = list()
    starting_index = 0
    for b in blocks:
        G1 = block_to_graph(b.IMAX,b.JMAX,b.KMAX,starting_index)
        starting_index += b.IMAX*b.JMAX*b.KMAX 
        graphs.append(G1)
    
    
    G = nx.compose_all(graphs)    

    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
                        "face_matches":face_matches, 
                        "outer_faces":outer_faces_to_keep,
                        "graph":G,
                        "block_sizes":block_sizes
                    },f)
    
    with open('connectivity-julia.json','w') as f:
        
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        json.dump({
                        "face_matches":face_matches, 
                        "outer_faces":outer_faces_to_keep,  
                        "block_sizes":block_sizes                      
                    },f,cls=NpEncoder)    
    write_plot3D('VSPT_Binary.xyz',blocks,binary=True)    # Writing plot3D to binary file
    
# # blocks = read_plot3D('VSPT_Binary.xyz',binary=True)
# with open('connectivity.pickle','rb') as f:
#     data = pickle.load(f)
#     face_matches = data['face_matches']
#     outer_faces = data['outer_faces']
#     G = data["graph"]
#     block_sizes = data["block_sizes"]
# print("done")

