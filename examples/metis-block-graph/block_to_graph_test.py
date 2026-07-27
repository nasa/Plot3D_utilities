#%% Import Scripts 
import sys, os
# sys.path.insert(1,'/mnt/c/GitHub/metis-python')
from plot3d import read_plot3D, write_plot3D, connectivity_fast, periodicity_fast, match_faces_dict_to_list, outer_face_dict_to_list
from plot3d.graph import block_to_graph, get_face_vertex_indices, add_connectivity_to_graph, block_connectivity_to_graph
import numpy as np
import networkx as nx
import pickle


def dump_data(data):
    with open('data.pickle','wb') as f:
        pickle.dump(data,f)

def read_data():
    with open('data.pickle','rb') as f:
        return pickle.load(f)
    
from IPython.display import Image, display

def view_pydot(pdot):
    plt = Image(pdot.create_png())
    display(plt)
    
    
if not os.path.exists('data.pickle'):
    
    blocks = read_plot3D('iso65_64blocks.xyz',binary=True,read_double=False)

    print('Finding connectivity')
    face_matches, outer_faces = connectivity_fast(blocks)
    [m.pop('match',None) for m in face_matches] # Remove the dataframe
    print('Organizing split and outerfaces')
    all_faces = match_faces_dict_to_list(blocks,face_matches)
    all_faces.extend(outer_face_dict_to_list(blocks,outer_faces))
    all_faces = [m.to_dict() for m in all_faces]

    block_sizes = [b.size for b in blocks]
    # This is where the edge weighting and block cell count weighting takes place
    G = block_connectivity_to_graph(face_matches,block_sizes)
    dump_data({'Graph':G})
#%%
import metis
num_partitions = 2              # Number of partitions
partition_weights = [0.5,0.5]   # This is how to weight each partition

G = read_data()['Graph']
G.graph['node_weight_attr'] = ['weight']
G.graph['edge_weight_attr'] = 'weight'
(edgecuts, parts) = metis.part_graph(G, 2,tpwgts=[0.5,0.5])

for i in range(num_partitions):
    print(f'Parition {i} has {parts.count(i)} blocks')

# colors = ['red','blue','green','magenta']
# for i, p in enumerate(parts):
#     G.nodes[i]['color'] = colors[p]
# nx.drawing.nx_pydot.write_dot(G, 'example.dot') # Requires pydot or pygraphviz
# pdot = nx.drawing.nx_pydot.to_pydot(G)

# red = parts.count(0)
# blue = parts.count(1)
# green = parts.count(2)
# magenta = parts.count(3)
# view_pydot(pdot)
# print('done')
