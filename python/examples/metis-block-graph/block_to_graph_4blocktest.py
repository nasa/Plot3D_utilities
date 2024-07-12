#%% Import Scripts 
import sys
sys.path.insert(0,'../../')
sys.path.insert(0,'~/miniconda3/envs/dev/lib/python3.10/site-packages/')
from plot3d.graph import block_connectivity_to_graph
import numpy as np
import networkx as nx
from read_glennht_to_conn import glennht_to_con
from IPython.display import Image, display

'''  _____ _____
    |     |     |
    |  0  |  1  |
    |_____|_____|
    |     |     |
    |  2  |  3  |
    |_____|_____|
    edge from 0 to 1:
'''
connections = [(0,1),(0,2),
                (2,3),(1,3)]
connection_weights = [10]      # 0 to 1
connection_weights.append(10) # 0 to 2
connection_weights.append(10) # 0 to 2
connection_weights.append(10) # 0 to 2

face_matches = []
for con,weight in zip(connections,connection_weights):
    face_matches.append({
        'block1':{'block_index':0,"IMAX":weight,"JMAX":1,"KMAX":1,
                  "IMIN":weight,"JMIN":1,"KMIN":1},
        "block2":{'block_index':1,"IMAX":weight,"JMAX":1,"KMAX":1,
                "IMIN":weight,"JMIN":1,"KMIN":1}})
nparts=2 # Number of splits 
block_sizes = [20,20,20,20] # Block size

G = block_connectivity_to_graph(face_matches,block_sizes)

#%% Split the Graph
import metis
# G.graph['node_weight_attr'] = ['weight']
G.graph['edge_weight_attr'] = 'weight'
nparts = 3
(edgecuts, parts) = metis.part_graph(G, nparts,tpwgts=None)
nodes_per_part = list()
# Print number of cells for each part
for i in range(nparts):
    indexes = [p for p, e in enumerate(parts) if e == i]
    nodes = 0 
    for idx in indexes:  
        nodes += block_sizes[idx-1]
    nodes_per_part.append(nodes)

for i in range(nparts):
    print(f'Parition {i} has {parts.count(i)} blocks')

# Determine work for each partition
communication_work = np.zeros((nparts,))
partition_edge_weights = np.zeros((nparts,))
for b in range(len(block_sizes)+1):
    partition_id = parts[b]
    for connected_block in G.adj[b]:
        connected_block_partition = parts[connected_block]
        edge_weight = G.adj[b][connected_block]['weight']
        if connected_block_partition != partition_id:
            communication_work[partition_id]+=1
            partition_edge_weights[partition_id] += edge_weight

for i in range(nparts):
    print(f'Parition {i} com_work {communication_work[i]} edge_work {partition_edge_weights[i]}')
