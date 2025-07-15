#%% Import Scripts 
from plot3d.graph import block_connectivity_to_graph
import numpy as np
import networkx as nx
from read_glennht_to_conn import glennht_to_con
from IPython.display import Image, display

def view_pydot(pdot):
    plt = Image(pdot.create_png())
    display(plt)
    
face_matches =  glennht_to_con('kenji_diced_1.p3d_conn')
max_block_index = [(f['block1']['block_index'], f['block2']['block_index']) for f in face_matches]
max_block_index = np.max(np.array(max_block_index).flatten())
block_sizes = [(np.random.randint(33,101)*np.random.randint(33,101)*np.random.randint(33,101)) for _ in range(max_block_index)]
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
for b in range(max_block_index+1):
    partition_id = parts[b]
    for connected_block in G.adj[b]:
        connected_block_partition = parts[connected_block]
        edge_weight = G.adj[b][connected_block]['weight']
        if connected_block_partition != partition_id:
            communication_work[partition_id]+=1
            partition_edge_weights[partition_id] += edge_weight

for i in range(nparts):
    print(f'Parition {i} com_work {communication_work[i]} edge_work {partition_edge_weights[i]}')

    
#%% Plotting
colors = ['red','blue','green','magenta']
for i, p in enumerate(parts):
    G.nodes[i]['color'] = colors[p]
nx.drawing.nx_pydot.write_dot(G, 'example.dot') # Requires pydot or pygraphviz
pdot = nx.drawing.nx_pydot.to_pydot(G)

red = parts.count(0)
blue = parts.count(1)
green = parts.count(2)
magenta = parts.count(3)
view_pydot(pdot)
print('done')
