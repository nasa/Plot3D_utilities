#%% Import Scripts 
import sys
sys.path.insert(0,'../../')
sys.path.insert(0,'~/miniconda3/envs/dev/lib/python3.10/site-packages/')
from plot3d.graph import block_connectivity_to_graph
import numpy as np
import networkx as nx
from read_glennht_to_conn import glennht_to_con
from IPython.display import Image, display

face_matches = []