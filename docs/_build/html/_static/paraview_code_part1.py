from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
sys.path.insert(1,"../Cascade")
from plot3d import write_plot3D, read_plot3D, find_periodicity
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

# Convert to binary because of size 
if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('[Path to your Plot3D File]', binary = False) # Paraview has issues with ASCII 
    write_plot3D('[Path to your Plot3D File]',blocks, binary = True)   # Converted to binary for easy reading 

    # Block 1 is the blade O-Mesh k=0
    face_matches, outer_faces_formatted = connectivity(blocks)
    with open('connectivity.pickle','wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)
