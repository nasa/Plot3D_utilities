from itertools import combinations
import os, sys
from copy import deepcopy
from math import *  
import numpy as np
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, rotated_periodicity, get_outer_faces
from glennht_con import export_to_glennht_conn
import pickle

blocks = read_plot3D('PahtCascade-ASCII.xyz', binary = False)
number_of_blades = 55
rotation_angle = 360/number_of_blades

# Lets find the connectivity and periodicity of the rotated blocks
copies = 3 # Lets use the example with 3 copies 

# Finding Periodicity between first and last block
outer_faces = list() 
for block in blocks:
    _, outer_face = get_outer_faces(block)
    outer_faces.append(outer_face.to_dict())

periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all = rotated_periodicity(blocks, outer_faces, rotation_angle=rotation_angle*4, rotation_axis = "x")

# Find Neighbor Connectivity 

# Rotate the Blocks 
rotated_blocks = list()
rotated_blocks.extend(blocks)

for i in range(1,copies):
    # Rotation matrix can be found on https://en.wikipedia.org/wiki/Rotation_matrix
    rotation_matrix = np.array([[1,0,0],
                            [0,cos(rotation_angle*i),-sin(rotation_angle*i)],
                            [0,sin(rotation_angle*i),cos(rotation_angle*i)]])
    blocks.append(deepcopy(blocks[i].rotate_block(rotation_matrix)))




# write_plot3D('finalmesh-paht_binary.xyz',blocks=blocks,binary=True)
# write_plot3D('finalmesh-paht_ascii.xyz',blocks=blocks,binary=False)