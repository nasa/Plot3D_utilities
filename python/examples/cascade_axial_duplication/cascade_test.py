from itertools import combinations
import os, sys
from copy import deepcopy
from math import *  
import numpy as np
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, rotated_periodicity, get_outer_faces,connectivity_fast
from glennht_con import export_to_glennht_conn
import pickle

blocks = read_plot3D('PahtCascade-ASCII.xyz', binary = False)
number_of_blades = 55
rotation_angle = 360.0/number_of_blades

# Lets find the connectivity and periodicity of the rotated blocks
copies = 3 # Lets use the example with 3 copies 

if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('PahtCascade-ASCII.xyz', binary = False)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity_fast(blocks)
    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

# Finding Periodicity between first and last block
# outer_faces_all = list() 
# for i in range(len(blocks)):
#     outer_faces,_ = get_outer_faces(blocks[i])
#     for o in outer_faces:
#         o.blockIndex=i
#     outer_faces_all.extend([o.to_dict() for o in outer_faces])


# Find Neighbor Connectivity / Interblock-to-block matching
periodic_faces, outer_faces_to_keep, _, _ = rotated_periodicity(blocks,face_matches, outer_faces, rotation_angle=rotation_angle, rotation_axis = "x")

# Append periodic surfaces to face_matches
face_matches.extend(periodic_faces)

with open('connectivity_periodic.pickle','wb') as f:
    # [m.pop('match',None) for m in face_matches] # Remove the dataframe
    pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_faces},f)

export_to_glennht_conn(face_matches,outer_faces_to_keep,'finalmesh')


# Note: 
#   If you rotate and copy the blocks by 3 periods (3 pie slices)
#   periodic_faces gets copied 3 times too. You have one for outer edges and 2 for adjacents

#   If you rotate and copy the blocks by 4 periods (4 pie slices)
#   periodic_faces gets copied 3 times too. You have one for outer edges and 3 for adjacents

 
# Rotate the Blocks 
rotated_blocks = list()
rotated_blocks.extend(blocks)

for i in range(1,copies):
    # Rotation matrix can be found on https://en.wikipedia.org/wiki/Rotation_matrix
    rotation_matrix = np.array([[1,0,0],
                            [0,cos(rotation_angle*i),-sin(rotation_angle*i)],
                            [0,sin(rotation_angle*i),cos(rotation_angle*i)]])
    blocks.append(deepcopy(blocks[i].rotate_block(rotation_matrix)))

write_plot3D('finalmesh_rotated_binary.xyz',blocks=blocks,binary=True)

# Modifying the Connectivity
