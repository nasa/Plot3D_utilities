import os, sys
from copy import deepcopy
from math import *  
import numpy as np
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, rotated_periodicity,connectivity_fast, rotate_block, create_rotation_matrix
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


# Find Neighbor Connectivity / Interblock-to-block matching
periodic_faces, outer_faces_to_keep, _, _ = rotated_periodicity(blocks,face_matches, outer_faces, rotation_angle=rotation_angle, rotation_axis = "x")

# Finding Periodicity between inner blocks
inner_periodicities = list()
for i in range(1,copies):
    temp = deepcopy(periodic_faces)
    for j in range(len(temp)):        
        temp[j]['block1']['block_index'] += int(i*len(blocks)) # Next set of blocks matches previous set
        temp[j]['block2']['block_index'] += int((i-1)*len(blocks)) # Right Block matches next set blocks
    inner_periodicities.extend(temp)

# Add outer periodicities
outer_periodicities = deepcopy(periodic_faces)
for ou in outer_periodicities:
    ou['block2']['block_index'] += int((copies-1)*len(blocks))

# Copy face matches
face_matches_all = list()
for i in range(1,copies):
    temp = deepcopy(face_matches)
    for ip in temp:
        ip['block1']['block_index'] += int(i*len(blocks))
        ip['block2']['block_index'] += int(i*len(blocks))
    face_matches_all.extend(temp)
# Append periodic surfaces to face_matches
inner_periodicities.extend(outer_periodicities)

face_matches_all.extend(face_matches)

with open('connectivity_periodic.pickle','wb') as f:
    # [m.pop('match',None) for m in face_matches] # Remove the dataframe
    pickle.dump({"face_matches":face_matches_all, "outer_faces":outer_faces_to_keep, "periodic_surfaces":inner_periodicities},f)

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
    rotation_matrix = create_rotation_matrix(radians(rotation_angle*i),'x')
    for i in range(len(blocks)):
        rotated_blocks.append(deepcopy(rotate_block(blocks[i],rotation_matrix)))

write_plot3D('finalmesh_rotated_binary.xyz',blocks=rotated_blocks,binary=True)

