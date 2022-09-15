# Imports
import sys
sys.path.insert(0,'C:\GitHub\Plot3D_utilities\python')
from typing import List
from plot3d import read_plot3D, connectivity_fast, rotated_periodicity, write_plot3D, Direction, split_blocks
import os, pickle
import numpy as np
import json 

mesh_filename = 'StageMesh.xyz'
print("Reading mesh")
# blocks = read_plot3D(mesh_filename,binary=True,big_endian=False)

# stator_blocks = [blocks[i] for i in range(0,4)]
# rotor_blocks = [blocks[i] for i in range(4,len(blocks))]

# write_plot3D('stator.xyz',stator_blocks,binary=True)
# stator_blocks_split = split_blocks(stator_blocks,500000, direction=Direction.i) # Splits the blocks while keeping the gcd. 380,000 is a rough number that it tries to match
# write_plot3D('stator_split.xyz',stator_blocks_split,binary=True)

# write_plot3D('rotor.xyz',rotor_blocks,binary=True)
# rotor_blocks_split = split_blocks(rotor_blocks,700000, direction=Direction.i)
# write_plot3D('rotor_split.xyz',rotor_blocks_split,binary=True)


def find_connectivity(filename:str,nblades:int):
    blocks = read_plot3D(f'{filename}.xyz',binary=True,big_endian=False)

    # Finds the connectivity 
    if not os.path.exists(f'{filename}_connectivity.pickle'):
        print('checking connectivity')
        face_matches, outer_faces_formatted = connectivity_fast(blocks)
        with open(f'{filename}_connectivity.pickle','wb') as f:
            [m.pop('match',None) for m in face_matches] # Remove the dataframe
            pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)


    # Finds the periodicity once connectivity is found   
    with open(f'{filename}_connectivity.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']
    
    print("Find periodicity")
    rotation_angle = 360.0/nblades 
    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = rotated_periodicity(blocks,face_matches,outer_faces,rotation_axis='x',rotation_angle=rotation_angle)

    with open(f'{filename}_connectivity_periodicity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
            "face_matches":face_matches,
            "periodic_faces":periodic_surfaces,
            "outer_faces":outer_faces_to_keep       
            },f)

find_connectivity(filename="stator_split",nblades=55)
find_connectivity(filename="rotor_split",nblades=60)