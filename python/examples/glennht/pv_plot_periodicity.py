'''
    plots the boundary conditions using paraview. 
    
    How to Run:
        Add the paraview executable to path
        run in terminal `paraview --script=pv_plot_bcs.py` this will plo
'''
import sys
import os
import json
from typing import List, Dict
sys.path.insert(0,os.getcwd()) # This allows you to select files locally
from pv_library import Load, ExtractBlocks, CreateSubset, plot_face_matches_for_block
from paraview.simple import *
import random


plot_face_matches = False
plot_periodicity = True
plot_outer_faces = False 
plot3d_filename = 'test.xyz'
with open('test.json','r') as f:
    data = json.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']
    periodic_faces = data['periodicity']
    volume_zones = data['volume_zones']
    
    # Boundary Conditions
    inlets = data['bc_group']['inlet']
    outlets = data['bc_group']['outlet']
    symm_slip = data['bc_group']['symm_slip']
    walls = data['bc_group']['wall']

blocks_to_extract = [f['block1']['block_index'] for f in periodic_faces]
blocks_to_extract.extend([f['block2']['block_index'] for f in periodic_faces])
blocks_to_extract.extend([f['block_index'] for f in outer_faces])
blocks_to_extract.extend([f['block1']['block_index'] for f in face_matches])
blocks_to_extract.extend([f['block2']['block_index'] for f in face_matches])
blocks_to_extract = list(set(blocks_to_extract))
blocks_to_extract.sort()
n = len(blocks_to_extract)

'''
Generate Random Colors 

This part generates random colors so that each color is associated with a face match. 
Doesn't matter what the block is a single match is assigned the same color. 
'''
rgb_periodic = list()
for i in range(len(periodic_faces)):
    rgb_periodic.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

# Load mesh
plot3D_source,plot3D_Display,View,LUT = Load(plot3d_filename,True)
print(f"Total number of blocks: {n}")

for b in blocks_to_extract: # Block indicies
    block_source,block_display,LUT = ExtractBlocks(plot3D_source,View,[b])
    vz = volume_zones[b]
    src_name = f'Block {b} {vz['zone_type']}-{vz['contiguous_id']}'
    RenameSource(src_name, block_source)
    block_source = FindSource(src_name)
    
    Hide(block_source, View)

    plot_face_matches_for_block(block_source,b,periodic_faces,"periodic-z",rgb_periodic)
   
     