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
from pv_library import Load, ExtractBlocks, CreateSubset, plot_faces_for_block
from paraview.simple import *
import random
from plot3d.gridpro import bc_faces_by_type


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
    inlets = bc_faces_by_type(data['bc_group'], 'inlet')
    outlets = bc_faces_by_type(data['bc_group'], 'outlet')
    symm_slip = bc_faces_by_type(data['bc_group'], 'symm_slip')
    walls = bc_faces_by_type(data['bc_group'], 'wall')

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


rgb_outer_faces = list()
for i in range(len(outer_faces)):
    rgb_outer_faces.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

rgb_inlet = [0,128,128] # teal
rgb_outlet = [1,141,161] # pink
rgb_wall = [255,165,0] # bright orange
rgb_symm_slip = [0, 158, 96] # shamrock 

rgb_inlet = [x/255 for x in rgb_inlet]
rgb_outlet = [x/255 for x in rgb_outlet]
rgb_wall = [x/255 for x in rgb_wall]
rgb_symm_slip = [x/255 for x in rgb_symm_slip]

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

    # Plot the inlets 
    plot_faces_for_block(block_source,b,inlets,"inlet",rgb_inlet) # type: ignore
    plot_faces_for_block(block_source,b,outlets,"outlet",rgb_outlet) # type: ignore
    plot_faces_for_block(block_source,b,walls,"wall",rgb_wall) # type: ignore
    plot_faces_for_block(block_source,b,symm_slip,"symm_slip",rgb_symm_slip) # type: ignore
    plot_faces_for_block(block_source,b,outer_faces,"outer_faces",rgb_outer_faces) # type: ignore

     
