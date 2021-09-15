import sys
import os
import pickle
import json
sys.path.insert(0,os.getcwd())
# from pv_library import Load, ExtractBlocks, ExtractSurface
# from paraview.simple import *
import random


with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

blocks_to_extract = [f['block1']['index'] for f in face_matches]
blocks_to_extract.extend([f['block2']['index'] for f in face_matches])
blocks_to_extract = list(set(blocks_to_extract))


# Generate Random Colors 
rgb_face_matches = list()
for i in range(len(face_matches)):
    rgb_face_matches.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

rgb_outer_faces = list()
for i in range(len(outer_faces)):
    rgb_outer_faces.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

# Load mesh
# plot3d_binary_filename = 'compressor_binary.xyz'
# plot3D_source,plot3D_Display,View,LUT = Load(plot3d_binary_filename)

# Plot the face matches 
for b in blocks_to_extract: # Block indicies 
    for match_indx, f in enumerate(face_matches):
        if f['block1']['index'] == b: 
            voi = [f['block1']['IMIN'], f['block1']['IMAX'], f['block1']['JMIN'], f['block1']['JMAX'],f['block1']['KMIN'], f['block1']['KMAX']]
            print(voi)
            
# print("done")