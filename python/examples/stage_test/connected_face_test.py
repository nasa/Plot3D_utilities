import sys, pickle, os
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast, rotated_periodicity, write_plot3D, Direction, split_blocks, find_connected_face

with open('stator_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

# Outerface[53] appears to be a face that is located on the stator surface
outer_faces.remove(53)
connected_faces = find_connected_face(outer_faces[53], outer_faces)
