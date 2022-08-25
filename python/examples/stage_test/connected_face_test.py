import sys, pickle, os
from venv import create
sys.path.insert(0,'../../')
from plot3d import Face, find_connected_face,create_face_from_diagonals,read_plot3D

blocks = read_plot3D('stator_split.xyz',binary=True)
with open('stator_split_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches_dict = data['face_matches']
    outer_faces_dict = data['outer_faces']

# Outerface[53] appears to be a face that is located on the stator surface
face = outer_faces_dict[51]
# outer_faces_dict.pop(53)

outer_faces = list()
for o in outer_faces_dict:
    outer_faces.append(create_face_from_diagonals(blocks[o['block_index']],o['IMIN'],o['JMIN'],o['KMIN'],o['IMAX'],o['JMAX'],o['KMAX']))
    outer_faces[-1].set_block_index(o['block_index'])

temp = create_face_from_diagonals(blocks[face['block_index']],face['IMIN'],face['JMIN'],face['KMIN'],face['IMAX'],face['JMAX'],face['KMAX'])
temp.set_block_index(face['block_index'])
face = temp
connected_faces = find_connected_face(face, outer_faces)