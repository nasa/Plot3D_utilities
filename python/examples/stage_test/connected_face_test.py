import sys, pickle, os
from venv import create
sys.path.insert(0,'../../')
import numpy as np 
from plot3d import Face, find_connected_face,create_face_from_diagonals,read_plot3D

blocks = read_plot3D('stator_split.xyz',binary=True)
with open('stator_split_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']


def find_face(block_index:int, indices:np.ndarray):
    outer_face_to_match = None
    for o in outer_faces:
        if o['block_index'] == block_index:
            a = np.array([o['IMIN'], o['JMIN'], o['KMIN'], o['IMAX'], o['JMAX'], o['KMAX']], dtype=int)
            if np.array_equal(a,indices):
                outer_face_to_match = create_face_from_diagonals(blocks[o['block_index']], o['IMIN'], o['JMIN'], o['KMIN'], o['IMAX'], o['JMAX'], o['KMAX'])
                outer_face_to_match.set_block_index(block_index)
    return outer_face_to_match

# Outerface[53] appears to be a face that is located on the stator surface
block_id = 14; indices = np.array([0,0,0,52,148,0], dtype=int) # this is the face we need to find matches for
outer_face_to_match = find_face(block_id, indices)

block_id = 15; indices = np.array([0,0,0,52,148,0], dtype=int) # this is the face that it should match
outer_face_should_match = find_face(block_id, indices)

outer_faces2 = list()
for o in outer_faces:
    outer_faces2.append(create_face_from_diagonals(blocks[o['block_index']],o['IMIN'],o['JMIN'],o['KMIN'],o['IMAX'],o['JMAX'],o['KMAX']))
    outer_faces2[-1].set_block_index(o['block_index'])

connected_faces = find_connected_face(outer_face_to_match, outer_faces2)
print('check')