import sys, pickle, os
from venv import create
import numpy as np 
from plot3d import Face, find_connected_face,find_face,read_plot3D

# Stator 
blocks = read_plot3D('stator_split.xyz',binary=True)
with open('stator_split_connectivity_periodicity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

# Stator body
block_id = 14; indices = np.array([0,0,0,52,148,0], dtype=int) # this is the face we need to find matches for
stator_face_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_faces,outer_faces = find_connected_face(blocks,stator_face_to_match, outer_faces)
stator_faces.append(stator_face_to_match.to_dict())

# Stator Hub
block_id = 0; indices = np.array([0,0,0,32,0,76],dtype=int)
stator_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_hub, outer_faces = find_connected_face(blocks,stator_hub_to_match, outer_faces)
stator_hub.append(stator_hub_to_match.to_dict())


# Stator Shroud
block_id = 0; indices = np.array([0,148,0,32,148,76],dtype=int)
stator_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_shroud, outer_faces = find_connected_face(blocks,stator_shroud_to_match, outer_faces)
stator_shroud.append(stator_shroud_to_match.to_dict())

# Mixing Plane
block_id = 9; indices = np.array([36, 0,0, 36,148,36], dtype=int)
mixing_plane_face_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_mixing_plane_faces, outer_faces = find_connected_face(blocks,stator_face_to_match, outer_faces)
stator_mixing_plane_faces.append(stator_face_to_match.to_dict())

data['outer_faces'] = outer_faces
data['mixing_plane'] = stator_mixing_plane_faces
data['stator_body'] = stator_faces
data['stator_shroud'] = stator_shroud
data['stator_hub'] = stator_hub

with open('stator_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

# Rotor
blocks = read_plot3D('rotor_split.xyz',binary=True)
with open('rotor_split_connectivity_periodicity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

# Mixing Plane
block_id = 0; indices = np.array([0,0,0,0,148,76], dtype=int) # this is the face we need to find matches for
rotor_mixing_plane_face_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_mixing_plane_faces, outer_faces = find_connected_face(blocks,rotor_mixing_plane_face_to_match, outer_faces)# Should match with block 10 [0,0,0,0,148,52]
rotor_mixing_plane_faces.append(rotor_mixing_plane_face_to_match.to_dict())  

# Rotor body
block_id = 15; indices = np.array([0,0,0,52,148,0],dtype=int)
rotor_face_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_faces, outer_faces = find_connected_face(blocks,rotor_mixing_plane_face_to_match, outer_faces)
rotor_faces.append(rotor_face_to_match.to_dict())

# Rotor Hub
block_id = 0; indices = np.array([0,148,0,32,148,76],dtype=int)
rotor_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_hub, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub.append(rotor_hub_to_match.to_dict())

# Rotor Shroud
block_id = 0; indices = np.array([0,0,0,32,0,76],dtype=int)
rotor_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_shroud, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_shroud.append(rotor_shroud_to_match.to_dict())

data['outer_faces'] = outer_faces
data['mixing_plane'] = rotor_mixing_plane_faces
data['rotor_body'] = rotor_faces
data['rotor_hub'] = rotor_hub
data['rotor_shroud'] = rotor_shroud

with open('rotor_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

