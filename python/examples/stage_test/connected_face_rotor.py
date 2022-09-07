import sys
sys.path.insert(0,'C:\GitHub\Plot3D_utilities\python')
from plot3d import Face, find_connected_face,find_face,read_plot3D
import pickle
import numpy as np 

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
block_id = 15; indices = np.array([0,0,0,68,148,0],dtype=int)
rotor_face_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_faces, outer_faces = find_connected_face(blocks,rotor_mixing_plane_face_to_match, outer_faces)
rotor_faces.append(rotor_face_to_match.to_dict())

# Rotor Hub
block_id = 0; indices = np.array([0,0,0,44,0,76],dtype=int)
rotor_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_hub, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub.append(rotor_hub_to_match.to_dict())

block_id = 10; indices = np.array([0,0,0,68,0,48],dtype=int) 
rotor_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_hub2, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub2.append(rotor_hub_to_match.to_dict())
rotor_hub.extend(rotor_hub2)

# Rotor Shroud
block_id = 0; indices = np.array([0,148,0,44,148,76],dtype=int)
rotor_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
rotor_shroud, outer_faces = find_connected_face(blocks,rotor_shroud_to_match, outer_faces)
rotor_shroud.append(rotor_shroud_to_match.to_dict())

block_id = 10; indices = np.array([0,148,0,68,148,48],dtype=int) 
rotor_shroud_to_match = find_face(blocks,block_id,indices,outer_faces)
rotor_shroud2, outer_faces = find_connected_face(blocks, rotor_shroud_to_match, outer_faces)
rotor_shroud2.append(rotor_shroud_to_match.to_dict())
rotor_shroud.extend(rotor_shroud2)

# Outlet
block_id = 23; indices = np.array([8,0,0,8,148,116],dtype=int) 
rotor_outlet = find_face(blocks,block_id, indices,outer_faces)


data['outer_faces'] = outer_faces
data['mixing_plane'] = rotor_mixing_plane_faces
data['rotor_body'] = rotor_faces
data['rotor_hub'] = rotor_hub
data['rotor_shroud'] = rotor_shroud
data['outlet'] = [rotor_outlet.to_dict()]

with open('rotor_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

