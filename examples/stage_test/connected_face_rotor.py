import sys
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
rotor_mixing_plane_face_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_mixing_plane_faces, outer_faces = find_connected_face(blocks,rotor_mixing_plane_face_to_match, outer_faces)# Should match with block 10 [0,0,0,0,148,52]
rotor_mixing_plane_faces.append(rotor_mixing_plane_face_to_match.to_dict())  

# Rotor body
block_id = 8; indices = np.array([0,0,0,96,148,0],dtype=int)
rotor_face_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_faces, outer_faces = find_connected_face(blocks,rotor_face_to_match, outer_faces)
rotor_faces.append(rotor_face_to_match.to_dict())

# Rotor Hub
block_id = 1; indices = np.array([0,0,0,60,0,76],dtype=int)
rotor_hub_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_hub, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub.append(rotor_hub_to_match.to_dict())

block_id = 7; indices = np.array([0,0,0,60,0,40],dtype=int) 
rotor_hub_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_hub2, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub2.append(rotor_hub_to_match.to_dict())
rotor_hub.extend(rotor_hub2)

block_id = 10; indices = np.array([0,0,0,96,0,48],dtype=int) 
rotor_hub_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_hub3, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub3.append(rotor_hub_to_match.to_dict())
rotor_hub.extend(rotor_hub3)

block_id = 14; indices = np.array([0,0,0,40,0,116],dtype=int) 
rotor_hub_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_hub4, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
rotor_hub4.append(rotor_hub_to_match.to_dict())
rotor_hub.extend(rotor_hub4)

# block_id = 13; indices = np.array([0,0,0,8,0,48],dtype=int) 
# rotor_hub_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
# rotor_hub5, outer_faces = find_connected_face(blocks,rotor_hub_to_match, outer_faces)
# rotor_hub5.append(rotor_hub_to_match.to_dict())
# rotor_hub.extend(rotor_hub5)

# Rotor Shroud
block_id = 1; indices = np.array([0,148,0,60,148,76],dtype=int)
rotor_shroud_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_shroud, outer_faces = find_connected_face(blocks,rotor_shroud_to_match, outer_faces)
rotor_shroud.append(rotor_shroud_to_match.to_dict())

block_id = 12; indices = np.array([0,148,0,96,148,48],dtype=int) 
rotor_shroud_to_match,indx = find_face(blocks,block_id,indices,outer_faces)
rotor_shroud2, outer_faces = find_connected_face(blocks, rotor_shroud_to_match, outer_faces)
rotor_shroud2.append(rotor_shroud_to_match.to_dict())
rotor_shroud.extend(rotor_shroud2)

block_id = 7; indices = np.array([0,148,0,60,148,40],dtype=int) 
rotor_shroud_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_shroud3, outer_faces = find_connected_face(blocks,rotor_shroud_to_match, outer_faces)
rotor_shroud3.append(rotor_hub_to_match.to_dict())
rotor_shroud.extend(rotor_shroud3)

block_id = 14; indices = np.array([0,148,0,40,148,116],dtype=int) 
rotor_shroud_to_match,indx = find_face(blocks,block_id, indices,outer_faces)
rotor_shroud4, outer_faces = find_connected_face(blocks,rotor_shroud_to_match, outer_faces)
rotor_shroud4.append(rotor_shroud_to_match.to_dict())
rotor_shroud.extend(rotor_shroud4)


data['outer_faces'] = outer_faces
data['mixing_plane'] = rotor_mixing_plane_faces
data['rotor_body'] = rotor_faces
data['rotor_hub'] = rotor_hub
data['rotor_shroud'] = rotor_shroud

with open('rotor_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

# Outlet
block_id = 17; indices = np.array([28,0,0,28,148,116],dtype=int)
rotor_outlet,indx = find_face(blocks,block_id, indices,outer_faces)
outer_faces.pop(indx)

# Add matching near leading edge
face_match_1,indx = find_face(blocks,10, np.array([12,0,48,52,148,48], dtype=int),outer_faces)
outer_faces.pop(indx)
face_match_2,indx = find_face(blocks,6, np.array([60,0,0,60,148,40], dtype=int),outer_faces)
outer_faces.pop(indx)
data['face_matches'].append({'block1':face_match_1.to_dict(),
                            'block2':face_match_2.to_dict()})

data['outlet'] = [rotor_outlet.to_dict()]
data['outer_faces'] = outer_faces

with open('rotor_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
