import sys
sys.path.insert(0,'C:\GitHub\Plot3D_utilities\python')
from plot3d import Face, find_connected_face,find_face,read_plot3D
import pickle
import numpy as np 
# Stator 
blocks = read_plot3D('stator_split.xyz',binary=True)
with open('stator_split_connectivity_periodicity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

# Stator body
block_id = 12; indices = np.array([0,0,0,68,148,0], dtype=int) # this is the face we need to find matches for
stator_face_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_faces,outer_faces = find_connected_face(blocks,stator_face_to_match, outer_faces)
stator_faces.append(stator_face_to_match.to_dict())

# Stator Hub
block_id = 0; indices = np.array([0,0,0,44,0,76],dtype=int)
stator_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_hub, outer_faces = find_connected_face(blocks,stator_hub_to_match, outer_faces)
stator_hub.append(stator_hub_to_match.to_dict())

block_id = 6; indices = np.array([0,0,0,60,0,36],dtype=int) 
stator_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_hub2, outer_faces = find_connected_face(blocks,stator_hub_to_match, outer_faces)
stator_hub2.append(stator_hub_to_match.to_dict())
stator_hub.extend(stator_hub2)

block_id = 9; indices = np.array([0,0,0,68,0,48],dtype=int) 
stator_hub_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_hub3, outer_faces = find_connected_face(blocks,stator_hub_to_match, outer_faces)
stator_hub3.append(stator_hub_to_match.to_dict())
stator_hub.extend(stator_hub3)

# Stator Shroud
block_id = 0; indices = np.array([0,148,0,44,148,76],dtype=int)
stator_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_shroud, outer_faces = find_connected_face(blocks,stator_shroud_to_match, outer_faces)
stator_shroud.append(stator_shroud_to_match.to_dict())

block_id = 6; indices = np.array([0,148,0,60,148,36],dtype=int) 
stator_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_shroud2, outer_faces = find_connected_face(blocks,stator_shroud_to_match, outer_faces)
stator_shroud2.append(stator_shroud_to_match.to_dict())
stator_shroud.extend(stator_shroud2)

block_id = 9; indices = np.array([0,148,0,68,148,48],dtype=int) 
stator_shroud_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_shroud3, outer_faces = find_connected_face(blocks,stator_shroud_to_match, outer_faces)
stator_shroud3.append(stator_shroud_to_match.to_dict())
stator_shroud.extend(stator_shroud3)

# Mixing Plane
block_id = 5; indices = np.array([36,0,0,36,148,76], dtype=int)
mixing_plane_face_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_mixing_plane_faces, outer_faces = find_connected_face(blocks,mixing_plane_face_to_match, outer_faces)
stator_mixing_plane_faces.append(mixing_plane_face_to_match.to_dict())

# Inlet
block_id = 0; indices = np.array([0,0,0,0,148,76], dtype=int)
inlet_face_to_match = find_face(blocks,block_id, indices,outer_faces)
stator_inlet_faces, outer_faces = find_connected_face(blocks,inlet_face_to_match, outer_faces)
stator_inlet_faces.append(inlet_face_to_match.to_dict())

data['outer_faces'] = outer_faces
data['mixing_plane'] = stator_mixing_plane_faces
data['stator_body'] = stator_faces
data['stator_shroud'] = stator_shroud
data['stator_hub'] = stator_hub
data['inlet'] = stator_inlet_faces

# This is useful for viewing the mesh
with open('stator_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("done")

# Additional Faces that weren't matched by plot3d code. Plot3D gets most things but some are missing
# Part of blade o-mesh
data['face_matches'].append({'block1':{'block_index':8,'IMIN':0,'JMIN':0,'KMIN':0,'IMAX':0,'JMAX':148,'KMAX':48},
                            'block2':{'block_index':13,'IMIN':52,'JMIN':0,'KMIN':0,'IMAX':52,'JMAX':148,'KMAX':48}})

# Part of blade TE (Periodic)
data['face_matches'].append({'block1':{'block_index':8,'IMIN':0,'JMIN':0,'KMIN':48,'IMAX':36,'JMAX':148,'KMAX':48},
                            'block2':{'block_index':7,'IMIN':0,'JMIN':0,'KMIN':0,'IMAX':0,'JMAX':148,'KMAX':36}})


with open('stator_split_connectivity_final.pickle','wb') as f:
    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print("done")
