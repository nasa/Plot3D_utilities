import sys, os, pickle 
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast, translational_periodicity, write_plot3D, Direction, split_blocks, block_connection_matrix, outer_face_dict_to_list, match_faces_dict_to_list

def dump_data(data):
    with open('CMC009_connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump(data
            ,f)

def read_data():
    with open('CMC009_connectivity.pickle','rb') as f:
        return pickle.load(f)
   
    
blocks = read_plot3D('CMC009_fine_binary.xyz',True)
if not os.path.exists(f'CMC009_connectivity.pickle'):    
    print('Finding connectivity')
    face_matches, outer_faces = connectivity_fast(blocks)
    print('Organizing split and outerfaces')
    matched_faces = match_faces_dict_to_list(blocks,face_matches)
    matched_faces.extend(outer_face_dict_to_list(blocks,outer_faces))
    matched_faces = [m.to_dict() for m in matched_faces]
    print('Creating block connection matrix')
    c,ic,jc,kc = block_connection_matrix(blocks,matched_faces)
    data = {
                "face_matches":face_matches, 
                "outer_faces":outer_faces,
                "connection_matrix":c,
                "connection_matrix_i":ic,
                "connection_matrix_j":jc,
                "connection_matrix_k":kc,
            }
    dump_data(data) 

data = read_data()    
face_matches = data['face_matches']
outer_faces = data['outer_faces']
connection_matrix = data['connection_matrix']
connection_matrix_i = data['connection_matrix_i']
connection_matrix_j = data['connection_matrix_j']
connection_matrix_k = data['connection_matrix_k']


# We need faces connected in I and J at KMIN
# lower_connected_faces, upper_connected_faces = translational_periodicity(blocks,connection_matrix, outer_faces,translational_direction='z')

# data['lower_connected_faces'] = lower_connected_faces
# data['upper_connected_faces'] = upper_connected_faces
# dump_data(data)

data = read_data()
lower_connected_faces = list(set(data['lower_connected_faces']))
upper_connected_faces = list(set(data['upper_connected_faces']))

# Shift lower connected face by the delta z and check for connectivity
dz = 1.5
lower_connected_faces 
print('done')

# periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=y_shift_distance,shift_direction='y')
# face_matches.extend(periodic_faces)



print('done')