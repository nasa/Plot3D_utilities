import sys, os, pickle 
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast, translational_periodicity, write_plot3D, Direction, split_blocks, block_connection_matrix

blocks = read_plot3D('CMC009_fine_binary.xyz',True)
if not os.path.exists(f'CMC009_connectivity.pickle'):
    c,ic,jc,kc = block_connection_matrix(blocks)
    face_matches, outer_faces = connectivity_fast(blocks)
    with open('CMC009_connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump(
            {
                "face_matches":face_matches, 
                "outer_faces":outer_faces,
                "connection_matrix":c,
                "connection_matrix_i":ic,
                "connection_matrix_j":jc,
                "connection_matrix_k":kc,
            },f)

    
with open('CMC009_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']
    connection_matrix = data['connection_matrix']
    connection_matrix_i = data['connection_matrix_i']
    connection_matrix_j = data['connection_matrix_j']
    connection_matrix_k = data['connection_matrix_k']

# We need faces connected in I and J at KMIN
periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=z_shift_distance, shift_direction='z')
face_matches.extend(periodic_faces)

periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=y_shift_distance,shift_direction='y')
face_matches.extend(periodic_faces)



print('done')