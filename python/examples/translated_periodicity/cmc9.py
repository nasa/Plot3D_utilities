import sys, os, pickle
import numpy as np
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast,translational_periodicity, translational_periodicity2, write_plot3D, Direction, split_blocks, block_connection_matrix, outer_face_dict_to_list, match_faces_dict_to_list

def dump_data(data):
    with open('cmc9_data.pickle','wb') as f:
        pickle.dump(data,f)

def read_data():
    with open('cmc9_data.pickle','rb') as f:
        return pickle.load(f)
   
blocks = read_plot3D('CMC009_fine_binary.xyz',True)

if not os.path.exists(f'cmc9_data.pickle'):    
    print('Finding connectivity')
    face_matches, outer_faces = connectivity_fast(blocks)
    [m.pop('match',None) for m in face_matches] # Remove the dataframe
    print('Organizing split and outerfaces')
    all_faces = match_faces_dict_to_list(blocks,face_matches)
    all_faces.extend(outer_face_dict_to_list(blocks,outer_faces))
    all_faces = [m.to_dict() for m in all_faces]
    data = {
                "face_matches":face_matches, 
                "outer_faces":outer_faces,
                "all_faces":all_faces
            }
    dump_data(data)
    print('Creating block connection matrix')
    c = block_connection_matrix(blocks,all_faces)
    data["connection_matrix"]=c
    dump_data(data)

data = read_data()    
all_faces = data['all_faces']
connection_matrix = data['connection_matrix']

# We need faces connected in I and J at KMIN
lower_connected_faces, upper_connected_faces = translational_periodicity(blocks,connection_matrix,all_faces,translational_direction='z')

data['lower_connected_faces'] = lower_connected_faces
data['upper_connected_faces'] = upper_connected_faces
dump_data(data)

data = read_data()
lower_connected_faces = data['lower_connected_faces']
upper_connected_faces = data['upper_connected_faces']
connection_matrix = data['connection_matrix']

# Shift lower connected face by the delta z and check for connectivity
# translational_periodicity2(blocks,lower_connected_faces,upper_connected_faces,direction="z")
# periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=y_shift_distance,shift_direction='y')
# face_matches.extend(periodic_faces)



print('done')