import sys, os, pickle
import numpy as np
sys.path.insert(0,'../../')
from plot3d import read_plot3D, connectivity_fast,translational_periodicity, write_plot3D, Direction, split_blocks, block_connection_matrix, find_bounding_faces
from plot3d import outer_face_dict_to_list, match_faces_dict_to_list

#%% Find connectivity 
def dump_data(data):
    with open('block_data.pickle','wb') as f:
        pickle.dump(data,f)

def read_data():
    with open('block_data.pickle','rb') as f:
        return pickle.load(f)

blocks = read_plot3D('iso65_64blocks.xyz',True)

if not os.path.exists(f'block_data.pickle'):    
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
    data["connectivity_matrix"]=c
    dump_data(data)

data = read_data()    
all_faces = data['all_faces']
connectivity_matrix = data['connectivity_matrix']

#%% Find bounding Faces
lower_bound, upper_bound,_,_ = find_bounding_faces(blocks,connectivity_matrix,all_faces,"z")
left_bound, right_bound,_,_ = find_bounding_faces(blocks,connectivity_matrix,all_faces,"y")
data['lower_bound'] = lower_bound
data['upper_bound'] = upper_bound
data['left_bound'] = left_bound
data['right_bound'] = right_bound
dump_data(data)

#%% Use bounding faces to find periodicity
data = read_data()
lower_bound = data['lower_bound']; upper_bound = data['upper_bound']
left_bound = data['left_bound']; right_bound = data['right_bound']
y_periodic_faces_export, periodic_faces = translational_periodicity(blocks,left_bound,right_bound,translational_direction='y')
z_periodic_faces_export, periodic_faces = translational_periodicity(blocks,lower_bound,upper_bound,translational_direction='z')
data['z_periodic'] = z_periodic_faces_export
data['y_periodic'] = y_periodic_faces_export
dump_data(data)

#%%  Lets check with faces are not periodic in the y-direction
data = read_data()
y_periodic_faces_export = data['y_periodic']
left_bound = data['left_bound']; right_bound = data['right_bound']

left_periodic_blocks_found = [p['block1']['block_index'] for p in y_periodic_faces_export]
left_faces_missing = [l for l in left_bound if l['block_index'] not in left_periodic_blocks_found]
print('Left faces missing')
[print(l) for l in left_faces_missing]

right_periodic_blocks_found = [p['block2']['block_index'] for p in y_periodic_faces_export]
right_faces_missing = [r for r in right_bound if r['block_index'] not in right_periodic_blocks_found]
print('Right faces missing')
[print(r) for r in right_faces_missing]

# Lets find out why it's missing