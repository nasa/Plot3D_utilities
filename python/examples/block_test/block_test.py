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
write_plot3D("dominick.xyz",[blocks[0]],False)
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
backward_bound,forward_bound,_,_ = find_bounding_faces(blocks,connectivity_matrix,all_faces,"x")
lower_bound, upper_bound,_,_ = find_bounding_faces(blocks,connectivity_matrix,all_faces,"z")
left_bound, right_bound,_,_ = find_bounding_faces(blocks,connectivity_matrix,all_faces,"y")
data['forward_bound'] = forward_bound
data['backward_bound'] = backward_bound
data['lower_bound'] = lower_bound
data['upper_bound'] = upper_bound
data['left_bound'] = left_bound
data['right_bound'] = right_bound
dump_data(data)

#%% Use bounding faces to find periodicity
data = read_data()
forward_bound = data['forward_bound']; backward_bound = data['backward_bound']
lower_bound = data['lower_bound']; upper_bound = data['upper_bound']
left_bound = data['left_bound']; right_bound = data['right_bound']
x_periodic_faces_export, periodic_faces = translational_periodicity(blocks,backward_bound,forward_bound,translational_direction='x')
y_periodic_faces_export, periodic_faces = translational_periodicity(blocks,left_bound,right_bound,translational_direction='y')
z_periodic_faces_export, periodic_faces = translational_periodicity(blocks,lower_bound,upper_bound,translational_direction='z')
data['x_periodic'] = x_periodic_faces_export
data['z_periodic'] = z_periodic_faces_export
data['y_periodic'] = y_periodic_faces_export
dump_data(data)

# Export to GlennHT Format
data = read_data()
matched_faces = data['x_periodic']
matched_faces.extend(data['y_periodic'])
matched_faces.extend(data['z_periodic'])
matched_faces.extend(data['face_matches'])

# Filter outer faces
outer_faces = data['outer_faces']
outer_faces = outer_face_dict_to_list(blocks,outer_faces)
matched_faces = list(set(match_faces_dict_to_list(blocks,matched_faces)))

outer_faces = [o.to_dict() for o in outer_faces if o not in matched_faces]
data['outer_faces'] = outer_faces
dump_data(data)

print(f"Number of outer_faces: {len(outer_faces)}")