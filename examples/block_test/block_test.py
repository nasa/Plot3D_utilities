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

blocks = read_plot3D('iso65_64blocks.xyz',binary=True,read_double=False)
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


#%% Use bounding faces to find periodicity
data = read_data()
z_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks,all_faces,translational_direction='z')
x_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks,outer_faces,translational_direction='x')
y_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks,outer_faces,translational_direction='y')
data['x_periodic'] = x_periodic_faces_export
data['z_periodic'] = z_periodic_faces_export
data['y_periodic'] = y_periodic_faces_export
data['outer_faces'] = outer_faces
dump_data(data)

# Export to GlennHT Format
import copy
data = read_data()
matched_faces = copy.deepcopy(data['x_periodic'])
matched_faces.extend(copy.deepcopy(data['y_periodic']))
matched_faces.extend(copy.deepcopy(data['z_periodic']))
matched_faces.extend(copy.deepcopy(data['face_matches']))

# Filter outer faces
outer_faces = data['outer_faces']
outer_faces = outer_face_dict_to_list(blocks,outer_faces)
matched_faces = list(set(match_faces_dict_to_list(blocks,matched_faces)))

outer_faces = [o.to_dict() for o in outer_faces if o not in matched_faces]
data['outer_faces'] = outer_faces
dump_data(data)

print(f"Number of outer_faces: {len(outer_faces)}")