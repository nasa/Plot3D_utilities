from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, periodicity, periodicity_fast
from plot3d import find_matching_blocks, get_outer_faces, connectivity_fast
from glennht_con import export_to_glennht_conn
import pickle

# Convert to binary because of size 

if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('vspt-ASCII.xyz', binary = False)
    write_plot3D('vspt.xyz', blocks,binary = True)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity_fast(blocks)
    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

blocks = read_plot3D('vspt.xyz', binary = True)

periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity_fast(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)


with open('connectivity_periodic.pickle','wb') as f:
    # [m.pop('match',None) for m in face_matches] # Remove the dataframe
    pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces},f)

# Append periodic surfaces to face_matches
face_matches.extend(periodic_surfaces)

export_to_glennht_conn(face_matches,outer_faces_to_keep,'vspt')


# write_plot3D('finalmesh-paht_binary.xyz',blocks=blocks,binary=True)
# write_plot3D('finalmesh-paht_ascii.xyz',blocks=blocks,binary=False)