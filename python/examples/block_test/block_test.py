import os, sys
sys.path.insert(0,'../../')
sys.path.insert(0,'../Cascade')
from plot3d import write_plot3D, read_plot3D, periodicity
from plot3d import find_matching_blocks, get_outer_faces, connectivity, connectivity_fast
from glennht_con import export_to_glennht_conn
import pickle


if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('connectivity_test.x', binary = False)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity_fast(blocks)
    # face_matches2, outer_faces_formatted2 = connectivity(blocks)
    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

export_to_glennht_conn(face_matches,outer_faces,'block_test')
