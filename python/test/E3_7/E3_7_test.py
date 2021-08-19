from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
sys.path.insert(1,'../Cascade')
from plot3d import write_plot3D, read_plot3D
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('../../../testfiles/E3_7_Assembly.p3d', binary = False)
    # write_plot3D('../../../testfiles/E3_7_Assembly.xyz',blocks=blocks,binary=True)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity(blocks)
    with open('connectivity.pickle','wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

export_to_glennht_conn(face_matches,outer_faces,'E3_7_Assembly')


# write_plot3D('finalmesh-paht_binary.xyz',blocks=blocks,binary=True)
# write_plot3D('finalmesh-paht_ascii.xyz',blocks=blocks,binary=False)