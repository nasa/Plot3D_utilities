from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, find_periodicity
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

# Convert to binary because of size 
# blocks = read_plot3D('../../../testfiles/finalmesh-ASCII.xyz', binary = False)
# write_plot3D('../../../testfiles/finalmesh.xyz',blocks,binary=True)
# blocks2 = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=True)

if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=False)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity(blocks)
    with open('connectivity.pickle','wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

blocks = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=False)
periodic_surfaces, outer_faces_to_keep = find_periodicity(blocks,outer_faces,periodic_direction='k')
with open('connectivity_periodic.pickle','wb') as f:
    pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces},f)

# Append periodic surfaces to face_matches
face_matches.extend(periodic_surfaces)
export_to_glennht_conn(face_matches,outer_faces_to_keep,'finalmesh')


# write_plot3D('finalmesh-paht_binary.xyz',blocks=blocks,binary=True)
# write_plot3D('finalmesh-paht_ascii.xyz',blocks=blocks,binary=False)