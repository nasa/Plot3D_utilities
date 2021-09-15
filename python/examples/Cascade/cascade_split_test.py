from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, find_periodicity, split_blocks, Direction
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

if not os.path.exists('connectivity-block-split.pickle'):
    blocks = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=False)
    blocks_split = split_blocks(blocks,400000, direction=Direction.i)
    write_plot3D('finalmesh_split.xyz',blocks_split,binary=True)
    # Note: Block splits may not be exactly matching with each other so we have to run the connecitvity code again 
    face_matches, outer_faces_formatted = connectivity(blocks_split)
    with open('connectivity-block-split.pickle','wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

with open('connectivity-block-split.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

blocks = read_plot3D('finalmesh_split.xyz', binary = True, big_endian=True)
periodic_surfaces, outer_faces_to_keep = find_periodicity(blocks,outer_faces,periodic_direction='k')
with open('connectivity-block-split.pickle','wb') as f:
    pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces},f)

# Append periodic surfaces to face_matches
face_matches.extend(periodic_surfaces)

export_to_glennht_conn(face_matches,outer_faces_to_keep,'finalmesh-block-split')

