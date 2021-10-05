from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, periodicity, split_blocks, Direction,get_face_intersection
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

'''
    # TODO: Connectivity Need to check if full facematch can help speed up matching by providing a starting I,J,K and ending I,J,K so that we find all datapoints
    # We also need to make sure matches that are not full face match have more than 4 points  

'''
def find_connectivity():
    if not os.path.exists('connectivity-block-split.pickle'):
        blocks = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=False)
        blocks_split = split_blocks(blocks,300000, direction=Direction.i)
        write_plot3D('finalmesh_split.xyz',blocks_split,binary=True)
        # Note: Block splits may not be exactly matching with each other so we have to run the connecitvity code again 
        face_matches, outer_faces_formatted = connectivity(blocks_split)
        with open('connectivity-block-split.pickle','wb') as f:
            pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

def find_periodicity():
    with open('connectivity-block-split.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

    blocks = read_plot3D('finalmesh_split.xyz', binary = True, big_endian=False)

    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)
    face_matches.extend(periodic_surfaces)

    with open('connectivity-block-split_v02.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
            "face_matches":face_matches, 
            "periodic_faces":periodic_surfaces,
            "outer_faces":outer_faces_to_keep       
            },f)

    export_to_glennht_conn(face_matches,outer_faces_to_keep,'finalmesh-block-split')

def debug():
    with open('connectivity-block-split_v02.pickle','rb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
            "face_matches":face_matches, 
            "periodic_faces":periodic_faces,
            "outer_faces":outer_faces_to_keep,
            "outer_faces_debug":outer_faces,
            },f)

    with open('connectivity-block-split_v02.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({
            "face_matches":face_matches, 
            "periodic_faces":periodic_faces,
            "outer_faces":outer_faces_to_keep,
            "outer_faces_debug":outer_faces,
            },f)

if __name__=="__main__":
    find_connectivity()
    find_periodicity()