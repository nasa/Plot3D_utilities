from itertools import combinations
import os, sys
sys.path.insert(0,'../../')
from plot3d import write_plot3D, read_plot3D, find_periodicity, split_blocks, Direction
from plot3d import find_matching_blocks, get_outer_faces, connectivity
from glennht_con import export_to_glennht_conn
import pickle

blocks = read_plot3D('../../../testfiles/finalmesh.xyz', binary = True, big_endian=False)
blocks_split = split_blocks(blocks,300000, direction=Direction.i)

print('done')
