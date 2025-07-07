from __future__ import absolute_import
from importlib import import_module     
import os, warnings

from .block import Block
from .blockfunctions import rotate_block, get_outer_bounds, block_connection_matrix,split_blocks, plot_blocks, reduce_blocks, find_matching_faces
from .block_merging_mixed_facepairs import combine_nxnxn_cubes_mixed_pairs
from .connectivity import find_matching_blocks, get_face_intersection, connectivity_fast, face_matches_to_dict
from .face import Face
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_connected_faces, find_bounding_faces,split_face,find_face_nearest_point,match_faces_dict_to_list,outer_face_dict_to_list,find_closest_block
from .read import read_plot3D, read_ap_nasa
from .write import write_plot3D
from .differencing import find_edges, find_face_edges
from .periodicity import periodicity, periodicity_fast, create_rotation_matrix, rotated_periodicity, translational_periodicity
from .point_match import point_match
from .split_block import split_blocks, Direction
from .listfunctions import unique_pairs

# Try importing metis
if os.getenv('METIS_DLL') is not None:
    if import_module('metis') is not None:
        import metis
        from .graph import block_to_graph,get_face_vertex_indices,get_starting_vertex,add_connectivity_to_graph, block_connectivity_to_graph
else:
    print("METIS_DLL is not set. metis may not be configured. plot3D will function without metis")