"""Plot3D multi-block structured grid toolkit.

Provides readers, writers, connectivity/periodicity matching, block splitting,
graph partitioning, and boundary-condition export for Plot3D (``.p3d``,
``.xyz``) mesh files.  Supports binary, ASCII, and Fortran-unformatted formats.

Subpackages:
    glennht: GlennHT boundary-condition and job-file export.
    gridpro: GridPro grid and connectivity importers.
    pointwise: Pointwise ``.inp`` / ``.fvbnd`` importers and GlennHT export.

Attributes:
    use_single_precision (bool): When ``True``, all :class:`Block` coordinate
        arrays are stored as ``float32`` instead of ``float64``.  Defaults to
        ``False``.
"""
from __future__ import absolute_import

use_single_precision = False

from .block import Block, compute_gcd, reduce_blocks
from .blockfunctions import rotate_block, get_outer_bounds, block_connection_matrix, plot_blocks, find_matching_faces
from .block_merging_mixed_facepairs import combine_nxnxn_cubes_mixed_pairs
from .connectivity import find_matching_blocks, get_face_intersection, connectivity, connectivity_fast, face_matches_to_dict, verify_connectivity
from .face import Face
from .face_pool import FacePool
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_bounding_faces, find_angular_bounding_faces, split_face, find_face_nearest_point, match_faces_dict_to_list, outer_face_dict_to_list, find_closest_block
from .read import read_plot3D, read_ap_nasa
from .write import write_plot3D
from .differencing import find_edges, find_face_edges
from .periodicity import create_rotation_matrix, rotated_periodicity, translational_periodicity, verify_periodicity, periodicity, periodicity_fast
from .point_match import point_match
from .split_block import split_blocks, Direction
from .listfunctions import unique_pairs

from .graph import write_ddcmp, build_weighted_graph_from_face_matches, csr_from_adj_and_weights, partition_from_face_matches
