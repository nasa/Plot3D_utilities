"""GridPro grid and connectivity import utilities.

Provides functions to read GridPro structured-grid text files into
:class:`~plot3d.block.Block` objects and to parse GridPro connectivity files
into ``plot3d.connectivity_fast``-compatible face-match dictionaries with
boundary-condition grouping.
"""
from __future__ import absolute_import
from .import_functions import (
    read_gridpro_to_blocks,
    read_gridpro_connectivity,
    bc_faces_by_type,
)
