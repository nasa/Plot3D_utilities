from __future__ import absolute_import
# Boundary condition related stuff 
from .class_definitions import BoundaryConditionType,InletBC_Subtype,InletBC_Direction,OutletSubtype,SymmetricSlipSubtype,WallSubtype,GIFCoordinate,GIFType,GIFOrder,TbModelType,BoundaryCondition,InletBC,OutletBC,SymmetricSlipBC,WallBC,GIF,BCGroup

# Job related stuff
from .class_definitions import JobFiles, JobControl, TurbModelInput, Plot3DParameters, InitialCond, TimeStpControl, SPDSchemeControl, RKSchemeControl, MGSchemeControl, GasPropertiesInput, ReferenceCondFull, Job

# export functions
from .export_functions import export_to_boundary_condition, export_to_job_file

# glennht-gpu connectivity.json export
from .gpu import (
    write_connectivity_json,
    tag_surfaces_from_diagonals,
    tag_surfaces_from_bc_codes,
    tag_surfaces_geometric,
    write_bc_codes_json,
    merge_connectivity_json,
    export_to_glennht_gpu,
)

# glennht-gpu boundary-condition dataclasses + YAML fragment writer
from .gpu_boundary_conditions import (
    GpuInletBC,
    GpuOutletBC,
    GpuWallBC,
    write_boundary_conditions_yaml,
)

# glennht-gpu single-file HDF5 bundle format (.graph_p3d). h5py is an
# OPTIONAL dependency -- it is only imported lazily inside graph_p3d.py's
# own functions, so this import stays safe even when h5py isn't installed.
from .graph_p3d import (
    write_graph_p3d,
    read_graph_p3d,
    graph_p3d_to_files,
)

# validation functions
from .validation import (
    check_handedness,
    check_corners,
    check_verify_connectivity,
    check_verify_periodicity,
    compute_pm,
    check_pm,
    check_face_coverage,
    check_negative_volumes,
    check_face_normals,
    run_all_checks,
)