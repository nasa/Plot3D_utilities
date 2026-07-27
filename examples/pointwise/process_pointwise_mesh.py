'''
Example of converting pointwise to GlennHT job and boundary condition files
'''

from pathlib import Path

from plot3d import read_plot3D, partition_from_face_matches, write_ddcmp
from plot3d.glennht.class_definitions import (
    BCGroup,
    BoundaryConditionType,
    InletBC,
    Job,
    OutletBC,
    SymmetricSlipBC,
    WallBC,
)
from plot3d.glennht.export_functions import (
    export_to_boundary_condition,
    export_to_glennht_conn,
    export_to_job_file,
)
from plot3d.gridpro import bc_faces_by_type
from plot3d.pointwise import build_pointwise_connectivity
import os.path as osp 
import pickle

SCRIPT_DIR = Path(__file__).resolve().parent

PLOT3D_FILE = SCRIPT_DIR / "[Insert filename for Plot3d ASCII]"
FVBND_FILE = SCRIPT_DIR / "[fvbnd filename]"
INP_FILE = SCRIPT_DIR / "[inp filename]"
SAVED_STATE = SCRIPT_DIR / "saved_state.pickle"

# Build connectivity + BC grouping from Pointwise inputs
if not osp.exists(SAVED_STATE): 
    blocks = read_plot3D(str(PLOT3D_FILE), binary=False)
    connectivity = build_pointwise_connectivity(
        blocks,
        FVBND_FILE,
        inp_path=INP_FILE,
    )
    pickle.dump({"blocks":blocks,"connectivity":connectivity}, open(SAVED_STATE,'wb'))


data = pickle.load(open(SAVED_STATE,'rb'))
blocks = data['blocks']
connectivity = data['connectivity']

def bc_list(bc_type: str):
    """Return faces belonging to a specific boundary condition type.

    Args:
        bc_type: Boundary condition type key expected by
            ``bc_faces_by_type`` (for example ``"inlet"`` or ``"wall"``).

    Returns:
        List of faces matching the requested boundary condition type with
        duplicate entries removed.
    """
    return bc_faces_by_type(connectivity["bc_group"], bc_type, unique=True)

# Define boundary conditions (one per unique face id)
inlets = [
    InletBC(
        BCType=BoundaryConditionType.Inlet,
        SurfaceID=face["id"],
        Name=f"{face['bc_name']}",
        inlet_subType=face["id"],
    )
    for face in bc_list("inlet")
]

outlets = [
    OutletBC(
        BCType=BoundaryConditionType.Outlet,
        SurfaceID=face["id"],
        Name=f"{face['bc_name']}",
        outlet_subType=face["id"],
    )
    for face in bc_list("outlet")
]

symmslips = [
    SymmetricSlipBC(
        BCType=BoundaryConditionType.SymmetryOrSlip,
        SurfaceID=face["id"],
        Name=f"{face['bc_name']}",
        slip_subType=face["id"],
    )
    for face in bc_list("symm_slip")
]

walls = [
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=face["id"],
        Name=f"{face['bc_name']}",
        wall_subType=face["id"],
    )
    for face in bc_list("wall")
]

bcg = BCGroup(Inlets=inlets, Outlets=outlets, SymmetricSlips=symmslips, Walls=walls)

# Job setup
job = Job()
job.JobFiles.ConnFILE = "connectivity.ght_conn"
job.JobFiles.DcmpFILE = "ddcmp.dat"
job.JobFiles.GridFile = PLOT3D_FILE.name
job.JobFiles.BCSpecFILE = "boundary_conditions.bcs"
job.ReferenceCondFull.reflen = 1.0

# Export GlennHT connectivity, BCs, job, and DDCMP
outer_faces = []
for entry in connectivity["bc_group"].values():
    faces = entry.get("faces", []) if isinstance(entry, dict) else entry
    outer_faces.extend(faces)

conn_path = SCRIPT_DIR / job.JobFiles.ConnFILE
bc_path = SCRIPT_DIR / job.JobFiles.BCSpecFILE
job_path = SCRIPT_DIR / "job"
ddcmp_path = SCRIPT_DIR / job.JobFiles.DcmpFILE

export_to_glennht_conn(
    connectivity["face_matches"],
    outer_faces,
    str(conn_path),
    connectivity.get("gif_pairs", []),
    connectivity["gif_faces"],
    connectivity["volume_zones"],
)

export_to_boundary_condition(
    file_path_to_write=str(bc_path),
    job_settings=job,
    bc_group=bcg,
    gif_pairs=connectivity["gif_faces"],
    volume_zones=connectivity["volume_zones"],
)

export_to_job_file(job, job_path, title="Pointwise Mesh")

nprocessors = 4
parts, adj_list, edge_w = partition_from_face_matches(
    connectivity["face_matches"], connectivity["blocksizes"], nprocessors, favor_blocksize=True
)
write_ddcmp(parts, connectivity["blocksizes"], adj_list, edge_w, str(ddcmp_path))

print("Wrote:", conn_path)
print("Wrote:", bc_path)
print("Wrote:", job_path)
print("Wrote:", ddcmp_path)
