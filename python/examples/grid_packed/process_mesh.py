from plot3d.glennht.export_functions import export_to_glennht_conn, export_to_boundary_condition, export_to_job_file
from plot3d.glennht.class_definitions import *
from plot3d.gridpro import read_gridpro_to_blocks, read_gridpro_connectivity, bc_faces_by_type
from plot3d import translational_periodicity, write_plot3D, read_plot3D, write_ddcmp,partition_from_face_matches
import glob
import os.path as osp
import pickle
import json
from pathlib import Path
import pandas as pd 

SCRIPT_DIR = Path(__file__).resolve().parent

INLET_IDS = {
    "Inlet": [1005],
    "Pressure Plenum Inlet": [1015],
    "Cavity Inflow": [1016],
    "Suction Plenum Inlet": [1017],
}
OUTLET_IDS = {"Outlet": [1006]}
WALL_IDS = {
    "Downstream Cavity": [1007, 1008],
    "Pressure Plenum Back": [2008],
    "Casing": [2007],
    "Blade": [3008],
    "Blade Top": [4008],
    "Downstream Hub": [5008],
    "Upstream Hub": [3007],
    "Upstream Cavity A": [4007],
    "Upstream Cavity B": [5007],
}
SYMM_SLIP_IDS = {
    "Pressure Plenum Top": [1004],
    "Suction Plenum Top": [2004],
}
#%% Read Connectivity File
connectivity_file = 'examples/grid_packed/grid_packed_binary.tmp.tmp.conn_n'
connectivity = read_gridpro_connectivity(
    connectivity_file,
    inlet_ids=INLET_IDS, # pyright: ignore[reportArgumentType]
    outlet_ids=OUTLET_IDS, # type: ignore
    wall_ids=WALL_IDS, # pyright: ignore[reportArgumentType]
    symm_slip_ids=SYMM_SLIP_IDS, # pyright: ignore[reportArgumentType]
)

#%% Define the boundary conditions 
def face_lookup_for(bc_type: str) -> dict[int, dict]:
    return {face['id']: face for face in bc_faces_by_type(connectivity['bc_group'], bc_type, unique=True)}

def get_face(face_lookup: dict[int, dict], face_id: int, label: str) -> dict:
    face = face_lookup.get(face_id)
    if face is None:
        raise KeyError(f"Missing {label} face id {face_id}")
    return face

inlet_face_lookup = face_lookup_for('inlet')
outlet_face_lookup = face_lookup_for('outlet')
wall_face_lookup = face_lookup_for('wall')
symm_face_lookup = face_lookup_for('symm_slip')

inlet_1005 = get_face(inlet_face_lookup, 1005, 'inlet')
inlet_1015 = get_face(inlet_face_lookup, 1015, 'inlet')
inlet_1016 = get_face(inlet_face_lookup, 1016, 'inlet')
inlet_1017 = get_face(inlet_face_lookup, 1017, 'inlet')

outlet_1006 = get_face(outlet_face_lookup, 1006, 'outlet')

wall_1007 = get_face(wall_face_lookup, 1007, 'wall')
wall_1008 = get_face(wall_face_lookup, 1008, 'wall')
wall_2008 = get_face(wall_face_lookup, 2008, 'wall')
wall_2007 = get_face(wall_face_lookup, 2007, 'wall')
wall_3008 = get_face(wall_face_lookup, 3008, 'wall')
wall_4008 = get_face(wall_face_lookup, 4008, 'wall')
wall_5008 = get_face(wall_face_lookup, 5008, 'wall')
wall_3007 = get_face(wall_face_lookup, 3007, 'wall')
wall_4007 = get_face(wall_face_lookup, 4007, 'wall')
wall_5007 = get_face(wall_face_lookup, 5007, 'wall')

symm_1004 = get_face(symm_face_lookup, 1004, 'symm_slip')
symm_2004 = get_face(symm_face_lookup, 2004, 'symm_slip')

inlets = [
    InletBC(
        BCType=BoundaryConditionType.Inlet,
        SurfaceID=inlet_1005['id'],
        Name=f"{inlet_1005['bc_name']}",
        P0_const=60.0,
        P0_const_unit="bar",
        T0_const=300.0,
        inlet_subType=inlet_1005['id'],
    ),
    InletBC(
        BCType=BoundaryConditionType.Inlet,
        SurfaceID=inlet_1015['id'],
        Name=f"{inlet_1015['bc_name']}",
        P0_const=60.0,
        P0_const_unit="bar",
        T0_const=300.0,
        inlet_subType=inlet_1015['id'],
    ),
    InletBC(
        BCType=BoundaryConditionType.Inlet,
        SurfaceID=inlet_1016['id'],
        Name=f"{inlet_1016['bc_name']}",
        P0_const=60.0,
        P0_const_unit="bar",
        T0_const=300.0,
        inlet_subType=inlet_1016['id'],
    ),
    InletBC(
        BCType=BoundaryConditionType.Inlet,
        SurfaceID=inlet_1017['id'],
        Name=f"{inlet_1017['bc_name']}",
        P0_const=60.0,
        P0_const_unit="bar",
        T0_const=300.0,
        inlet_subType=inlet_1017['id'],
    ),
]

outlets = [
    OutletBC(
        BCType=BoundaryConditionType.Outlet,
        SurfaceID=outlet_1006['id'],
        Name=f"{outlet_1006['bc_name']}",
        Pback_const=1.0,
        Pback_const_unit="bar",
        outlet_subType=outlet_1006['id'],
    ),
]

walls = [
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_1007['id'],
        Name=f"{wall_1007['bc_name']}",
        wall_subType=wall_1007['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_1008['id'],
        Name=f"{wall_1008['bc_name']}",
        wall_subType=wall_1008['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_2008['id'],
        Name=f"{wall_2008['bc_name']}",
        wall_subType=wall_2008['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_2007['id'],
        Name=f"{wall_2007['bc_name']}",
        wall_subType=wall_2007['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_3008['id'],
        Name=f"{wall_3008['bc_name']}",
        wall_subType=wall_3008['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_4008['id'],
        Name=f"{wall_4008['bc_name']}",
        wall_subType=wall_4008['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_5008['id'],
        Name=f"{wall_5008['bc_name']}",
        wall_subType=wall_5008['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_3007['id'],
        Name=f"{wall_3007['bc_name']}",
        wall_subType=wall_3007['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_4007['id'],
        Name=f"{wall_4007['bc_name']}",
        wall_subType=wall_4007['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=wall_5007['id'],
        Name=f"{wall_5007['bc_name']}",
        wall_subType=wall_5007['id'],
    ),
]

symmslips = [
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=symm_1004['id'],
        Name=f"{symm_1004['bc_name']}",
        wall_subType=symm_1004['id'],
    ),
    WallBC(
        BCType=BoundaryConditionType.Wall,
        SurfaceID=symm_2004['id'],
        Name=f"{symm_2004['bc_name']}",
        wall_subType=symm_2004['id'],
    ),
]

bcg = BCGroup(Inlets=inlets, Outlets=outlets, SymmetricSlips=symmslips, Walls=walls)

job = Job()
job.JobFiles.ConnFILE = "connectivity.ght_conn" # Name of connectivity file
job.JobFiles.DcmpFILE = "ddcmp.dat"             # Name of ddcmp file
job.JobFiles.GridFile = "blk.xyz"               # Name of grid pro mesh file to plot3d mesh file
job.JobFiles.BCSpecFILE = 'boundary_conditions.bcs' # Name of boundary condition file
job.ReferenceCondFull.reflen = 5.119                # inches, used to calculate the reynolds number should be same scale as mesh

outer_faces = []
for entry in connectivity['bc_group'].values():
    faces = entry.get('faces', []) if isinstance(entry, dict) else entry
    outer_faces.extend(faces)        # Export GlennHT connectivity file

conn_path = SCRIPT_DIR / job.JobFiles.ConnFILE
bc_path = SCRIPT_DIR / job.JobFiles.BCSpecFILE
job_path = SCRIPT_DIR / 'job'
ddcmp_path = SCRIPT_DIR / job.JobFiles.DcmpFILE

export_to_glennht_conn(
    connectivity['face_matches'],
    outer_faces,
    str(conn_path),
    [],
    connectivity['gif_faces'],
    connectivity['volume_zones'],
) # type: ignore

# Export GlennHT boundary conditions file
export_to_boundary_condition(file_path_to_write=str(bc_path),
                                job_settings=job, 
                                bc_group=bcg, 
                                gif_pairs=connectivity['gif_faces'], # type: ignore
                                volume_zones=connectivity['volume_zones']) # type: ignore

# Export job file
export_to_job_file(job, job_path, title="Grid Packed")

# Export DDCMP
nprocessors = 10
parts, adj_list, edge_w = partition_from_face_matches(connectivity['face_matches'],connectivity['blocksizes'],nprocessors,favor_blocksize=True) # type: ignore
write_ddcmp(parts, connectivity['blocksizes'], adj_list, edge_w, str(ddcmp_path))

#e  Export the plot3d
# write_plot3D(job.JobFiles.GridFile,blocks,binary=False)
        
