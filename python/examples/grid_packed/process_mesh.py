from plot3d.glennht.export_functions import export_to_glennht_conn, export_to_boundary_condition, export_to_job_file
from plot3d.glennht.class_definitions import *
from plot3d.gridpro import read_gridpro_to_blocks, read_gridpro_connectivity, bc_faces_by_type
import glob
import os.path as osp
from plot3d import translational_periodicity, write_plot3D, read_plot3D, write_ddcmp,partition_from_face_matches
import pickle
import json
from pathlib import Path
import pandas as pd 

connectivity_file = 'examples/grid_packed/grid_packed_binary.tmp.tmp.conn_n'
connectivity = read_gridpro_connectivity(connectivity_file,
                                         inlet_ids={
                                             "Inlet":[1005],
                                             "Pressure Plenum Inlet":[1015],
                                             "Cavity Inflow":[1016],
                                             "Suction Plenum Inlet":[1017]}, # pyright: ignore[reportArgumentType]
                                         outlet_ids={"Outlet":[1006]}, # type: ignore
                                         wall_ids={
                                             "Downstream Cavity":[1007,1008],
                                             "Pressure Plenum Back":[2008],
                                             "Casing":[2007],
                                             "Blade":[3008],
                                             "Blade Top":[4008],
                                             "Downstream Hub":[5008],
                                             "Upstream Hub":[3007],
                                             "Upstream Cavity A":[4007],
                                             "Upstream Cavity B":[5007]}, # pyright: ignore[reportArgumentType]
                                         symm_slip_ids={
                                             "Pressure Plenum Top":[1004],
                                             "Suction Plenum Top":[2004]
                                             }) # pyright: ignore[reportArgumentType]
# grid_file = "grid_packed_binary.tmp.tmp.p3d"
# saved_state = "saved_state.pkl"
# if not osp.exists(saved_state):
#     # Grid Grid Pro to Blocks 
#     blocks = read_plot3D(grid_file,binary=True)
#     connectivity = read_gridpro_connectivity(connectivity_file)
#     write_plot3D("grid_packed-ascii.p3d",blocks,binary=False)
# else:
#     data = pickle.load(open(saved_state,'rb'))
#     blocks = data['mesh']
#     connectivity = data['connectivity']

# outer_faces = connectivity['bc_group']['symm_slip']['faces'] # type: ignore

# Add Translational Periodicity: symmetric slip faces should be periodic to each other now. 
# z_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks,outer_faces,translational_direction='z') # type: ignore

# connectivity['bc_group']['symm_slip'] = {
#     "type": "symm_slip",
#     "pty_ids": [],
#     "faces": outer_faces,
# } # type: ignore # Add unmatched faces to symmetric slip


# Define the boundary conditions 
def unique_faces_by_type(bc_type: str):
    return bc_faces_by_type(connectivity['bc_group'], bc_type, unique=True)

inlets = []
for inlet in unique_faces_by_type('inlet'):
    inlets.append(InletBC(
        BCType=BoundaryConditionType.Inlet, SurfaceID=inlet['id'], Name=f"{inlet['bc_name']}",
        P0_const=60.0, P0_const_unit="bar", T0_const=300.0, inlet_subType=inlet['id'],
    ))          # Use absolute conditions, code will automatically normalize
outlets = []
for outlet in unique_faces_by_type('outlet'):
    outlets.append(OutletBC(
        BCType=BoundaryConditionType.Outlet, SurfaceID=outlet['id'], Name=f"{outlet['bc_name']}",
        Pback_const=1.0, Pback_const_unit="bar", outlet_subType=outlet['id']
    ))          # Use absolute conditions, code will automatically normalize
walls = [] 
for wall in unique_faces_by_type('wall'):
    walls.append(WallBC(
        BCType=BoundaryConditionType.Wall, SurfaceID=wall['id'], Name=f"{wall['bc_name']}", wall_subType=wall['id']
    ))

symmslips = [] 
for slip in unique_faces_by_type('symm_slip'):
    symmslips.append(WallBC(
        BCType=BoundaryConditionType.Wall, SurfaceID=slip['id'], Name=f"{slip['bc_name']}", wall_subType=slip['id']
    ))
bcg = BCGroup(Inlets=inlets, Outlets=outlets, SymmetricSlips=[], Walls=walls)

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
export_to_glennht_conn(connectivity['face_matches'], outer_faces, job.JobFiles.ConnFILE,[],connectivity['gif_faces'],connectivity['volume_zones']) # type: ignore

# Export GlennHT boundary conditions file
export_to_boundary_condition(file_path_to_write=job.JobFiles.BCSpecFILE,
                                job_settings= job, 
                                bc_group=bcg, 
                                gif_pairs=connectivity['gif_faces'], # type: ignore
                                volume_zones=connectivity['volume_zones']) # type: ignore

# Export job file
export_to_job_file(job,'job',title="Grid Packed")

# Export DDCMP
nprocessors = 10
parts, adj_list, edge_w = partition_from_face_matches(connectivity['face_matches'],connectivity['blocksizes'],nprocessors,favor_blocksize=True) # type: ignore
write_ddcmp(parts, connectivity['blocksizes'], adj_list, edge_w, job.JobFiles.DcmpFILE)

#e  Export the plot3d
# write_plot3D(job.JobFiles.GridFile,blocks,binary=False)
        
