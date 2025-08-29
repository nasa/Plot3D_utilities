from plot3d.glennht.export_functions import export_to_glennht_conn, export_to_boundary_condition, export_to_job_file
from plot3d.glennht.class_definitions import *
from plot3d.gridpro import read_gridpro_to_blocks, read_gridpro_connectivity
import glob
import os.path as osp
from plot3d import translational_periodicity, write_plot3D, write_ddcmp,partition_from_face_matches
import pickle
import json
from pathlib import Path

if __name__ == "__main__":      
    root_path = "[path to gridpro export]"  
    folders = [ 
                f"{root_path}/Quartermillion",
                f"{root_path}/Halfmillion",
                f"{root_path}/OneMillion",
            ]
    # Loop through the folders 
    for folder in folders:
        connectivity_file = osp.join(folder,'blk.tmp.tmp.tmp.conn_n')
        gridpro_file = osp.join(folder,'blk.tmp.tmp.tmp')
        
        # Grid Grid Pro to Blocks 
        blocks = read_gridpro_to_blocks(gridpro_file)
        connectivity = read_gridpro_connectivity(connectivity_file)
        outer_faces = connectivity['bc_group']['symm_slip'] # type: ignore
        
        # Add Translational Periodicity: symmetric slip faces should be periodic to each other now. 
        z_periodic_faces_export, periodic_faces, outer_faces = translational_periodicity(blocks,outer_faces,translational_direction='z') # type: ignore
        connectivity['bc_group']['symm_slip'] = outer_faces # type: ignore # Add unmatched faces to symmetric slip

        # Write connectivity.json for debugging purposes
        with open(osp.join(folder,"connectivity.json"), "w") as json_file:
            connectivity['periodicity'] = z_periodic_faces_export
            json.dump(connectivity, json_file, indent=4) # type: ignore
            
        for d in z_periodic_faces_export:
            del d['mapping']
            del d['mode']
        connectivity['face_matches'].extend(z_periodic_faces_export) # type: ignore
        
        # Define the boundary conditions 
        inlets = []
        for i,inlet in enumerate(connectivity['bc_group']['inlet']): # type: ignore
            inlets.append(InletBC(
                BCType=BoundaryConditionType.Inlet, SurfaceID=inlet['id'], Name=f"Inlet-{i}",
                P0_const=60.0, P0_const_unit="bar", T0_const=300.0, inlet_subType=inlet['id'],
            ))          # Use absolute conditions, code will automatically normalize
        outlets = []
        for i,outlet in enumerate(connectivity['bc_group']['outlet']): # type: ignore
            outlets.append(OutletBC(
                BCType=BoundaryConditionType.Outlet, SurfaceID=outlet['id'], Name=f"Outlet-{i}",
                Pback_const=1.0, Pback_const_unit="bar", outlet_subType=outlet['id']
            ))          # Use absolute conditions, code will automatically normalize
        walls = [] 
        for i,wall in enumerate(connectivity['bc_group']['wall']): # type: ignore
            walls.append(WallBC(
                BCType=BoundaryConditionType.Wall, SurfaceID=wall['id'], Name=f"Wall-{i}", wall_subType=wall['id']
            ))
            
        bcg = BCGroup(Inlets=inlets, Outlets=outlets, SymmetricSlips=[], Walls=walls)

        job = Job()
        job.JobFiles.ConnFILE = "connectivity.ght_conn" # Name of connectivity file
        job.JobFiles.DcmpFILE = "ddcmp.dat"             # Name of ddcmp file
        job.JobFiles.GridFile = "blk.xyz"               # Name of grid pro mesh file to plot3d mesh file
        job.JobFiles.BCSpecFILE = 'boundary_conditions.bcs' # Name of boundary condition file
        job.ReferenceCondFull.reflen = 5.119                # inches, used to calculate the reynolds number should be same scale as mesh
        
        outer_faces = []
        for k,v in connectivity['bc_group'].items():
            outer_faces.extend(v)        # Export GlennHT connectivity file
        export_to_glennht_conn(connectivity['face_matches'], outer_faces, osp.join(folder, job.JobFiles.ConnFILE),connectivity['gif_faces'],connectivity['volume_zones']) # type: ignore
        
        # Export GlennHT boundary conditions file
        export_to_boundary_condition(file_path_to_write=osp.join(folder, job.JobFiles.BCSpecFILE),
                                     job_settings= job, 
                                     bc_group=bcg, 
                                     gifs=connectivity['gif_faces'], # type: ignore
                                     volume_zones=connectivity['volume_zones']) # type: ignore
        
        # Export job file
        p = Path(folder)
        last_folder = p.name
        export_to_job_file(job,osp.join(folder,'job'),title=last_folder)
        
        # Export DDCMP
        nprocessors = 10
        parts, adj_list, edge_w = partition_from_face_matches(connectivity['face_matches'],blocks,nprocessors,favor_blocksize=True) # type: ignore
        write_ddcmp(parts,blocks,adj_list,edge_w,osp.join(folder,job.JobFiles.DcmpFILE))
        
        # Export the plot3d
        write_plot3D(osp.join(folder,job.JobFiles.GridFile),blocks,binary=False)
        