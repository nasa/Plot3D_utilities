import sys
import os
import pickle
from typing import List, Dict
sys.path.insert(0,os.getcwd()) # This allows you to select files locally
from pv_library import Load, ExtractBlocks, CreateSubset
from paraview.simple import *
import random

def CheckDictionary(data:Dict[str,List],name:str):
    """Checks for key in dictionary and returns empty list if key is not found 

    Args:
        data (Dict[str,List]): _description_
        name (str): name of key in dictionary 

    Returns:
        _type_: _description_
    """
    if name in data:
        print(f'{name} found')
        return data[name]
    else:
        return list() 


'''
Main Code
'''

if __name__=="__main__":
    '''
    Read the connectivity file
    '''
    plot3d_filename = 'iso65_64blocks.xyz'
    with open('block_data.pickle','rb') as f:
        data = pickle.load(f)
        outer_faces = CheckDictionary(data,'outer_faces')
        x_periodic_faces = CheckDictionary(data,'x_periodic')
        y_periodic_faces = CheckDictionary(data,'y_periodic')
        z_periodic_faces = CheckDictionary(data,'z_periodic')
        face_matches = CheckDictionary(data,'face_matches')


    blocks_to_extract = [o['block_index'] for o in outer_faces]
    blocks_to_extract.extend([f['block1']['block_index'] for f in x_periodic_faces])
    blocks_to_extract.extend([f['block2']['block_index'] for f in x_periodic_faces])

    blocks_to_extract.extend([f['block1']['block_index'] for f in y_periodic_faces])
    blocks_to_extract.extend([f['block2']['block_index'] for f in y_periodic_faces])

    blocks_to_extract.extend([f['block1']['block_index'] for f in z_periodic_faces])
    blocks_to_extract.extend([f['block2']['block_index'] for f in z_periodic_faces])
    blocks_to_extract = list(set(blocks_to_extract))
    blocks_to_extract.sort()
    n = len(blocks_to_extract)

    
    '''
    Generate Random Colors 
    
    This part generates random colors so that each color is associated with a face match. 
    Doesn't matter what the block is a single match is assigned the same color. 
    '''
    rgb_outer_faces = list()
    for i in range(len(outer_faces)):
        rgb_outer_faces.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    rgb_periodic_x = list()
    for i in range(len(x_periodic_faces)):
        rgb_periodic_x.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    rgb_periodic_y = list()
    for i in range(len(y_periodic_faces)):
        rgb_periodic_y.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    rgb_periodic_z = list()
    for i in range(len(z_periodic_faces)):
        rgb_periodic_z.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    rgb_face_matches = list()
    for i in range(len(face_matches)):
        rgb_face_matches.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    # Load mesh
    plot3D_source,plot3D_Display,View,LUT = Load(plot3d_filename,True)
    print(f"Total number of blocks: {n}")
    def check_and_swap(ijkmin, ijkmax):
        if (ijkmin> ijkmax):
            temp = ijkmax
            ijkmax = ijkmin
            ijkmin = temp
        return ijkmin, ijkmax
    
    for b in range(10):# blocks_to_extract: # Block indicies 
        block_source,block_display,LUT = ExtractBlocks(plot3D_source,View,[b])
        RenameSource('Block '+str(b), block_source)
        block_source = FindSource('Block '+str(b))
        
        for surface_indx, o in enumerate(outer_faces):
            # Add Plots for Outer Faces
            if o['block_index'] == b:
                voi = [o['IMIN'], o['IMAX'], o['JMIN'], o['JMAX'],o['KMIN'], o['KMAX']]
                CreateSubset(block_source, voi, name='outer_face '+str(surface_indx),opacity=1)
        
        # Plot the periodic faces  
        for periodic_indx, p in enumerate(x_periodic_faces):
            print("Plotting x-periodic")
            # Add Plots for Outer Faces
            if p['block1']['block_index'] == b and p['block2']['block_index'] == b: # Periodicity within the block 
                p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                CreateSubset(block_source, voi, name='x-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_x[periodic_indx])

                p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='x-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_x[periodic_indx])

            elif p['block1']['block_index'] == b or p['block2']['block_index'] == b: # Periodicity from block to block 
                if p['block1']['block_index'] == b:
                    p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                    p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                    p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                    voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                else:
                    p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                    p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                    p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                    voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='x-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_x[periodic_indx])
        
        # Plot the periodic faces  
        for periodic_indx, p in enumerate(y_periodic_faces):
            # Add Plots for Outer Faces
            if p['block1']['block_index'] == b and p['block2']['block_index'] == b: # Periodicity within the block 
                p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                CreateSubset(block_source, voi, name='y-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_y[periodic_indx])

                p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='y-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_y[periodic_indx])

            elif p['block1']['block_index'] == b or p['block2']['block_index'] == b: # Periodicity from block to block 
                if p['block1']['block_index'] == b:
                    p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                    p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                    p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                    voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                else:
                    p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                    p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                    p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                    voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='y-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_y[periodic_indx])

        # Plot the periodic faces  
        for periodic_indx, p in enumerate(z_periodic_faces):
            # Add Plots for Outer Faces
            if p['block1']['block_index'] == b and p['block2']['block_index'] == b: # Periodicity within the block 
                p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                CreateSubset(block_source, voi, name='z-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_z[periodic_indx])

                p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='z-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_z[periodic_indx])

            elif p['block1']['block_index'] == b or p['block2']['block_index'] == b: # Periodicity from block to block 
                if p['block1']['block_index'] == b:
                    p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                    p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                    p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                    voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                else:
                    p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                    p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                    p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                    voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='z-periodic '+str(periodic_indx),rgb_face_matches=rgb_periodic_z[periodic_indx])
        
        # Plot the periodic faces  
        for idx, p in enumerate(face_matches):
            # Add Plots for Outer Faces
            if p['block1']['block_index'] == b and p['block2']['block_index'] == b: # Periodicity within the block 
                p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                CreateSubset(block_source, voi, name='face_match '+str(idx),rgb_face_matches=rgb_face_matches[idx])

                p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='face_match '+str(idx),rgb_face_matches=rgb_face_matches[idx])

            elif p['block1']['block_index'] == b or p['block2']['block_index'] == b: # Periodicity from block to block 
                if p['block1']['block_index'] == b:
                    p['block1']['IMIN'], p['block1']['IMAX'] = check_and_swap(p['block1']['IMIN'], p['block1']['IMAX'])
                    p['block1']['JMIN'], p['block1']['JMAX'] = check_and_swap(p['block1']['JMIN'], p['block1']['JMAX'])
                    p['block1']['KMIN'], p['block1']['KMAX'] = check_and_swap(p['block1']['KMIN'], p['block1']['KMAX'])
                    voi = [p['block1']['IMIN'], p['block1']['IMAX'], p['block1']['JMIN'], p['block1']['JMAX'],p['block1']['KMIN'], p['block1']['KMAX']]
                else:
                    p['block2']['IMIN'], p['block2']['IMAX'] = check_and_swap(p['block2']['IMIN'], p['block2']['IMAX'])
                    p['block2']['JMIN'], p['block2']['JMAX'] = check_and_swap(p['block2']['JMIN'], p['block2']['JMAX'])
                    p['block2']['KMIN'], p['block2']['KMAX'] = check_and_swap(p['block2']['KMIN'], p['block2']['KMAX'])
                    voi = [p['block2']['IMIN'], p['block2']['IMAX'], p['block2']['JMIN'], p['block2']['JMAX'],p['block2']['KMIN'], p['block2']['KMAX']]
                CreateSubset(block_source, voi, name='face_match '+str(idx),rgb_face_matches=rgb_face_matches[idx])