'''
    Combines the connectivty and constructs the glennht boundary conditions file 
'''
from typing import Dict, List
import pickle
from glennht_library import CheckDictionary
from plot3d import Face, find_connected_face,find_face,read_plot3D

def read_and_shift_block_index(filename:str, block_starting:int=0,blade_name:str='stator'):
    def advance_block_index1(data:List[Dict[str,int]]):
        for i in range(len(data)):
            data[i]['block_index'] += block_starting
        return data

    def advance_block_index2(matches:List[Dict[str,int]]):
        for i in range(len(matches)):
            matches[i]['block1']['block_index']+=block_starting
        return matches

    with open(filename,'rb') as f:
        data = pickle.load(f)
        matches = data['face_matches']
        advance_block_index2(matches)
        
        shroud = CheckDictionary(data,f'{blade_name}_shroud')
        advance_block_index1(shroud)

        hub = CheckDictionary(data,f'{blade_name}_hub')
        advance_block_index1(hub)

        body = CheckDictionary(data,f'{blade_name}_body')
        advance_block_index1(body)

        outer_faces = CheckDictionary(data,'outer_faces')
        advance_block_index1(outer_faces)

        periodic = CheckDictionary(data,'periodic_faces')
        advance_block_index2(periodic)

        mixing_plane = CheckDictionary(data,'mixing_plane')
        advance_block_index1(mixing_plane)
        
        inlet = CheckDictionary(data,'inlet')
        advance_block_index1(inlet)
        
        outlet = CheckDictionary(data,'outlet')
        advance_block_index1(outlet)
        
    return matches,periodic,shroud,hub,body,outer_faces,mixing_plane,inlet,outlet

def print_connectivity(filename:str,matches:List[Dict[str,int]],faces_and_types:Dict[int,List[Dict[str,int]]],gifs:Dict[str,int],zones:Dict[str,List[int]]):
    """Writes the Connectivity File

    Args:
        filename (str): _description_
        matches (List[Dict[str,int]]): _description_
        faces_and_types (Dict[int,List[Dict[str,int]]]): _description_

    Returns:
        _type_: _description_
    """
    filename=filename.split('.')[0]
    def print_matches(matches):
        lines = list() 
        match_keys = ['block1','block2'] # block1 and block2 are arbitrary names, the key is the block index 
        nMatches = len(matches)
        lines.append(f'{nMatches}\n') # Print number of matches 
        for match in matches:                        
            for block in match_keys:
                block_indx = match[block]['block_index']+1 
                block_IMIN = match[block]['IMIN']+1
                block_JMIN = match[block]['JMIN']+1
                block_KMIN = match[block]['KMIN']+1

                block_IMAX = match[block]['IMAX']+1
                block_JMAX = match[block]['JMAX']+1
                block_KMAX = match[block]['KMAX']+1

                lines.append(f"{block_indx:3d}\t{block_IMIN:5d} {block_JMIN:5d} {block_KMIN:5d}\t{block_IMAX:5d} {block_JMAX:5d} {block_KMAX:5d}\n")
        return lines

    def print_face_group(id:int,faces_and_types:List[Dict[str,int]]):
        lines = list()
        for f in faces_and_types['faces']:
            block_index = f["block_index"]
            IMIN = f['IMIN']+1
            JMIN = f['JMIN']+1
            KMIN = f['KMIN']+1
            
            IMAX = f['IMAX']+1
            JMAX = f['JMAX']+1
            KMAX = f['KMAX']+1
            lines.append(f"{block_index:3d}\t{IMIN:5d} {JMIN:5d} {KMIN:5d}\t{IMAX:5d} {JMAX:5d} {KMAX:5d}\t{id:4d} \n")
        return lines 

    matches = print_matches(matches)
    surfaces = list()
    for k,v in faces_and_types.items():
        surfaces.extend(print_face_group(k,v))
        
    
    with open(f'{filename}.ght_conn','w') as fp:
        # Print matches
        fp.writelines(matches)
        # Print number of surfaces
        fp.write(f'{len(surfaces)}\n')
        fp.writelines(surfaces)

        # Lets write the number of gifs (mixing planes, solid|fluid)
        fp.write(f'{len(gifs)}\n')
        for gif in gifs:
            surface_pairs = gif['surface_pairs']
            gif_type = gif['gif_type']
            fp.write(f'{surface_pairs[0]} {surface_pairs[1]} {gif_type} 1\n') # 2 = mixing plane, -2 = use polar coordinates, 1 is for linear interpolation
        
        # Write the zones
        num_zones = len(zones['fluid'])>0 + len(zones['solid'])>0
        fp.write(f'{num_zones}\n')
        for f in zones['fluid']: # loop through list
            fp.write(f'{f}')
        fp.write('\n')
        
def modify_bc_template():
    pass
'''
    We shouldn't have any outerfaces in this step. All outerfaces should either be a hub, shroud, or body
'''
stator_blocks = read_plot3D('stator_split.xyz',binary=True)
rotor_blocks = read_plot3D('rotor_split.xyz',binary=True)
# Read in all the connectivity 
stator_matches,stator_periodic,stator_shroud,stator_hub,stator_body,_,stator_mixing_plane,stator_inlet,_ = read_and_shift_block_index('stator_split_connectivity_final.pickle')
rotor_matches,rotor_periodic,rotor_shroud,rotor_hub,rotor_body,_,rotor_mixing_plane,_,rotor_outlet = read_and_shift_block_index('rotor_split_connectivity_final.pickle',block_starting=len(stator_blocks),blade_name='rotor')

# Lets combine these connectivities
stator_matches.extend(rotor_matches)
stator_matches.extend(stator_periodic)
stator_matches.extend(rotor_periodic)

faces_and_types = {
                    0:{
                        "faces":stator_inlet,
                        "type":"inlet"
                    },
                    1:{
                        "faces":rotor_outlet,
                        "type":"outlet"
                    },
                    2:{
                        "faces":stator_shroud,
                        "type":"wall"
                    },
                    3:{
                        "faces":stator_hub,
                        "type":"wall"
                    },
                    4:{
                        "faces":stator_body,
                        "type":"wall"
                    },
                    5:{
                        "faces":rotor_shroud,
                        "type":"wall"
                    },
                    6:{
                        "faces":rotor_hub,
                        "type":"wall"
                    },
                    7:{
                        "faces":rotor_body,
                        "type":"wall"
                    },
                    8:{
                        "faces":stator_mixing_plane,
                        "type":"gif"
                    },
                    9:{
                        "faces":rotor_mixing_plane,
                        "type":"gif"
                    }
                }
gifs = [{
    'surface_pairs': [8,9],
    'gif_type': -2
}]
zones = {'fluid':[], 'solid':[]} # which block indicies are fluid or solid 

# all blocks are fluid
stator_blocks.extend(rotor_blocks)
for i in range(len(stator_blocks)):
    zones['fluid'].append(1)

print_connectivity('stage_connectivity.ght_conn',stator_matches,faces_and_types, gifs, zones)

print('check')