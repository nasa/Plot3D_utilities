'''
    Combines the connectivty and constructs the glennht boundary conditions file 
'''
from typing import Dict, List
import pickle
from glennht_library import CheckDictionary, print_connectivity
from plot3d import Face, Block,read_plot3D

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

# Scale mesh to metric units
stator_blocks = 
print('check')