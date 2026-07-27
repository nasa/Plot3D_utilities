'''
    Combines the connectivty and constructs the glennht boundary conditions file 
'''
import sys
sys.path.insert(0,'C:\GitHub\Plot3D_utilities\python')
from typing import Dict, List
import pickle, json, math, copy
from glennht_library import CheckDictionary, modify_bc_template_file, print_connectivity, sutherland
from plot3d import Block,read_plot3D, write_plot3D
from thermo import Mixture


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
    We shouldn't have any outerfaces in this step. All outerfaces should either be a hub, shroud, or body.
'''
stator_blocks = read_plot3D('stator_split.xyz',binary=True)
rotor_blocks = read_plot3D('rotor_split.xyz',binary=True)

# Read in all the connectivity 
stator_matches,stator_periodic,stator_shroud,stator_hub,stator_body,_,stator_mixing_plane,stator_inlet,_ = read_and_shift_block_index('stator_split_connectivity_final.pickle')
rotor_matches,rotor_periodic,rotor_shroud,rotor_hub,rotor_body,_,rotor_mixing_plane,_,rotor_outlet = read_and_shift_block_index('rotor_split_connectivity_final.pickle', block_starting=len(stator_blocks),blade_name='rotor')

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

# All blocks are fluid
stator_blocks.extend(rotor_blocks)
for i in range(len(stator_blocks)):
    zones['fluid'].append(1)

print_connectivity('stage_connectivity.ght_conn',stator_matches,faces_and_types, gifs, zones)

# Scale mesh in inches to meters
[b.scale(0.0254) for b in stator_blocks] # convert to meters
write_plot3D('StageMesh_metric.xyz',stator_blocks,binary=True)
# Normalize the boundary conditions
with open('settings.json','r') as f:
    data = json.load(f)
    bcs = data['boundary_conditions']
    # Perform calculations
    air = Mixture('air', T=bcs['T0_Inlet'], P=bcs['P0_Inlet'])
    gamma = air.Cp/air.Cvg
    mu,k = sutherland(bcs['T0_Inlet'])
    Pr = air.Cp * mu / k
    Ts = 1/(1+(gamma-1)/2 * bcs['Mach_Inlet']*bcs['Mach_Inlet']) * bcs['T0_Inlet']
    V_rel = bcs['Mach_Inlet'] * math.sqrt(gamma*287*Ts)

    # Non dimensionalize
    bcs_non_dim = copy.deepcopy(bcs)
    bcs_non_dim['Ps_outlet'] /= bcs['P0_Inlet']
    bcs_non_dim['Pref'] = bcs['P0_Inlet']
    bcs_non_dim['Tref'] = bcs['T0_Inlet']
    bcs_non_dim['refPr'] = Pr
    bcs_non_dim['k'] = k
    bcs_non_dim['mu'] = mu
    bcs_non_dim['rho'] = air.rho

    bcs_non_dim['gamma'] = gamma
    bcs_non_dim['Cp'] = air.Cp
    bcs_non_dim['V_rel'] = V_rel
    bcs_non_dim['Re'] = air.rho * V_rel * bcs_non_dim['refLen'] / mu
    
    modify_bc_template_file('boundary_conditions_template.bcs', bcs_non_dim)