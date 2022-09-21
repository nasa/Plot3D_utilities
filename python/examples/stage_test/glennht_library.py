'''
    Reads in the boundary conditions as a JSON 
'''
import sys
from typing import Dict, List
import numpy as np 
import matplotlib.pyplot as plt
import math 
from thermo import Mixture
import copy
import argparse
import json 

# Boundary conditions in Dimensional Quantities 
# * Note: hbl, tbl = boundary layer thickness % span 
# * tfree = freestream turbulence level
# * lfree = freestream turbulence length
n = 100 # Number of profile points

def get_args_parser():
    parser = argparse.ArgumentParser(description='Self-Attention GAN trainer')
    parser.add_argument('--input_file', metavar='PATH', type=str, help='JSON file containing boundary conditions')
    parser.add_argument('--geom_file', metavar='PATH', type=str, help='JSON file the geometry information')
    parser.add_argument('--output_file', metavar='PATH', type=str, help='Path to save non-dimensional boundary condition data')
    parser.add_argument('--english_units', default=True, type=bool, help='Are you using english units of (psi) and (Rankine) for pressure and temperature')
    return parser

def sutherland(T:float):
    """Sutherland's law calculation for mu and k

    https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html

    Args:
        T (float): Temperature to calculate viscosity and k

    (tuple): tuple containing:

        - **mu** (float): dynamic viscosity (Pa*s)
        - **k** (float): air thermal conductivity W/(m*K)
    """
    Tref = 273
    mu_ref = 1.716E-5
    S_mu = 111
    mu = mu_ref * math.pow(T/Tref,3/2)*(Tref+S_mu)/(T+S_mu)

    k_ref = 0.0241
    S_k = 194
    k = k_ref * math.pow(T/Tref,3/2)*(Tref+S_k)/(T+S_k)

    return mu,k

def two_point_boundary(h:float,beta:float=1.2,xi:float=None,num_points:int=100):
    """Spaces points along a grid by concentrating more points near the ends


    Args:
        h (float): max value 0<=x<=h
        beta (float, optional): Expansion ratio. Defaults to 1.2.
        xi (float, optional): Spacing in the time domain (constant spacing). Defaults to None.
        num_points (int, optional): number of points. Defaults to 100.

    Returns:
        [type]: [description]
    """
    if (not xi):
        xi = np.linspace(0,1,num_points)
        
    num  = h*(beta+1)* np.power((beta+1)/(beta-1), 2*xi-1)-beta+1
    den = 2*(1+np.power((beta+1)/(beta-1), 2*xi-1))
    x = num/den 
    return x

def create_profile(mean_val:float,bl_height:float,pspan:np.ndarray):
    """Creates a profile to describe inlet P0,T0, angle etc

    Args:
        mean_val ([type]): [description]
        bl_height ([type]): [description]

    Returns:
        [type]: [description]
    """
    quantity=pspan*0
    for i in range(len(pspan)):
        if pspan[i]<bl_height:
            quantity[i] = mean_val*np.power(pspan[i]/bl_height,1/7)
        elif pspan[i]>(max(pspan)-bl_height):
            quantity[i] = mean_val*np.power((pspan[-1]-pspan[i])/bl_height,1/7)
        else:
            quantity[i] = mean_val
    return quantity

def write_inlet_profile(r:np.ndarray,P0:np.ndarray,T0:np.ndarray,alpha:np.ndarray,beta:np.ndarray):
    """Creates inlet_profile.dat file with headers 

    Args:
        P0 (np.ndarray): Total Pressure Array
        T0 (np.ndarray): Total Temperature Array
        alpha (np.ndarray): inlet flow angle in x,y direction
        beta (np.ndarray): inlet flow angle in radial direction
    """
    with open('inlet_profile.dat','w') as f:
        for i in range(len(P0)):
            f.write(f'{r[i]}\t{P0[i]}\t{T0[i]}\t{alpha[i]}\t{beta[i]}\n')

def modify_job_template_file(filename:str,bcs:dict,job_settings:dict, output_filename:str):
    with open(filename, "r") as f:
        text = f.read()

        # Job .ght related 
        text = text.replace('{initP0}','{0:.4f}'.format(bcs['P0_Inlet']))
        text = text.replace('{initT0}','{0:.4f}'.format(bcs['T0_Inlet']))
        text = text.replace('{initMach}','{0:.4f}'.format(bcs['Mach_Inlet']))
        text = text = text.replace('{initAlpha}','{0:.4f}'.format(bcs['alpha']))
        text = text.replace('{initTu}','{0:.4f}'.format(bcs['Tu_intensity_init']))
        text = text.replace('{initTs}','{0:.4f}'.format(bcs['Tu_length_scale_init']))        

        text = text.replace('{refLen}','{0:.4f}'.format(bcs['refLen']))
        text = text.replace('{refP0}','{0:.4f}'.format(bcs['Pref']))
        text = text.replace('{refT0}','{0:.4f}'.format(bcs['Tref']))
        text = text.replace('{refVisc}', '{0:5e}'.format(bcs['mu']))
        text = text.replace('{refGamma}','{0:.4f}'.format(bcs['gamma']))
        text = text.replace('{refCp}','{0:.4f}'.format(bcs['Cp']))
        text = text.replace('{refVel}','{0:.4f}'.format(bcs['V_rel']))
        text = text.replace('{refRho}','{0:.4f}'.format(bcs['rho']))
        text = text.replace('{refCond}','{0:.4f}'.format(bcs['k']))
        text = text.replace('{Tref_fluid}','{0:.4f}'.format(bcs['Tref']))
        text = text.replace('{Prref_fluid}','{0:.4f}'.format(bcs['refPr']))
        text = text.replace('{refPr}','{0:.4f}'.format(bcs['refPr']))
        text = text.replace('{refPrandtl}','{0:.4f}'.format(bcs['refPr']))
        text = text.replace('{refReynolds}','{0:.4f}'.format(bcs['Re']))

        text = text.replace('{RestartSoln}','{0}'.format(job_settings['RestartSoln']))
        text = text.replace('{mRunLevel}','{0:d}'.format(job_settings['mRunLevel']))
        text = text.replace('{nTimeSteps}','{0:4d}'.format(job_settings['nTimeSteps']))
        text = text.replace('{maxPseudoSteps}','{0:4d}'.format(job_settings['maxPseudoSteps']))
        text = text.replace('{SolnInFile}','"{0}"'.format(job_settings['SolnInFile']))

        text = text.replace('{CFLn}','{0}'.format(job_settings['CFLn']))
        text = text.replace('{CFLr}','{0}'.format(job_settings['CFLr']))
        text = text.replace('{ddcmp}','{0}'.format(job_settings['ddcmp']))

        text = text.replace('{mesh_file}','{0}'.format(job_settings['mesh_file']))
        text = text.replace('{connectivity}','{0}'.format(job_settings['connectivity']))
        
        text = text.replace('{grid_coarsest}','{0}'.format(job_settings['grid_coarsest']))
        text = text.replace('{pre_mg_sweeps}','{0}'.format(job_settings['pre_mg_sweeps']))
        text = text.replace('{mg_sweeps}','{0}'.format(job_settings['mg_sweeps']))
    with open(output_filename,'w') as ff:
        ff.write(text)
    
def modify_bc_template_file(filename:str,bcs:dict):
    """Modifies the template boundary condition file 

    Args:
        filename (str): name of the template file
        bcs (dict): boundary conditions 
    """

    with open(filename, "r") as f:
        text = f.read()
        # Boundary condition related
        text = text.replace('{Mach_Inlet}','{0:.4f}'.format(bcs['Mach_Inlet']))
        text = text.replace('{bl_height}','{0:.4f}'.format(bcs['bl_height'])) # ! Check to see if this is 0.1336
        text = text.replace('{NBlades_inlet}','{0:d}'.format(bcs['NBlades'][0]))
        text = text.replace('{T0_Inlet}','{0:.4f}'.format(bcs['T0_Inlet']))
        text = text.replace('{P0_inlet}','{0:.4f}'.format(bcs['P0_Inlet']))
        text = text.replace('{Ps_outlet}','{0:.4f}'.format(bcs['Ps_outlet']))
        text = text.replace('{NBlades_outlet}','{0:d}'.format(bcs['NBlades'][-1]))
        text = text.replace('{Tu_const}','{0:.4f}'.format(bcs['Tu_intensity_init']))
        text = text.replace('{Tu_length_scale_init}','{0:.4f}'.format(bcs['Tu_length_scale_init']))
        text = text.replace('{initAlpha}','{0:.4f}'.format(bcs['alpha']))
        text = text.replace('{Tref_fluid}','{0:.4f}'.format(bcs['Tref']))
        text = text.replace('{Prref_fluid}','{0:.4f}'.format(bcs['refPr']))

        with open(filename.replace('_template',''),'w') as ff:
            ff.write(text)

def print_inlet_profile(P0:np.ndarray,T0:np.ndarray,alpin:np.ndarray,pspan:np.ndarray):
    """Prints the inlet profile by plotting [P0, T0, alpin] vs pspan

    Args:
        P0 (np.ndarray): Total Pressure as an array
        T0 (np.ndarray): Total Temperature an as array 
        alpin (np.ndarray): Inlet angle 
        pspan (np.ndarray): percentage span 
    """
    # * Plot the boundary conditions
    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(8,6), dpi=150)
    ax[0,0].plot(P0,pspan*100,label='Pt-in')
    ax[0,0].set_xlabel('Inlet Total Pressure (psia)',fontsize=14)
    ax[0,0].set_ylabel('percent span',fontsize=14)
    ax[0,0].set_ylim(0,100)

    ax[0,1].plot(P0,pspan*100,label='Pt-in')
    ax[0,1].set_xlabel('Inlet Total Pressure (psia)',fontsize=14)
    ax[0,1].set_ylabel('percent span',fontsize=14)
    ax[0,1].set_ylim(0,100)

    ax[1,0].plot(T0,pspan*100,label='Tt-in')
    ax[1,0].set_xlabel('Inlet Total Temperature',fontsize=14)
    ax[1,0].set_ylabel('percent span',fontsize=14)
    ax[1,0].set_ylim(0,100)

    ax[1,1].plot(alpin,pspan*100,label='alpha-in')
    ax[1,1].set_xlabel('Inlet flow angle - alpha',fontsize=14)
    ax[1,1].set_ylabel('percent span',fontsize=14)
    ax[1,1].set_ylim(0,100)
    fig.tight_layout(pad=1.0)
    plt.savefig('boundary conditions.png',dpi=150)


def modify_boundary_conditions_file(settings_file:str,geom_file:str):
    print(args)
    
    with open(settings_file,'r') as json_file:
        settings = json.load(json_file)
        bcs = settings['boundary_conditions']

    with open(geom_file,'r') as geom_file:
        geom = json.load(geom_file)
        xhub = geom['channel']['hub']['x']
        rhub = geom['channel']['hub']['r']
        xshroud = geom['channel']['shroud']['x']
        rshroud = geom['channel']['shroud']['r']

    # Conversion from english to kelvin 
    if (args.english_units):
        bcs['T0_Inlet'] *= 0.555556     # Convert to Kelvin
        bcs['P0_Inlet'] *= 6894.76      # Convert to Pascal 
        bcs['Ps_outlet'] *= 6894.76

    bcs['alpha'] = bcs['alpha']
    air = Mixture('air', T=bcs['T0_Inlet'], P=bcs['P0_Inlet'])
    gamma = air.Cp/air.Cvg
    pspan = two_point_boundary(1,beta=1.1,num_points=100)
    rspan = rhub[0] + (rshroud[0] - rhub[0]) * pspan 

    # Inlet profile Radius, P0, T0, alpha, beta
    # Make it flat at 0.6
    # Pt at hub and case endwall equals static pressure at inlet 
    Ps_inlet = math.pow( 1+(gamma-1)/2*bcs['Mach_Inlet']*bcs['Mach_Inlet'], -gamma/(gamma-1)) * bcs['P0_Inlet']
    P0 = Ps_inlet + create_profile(mean_val=bcs['P0_Inlet']-Ps_inlet,bl_height=bcs['bl_height'], pspan=pspan)     # psia
    T0 = bcs['T0_Inlet']*np.ones(len(P0),)     # Rankine
    alpin = bcs['alpha']*np.ones(len(P0),)
    beta = 0*np.zeros(len(P0),)
    
    # Additional calculations 
    mu,k = sutherland(bcs['T0_Inlet'])
    air = Mixture('air', T=bcs['T0_Inlet'], P=bcs['P0_Inlet'])
    gamma = air.Cp/air.Cvg
    
    Ts = 1/(1+(gamma-1)/2 * bcs['Mach_Inlet']*bcs['Mach_Inlet']) * bcs['T0_Inlet']
    V_rel = bcs['Mach_Inlet'] * math.sqrt(gamma*287*Ts)
    Pr = air.Cp * mu / k

    # Non dimensionalize
    bcs_non_dim = copy.deepcopy(bcs)
    bcs_non_dim['T0_Inlet'] = 1
    bcs_non_dim['P0_Inlet'] = 1
    bcs_non_dim['Ts_Inlet'] = Ts / bcs['T0_Inlet']
    bcs_non_dim['Ps_outlet'] = bcs['Ps_outlet']/bcs['P0_Inlet']
    bcs_non_dim['Pref'] = bcs['P0_Inlet']
    bcs_non_dim['Tref'] = bcs['T0_Inlet']
    bcs_non_dim['refPr'] = Pr
    bcs_non_dim['k'] = k
    bcs_non_dim['mu'] = mu
    bcs_non_dim['rho'] = air.rho

    # P = rho R T; P/P0 = rho * R * T/T0; P/P0 * T0/T * 1/R = rho
    bcs_non_dim['gamma'] = gamma
    bcs_non_dim['Cp'] = air.Cp
    bcs_non_dim['V_rel'] = V_rel
    bcs_non_dim['Re'] = air.rho * V_rel * bcs_non_dim['refLen'] / mu 

    
    for job in settings['job_settings']:
        name = job["name"]
        modify_job_template_file('job_template.ght',bcs_non_dim,job,output_filename=f"job_{name}.ght") # Write out the modified job file
        
    modify_bc_template_file('boundary_conditions_template.bcs', bcs_non_dim) # Write out the modified job file
    write_inlet_profile(rspan,P0/bcs['P0_Inlet'],T0/bcs['T0_Inlet'],alpin,beta)


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


def print_connectivity(filename:str,matches:List[Dict[str,int]],faces_and_types:Dict[int,List[Dict[str,int]]],gifs:Dict[str,int],zones:Dict[str,List[int]]):
    """Writes the Connectivity File

    Args:
        filename (str): Name of connectivity file
        matches (List[Dict[str,int]]): _description_
        faces_and_types (Dict[int,List[Dict[str,int]]]): _description_

    Returns:
        None: outputs a connectivity file
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
        
