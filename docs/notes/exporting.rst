Exporting to GlennHT Connectivity Files
###############################################
This is an example of how I read a mesh, convert it to ascii, save it, find connectivity, find periodicty and export to glennht format. 
In this example we will use the file `PahtCascade-ASCII <https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/PahtCascade-ASCII.xyz>`_

.. code-block:: python 
    :linenos: 

    from itertools import combinations
    import os, sys
    sys.path.insert(0,'../../')
    from plot3d import write_plot3D, read_plot3D, find_periodicity
    from plot3d import find_matching_blocks, get_outer_faces, connectivity
    from glennht_con import export_to_glennht_conn
    import pickle

    # Convert to binary because of size 
    blocks = read_plot3D('PahtCascade-ASCII.xyz', binary = False)
    write_plot3D('PahtCascade.xyz',blocks,binary=True)

    blocks = read_plot3D('PahtCascade.xyz', binary = True, big_endian=True)

    if not os.path.exists('connectivity.pickle'):
        blocks = read_plot3D('PahtCascade.xyz', binary = True, big_endian=True)
        # Block 1 is the blade O-Mesh k=0
        # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
        face_matches, outer_faces_formatted = connectivity(blocks)
        with open('connectivity.pickle','wb') as f:
            pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

    with open('connectivity.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

    blocks = read_plot3D('PahtCascade.xyz', binary = True, big_endian=True)
    periodic_surfaces, outer_faces_to_keep = find_periodicity(blocks,outer_faces,periodic_direction='k')
    # Append periodic surfaces to face_matches
    face_matches.extend(periodic_surfaces)

    export_to_glennht_conn(face_matches,outer_faces_to_keep,'finalmesh')


.. code-block:: python 
    :linenos:

    from typing import List
    def export_to_glennht_conn(nblocks:int,matches:List,block_surfaces:dict,filename:str,other_settings:dict=None):
        """Exports the connectivity to GlennHT format 

        Args:
            matches (dict): Any matching faces between blocks 
            block_surfaces (dict): Non matching faces of all blocks or surfaces to consider 
            filename (str): filename to write to 
            other_settings (dict, optional): contains general interfaces, zones, blocks that are in the zone. Defaults to None.
        """

        blocks = ['block1','block2'] # Block 1 and Block 2 are arbitrary names. Their index matters
        with open(filename + '.ght_conn','w') as fp:
            # Print matches
            nMatches = len(matches)
            fp.write(f'{nMatches}\n') # Print number of matches 
            for match in matches:                        
                for block in blocks:
                    block_indx = match[block]['index']+1 # block1 and block2 are arbitrary names, the key is the block index 
                    block_IMIN = match[block]['IMIN']+1
                    block_JMIN = match[block]['JMIN']+1
                    block_KMIN = match[block]['KMIN']+1

                    block_IMAX = match[block]['IMAX']+1
                    block_JMAX = match[block]['JMAX']+1
                    block_KMAX = match[block]['KMAX']+1

                    fp.write(f"{block_indx:3d}\t{block_IMIN:5d} {block_JMIN:5d} {block_KMIN:5d}\t{block_IMAX:5d} {block_JMAX:5d} {block_KMAX:5d}\n")
            # Print Surfaces 
            # Get total number of surfaces 
            id = 1
            lines = list()
            for block in block_surfaces:
                block_indx = block['index']+1
                for surface in block['surfaces']:                
                    IMIN = surface['IMIN']+1
                    JMIN = surface['JMIN']+1
                    KMIN = surface['KMIN']+1
                    
                    IMAX = surface['IMAX']+1
                    JMAX = surface['JMAX']+1
                    KMAX = surface['KMAX']+1
                    lines.append(f"{block_indx:3d}\t{IMIN:5d} {JMIN:5d} {KMIN:5d}\t{IMAX:5d} {JMAX:5d} {KMAX:5d}\t{id:4d}\n")
                    id+=1
                    
            fp.write(f'{len(lines)}\n')
            [fp.write(line) for line in lines]

            # Write general interfaces
            n_gif = len(other_settings['general_interfaces'])
            fp.write(f'{n_gif:d}\n')
            for gif in other_settings['general_interfaces']:
                surf1 = gif['surface_pairs'][0]
                surf2 = gif['surface_pairs'][1]
                if gif['gif_type'].lower() == 'mixing_plane':
                    gif_kind = 2 
                else:
                    gif_kind = 1 # for conjugate
                if gif['is_polar']:
                    gif_kind*=-1
                fp.write(f'{surf1} {surf2} {gif_kind} 1 ') # Assume no new lines and just write on a single line. 
            if len(other_settings['general_interfaces'])>1:
                fp.write('\n')

            def zonetype_to_glennht(type:str):
                if type == 'fluid':
                    return '1'
                else:                   # Second zone is solid
                    return '2'

            def lookup_zone_index(block_indx:int) -> int:
                """Searches through all the zones for the block_index 

                Args:
                    block_indx (int): index of block

                Returns:
                    int: the zone index 
                """
                zone_indx = 1 
                for zone in other_settings['zones']:
                    if block_indx in zone['blocks']:
                        return zone_indx
                    zone_index += 1  
            
            # Write the Zones 
            n_zones = len(other_settings['zones'])
            fp.write(f'{n_zones:d}\n')
            for zone in other_settings['zones']:    # Write which zone is fluid and which is solid
                fp.write(zonetype_to_glennht(zone['type']) + ' ')
            if len(other_settings['zones'])>0:
                fp.write('\n')
            
            # Write out which block belongs to what zone index
            for block_indx in range(nblocks):   
                zone_index = lookup_zone_index(block_indx) 
                if (block_indx+1) % 6 == 0:
                    fp.write(f'{zone_index:d}')
                    fp.write('\n')
                else:
                    fp.write(f'{zone_index:d}' + ' ')


