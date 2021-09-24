def export_to_glennht_conn(matches:dict,block_surfaces:dict,filename:str):
    """Exports the connectivity to GlennHT format 

    Args:
        matches (dict): Any matching faces between blocks 
        block_surfaces (dict): Non matching faces of all blocks or surfaces to consider 
    """

    blocks = ['block1','block2']
    with open(filename + '.ght_conn','w') as fp:
        # Print matches
        nMatches = len(matches)
        fp.write(f'{nMatches}\n') # Print number of matches 
        for match in matches:                        
            for block in blocks:
                block_indx = match[block]['block_index']+1 # block1 and block2 are arbitrary names, the key is the block index 
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
        for surface in block_surfaces:
            block_indx = surface['block_index']+1
            IMIN = surface['IMIN']+1
            JMIN = surface['JMIN']+1
            KMIN = surface['KMIN']+1
            
            IMAX = surface['IMAX']+1
            JMAX = surface['JMAX']+1
            KMAX = surface['KMAX']+1
            lines.append(f"{block_indx:3d}\t{IMIN:5d} {JMIN:5d} {KMIN:5d}\t{IMAX:5d} {JMAX:5d} {KMAX:5d}\t{id:4d}\n")
            id += 1
             
        fp.write(f'{len(lines)}\n')
        [fp.write(line) for line in lines]