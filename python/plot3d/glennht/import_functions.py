from typing import Any, Dict, List
import numpy as np 

def read_ght_conn(filename:str):
    """Reads in the glennht connectivity file and splits it into partitions using metis 

    Args:
        filename (str): _description_

    Returns:
        
        Tuple containing: 
            **max_block_index** (int): Maximum block index
            **face_matches** (List[Dict[str,int]]): Face matches containing information which block is connected to which
            **num_connections_foreach_block** (np.ndarray): Number of connections for each block 
            **nZones** (int): Number of zones
            **zone_types** (List[int]): Zone types
            **Zones** (List[int]): Zones
            **GIFs** (List[Dict[str,int]]): General Interfaces
            
    """
    IMIN = 0; JMIN = 0; KMIN = 0; IMAX = 0; JMAX = 0; KMAX = 0  # Preallocate  # noqa: E702
    # Read in the glennht connectivity file
    with open(filename, "r") as f:        
        line = f.readline()
        pairs = int(line.lstrip().rstrip().split(" ")[0])
        idx = 0
        blk_id1 = 0
        block_to_block = list()
        # Read matches 
        for _ in range(pairs*2): 
            line = f.readline()
            idx += 1
            temp = line.lstrip().rstrip().split(" ")
            temp = [int(t)-1 for t in temp if t.strip() != '']
            if (len(temp) == 7):
                if idx % 2 == 0: # When there is a pair of blockd, write to the list
                    block_to_block.append({
                        'block1':{'block_index':blk_id1,
                                  "IMIN":IMIN,"JMIN":JMIN,"KMIN":KMIN,
                                  "IMAX":IMAX,"JMAX":JMAX,"KMAX":KMAX},
                        'block2':{'block_index':temp[0],
                                  "IMIN":temp[1],"JMIN":temp[2],"KMIN":temp[3],
                                  "IMAX":temp[4],"JMAX":temp[5],"KMAX":temp[6]},
                        
                    })                   
                else:
                    blk_id1 = temp[0]
                    IMIN = temp[1]; JMIN = temp[2]; KMIN = temp[3]  # noqa: E702
                    IMAX = temp[4]; JMAX = temp[5]; KMAX = temp[6]  # noqa: E702
        line = f.readline()
        # Read OuterFaces
        outer_faces = []
        nFaces = int(line.lstrip().rstrip().split(" ")[0])
        for id in range(int(nFaces)): 
            line = f.readline() # Skip outer faces
            temp  = line.lstrip().rstrip().split(" ") 
            temp = [int(x)-1 for x in temp if x]
            outer_faces.append({
                "block_index":temp[0],
                "IMIN":temp[1],"JMIN":temp[2],"KMIN":temp[3],
                "IMAX":temp[4],"JMAX":temp[5],"KMAX":temp[6],
                "id":temp[7]
                })
        
        # Read GIF
        line = f.readline() 
        nGIF = int(line.lstrip().rstrip().split(" ")[0])
        GIFs = list()
        for _ in range(int(nGIF)):
            line = f.readline()
            temp = line.lstrip().rstrip().split(" ")
            temp = [int(t) for t in temp if t.strip() != '']
            GIFs.append({'S1':temp[0],'S2':temp[1],'GIF_TYPE':temp[2],'GIF_ORDER':temp[3]})
        
        # Read Zones
        line = f.readline() 
        nZones = int(line.lstrip().rstrip().split(" ")[0]) # 4 - number of zones
        temp = line.lstrip().rstrip().split(" ")
        zone_types = [int(t) for t in temp if t.strip() != ''] # 1 1 1 2 -> fluid, fluid, fluid, solid; 4 zones
        Zones = list()
        while True:
            line = f.readline()
            if not line:
                break
            else:
                temp = line.lstrip().rstrip().split(" ")
                temp = [int(t) for t in temp if t.strip() != '']
                Zones.extend(temp)
    
        # Find out how many connecting nodes each block has 
        l1 = [(row['block1']['block_index'],
              (row['block1']['IMAX']-row['block1']['IMIN'])*
              (row['block1']['JMAX']-row['block1']['IMIN'])*
              (row['block1']['KMAX']-row['block1']['KMIN']))
              for row in block_to_block]
        l2= [(row['block2']['block_index'],
              (row['block2']['IMAX']-row['block2']['IMIN'])*  
              (row['block2']['JMAX']-row['block2']['IMIN'])*
              (row['block2']['KMAX']-row['block2']['KMIN']))
              for row in block_to_block]
        
        l1.extend(l2)
        l1 = np.array(l1)
        total_blocks = l1[:,0].max()
        _,indices = np.unique(l1[:,0],return_index=True)
        num_connections_foreach_block = l1[indices,:]
        num_connections_foreach_block.sort(axis=0)
        face_matches = block_to_block
        return total_blocks, face_matches,outer_faces,num_connections_foreach_block,nZones,zone_types,Zones,GIFs
