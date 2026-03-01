"""Import functions for reading GlennHT connectivity files.

This module provides utilities for parsing ``.ght_conn`` files produced by
GlennHT or by the export helpers in
:mod:`~plot3d.glennht.export_functions`.  Block indices read from the file
are converted from 1-based (GlennHT convention) to 0-based (Python
convention) before being returned.
"""
from typing import Any, Dict, List
import numpy as np


def read_ght_conn(filename: str):
    """Read a GlennHT connectivity file and return its parsed sections.

    All block and face indices are returned 0-based (the file stores them
    1-based; this function subtracts 1 on read).

    Args:
        filename: Path to the ``.ght_conn`` connectivity file.

    Returns:
        A tuple of eight elements:

        - **max_block_index** (*int*): Maximum 0-based block index found in
          the face-match section.
        - **face_matches** (*List[Dict[str, dict]]*): Block-to-block face
          match records.  Each entry is a dict with keys ``"block1"`` and
          ``"block2"``, where each value contains ``"block_index"``,
          ``"lb"`` (tuple of 3 ints), and ``"ub"`` (tuple of 3 ints)
          (all 0-based).
        - **outer_faces** (*List[Dict[str, int]]*): Boundary surface records
          with keys ``"block_index"``, ``"lb"`` (tuple of 3 ints),
          ``"ub"`` (tuple of 3 ints) (0-based), and ``"id"``
          (0-based surface ID).
        - **num_connections_foreach_block** (*np.ndarray*): Array of shape
          ``(N, 2)`` where column 0 is the 0-based block index and column 1
          is the number of shared interface nodes for that block.
        - **nZones** (*int*): Total number of volume zones declared in the
          file.
        - **zone_types** (*List[int]*): Zone type codes (e.g. 1 = fluid,
          2 = solid) read from the zone header line.
        - **Zones** (*List[int]*): Flat list of contiguous zone indices
          assigned to each block.
        - **GIFs** (*List[Dict[str, int]]*): General Interface records, each
          with keys ``"S1"``, ``"S2"``, ``"GIF_TYPE"``, and ``"GIF_ORDER"``.
    """
    lb1 = (0, 0, 0)  # Preallocate lower bound tuple
    ub1 = (0, 0, 0)  # Preallocate upper bound tuple
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
                if idx % 2 == 0: # When there is a pair of blocks, write to the list
                    block_to_block.append({
                        'block1':{'block_index':blk_id1,
                                  'lb':lb1,
                                  'ub':ub1},
                        'block2':{'block_index':temp[0],
                                  'lb':(temp[1],temp[2],temp[3]),
                                  'ub':(temp[4],temp[5],temp[6])},
                    })
                else:
                    blk_id1 = temp[0]
                    lb1 = (temp[1], temp[2], temp[3])
                    ub1 = (temp[4], temp[5], temp[6])
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
                "lb":(temp[1],temp[2],temp[3]),
                "ub":(temp[4],temp[5],temp[6]),
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
              (row['block1']['ub'][0]-row['block1']['lb'][0])*
              (row['block1']['ub'][1]-row['block1']['lb'][1])*
              (row['block1']['ub'][2]-row['block1']['lb'][2]))
              for row in block_to_block]
        l2= [(row['block2']['block_index'],
              (row['block2']['ub'][0]-row['block2']['lb'][0])*
              (row['block2']['ub'][1]-row['block2']['lb'][1])*
              (row['block2']['ub'][2]-row['block2']['lb'][2]))
              for row in block_to_block]
        
        l1.extend(l2)
        l1 = np.array(l1)
        total_blocks = l1[:,0].max()
        _,indices = np.unique(l1[:,0],return_index=True)
        num_connections_foreach_block = l1[indices,:]
        num_connections_foreach_block.sort(axis=0)
        face_matches = block_to_block
        return total_blocks, face_matches,outer_faces,num_connections_foreach_block,nZones,zone_types,Zones,GIFs
