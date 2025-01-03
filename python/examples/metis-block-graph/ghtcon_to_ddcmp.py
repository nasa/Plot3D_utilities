import os
from typing import List, Optional
import numpy as np 
import metis 
from plot3d.graph import block_connectivity_to_graph
import pickle

def dump_data(data,filename:str):
    with open(filename,'wb') as f:
        pickle.dump(data,f)

def read_data(filename:str):
    with open(filename,'rb') as f:
        return pickle.load(f)
        
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
            temp = [int(t) for t in temp if t.strip() != '']
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
        # Read Faces
        nFaces = int(line.lstrip().rstrip().split(" ")[0])
        [f.readline() for _ in range(int(nFaces))] # Skip outer faces
        
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
        return total_blocks, face_matches,num_connections_foreach_block,nZones,zone_types,Zones,GIFs

def create_ddcmp(ght_conn:str,ddcmp_file:str='ddcmp.dat',nparts:int=3,pkl_file_blocksizes:Optional[str]=None):
    """_summary_

    Args:
        ght_conn (str): GlennHT Connectivity File
        ddcmp_file (str, optional): Name of the file to write to. Defaults to 'ddcmp.dat'.
        nparts (int, optional): Number of partitions. Defaults to 3.
        pkl_file_blocksizes (Optional[str], optional): Pickle file containing blocksizes. Defaults to None.
    """
    def write_ddcmp(parts:List[int],nZones:int,filename:str='ddcmp.dat'):
        """Reads in the partitions from metis and exports to ddcmp file for solver to split the work 

        Args:
            parts (List[int]): len = number of blocks, each element is the partition number
            nZones (int): Number of zones
            filename (str): Name of the file to write to

        Returns:
            np.ndarray: Number of connections between partitions
        """
        nProc = max(parts)+1 # Fortran starts at 1
        nISP = nZones
        nBlocks = len(parts)
        with open(filename,'w') as f:
            f.write(f'{nProc:d}\n')
            f.write(f'{nISP:d}\n') # ISP or number of zones
            f.write(f'{nBlocks:d}\n')
            # Write which block goes to which ISP 
            for isp in range(nISP):
                for bIdx in range(nBlocks):
                    f.write(f'{bIdx+1:d} {isp+1:d}\n')
            # Write which block goes to which Processor 
            for block_indx,p in enumerate(parts):
                f.write(f'{block_indx+1:d} {p+1:d}\n')
           
    max_block_index, face_matches,connections,nZones,zone_types,Zones,GIFs = read_ght_conn(ght_conn)
    G = block_connectivity_to_graph(face_matches,connections[:,1])
    
    # G.graph['node_weight_attr'] = ['weight']
    G.graph['edge_weight_attr'] = 'weight'
    (edgecuts, parts) = metis.part_graph(G, nparts,tpwgts=None)
    
    if not pkl_file_blocksizes:
        blocksizes = read_data(pkl_file_blocksizes)         # Useful for knowing how many nodes per partition
        nodes_per_part = list()
        # Print number of cells for each part
        for i in range(nparts):
            indexes = [p for p, e in enumerate(parts) if e == i]
            nodes = 0 
            for idx in indexes:  
                nodes += blocksizes[idx-1]
            nodes_per_part.append(nodes)

    for i in range(nparts):
        print(f'Parition {i} has {parts.count(i)} blocks')

    # Determine work for each partition
    communication_work = np.zeros((nparts,))
    partition_edge_weights = np.zeros((nparts,))
    for b in range(1,max_block_index+1):
        partition_id = parts[b]
        for connected_block in G.adj[b]:
            connected_block_partition = parts[connected_block]
            edge_weight = G.adj[b][connected_block]['weight']
            if connected_block_partition != partition_id:
                communication_work[partition_id]+=1
                partition_edge_weights[partition_id] += edge_weight

    with open('partition.out','w') as f:
        f.write(f'Number of partitions {nparts}\n')
        for i in range(nparts):
            f.write(f'Parition {i} has communication work {communication_work[i]} edge_work {partition_edge_weights[i]}\n')
    
    write_ddcmp(parts[1:],nZones=nZones,filename=ddcmp_file) # Fortran starts at 1 

if __name__=="__main__":
    create_ddcmp("python/examples/metis-block-graph/CMC9/conn.ght_conn",ddcmp_file='python/examples/metis-block-graph/CMC9/ddcmp.dat',nparts=3,pkl_file_blocksizes='python/examples/metis-block-graph/CMC9/blocksizes.pickle')
    create_ddcmp("python/examples/metis-block-graph/EEE-Stator/conn.ght_conn",ddcmp_file='python/examples/metis-block-graph/EEE-Stator/ddcmp.dat',nparts=3,pkl_file_blocksizes='python/examples/metis-block-graph/EEE-Stator/blocksizes.pickle')




       

