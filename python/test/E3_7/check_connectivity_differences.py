from typing import List
import numpy as np
import sys

def read_connectivity_file(filename:str):
    connections = list()
    with open(filename,'r') as f:
        nConnectionPairs = int(f.readline())        
        for i in range(nConnectionPairs):
            pair1 = [int(w.replace('\n','')) for w in f.readline().split(' ') if w]
            pair2 = [int(w.replace('\n','')) for w in f.readline().split(' ') if w]
            connections.append({ "block1":{ 'index':pair1[0], 
                                    'IMIN':pair1[1],'JMIN':pair1[2],'KMIN':pair1[3],
                                    'IMAX':pair1[4],'JMAX':pair1[5],'KMAX':pair1[6] },
                                  "block2":{ 'index':pair2[0], 
                                    'IMIN':pair2[1],'JMIN':pair2[2],'KMIN':pair2[3],
                                    'IMAX':pair2[4],'JMAX':pair2[5],'KMAX':pair2[6] }
                                })
    return connections

def compare_connectivity(con1:List,con2:List):
    """
        returns stuff in list 1 missing from list 2 
    """

    # check what connections in 1 are present in 2
    matching_con = list() 
    nonmatching_con = list()
    for p in range(len(con1)):
        bMatching = False
        pairs_1 = con1[p]    
        p1_block_1_index = pairs_1['block1']['index']
        p1_block_1_corner = np.array([pairs_1['block1']['IMIN'], pairs_1['block1']['JMIN'], pairs_1['block1']['KMIN'],
                                        pairs_1['block1']['IMAX'], pairs_1['block1']['JMAX'], pairs_1['block1']['KMAX']])
        
        p1_block_2_index = pairs_1['block2']['index']
        p1_block_2_corner = np.array([pairs_1['block2']['IMIN'], pairs_1['block2']['JMIN'], pairs_1['block2']['KMIN'],
                                        pairs_1['block2']['IMAX'], pairs_1['block2']['JMAX'], pairs_1['block2']['KMAX']])
        
        
        for q in range(len(con2)):
            pairs_2 = con2[q]  
            p2_block_1_index = pairs_2['block1']['index']
            p2_block_1_corner = np.array([pairs_2['block1']['IMIN'], pairs_2['block1']['JMIN'], pairs_2['block1']['KMIN'],
                                            pairs_2['block2']['IMAX'], pairs_2['block1']['JMAX'], pairs_2['block1']['KMAX']])
            
            p2_block_2_index = pairs_2['block2']['index']
            p2_block_2_corner = np.array([pairs_2['block2']['IMIN'], pairs_2['block2']['JMIN'], pairs_2['block2']['KMIN'],
                                            pairs_2['block2']['IMAX'], pairs_2['block2']['JMAX'], pairs_2['block2']['KMAX']])

            # Check for matches 
            ## Check if the match is diagonal 
            if (p1_block_1_index == p2_block_1_index and p1_block_2_index==p2_block_2_index):
                
                if (p1_block_1_corner-p2_block_1_corner).sum() == 0 and (p1_block_2_corner-p2_block_2_corner).sum() == 0:
                    matching_con.append((p,q)) 
                    bMatching = True

            elif (p1_block_1_index == p2_block_2_index and p1_block_2_index==p2_block_1_index):
                if (p1_block_1_corner-p2_block_2_corner).sum() == 0 and (p2_block_2_corner-p1_block_1_corner).sum() == 0:
                    matching_con.append((p,q)) 
                    bMatching = True

        if not bMatching:
            nonmatching_con.append(p)

    return matching_con,nonmatching_con

def write_connectivity(matches):
    blocks = ['block1','block2']
    with open('test' + '.ght_conn','w') as f:
        # Print matches
        nMatches = len(matches)
        f.write(f'{nMatches}\n') # Print number of matches 
        for match in matches:                        
            for block in blocks:
                block_indx = match[block]['index'] # block1 and block2 are arbitrary names, the key is the block index 
                block_IMIN = match[block]['IMIN']
                block_JMIN = match[block]['JMIN']
                block_KMIN = match[block]['KMIN']

                block_IMAX = match[block]['IMAX']
                block_JMAX = match[block]['JMAX']
                block_KMAX = match[block]['KMAX']

                f.write(f"{block_indx:3d}\t{block_IMIN:5d} {block_JMIN:5d} {block_KMIN:5d}\t{block_IMAX:5d} {block_JMAX:5d} {block_KMAX:5d}\n")        

ali = read_connectivity_file('ali.ght_conn')
paht = read_connectivity_file('E3_7_Assembly_6_blocks.ght_conn')

matching,non_matching = compare_connectivity(ali,paht)
print_non_matching = [ali[i] for i in non_matching]
write_connectivity(print_non_matching)

print('done')