import os
import numpy as np 

def glennht_to_con(filename:str):
    with open(filename, "r") as f:        
        line = f.readline()
        pairs = int(line.lstrip().rstrip().split(" ")[0])
        
        idx = 0; blk_id1 = 0
        block_to_block = list() 
        for line in f.readlines():    
            idx += 1
            temp = line.lstrip().rstrip().split(" ")
            temp = [int(t) for t in temp if t.strip() != '']
            if (len(temp) == 7):
                if idx % 2 == 0:
                    block_to_block.append({
                        'block1':{'block_index':blk_id1},
                        'block2':{'block_index':temp[0]},
                    })
                                            
                else:
                    blk_id1 = temp[0]
                
            if idx>pairs*2:
                break
        return block_to_block
    return None
            
                


       

