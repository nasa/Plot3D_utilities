from os import write
import numpy as np 
import os.path as osp
import struct
from typing import List

from pandas.core.indexing import need_slice
from .block import Block

def __write_plot3D_block_binary(f,B:Block):
    """Write binary plot3D block which contains X,Y,Z
        default format is Big-Endian

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
    """
    '''
        https://docs.python.org/3/library/struct.html
    '''
    def write_var(V:np.ndarray):
        for k in range(B.KMAX):
            for j in range(B.JMAX):
                for i in range(B.IMAX):
                    f.write(struct.pack('f',V[i,j,k]))
    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)


def __write_plot3D_block_ASCII(f,B:Block,columns:int=6):
    """Write plot3D block in ascii format 

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
        columns (int, optional): Number of columns in the file. Defaults to 6.
    """
    def write_var(V:np.ndarray):
        bNewLine = False
        indx = 0
        for k in range(B.KMAX):
            for j in range(B.JMAX):
                for i in range(B.IMAX):
                    f.write('{0:8.8f} '.format(V[i,j,k]))
                    bNewLine=False
                    indx+=1
                    if (indx % columns) == 0:
                        f.write('\n')
                        bNewLine=True
                    
        if not bNewLine:
            f.write('\n')
    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)

def write_plot3D(filename:str,blocks:List[Block],binary:bool=True):
    """Writes blocks to a Plot3D file

    Args:
        filename (str): name of the file to create 
        blocks (List[Block]): List containing all the blocks to write
        binary (bool, optional): Binary big endian. Defaults to True.
    """
    if binary:
        with open(filename,'wb') as f:
            f.write(struct.pack('I',len(blocks)))
            for b in blocks:
                IMAX,JMAX,KMAX = b.X.shape
                f.write(struct.pack('I',IMAX))
                f.write(struct.pack('I',JMAX))
                f.write(struct.pack('I',KMAX))
            for b in blocks:
                __write_plot3D_block_binary(f,b)
    else:
        with open(filename,'w') as f:
            f.write('{0:d}\n'.format(len(blocks)))
            for b in blocks:
                IMAX,JMAX,KMAX = b.X.shape
                f.write('{0:d} {1:d} {2:d}\n'.format(IMAX,JMAX,KMAX))            
            for b in blocks:
                __write_plot3D_block_ASCII(f,b)