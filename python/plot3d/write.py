from os import write
import numpy as np 
import os.path as osp
import struct
from typing import List
from tqdm import tqdm
from .block import Block

def __write_plot3D_block_binary(f,B:Block,double_precision:bool=True,batch_size:int=100):
    """Write binary plot3D block which contains X,Y,Z
        default format is Big-Endian

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
        double_precision (bool): writes to binary using double precision
        batch_size (int, optional): number of packed values to buffer before writing. Defaults to 100.
    """
    '''
        https://docs.python.org/3/library/struct.html
    '''
    def write_var(V:np.ndarray):
        fmt = '<d' if double_precision else '<f'
        value_size = struct.calcsize(fmt)
        buffer = bytearray()
        buffer_limit = value_size * batch_size

        for k in range(B.KMAX):
            for j in range(B.JMAX):
                for i in range(B.IMAX):
                    buffer.extend(struct.pack(fmt,V[i,j,k]))
                    if len(buffer) >= buffer_limit:
                        f.write(buffer)
                        buffer.clear()

        if buffer:
            f.write(buffer)
    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)


def __write_plot3D_block_ASCII(f,B:Block,columns:int=6,batch_size:int=100):
    """Write plot3D block in ascii format 

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
        columns (int, optional): Number of columns in the file. Defaults to 6.
        batch_size (int, optional): number of lines to buffer before writing. Defaults to 100.
    """
    def write_var(V:np.ndarray):
        line_entries = []
        line_batch = []

        def flush_batch():
            if line_batch:
                f.writelines(line_batch)
                line_batch.clear()

        for k in range(B.KMAX):
            for j in range(B.JMAX):
                for i in range(B.IMAX):
                    line_entries.append('{0:8.8f}'.format(V[i,j,k]))
                    if len(line_entries) == columns:
                        line_batch.append(' '.join(line_entries) + '\n')
                        line_entries.clear()
                        if len(line_batch) == batch_size:
                            flush_batch()

        if line_entries:
            line_batch.append(' '.join(line_entries) + '\n')

        flush_batch()
    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)

def write_plot3D(filename:str,blocks:List[Block],binary:bool=True,double_precision:bool=True,batch_size:int=100):
    """Writes blocks to a Plot3D file

    Args:
        filename (str): name of the file to create 
        blocks (List[Block]): List containing all the blocks to write
        binary (bool, optional): Binary big endian. Defaults to True.
        double_precision (bool, optional). Writes to binary file using double precision. Defaults to True
        batch_size (int, optional): number of items (binary) or lines (ASCII) to buffer before writing. Defaults to 100.
    """
    if binary:
        with open(filename,'wb') as f:
            f.write(struct.pack('I',len(blocks)))
            for b in blocks:
                IMAX,JMAX,KMAX = b.X.shape
                f.write(struct.pack('I',IMAX))
                f.write(struct.pack('I',JMAX))
                f.write(struct.pack('I',KMAX))
            for b in tqdm(blocks, desc="Writing binary blocks", unit="block"):
                __write_plot3D_block_binary(f,b,double_precision,batch_size)
    else:
        with open(filename,'w') as f:
            f.write('{0:d}\n'.format(len(blocks)))
            for b in blocks:
                IMAX,JMAX,KMAX = b.X.shape
                f.write('{0:d} {1:d} {2:d}\n'.format(IMAX,JMAX,KMAX))            
            for b in tqdm(blocks, desc="Writing ASCII blocks", unit="block"):
                __write_plot3D_block_ASCII(f,b,batch_size=batch_size)
