import numpy as np
import os.path as osp
import struct
from typing import List
from tqdm import tqdm
from scipy.io import FortranFile
from .block import Block

def __write_plot3D_block_binary(f, B:Block, big_endian:bool=False, double_precision:bool=True, batch_size:int=100):
    """Write binary plot3D block which contains X,Y,Z.

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
        big_endian (bool): use big-endian byte order. Defaults to False (little-endian).
        double_precision (bool): writes to binary using double precision. Defaults to True.
        batch_size (int, optional): number of packed values to buffer before writing. Defaults to 100.

    See: https://docs.python.org/3/library/struct.html
    """
    def write_var(V:np.ndarray):
        endian = '>' if big_endian else '<'
        fmt = f'{endian}d' if double_precision else f'{endian}f'
        value_size = struct.calcsize(fmt)
        buffer = bytearray()
        buffer_limit = value_size * batch_size

        for k in range(B.KMAX):
            for j in range(B.JMAX):
                for i in range(B.IMAX):
                    buffer.extend(struct.pack(fmt, V[i,j,k]))
                    if len(buffer) >= buffer_limit:
                        f.write(buffer)
                        buffer.clear()

        if buffer:
            f.write(buffer)
    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)


def __write_plot3D_block_ASCII(f, B:Block, double_precision:bool=True, columns:int=6, batch_size:int=100):
    """Write plot3D block in ASCII format using scientific notation.

    Args:
        f (IO): file handle
        B (Block): writes a single block to a file
        double_precision (bool, optional): Use double precision format (15 decimals). Defaults to True.
        columns (int, optional): Number of columns in the file. Defaults to 6.
        batch_size (int, optional): number of lines to buffer before writing. Defaults to 100.
    """
    # Scientific notation format: width.precision E
    fmt = '{0:23.15E}' if double_precision else '{0:15.8E}'

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
                    line_entries.append(fmt.format(V[i,j,k]))
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

def write_plot3D(filename:str, blocks:List[Block], binary:bool=True, big_endian:bool=False,
                 double_precision:bool=True, fortran:bool=False, batch_size:int=100):
    """Writes blocks to a Plot3D file.

    Args:
        filename (str): name of the file to create
        blocks (List[Block]): List containing all the blocks to write
        binary (bool, optional): Write in binary format. Defaults to True.
        big_endian (bool, optional): Use big-endian byte order for binary. Defaults to False (little-endian).
        double_precision (bool, optional): Write using double precision (8-byte). Defaults to True.
        fortran (bool, optional): Write Fortran unformatted binary with record markers. Defaults to False.
        batch_size (int, optional): number of items (binary) or lines (ASCII) to buffer before writing. Defaults to 100.
    """
    if fortran:
        # Fortran unformatted binary with record markers
        dtype = np.float64 if double_precision else np.float32
        with FortranFile(filename, 'w') as f:
            # Write nblocks as one record
            f.write_record(np.array([len(blocks)], dtype=np.int32))
            # Write all dimensions in one record
            dims = np.array([[b.X.shape[0], b.X.shape[1], b.X.shape[2]]
                             for b in blocks], dtype=np.int32).flatten()
            f.write_record(dims)
            # Write each coordinate array as a record (Fortran column-major order)
            for b in tqdm(blocks, desc="Writing Fortran blocks", unit="block"):
                f.write_record(b.X.flatten(order='F').astype(dtype))
                f.write_record(b.Y.flatten(order='F').astype(dtype))
                f.write_record(b.Z.flatten(order='F').astype(dtype))
    elif binary:
        endian = '>' if big_endian else '<'
        with open(filename, 'wb') as f:
            f.write(struct.pack(f'{endian}I', len(blocks)))
            for b in blocks:
                IMAX, JMAX, KMAX = b.X.shape
                f.write(struct.pack(f'{endian}I', IMAX))
                f.write(struct.pack(f'{endian}I', JMAX))
                f.write(struct.pack(f'{endian}I', KMAX))
            for b in tqdm(blocks, desc="Writing binary blocks", unit="block"):
                __write_plot3D_block_binary(f, b, big_endian, double_precision, batch_size)
    else:
        with open(filename, 'w') as f:
            f.write('{0:d}\n'.format(len(blocks)))
            for b in blocks:
                IMAX, JMAX, KMAX = b.X.shape
                f.write('{0:d} {1:d} {2:d}\n'.format(IMAX, JMAX, KMAX))
            for b in tqdm(blocks, desc="Writing ASCII blocks", unit="block"):
                __write_plot3D_block_ASCII(f, b, double_precision=double_precision, batch_size=batch_size)
