import numpy as np
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
        batch_size (int, optional): unused, kept for API compatibility.

    See: https://docs.python.org/3/library/struct.html
    """
    byte_order = '>' if big_endian else '<'
    dtype = np.dtype(f'{byte_order}f8' if double_precision else f'{byte_order}f4')
    def write_var(V: np.ndarray):
        f.write(V.transpose(2, 1, 0).astype(dtype).tobytes())
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
        batch_size (int, optional): unused, kept for API compatibility.
    """
    fmt = '%.15f' if double_precision else '%.8f'

    def write_var(V: np.ndarray):
        flat = V.transpose(2, 1, 0).ravel()
        lines = []
        for start in range(0, len(flat), columns):
            chunk = flat[start:start + columns]
            lines.append(' '.join(fmt % v for v in chunk) + '\n')
        f.writelines(lines)
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
        # Fortran unformatted binary with record markers, standard PLOT3D
        # layout (one record per block containing X, Y, Z concatenated).
        endian = '>' if big_endian else '<'
        header_dtype = np.dtype(f'{endian}u4')
        int_dtype = np.dtype(f'{endian}i4')
        real_dtype = np.dtype(f'{endian}f8' if double_precision else f'{endian}f4')
        with FortranFile(filename, 'w', header_dtype) as f:
            f.write_record(np.array([len(blocks)], dtype=int_dtype))
            dims = np.array([[b.X.shape[0], b.X.shape[1], b.X.shape[2]]
                             for b in blocks], dtype=int_dtype).flatten()
            f.write_record(dims)
            for b in tqdm(blocks, desc="Writing Fortran blocks", unit="block"):
                xyz = np.concatenate([
                    b.X.flatten(order='F'),
                    b.Y.flatten(order='F'),
                    b.Z.flatten(order='F'),
                ]).astype(real_dtype)
                f.write_record(xyz)
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
