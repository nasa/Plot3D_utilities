"""Write Plot3D multi-block structured grid files.

Supports three file formats:

- **Binary** (default) -- raw byte streams with optional big-endian and
  single/double precision control.
- **ASCII** -- whitespace-delimited scientific-notation text files.
- **Fortran unformatted** -- binary files with Fortran record markers,
  written via :mod:`scipy.io.FortranFile`.

See also :mod:`plot3d.read` for the corresponding readers.
"""
import numpy as np
import os.path as osp
import struct
from typing import List
from tqdm import tqdm
from scipy.io import FortranFile
from .block import Block


def __write_plot3D_block_binary(f, B: Block, big_endian: bool = False,
                                double_precision: bool = True, batch_size: int = 100):
    """Write a single block's X, Y, Z arrays to a binary Plot3D file.

    Uses vectorized NumPy writes for fast I/O. Coordinates are transposed
    from C order ``(I, J, K)`` to Fortran order ``(K, J, I)`` before writing.

    Args:
        f: Open binary file handle.
        B: Block to write.
        big_endian: If ``True``, write in big-endian byte order.
            Defaults to ``False`` (little-endian).
        double_precision: If ``True``, write 8-byte doubles; otherwise
            4-byte floats. Defaults to ``True``.
        batch_size: Unused (kept for API compatibility). Defaults to 100.
    """
    endian = '>' if big_endian else '<'
    dtype_char = 'f8' if double_precision else 'f4'
    dtype = np.dtype(f'{endian}{dtype_char}')

    def write_var(V: np.ndarray):
        # Plot3D stores in Fortran order: transpose (I,J,K) -> (K,J,I) then flatten
        data = V.transpose(2, 1, 0).astype(dtype).tobytes()
        f.write(data)

    write_var(B.X)
    write_var(B.Y)
    write_var(B.Z)


def __write_plot3D_block_ASCII(f, B: Block, double_precision: bool = True,
                               columns: int = 6, batch_size: int = 100):
    """Write a single block's X, Y, Z arrays to an ASCII Plot3D file.

    Values are written in scientific notation, with configurable column
    count and buffered line output for performance.

    Args:
        f: Open text file handle.
        B: Block to write.
        double_precision: If ``True``, use 15-decimal precision; otherwise
            8-decimal. Defaults to ``True``.
        columns: Number of values per line. Defaults to 6.
        batch_size: Number of lines to buffer before flushing.
            Defaults to 100.
    """
    fmt = '{0:23.15f}' if double_precision else '{0:15.8f}'

    def write_var(V: np.ndarray):
        # Flatten in Fortran order (K varies slowest, I fastest)
        flat = V.transpose(2, 1, 0).ravel()
        line_entries = []
        line_batch = []

        def flush_batch():
            if line_batch:
                f.writelines(line_batch)
                line_batch.clear()

        for val in flat:
            line_entries.append(fmt.format(float(val)))
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


def write_plot3D(filename: str, blocks: List[Block], binary: bool = True,
                 big_endian: bool = False, double_precision: bool = True,
                 fortran: bool = False, batch_size: int = 100):
    """Write a list of blocks to a multi-block Plot3D grid file.

    The file header contains the number of blocks followed by the
    ``(IMAX, JMAX, KMAX)`` dimensions for each block, then the X, Y, Z
    coordinate arrays in Fortran (column-major) order.

    Args:
        filename: Output file path (``.p3d``, ``.xyz``, or ``.plot3d``).
        blocks: List of blocks to write.
        binary: If ``True``, write as raw binary. Defaults to ``True``.
        big_endian: Use big-endian byte order for binary writes.
            Defaults to ``False`` (little-endian).
        double_precision: Write 8-byte doubles (``True``) or 4-byte floats
            (``False``). Defaults to ``True``.
        fortran: Write Fortran unformatted binary with record markers.
            Defaults to ``False``.
        batch_size: Buffer size for batched writing. Defaults to 100.
    """
    if fortran:
        dtype = np.float64 if double_precision else np.float32
        with FortranFile(filename, 'w') as f:
            f.write_record(np.array([len(blocks)], dtype=np.int32))
            dims = np.array([[b.X.shape[0], b.X.shape[1], b.X.shape[2]]
                             for b in blocks], dtype=np.int32).flatten()
            f.write_record(dims)
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
                f.write(struct.pack(f'{endian}III', IMAX, JMAX, KMAX))
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
