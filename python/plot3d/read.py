"""Read Plot3D multi-block structured grid files.

Supports three file formats:

- **Binary** (default) -- raw byte streams with optional big-endian and
  single/double precision control.
- **ASCII** -- whitespace-delimited scientific-notation text files.
- **Fortran unformatted** -- binary files with Fortran record markers,
  read via :mod:`scipy.io.FortranFile`.

All readers return a list of :class:`~plot3d.block.Block` objects.
"""
import numpy as np
import os.path as osp
import struct
from typing import List
from .block import Block
from scipy.io import FortranFile
from tqdm import tqdm


def __read_plot3D_chunk_binary(f, IMAX: int, JMAX: int, KMAX: int,
                               big_endian: bool = False, read_double: bool = True):
    """Read a single coordinate variable from a binary Plot3D file.

    Uses vectorized NumPy reads for fast I/O (10-100x faster than
    element-by-element ``struct.unpack`` on large meshes).

    Args:
        f: Open binary file handle positioned at the start of the chunk.
        IMAX: Number of nodes in the I direction.
        JMAX: Number of nodes in the J direction.
        KMAX: Number of nodes in the K direction.
        big_endian: If ``True``, read in big-endian byte order.
            Defaults to ``False`` (little-endian).
        read_double: If ``True``, read 8-byte doubles; otherwise 4-byte
            floats. Defaults to ``True``.

    Returns:
        numpy.ndarray: 3-D coordinate array with shape ``(IMAX, JMAX, KMAX)``
        in C order (transposed from the Fortran-order storage).
    """
    n = IMAX * JMAX * KMAX
    endian = '>' if big_endian else '<'
    dtype_char = 'f8' if read_double else 'f4'
    dtype = np.dtype(f'{endian}{dtype_char}')
    byte_size = dtype.itemsize * n
    buf = f.read(byte_size)
    A = np.frombuffer(buf, dtype=dtype, count=n)
    # Plot3D stores data in Fortran order: K varies slowest, I varies fastest
    # After reading flat array, reshape as (K,J,I) then transpose to (I,J,K)
    A = A.reshape((KMAX, JMAX, IMAX)).transpose(2, 1, 0).astype(np.float64)
    return A


def read_word(f):
    """Yield individual float values from an ASCII file, line by line.

    Splits each line on whitespace and yields each token as a float.
    Used internally by :func:`__read_plot3D_chunk_ASCII`.

    Args:
        f: Open text file handle.

    Yields:
        float: The next numeric value from the file.
    """
    for line in f:
        line = line.strip().replace('\n', '').split(' ')
        tokenArray = [float(entry) for entry in line if entry]
        for token in tokenArray:
            yield token


def __read_plot3D_chunk_ASCII(f, IMAX: int, JMAX: int, KMAX: int):
    """Read a single coordinate variable from an ASCII Plot3D file.

    Args:
        f: Open text file handle positioned after the header.
        IMAX: Number of nodes in the I direction.
        JMAX: Number of nodes in the J direction.
        KMAX: Number of nodes in the K direction.

    Returns:
        numpy.ndarray: 3-D coordinate array with shape ``(IMAX, JMAX, KMAX)``.
    """
    n = IMAX * JMAX * KMAX
    values = []
    for w in read_word(f):
        values.append(w)
        if len(values) >= n:
            break

    A = np.array(values, dtype=np.float64).reshape((KMAX, JMAX, IMAX))
    A = np.transpose(A, [2, 1, 0])
    return A


def read_ap_nasa(filename: str):
    """Read an AP NASA file and convert it to a Plot3D :class:`~plot3d.block.Block`.

    The AP NASA format stores a single-block blade grid in cylindrical
    coordinates ``(x, r, theta)``. The first record contains seven integers:
    ``il, jl, kl, ile, ite, jtip, nbld``.

    Args:
        filename: Path to the ``.ap`` file.

    Returns:
        tuple: A 2-element tuple containing:

            - **block** (:class:`~plot3d.block.Block`): The grid converted
              to Cartesian coordinates.
            - **nbld** (int): Number of blades.
    """
    f = FortranFile(filename, 'r')

    ints = f.read_ints(np.int32)
    idim = np.array([ints[0], ints[1], ints[2]])
    mdim = np.array([3, ints[0] * ints[2]])
    il = ints[0]
    jl = ints[1]
    kl = ints[2]
    jdim = jl

    ile = ints[3]
    ite = ints[4]
    jtip = ints[5]
    nbld = ints[6]

    for j in range(0, jdim):
        jmeshxrt = f.read_reals(dtype='f4').reshape(mdim)
        meshi = np.array(jmeshxrt[0, :])
        meshj = np.array(jmeshxrt[1, :])
        meshk = np.array(jmeshxrt[2, :])
        if j == 0:
            meshx = meshi
            meshr = meshj
            mesht = meshk
        else:
            meshx = np.append(meshx, meshi)
            meshr = np.append(meshr, meshj)
            mesht = np.append(mesht, meshk)

    meshx = meshx.reshape(ints[1], ints[2], ints[0])
    meshr = meshr.reshape(ints[1], ints[2], ints[0])
    mesht = mesht.reshape(ints[1], ints[2], ints[0])

    # Convert from x,r,theta to x,y,z
    z = meshr * np.sin(mesht)
    y = meshr * np.cos(mesht)

    return Block(X=meshx, Y=y, Z=z), nbld


def read_plot3D(filename: str, binary: bool = True, big_endian: bool = False,
                read_double: bool = True, fortran: bool = False):
    """Read a multi-block Plot3D grid file.

    Supports binary, ASCII, and Fortran unformatted formats. The file
    header contains the number of blocks followed by the ``(IMAX, JMAX,
    KMAX)`` dimensions for each block, then the X, Y, Z coordinate arrays
    in Fortran (column-major) order.

    Args:
        filename: Path to the file (``.p3d``, ``.xyz``, or ``.plot3d``).
        binary: If ``True``, read as raw binary. Defaults to ``True``.
        big_endian: Use big-endian byte order for binary reads.
            Defaults to ``False``.
        read_double: Read 8-byte doubles (``True``) or 4-byte floats
            (``False``). Defaults to ``True``.
        fortran: Read Fortran unformatted binary with record markers.
            Defaults to ``False``.

    Returns:
        List[Block]: List of blocks read from the file. Returns an empty
        list if the file does not exist.
    """
    blocks = list()
    if osp.isfile(filename):
        if fortran:
            # Fortran unformatted binary with record markers
            dtype = 'f8' if read_double else 'f4'
            with FortranFile(filename, 'r') as f:
                nblocks = f.read_ints('i4')[0]
                dims = f.read_ints('i4')
                IMAX = dims[0::3]
                JMAX = dims[1::3]
                KMAX = dims[2::3]
                for b in tqdm(range(nblocks), desc="Reading Fortran blocks", unit="block"):
                    X = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    Y = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    Z = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    blocks.append(Block(X, Y, Z))
        elif binary:
            with open(filename, 'rb') as f:
                endian = '>' if big_endian else '<'
                nblocks = struct.unpack(f'{endian}I', f.read(4))[0]
                dim_data = f.read(nblocks * 3 * 4)
                dims = struct.unpack(f'{endian}{nblocks * 3}I', dim_data)
                IMAX = [dims[b * 3] for b in range(nblocks)]
                JMAX = [dims[b * 3 + 1] for b in range(nblocks)]
                KMAX = [dims[b * 3 + 2] for b in range(nblocks)]

                for b in tqdm(range(nblocks), desc="Reading binary blocks", unit="block"):
                    X = __read_plot3D_chunk_binary(f, IMAX[b], JMAX[b], KMAX[b], big_endian, read_double)
                    Y = __read_plot3D_chunk_binary(f, IMAX[b], JMAX[b], KMAX[b], big_endian, read_double)
                    Z = __read_plot3D_chunk_binary(f, IMAX[b], JMAX[b], KMAX[b], big_endian, read_double)
                    blocks.append(Block(X, Y, Z))
        else:
            with open(filename, 'r') as f:
                nblocks = int(f.readline())
                IMAX = list(); JMAX = list(); KMAX = list()

                for b in range(nblocks):
                    IJK = f.readline().replace('\n', '').split(' ')
                    tokens = [int(w) for w in IJK if w]
                    IMAX.append(tokens[0])
                    JMAX.append(tokens[1])
                    KMAX.append(tokens[2])

                for b in tqdm(range(nblocks), desc="Reading ASCII blocks", unit="block"):
                    X = __read_plot3D_chunk_ASCII(f, IMAX[b], JMAX[b], KMAX[b])
                    Y = __read_plot3D_chunk_ASCII(f, IMAX[b], JMAX[b], KMAX[b])
                    Z = __read_plot3D_chunk_ASCII(f, IMAX[b], JMAX[b], KMAX[b])
                    b_temp = Block(X, Y, Z)
                    blocks.append(b_temp)
    return blocks
