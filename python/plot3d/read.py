import numpy as np
import os.path as osp
import struct
from typing import List
from .block import Block
from scipy.io import FortranFile
from tqdm import tqdm


def __read_plot3D_chunk_binary(f, IMAX: int, JMAX: int, KMAX: int,
                               big_endian: bool = False, read_double: bool = True):
    """Reads and formats a binary chunk of data into a plot3D block.

    Uses vectorized numpy reads instead of element-by-element struct.unpack
    for dramatically faster I/O (10-100x speedup on large meshes).

    Args:
        f (io): file handle
        IMAX (int): maximum I index
        JMAX (int): maximum J index
        KMAX (int): maximum K index
        big_endian (bool, optional): Use big endian format for reading binary files. Defaults False.
        read_double (bool, optional): When ``True`` read 8-byte doubles, otherwise read 4-byte floats.

    Returns:
        numpy.ndarray: Plot3D variable either X,Y, or Z
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
    """Continously read a word from an ascii file

    Args:
        f (io): file handle

    Yields:
        float: value from ascii file
    """
    for line in f:
        line = line.strip().replace('\n','').split(' ')
        tokenArray = [float(entry) for entry in line if entry]
        for token in tokenArray:
            yield token

def __read_plot3D_chunk_ASCII(f,IMAX:int,JMAX:int,KMAX:int):
    """Reads and formats an ASCII chunk of data into a plot3D block.

    Args:
        f (io): file handle
        IMAX (int): maximum I index
        JMAX (int): maximum J index
        KMAX (int): maximum K index

    Returns:
        numpy.ndarray: Plot3D variable either X,Y, or Z
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

def read_ap_nasa(filename:str):
    """Reads an AP NASA File and converts it to Block format which can be exported to a plot3d file
        AP NASA file represents a single block. The first 7 integers are il,jl,kl,ile,ite,jtip,nbld

    Args:
        filename (str): location of the .ap file

    Returns:
        Tuple containing:

            *block* (Block): file in block format
            *nbld* (int): Number of blades
    """

    f = FortranFile(filename, 'r')

    ints = f.read_ints(np.int32)
    idim = np.array([ints[0],ints[1],ints[2]])
    mdim = np.array([3,ints[0]*ints[2]])
    il = ints[0]
    jl = ints[1]
    kl = ints[2]
    jdim = jl

    ile = ints[3]
    ite = ints[4]
    jtip = ints[5]
    nbld = ints[6]

    for j in range(0,jdim):
        jmeshxrt = f.read_reals(dtype='f4').reshape(mdim)
        meshi    = np.array(jmeshxrt[0,:])
        meshj    = np.array(jmeshxrt[1,:])
        meshk    = np.array(jmeshxrt[2,:])
        if j == 0:
            meshx   = meshi
            meshr   = meshj
            mesht   = meshk
        else:
            meshx = np.append(meshx,meshi)
            meshr = np.append(meshr,meshj)
            mesht = np.append(mesht,meshk)

    meshx = meshx.reshape(ints[1],ints[2],ints[0])
    meshr = meshr.reshape(ints[1],ints[2],ints[0])
    mesht = mesht.reshape(ints[1],ints[2],ints[0])

    # Convert from x,r,theta to x,y,z
    z = meshr*np.sin(mesht)
    y = meshr*np.cos(mesht)

    return Block(X=meshx,Y=y,Z=z), nbld


def read_plot3D(filename:str, binary:bool=True, big_endian:bool=False, read_double:bool=True, fortran:bool=False):
    """Reads a Plot3D file and returns blocks.

    Args:
        filename (str): Name of the file to read, e.g. ``.p3d``, ``.xyz`` or ``.plot3d``.
        binary (bool, optional): Indicates if the file is binary. Defaults to True.
        big_endian (bool, optional): Use big endian format when reading binary files. Defaults to False.
        read_double (bool, optional): Read 8-byte doubles when ``True`` and 4-byte floats otherwise.
        fortran (bool, optional): Read Fortran unformatted binary with record markers. Defaults to False.

    Returns:
        List[Block]: List of blocks inside the Plot3D file.
    """

    blocks = list()
    if osp.isfile(filename):
        if fortran:
            # Fortran unformatted binary with record markers
            dtype = 'f8' if read_double else 'f4'
            with FortranFile(filename, 'r') as f:
                # Read nblocks
                nblocks = f.read_ints('i4')[0]
                # Read all dimensions
                dims = f.read_ints('i4')
                IMAX = dims[0::3]  # Every 3rd starting at 0
                JMAX = dims[1::3]  # Every 3rd starting at 1
                KMAX = dims[2::3]  # Every 3rd starting at 2
                # Read coordinate arrays
                for b in tqdm(range(nblocks), desc="Reading Fortran blocks", unit="block"):
                    X = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    Y = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    Z = f.read_reals(dtype).reshape((IMAX[b], JMAX[b], KMAX[b]), order='F')
                    blocks.append(Block(X, Y, Z))
        elif binary:
            with open(filename, 'rb') as f:
                # Read nblocks
                endian = '>' if big_endian else '<'
                nblocks = struct.unpack(f'{endian}I', f.read(4))[0]
                # Read all dimensions at once for efficiency
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
            with open(filename,'r') as f:
                nblocks = int(f.readline())
                IMAX = list(); JMAX = list(); KMAX = list()

                for b in range(nblocks):
                    IJK = f.readline().replace('\n','').split(' ')
                    tokens = [int(w) for w in IJK if w]
                    IMAX.append(tokens[0])
                    JMAX.append(tokens[1])
                    KMAX.append(tokens[2])

                for b in tqdm(range(nblocks), desc="Reading ASCII blocks", unit="block"):
                    X = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Y = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    Z = __read_plot3D_chunk_ASCII(f,IMAX[b],JMAX[b],KMAX[b])
                    b_temp = Block(X,Y,Z)
                    blocks.append(b_temp)
    return blocks
