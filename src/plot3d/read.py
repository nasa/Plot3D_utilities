import numpy as np 
import os.path as osp
import struct
from typing import List, Optional
from .block import Block
from scipy.io import FortranFile
from tqdm import tqdm

def __read_plot3D_chunk_binary(f,IMAX:int,JMAX:int,KMAX:int, big_endian:bool=False,read_double:bool=True):
    """Reads and formats a binary chunk of data into a plot3D block.

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
    byte_order = '>' if big_endian else '<'
    dtype = np.dtype(f'{byte_order}f8' if read_double else f'{byte_order}f4')
    data = np.frombuffer(f.read(n * dtype.itemsize), dtype=dtype)
    return data.reshape((KMAX, JMAX, IMAX)).transpose(2, 1, 0).copy()

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
    while len(values) < n:
        line = f.readline()
        if not line:
            break
        values.extend(line.split())
    A = np.array(values[:n], dtype=np.float64)
    return A.reshape((KMAX, JMAX, IMAX)).transpose(2, 1, 0).copy()

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

    meshx = meshx.reshape(ints[1],ints[2],ints[0]).transpose(2,0,1)
    meshr = meshr.reshape(ints[1],ints[2],ints[0]).transpose(2,0,1)
    mesht = mesht.reshape(ints[1],ints[2],ints[0]).transpose(2,0,1)

    # Convert from x,r,theta to x,y,z
    z = meshr*np.sin(mesht)
    y = meshr*np.cos(mesht)

    return Block(X=meshx,Y=y,Z=z), nbld


def _detect_plot3d_format(filename: str):
    """Auto-detect Plot3D file format: ASCII/binary, endianness, precision, Fortran.

    Returns:
        dict with keys: binary, big_endian, read_double, fortran
    """
    file_size = osp.getsize(filename)

    # ── Try ASCII first ──
    try:
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            nblocks = int(first_line)
            if 1 <= nblocks <= 10000:
                # Read dimensions to confirm
                dims_line = f.readline()
                if dims_line:
                    return {'binary': False, 'big_endian': False, 'read_double': True, 'fortran': False}
    except (ValueError, UnicodeDecodeError):
        pass

    with open(filename, 'rb') as f:
        header = f.read(16)

    if len(header) < 8:
        # Fallback
        return {'binary': True, 'big_endian': False, 'read_double': True, 'fortran': False}

    def _try_binary(big_endian, fortran):
        """Try reading header as binary with given settings. Returns (nblocks, dims) or None."""
        fmt = '>' if big_endian else '<'
        try:
            with open(filename, 'rb') as f:
                if fortran:
                    # Fortran record: [rec_len][nblocks][rec_len]
                    rec_len = struct.unpack(f'{fmt}i', f.read(4))[0]
                    if rec_len != 4:
                        return None
                    nblocks = struct.unpack(f'{fmt}i', f.read(4))[0]
                    rec_end = struct.unpack(f'{fmt}i', f.read(4))[0]
                    if rec_end != 4 or not (1 <= nblocks <= 10000):
                        return None
                    # Read dims record
                    rec_len2 = struct.unpack(f'{fmt}i', f.read(4))[0]
                    if rec_len2 != nblocks * 3 * 4:
                        return None
                    dims = struct.unpack(f'{fmt}{nblocks*3}i', f.read(nblocks * 3 * 4))
                    rec_end2 = struct.unpack(f'{fmt}i', f.read(4))[0]
                    if rec_end2 != rec_len2:
                        return None
                else:
                    nblocks = struct.unpack(f'{fmt}I', f.read(4))[0]
                    if not (1 <= nblocks <= 10000):
                        return None
                    dims = struct.unpack(f'{fmt}{nblocks*3}I', f.read(nblocks * 3 * 4))

                # Validate dimensions
                for d in dims:
                    if not (1 <= d <= 100000):
                        return None
                return nblocks, dims
        except (struct.error, EOFError):
            return None

    def _expected_data_size(nblocks, dims, precision_bytes, fortran):
        """Compute expected file size for given precision."""
        header_size = 4 + nblocks * 3 * 4  # nblocks int + dimension ints
        data_size = 0
        for b in range(nblocks):
            ni, nj, nk = dims[b*3], dims[b*3+1], dims[b*3+2]
            data_size += 3 * ni * nj * nk * precision_bytes
        if fortran:
            # Record markers: 2*4 bytes per record
            # 1 record for nblocks, 1 for dims, nblocks records for data (one per block)
            n_records = 2 + nblocks  # nblocks record + dims record + one data record per block
            header_size += n_records * 2 * 4
        return header_size + data_size

    # ── Try each combination ──
    for fortran in [False, True]:
        for big_endian in [False, True]:
            result = _try_binary(big_endian, fortran)
            if result is None:
                continue
            nblocks, dims = result

            size_single = _expected_data_size(nblocks, dims, 4, fortran)
            size_double = _expected_data_size(nblocks, dims, 8, fortran)

            if file_size == size_double:
                return {'binary': True, 'big_endian': big_endian, 'read_double': True, 'fortran': fortran}
            if file_size == size_single:
                return {'binary': True, 'big_endian': big_endian, 'read_double': False, 'fortran': fortran}

    # Fallback
    return {'binary': True, 'big_endian': False, 'read_double': True, 'fortran': False}


def read_plot3D(filename:str, binary:Optional[bool]=None, big_endian:Optional[bool]=None, read_double:Optional[bool]=None, fortran:Optional[bool]=None):
    """Reads a Plot3D file and returns blocks.

    When any format parameter is ``None`` (the default), the file format is
    auto-detected by probing the header and comparing file size against
    expected sizes for single/double precision.

    Args:
        filename (str): Name of the file to read, e.g. ``.p3d``, ``.xyz`` or ``.plot3d``.
        binary (bool, optional): ``True`` for binary, ``False`` for ASCII. Auto-detected when ``None``.
        big_endian (bool, optional): Use big endian format for binary files. Auto-detected when ``None``.
        read_double (bool, optional): Read 8-byte doubles (``True``) or 4-byte floats (``False``). Auto-detected when ``None``.
        fortran (bool, optional): Read Fortran unformatted binary with record markers. Auto-detected when ``None``.

    Returns:
        List[Block]: List of blocks inside the Plot3D file.
    """
    # Auto-detect any unspecified parameters
    if any(p is None for p in [binary, big_endian, read_double, fortran]):
        detected = _detect_plot3d_format(filename)
        if binary is None:
            binary = detected['binary']
        if big_endian is None:
            big_endian = detected['big_endian']
        if read_double is None:
            read_double = detected['read_double']
        if fortran is None:
            fortran = detected['fortran']
    assert binary is not None and big_endian is not None
    assert read_double is not None and fortran is not None

    blocks = list()
    if osp.isfile(filename):
        if fortran:
            # Fortran unformatted binary with record markers
            endian = '>' if big_endian else '<'
            header_dtype = np.dtype(f'{endian}u4')
            int_dtype = np.dtype(f'{endian}i4')
            real_dtype = np.dtype(f'{endian}f8') if read_double else np.dtype(f'{endian}f4')

            with FortranFile(filename, 'r', header_dtype) as f:
                # Read nblocks
                nblocks = f.read_ints(int_dtype)[0]
                # Read all dimensions
                dims = f.read_ints(int_dtype)
                IMAX = dims[0::3]  # Every 3rd starting at 0
                JMAX = dims[1::3]  # Every 3rd starting at 1
                KMAX = dims[2::3]  # Every 3rd starting at 2
                # Read coordinate arrays. The standard PLOT3D convention writes
                # X, Y, Z as a single record per block; older files produced by
                # this library wrote them as three separate records. Detect
                # which layout is present from the size of the first record.
                for b in tqdm(range(nblocks), desc="Reading Fortran blocks", unit="block"):
                    npts = IMAX[b] * JMAX[b] * KMAX[b]
                    arr = f.read_reals(real_dtype)
                    shape = (IMAX[b], JMAX[b], KMAX[b])

                    if arr.size == 3 * npts:
                        X = arr[:npts].reshape(shape, order='F')
                        Y = arr[npts:2 * npts].reshape(shape, order='F')
                        Z = arr[2 * npts:].reshape(shape, order='F')
                    elif arr.size == npts:
                        X = arr.reshape(shape, order='F')
                        Y = f.read_reals(real_dtype).reshape(shape, order='F')
                        Z = f.read_reals(real_dtype).reshape(shape, order='F')
                    else:
                        raise ValueError(
                            f"Unexpected Fortran record size for block {b}: "
                            f"got {arr.size} reals, expected {npts} or {3 * npts}"
                        )

                    blocks.append(Block(X, Y, Z))
        elif binary:
            with open(filename,'rb') as f:
                nblocks = struct.unpack(">I",f.read(4))[0] if big_endian else struct.unpack("I",f.read(4))[0] # Read bytes
                IMAX = list(); JMAX = list(); KMAX = list()
                for b in range(nblocks):
                    if big_endian:
                        IMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        KMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                    else:
                        IMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                        JMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                        KMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes

                for b in tqdm(range(nblocks), desc="Reading binary blocks", unit="block"):
                    X = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Y = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    Z = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian,read_double)
                    b_temp = Block(X,Y,Z)
                    blocks.append(b_temp)
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
