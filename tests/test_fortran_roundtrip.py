"""Round-trip and legacy-format coverage for the Fortran unformatted
PLOT3D I/O path (``fortran=True``).

The standard PLOT3D / FUN3D layout writes X, Y, Z as a single Fortran
record per block. Earlier versions of this library wrote them as three
separate records; the reader must continue to handle those existing
files.
"""

import os
import tempfile

import numpy as np
import pytest
from scipy.io import FortranFile

from plot3d import Block, read_plot3D, write_plot3D


def _make_blocks():
    rng = np.random.default_rng(0)
    shapes = [(3, 4, 2), (5, 2, 3)]
    blocks = []
    for shape in shapes:
        X = rng.random(shape).astype(np.float64)
        Y = rng.random(shape).astype(np.float64)
        Z = rng.random(shape).astype(np.float64)
        blocks.append(Block(X, Y, Z))
    return blocks


def _assert_blocks_equal(a, b, atol):
    assert len(a) == len(b)
    for ba, bb in zip(a, b):
        np.testing.assert_allclose(ba.X, bb.X, atol=atol)
        np.testing.assert_allclose(ba.Y, bb.Y, atol=atol)
        np.testing.assert_allclose(ba.Z, bb.Z, atol=atol)


@pytest.mark.parametrize("big_endian", [False, True])
@pytest.mark.parametrize("double_precision", [False, True])
def test_fortran_roundtrip(big_endian, double_precision):
    blocks = _make_blocks()
    atol = 0.0 if double_precision else 1e-6

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "mesh.xyz")
        write_plot3D(path, blocks, binary=True, big_endian=big_endian,
                     double_precision=double_precision, fortran=True)
        read_back = read_plot3D(path, binary=True, big_endian=big_endian,
                                read_double=double_precision, fortran=True)

    _assert_blocks_equal(blocks, read_back, atol)


@pytest.mark.parametrize("big_endian", [False, True])
def test_fortran_reader_accepts_legacy_three_record_layout(big_endian):
    """Files written by older versions of this library used three Fortran
    records per block (one each for X, Y, Z). The reader must still
    handle them."""
    blocks = _make_blocks()
    endian = '>' if big_endian else '<'
    header_dtype = np.dtype(f'{endian}u4')
    int_dtype = np.dtype(f'{endian}i4')
    real_dtype = np.dtype(f'{endian}f8')

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "legacy.xyz")
        with FortranFile(path, 'w', header_dtype) as f:
            f.write_record(np.array([len(blocks)], dtype=int_dtype))
            dims = np.array(
                [[b.X.shape[0], b.X.shape[1], b.X.shape[2]] for b in blocks],
                dtype=int_dtype,
            ).flatten()
            f.write_record(dims)
            for b in blocks:
                f.write_record(b.X.flatten(order='F').astype(real_dtype))
                f.write_record(b.Y.flatten(order='F').astype(real_dtype))
                f.write_record(b.Z.flatten(order='F').astype(real_dtype))

        read_back = read_plot3D(path, binary=True, big_endian=big_endian,
                                read_double=True, fortran=True)

    _assert_blocks_equal(blocks, read_back, atol=0.0)


@pytest.mark.parametrize("big_endian", [False, True])
@pytest.mark.parametrize("double_precision", [False, True])
def test_fortran_autodetect(big_endian, double_precision):
    """read_plot3D with no format args should auto-detect that the file
    is Fortran unformatted, plus the correct endian and precision."""
    blocks = _make_blocks()
    atol = 0.0 if double_precision else 1e-6

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "mesh.xyz")
        write_plot3D(path, blocks, binary=True, big_endian=big_endian,
                     double_precision=double_precision, fortran=True)
        read_back = read_plot3D(path)

    _assert_blocks_equal(blocks, read_back, atol)


def test_fortran_reader_raises_on_unexpected_record_size():
    """A Fortran record whose size doesn't match either layout should
    raise a clear ValueError rather than silently producing wrong
    geometry."""
    big_endian = False
    endian = '<'
    header_dtype = np.dtype(f'{endian}u4')
    int_dtype = np.dtype(f'{endian}i4')
    real_dtype = np.dtype(f'{endian}f8')

    shape = (3, 4, 2)
    npts = shape[0] * shape[1] * shape[2]

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "corrupt.xyz")
        with FortranFile(path, 'w', header_dtype) as f:
            f.write_record(np.array([1], dtype=int_dtype))
            f.write_record(np.array(shape, dtype=int_dtype))
            # Write a record that is neither npts nor 3*npts in size.
            f.write_record(np.zeros(2 * npts, dtype=real_dtype))

        with pytest.raises(ValueError, match="Unexpected Fortran record size"):
            read_plot3D(path, binary=True, big_endian=big_endian,
                        read_double=True, fortran=True)
