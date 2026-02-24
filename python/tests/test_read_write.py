"""Tests for read/write round-trip correctness."""
import numpy as np
import os
import tempfile
import pytest
from plot3d import Block, read_plot3D, write_plot3D


def _make_blocks():
    """Create a list of 2 blocks with known data for round-trip testing."""
    x1 = np.linspace(0, 1, 5)
    y1 = np.linspace(0, 1, 3)
    z1 = np.linspace(0, 1, 7)
    X1, Y1, Z1 = np.meshgrid(x1, y1, z1, indexing='ij')

    x2 = np.linspace(1, 2, 9)
    y2 = np.linspace(0, 1, 5)
    z2 = np.linspace(0, 1, 5)
    X2, Y2, Z2 = np.meshgrid(x2, y2, z2, indexing='ij')

    return [Block(X1, Y1, Z1), Block(X2, Y2, Z2)]


def _assert_blocks_equal(original, loaded, tol=1e-12):
    """Assert two block lists have matching data."""
    assert len(original) == len(loaded)
    for bo, bl in zip(original, loaded):
        assert bo.IMAX == bl.IMAX
        assert bo.JMAX == bl.JMAX
        assert bo.KMAX == bl.KMAX
        np.testing.assert_allclose(bl.X, bo.X, atol=tol)
        np.testing.assert_allclose(bl.Y, bo.Y, atol=tol)
        np.testing.assert_allclose(bl.Z, bo.Z, atol=tol)


class TestBinaryRoundTrip:
    def test_binary_double_little_endian(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, binary=True, big_endian=False, double_precision=True)
        loaded = read_plot3D(path, binary=True, big_endian=False, read_double=True)
        _assert_blocks_equal(blocks, loaded)

    def test_binary_double_big_endian(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, binary=True, big_endian=True, double_precision=True)
        loaded = read_plot3D(path, binary=True, big_endian=True, read_double=True)
        _assert_blocks_equal(blocks, loaded)

    def test_binary_single_little_endian(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, binary=True, big_endian=False, double_precision=False)
        loaded = read_plot3D(path, binary=True, big_endian=False, read_double=False)
        _assert_blocks_equal(blocks, loaded, tol=1e-6)

    def test_binary_single_big_endian(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, binary=True, big_endian=True, double_precision=False)
        loaded = read_plot3D(path, binary=True, big_endian=True, read_double=False)
        _assert_blocks_equal(blocks, loaded, tol=1e-6)


class TestASCIIRoundTrip:
    def test_ascii_double(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.xyz")
        write_plot3D(path, blocks, binary=False, double_precision=True)
        loaded = read_plot3D(path, binary=False)
        _assert_blocks_equal(blocks, loaded, tol=1e-12)

    def test_ascii_single(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.xyz")
        write_plot3D(path, blocks, binary=False, double_precision=False)
        loaded = read_plot3D(path, binary=False)
        _assert_blocks_equal(blocks, loaded, tol=1e-6)


class TestFortranRoundTrip:
    def test_fortran_double(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, fortran=True, double_precision=True)
        loaded = read_plot3D(path, fortran=True, read_double=True)
        _assert_blocks_equal(blocks, loaded)

    def test_fortran_single(self, tmp_path):
        blocks = _make_blocks()
        path = str(tmp_path / "test.p3d")
        write_plot3D(path, blocks, fortran=True, double_precision=False)
        loaded = read_plot3D(path, fortran=True, read_double=False)
        _assert_blocks_equal(blocks, loaded, tol=1e-6)


class TestSingleBlock:
    def test_single_block_binary(self, tmp_path, simple_block):
        path = str(tmp_path / "single.p3d")
        write_plot3D(path, [simple_block], binary=True)
        loaded = read_plot3D(path, binary=True)
        assert len(loaded) == 1
        np.testing.assert_allclose(loaded[0].X, simple_block.X, atol=1e-12)

    def test_nonexistent_file(self):
        result = read_plot3D("/tmp/nonexistent_file_xyz.p3d")
        assert result == []
