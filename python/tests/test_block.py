"""Tests for Block class operations."""
import numpy as np
import pytest
from plot3d import Block
from plot3d.block import compute_gcd, reduce_blocks


class TestBlock:
    def test_create_block(self, simple_block):
        assert simple_block.IMAX == 5
        assert simple_block.JMAX == 5
        assert simple_block.KMAX == 5

    def test_size(self, simple_block):
        assert simple_block.size == 125

    def test_centroid(self, simple_block):
        assert abs(simple_block.cx - 0.5) < 1e-12
        assert abs(simple_block.cy - 0.5) < 1e-12
        assert abs(simple_block.cz - 0.5) < 1e-12

    def test_copy(self, simple_block):
        b2 = simple_block.copy()
        # Separate objects
        assert b2 is not simple_block
        assert b2.X is not simple_block.X
        # Same data
        np.testing.assert_array_equal(b2.X, simple_block.X)
        np.testing.assert_array_equal(b2.Y, simple_block.Y)
        np.testing.assert_array_equal(b2.Z, simple_block.Z)
        # Mutation doesn't affect original
        b2.X[0, 0, 0] = 999.0
        assert simple_block.X[0, 0, 0] != 999.0

    def test_scale(self, simple_block):
        simple_block.scale(2.0)
        assert abs(simple_block.X.max() - 2.0) < 1e-12
        assert abs(simple_block.Y.max() - 2.0) < 1e-12
        assert abs(simple_block.Z.max() - 2.0) < 1e-12

    def test_shift_z(self, simple_block):
        simple_block.shift(5.0, 'z')
        assert abs(simple_block.Z.min() - 5.0) < 1e-12
        assert abs(simple_block.Z.max() - 6.0) < 1e-12

    def test_shift_x(self, simple_block):
        simple_block.shift(-1.0, 'x')
        assert abs(simple_block.X.min() - (-1.0)) < 1e-12

    def test_repr(self, simple_block):
        assert repr(simple_block) == "(5,5,5)"

    def test_get_faces(self, simple_block):
        faces = simple_block.get_faces()
        assert set(faces.keys()) == {'imin', 'imax', 'jmin', 'jmax', 'kmin', 'kmax'}
        # imin face should have X=0 everywhere
        xf, yf, zf = faces['imin']
        np.testing.assert_allclose(xf, 0.0)
        # imax face should have X=1 everywhere
        xf, yf, zf = faces['imax']
        np.testing.assert_allclose(xf, 1.0)


class TestComputeGCD:
    def test_uniform_blocks(self):
        """All blocks with dimensions divisible by 4."""
        x = np.linspace(0, 1, 5)
        X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
        blocks = [Block(X, Y, Z) for _ in range(3)]
        assert compute_gcd(blocks) == 4

    def test_mixed_dimensions(self):
        """Blocks with different per-block GCDs -> min per-block GCD is taken."""
        x4 = np.linspace(0, 1, 5)  # dim-1=4, GCD(4,4,4)=4
        x6 = np.linspace(0, 1, 7)  # dim-1=6, GCD(6,6,6)=6
        X4, Y4, Z4 = np.meshgrid(x4, x4, x4, indexing='ij')
        X6, Y6, Z6 = np.meshgrid(x6, x6, x6, indexing='ij')
        blocks = [Block(X4, Y4, Z4), Block(X6, Y6, Z6)]
        # compute_gcd returns min of per-block GCDs: min(4, 6) = 4
        assert compute_gcd(blocks) == 4


class TestReduceBlocks:
    def test_reduce_by_factor(self):
        x = np.linspace(0, 1, 9)
        X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
        blocks = [Block(X, Y, Z)]
        reduced = reduce_blocks([b.copy() for b in blocks], 2)
        assert reduced[0].IMAX == 5
        assert reduced[0].JMAX == 5
        assert reduced[0].KMAX == 5
        # Corner values preserved
        np.testing.assert_allclose(reduced[0].X[0, 0, 0], 0.0)
        np.testing.assert_allclose(reduced[0].X[-1, -1, -1], 1.0)
