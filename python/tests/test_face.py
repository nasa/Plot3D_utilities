"""Tests for Face class operations."""
import numpy as np
import pytest
from plot3d import Block
from plot3d.face import Face
from plot3d.facefunctions import get_outer_faces, create_face_from_diagonals, split_face


class TestFace:
    def test_create_face(self):
        f = Face()
        assert f.nvertex == 0
        assert f.blockIndex == 0

    def test_copy(self):
        f = Face()
        f.x = np.array([0.0, 1.0, 1.0, 0.0])
        f.y = np.array([0.0, 0.0, 1.0, 1.0])
        f.z = np.zeros(4)
        f.I = np.array([0, 4, 4, 0])
        f.J = np.array([0, 0, 4, 4])
        f.K = np.zeros(4, dtype=np.int64)
        f.nvertex = 4
        f.blockIndex = 3
        f.cx = 0.5
        f.cy = 0.5

        g = f.copy()
        assert g is not f
        assert g.blockIndex == 3
        assert g.nvertex == 4
        np.testing.assert_array_equal(g.x, f.x)
        np.testing.assert_array_equal(g.I, f.I)
        # Mutation doesn't affect original
        g.x[0] = 999.0
        assert f.x[0] == 0.0

    def test_to_dict(self):
        f = Face()
        f.I = np.array([0, 4, 4, 0])
        f.J = np.array([0, 0, 4, 4])
        f.K = np.array([0, 0, 0, 0])
        f.blockIndex = 2
        f.id = 7
        d = f.to_dict()
        assert d['IMIN'] == 0
        assert d['IMAX'] == 4
        assert d['JMIN'] == 0
        assert d['JMAX'] == 4
        assert d['KMIN'] == 0
        assert d['KMAX'] == 0
        assert d['block_index'] == 2
        assert d['id'] == 7


class TestGetOuterFaces:
    def test_unit_cube_has_6_faces(self, simple_block):
        faces, _ = get_outer_faces(simple_block)
        assert len(faces) == 6

    def test_face_coordinates(self, simple_block):
        """Outer faces should have correct coordinates at block boundaries."""
        faces, _ = get_outer_faces(simple_block)
        # Find the imin face (all I indices should be 0)
        imin_faces = [f for f in faces if f.I.max() == 0]
        assert len(imin_faces) == 1
        f = imin_faces[0]
        # All x-coords should be 0
        assert abs(f.x.max()) < 1e-12


class TestCreateFaceFromDiagonals:
    def test_create_imin_face(self, simple_block):
        f = create_face_from_diagonals(simple_block, 0, 0, 0, 0, 4, 4)
        assert f.nvertex == 4
        # This is an i-constant face at i=0
        assert int(f.I.min()) == 0
        assert int(f.I.max()) == 0
        assert int(f.J.min()) == 0
        assert int(f.J.max()) == 4
        assert int(f.K.min()) == 0
        assert int(f.K.max()) == 4


class TestSplitFace:
    def test_split_creates_valid_remnants(self, simple_block):
        """Splitting a face should produce remnants that cover remaining area."""
        full_face = create_face_from_diagonals(simple_block, 0, 0, 0, 0, 4, 4)
        # split_face(face, block, imin, jmin, kmin, imax, jmax, kmax)
        remnants = split_face(full_face, simple_block, 0, 0, 0, 0, 2, 4)
        # Should produce at least 1 remnant covering the j=2..4 region
        assert len(remnants) >= 1
        for r in remnants:
            assert r.nvertex == 4
