"""Tests for FacePool theta-bucketed spatial index."""
import numpy as np
import math
import pytest
from plot3d import Block, FacePool
from plot3d.face import Face
from plot3d.facefunctions import get_outer_faces


def _make_cylindrical_block(x_range, r_range, theta_range, nx=5, nr=5, nt=9):
    """Create a block in cylindrical coordinates around x-axis."""
    x = np.linspace(*x_range, nx)
    r = np.linspace(*r_range, nr)
    theta = np.linspace(*theta_range, nt)

    X = np.zeros((nx, nr, nt))
    Y = np.zeros((nx, nr, nt))
    Z = np.zeros((nx, nr, nt))
    for i in range(nx):
        for j in range(nr):
            for k in range(nt):
                X[i, j, k] = x[i]
                Y[i, j, k] = r[j] * np.cos(theta[k])
                Z[i, j, k] = r[j] * np.sin(theta[k])
    return Block(X, Y, Z)


class TestFacePoolConstruction:
    def test_build_face_pool(self, cylindrical_block_pair):
        blocks = cylindrical_block_pair
        faces = []
        for bi, b in enumerate(blocks):
            outs, _ = get_outer_faces(b)
            for f in outs:
                f.blockIndex = bi
            faces.extend(outs)

        pool = FacePool(faces, blocks, rotation_axis='x')
        assert len(pool.faces) == len(faces)
        assert len(pool.cyl_info) == len(faces)
        assert len(pool._theta_sorted) == len(faces)

    def test_consume_and_active(self, cylindrical_block_pair):
        blocks = cylindrical_block_pair
        faces = []
        for bi, b in enumerate(blocks):
            outs, _ = get_outer_faces(b)
            for f in outs:
                f.blockIndex = bi
            faces.extend(outs)

        pool = FacePool(faces, blocks, rotation_axis='x')
        n = len(faces)
        assert len(pool.active_indices()) == n

        pool.consume(0)
        assert pool.is_consumed(0)
        assert len(pool.active_indices()) == n - 1


class TestFacePoolCandidates:
    def test_find_candidates_at_same_theta(self):
        """Two blocks with overlapping axial/radial ranges and similar theta should find each other."""
        b1 = _make_cylindrical_block((0, 1), (0.5, 1.0), (0, np.pi / 4))
        b2 = _make_cylindrical_block((0, 1), (0.5, 1.0), (np.pi / 4, np.pi / 2))
        blocks = [b1, b2]

        faces = []
        for bi, b in enumerate(blocks):
            outs, _ = get_outer_faces(b)
            for f in outs:
                f.blockIndex = bi
            faces.extend(outs)

        pool = FacePool(faces, blocks, rotation_axis='x')

        # Look for faces from block 1 near the shared boundary (theta = pi/4)
        # using a broad enough theta tolerance
        found_block1 = False
        for i, f in enumerate(pool.faces):
            if f.blockIndex == 0:
                info = pool.cyl_info[i]
                # Search at the shared boundary theta
                cands = pool.find_candidates(
                    np.pi / 4,
                    (info.axial_min, info.axial_max),
                    (info.radial_min, info.radial_max),
                    theta_tol=0.5,
                )
                cand_blocks = {pool.faces[c].blockIndex for c in cands}
                if 1 in cand_blocks:
                    found_block1 = True
                    break
        assert found_block1, "Should find block 1 faces near theta=pi/4"

    def test_add_face(self, cylindrical_block_pair):
        blocks = cylindrical_block_pair
        outs0, _ = get_outer_faces(blocks[0])
        for f in outs0:
            f.blockIndex = 0

        pool = FacePool(outs0, blocks, rotation_axis='x')
        n_before = len(pool.faces)

        outs1, _ = get_outer_faces(blocks[1])
        outs1[0].blockIndex = 1
        pool.add_face(outs1[0])
        assert len(pool.faces) == n_before + 1
