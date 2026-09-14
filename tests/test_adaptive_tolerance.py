"""Adaptive node-matching tolerance for connectivity detection.

Mirrors plot3d-rs ``connectivity::adaptive_tolerance``/``TOL_FLOOR``
(commit c7e9cf6, ``src/connectivity.rs``). A fixed ``1e-6`` tolerance only
makes sense for coordinates of order 1; coordinate storage loses precision
proportionally to magnitude, so the tolerance is scaled by the mesh's own
coordinate magnitude, floored at the historical ``1e-6`` (so nothing that
matched before can stop matching), and capped at a quarter of the mesh's
finest cell-corner-to-corner spacing (so it can never start matching
genuinely separate faces or cross-pair a node with the wrong close
neighbour on a sheared cell).
"""

import numpy as np
import pytest

from plot3d import Block
from plot3d.connectivity import (
    TOL_FLOOR,
    adaptive_tolerance,
    connectivity,
    connectivity_fast,
)


def _grid_block(n: int, spacing: float, origin: float = 0.0) -> Block:
    """An ``n``x``n``x``n`` axis-aligned grid block with uniform ``spacing``,
    starting at ``origin`` along every axis."""
    coords = origin + spacing * np.arange(n, dtype=np.float64)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing="ij")
    return Block(X.copy(), Y.copy(), Z.copy())


def _grid_block_x_offset(n: int, spacing: float, x_origin: float = 0.0) -> Block:
    """Like :func:`_grid_block`, but only the X axis is shifted by
    ``x_origin`` -- Y and Z stay at ``[0, (n - 1) * spacing]``. Used to build
    two blocks that share a face in the Y-Z plane."""
    x_coords = x_origin + spacing * np.arange(n, dtype=np.float64)
    yz_coords = spacing * np.arange(n, dtype=np.float64)
    X, Y, Z = np.meshgrid(x_coords, yz_coords, yz_coords, indexing="ij")
    return Block(X.copy(), Y.copy(), Z.copy())


class TestSmallCoordinateMesh:
    """max|coord| <= 1 -> adaptive_tolerance returns exactly TOL_FLOOR.

    This is the "bit-identical to the historical constant" invariant: at
    coordinate magnitude 1 the adaptive formula agrees with the old fixed
    tolerance bit-for-bit, and the (more expensive) spacing pass is skipped
    entirely.
    """

    def test_unit_cube_returns_exactly_tol_floor(self):
        block = _grid_block(n=5, spacing=0.25)  # coords in [0, 1]
        assert adaptive_tolerance([block]) == TOL_FLOOR

    def test_sub_unit_cube_returns_exactly_tol_floor(self):
        block = _grid_block(n=4, spacing=0.1)  # coords in [0, 0.3]
        assert adaptive_tolerance([block]) == TOL_FLOOR

    def test_empty_block_list_returns_tol_floor(self):
        assert adaptive_tolerance([]) == TOL_FLOOR


class TestLargeCoordinateMesh:
    """Large coordinates + coarse spacing -> noise-scaled tolerance, above
    TOL_FLOOR, unconstrained by the spacing ceiling."""

    def test_large_coarse_mesh_is_noise_scaled(self):
        # 5 points spaced 1e6 apart: coords in [0, 4e6], edge length 1e6.
        block = _grid_block(n=5, spacing=1.0e6)
        tol = adaptive_tolerance([block])

        scale = 4.0e6
        expected_noise = 1e-6 * scale  # = 4.0
        # Ceiling (0.25 * edge length = 2.5e5) is far above the noise
        # estimate here, so it should not engage.
        assert tol == pytest.approx(expected_noise)
        assert tol > TOL_FLOOR

    def test_moderately_large_mesh_exceeds_floor(self):
        block = _grid_block(n=3, spacing=10.0)  # coords in [0, 20]
        tol = adaptive_tolerance([block])
        assert tol > TOL_FLOOR
        assert tol == pytest.approx(1e-6 * 20.0)


class TestFineCellClamp:
    """A very fine/thin cell somewhere in the mesh clamps the tolerance
    below the pure magnitude-noise estimate, via the spacing ceiling."""

    def test_fine_cell_clamps_below_pure_noise_estimate(self):
        # Coarse, large-magnitude block: drives up the noise estimate.
        coarse = _grid_block(n=3, spacing=1.0e6)  # coords in [0, 2e6]

        # A separate block with one very fine cell (spacing 1e-2), well
        # within the coordinate magnitude of `coarse` so it does not change
        # `scale`, but which sets the mesh-wide minimum corner spacing.
        fine = _grid_block(n=2, spacing=1.0e-2)  # coords in [0, 0.01]

        pure_noise = 1e-6 * 2.0e6  # = 2.0, what adaptive_tolerance would
        # return if only `coarse` existed (no spacing ceiling to apply).
        assert adaptive_tolerance([coarse]) == pytest.approx(pure_noise)

        tol = adaptive_tolerance([coarse, fine])
        expected_ceiling = 0.25 * 1.0e-2
        assert tol == pytest.approx(expected_ceiling)
        assert tol < pure_noise
        assert tol > TOL_FLOOR

    def test_single_node_blocks_have_no_spacing_and_return_noise_or_floor(self):
        # A block with a single node has no corner pairs at all; per Rust,
        # min_cell_corner_spacing returns None and the historical value is
        # kept even though coordinate magnitude is large.
        X = np.array([[[5.0e6]]])
        Y = np.array([[[0.0]]])
        Z = np.array([[[0.0]]])
        block = Block(X, Y, Z)
        assert adaptive_tolerance([block]) == TOL_FLOOR


class TestConnectivityEndToEnd:
    """connectivity()/connectivity_fast() still work with tol=None
    (default), reproducing pre-existing behaviour on an ordinary
    order-one mesh -- no regression for existing callers."""

    @staticmethod
    def _two_unit_cubes():
        # Two 3x3x3 unit-spacing blocks sharing the x=2 face (Y, Z ranges
        # match; only X is offset).
        block1 = _grid_block_x_offset(n=3, spacing=1.0, x_origin=0.0)
        block2 = _grid_block_x_offset(n=3, spacing=1.0, x_origin=2.0)
        return [block1, block2]

    def test_connectivity_default_tol_finds_shared_face(self):
        blocks = self._two_unit_cubes()
        face_matches, outer_faces = connectivity(blocks)
        assert len(face_matches) == 1
        assert len(outer_faces) == 10  # 6 + 6 - 2 matched faces

    def test_connectivity_fast_default_tol_finds_shared_face(self):
        blocks = self._two_unit_cubes()
        face_matches, outer_faces = connectivity_fast(blocks)
        assert len(face_matches) == 1
        assert len(outer_faces) == 10

    def test_connectivity_explicit_tol_matches_default(self):
        blocks_a = self._two_unit_cubes()
        blocks_b = self._two_unit_cubes()
        default_matches, default_outer = connectivity(blocks_a)
        explicit_matches, explicit_outer = connectivity(blocks_b, tol=TOL_FLOOR)
        assert len(default_matches) == len(explicit_matches) == 1
        assert len(default_outer) == len(explicit_outer) == 10
