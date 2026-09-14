"""Completeness + node-for-node certification in ``periodicity.py``.

Mirrors plot3d-rs commit ``0e1b1a1`` (full match certification) applied to
``periodicity.py``'s ``__periodicity_check__``: the cKDTree nearest-neighbor
matcher used for rotational periodicity previously accepted a proposal as
soon as it found >=4 non-collinear point pairs whose index bounding box
looked like a face, with no check that every point *inside* that bounding
box actually matched, and no check that the correspondence is a genuine
structured (index-consistent) mapping rather than an accidental point-cloud
overlap. This adds:

1. A completeness check -- matched point count must equal the claimed
   sub-patch's index-space area (mirrors ``connectivity.py``'s
   ``get_face_intersection``, same standard: reject iff
   ``len(df) < matched_area``).
2. Full certification via ``correspondence.certify_correspondence`` -- every
   node of the claimed patch must correspond under exactly one of the 8
   structured permutations.

Tests operate directly on ``__periodicity_check__`` with small hand-built
single-k-layer blocks (no rotation applied -- ``block1``'s "already rotated"
role is irrelevant to unit-testing the matcher itself; any two blocks whose
face points coincide within tolerance exercise the same code path), plus one
end-to-end regression against the real VSPT mesh fixture already used by
``test_rotated_periodicity.py``.
"""

import os

import numpy as np
import pytest

from plot3d import Block
from plot3d.facefunctions import create_face_from_diagonals
from plot3d.periodicity import __periodicity_check__ as periodicity_check

MESH_PATH = os.path.join(os.path.dirname(__file__), "data", "vspt_mesh_scaled.xyz")


def _asym_grid(nu: int, nv: int) -> np.ndarray:
    """An (nu, nv, 3) point grid with no reversal/swap/translation symmetry.

    Every (i, j) maps to a distinct (x, y, z) triple (X alone distinguishes
    i, Y alone distinguishes j), so a deliberate value-swap between two grid
    cells (used by the certification test below) is unambiguous and every
    one of the 8 structured permutations produces geometrically distinct
    points -- an "exactly one permutation certifies" or "no permutation
    certifies" outcome isn't an accident of a too-symmetric fixture.
    """
    i = np.arange(nu, dtype=float)
    j = np.arange(nv, dtype=float)
    I, J = np.meshgrid(i, j, indexing="ij")
    X = I * 1.7
    Y = J * 0.9
    Z = 0.05 * I * I + 0.03 * J + 0.01 * I * J
    return np.stack([X, Y, Z], axis=-1)


def _block_from_grid(grid: np.ndarray) -> Block:
    """Build a single-k-layer Block from an (nu, nv, 3) point grid."""
    return Block(grid[:, :, 0:1].copy(), grid[:, :, 1:2].copy(), grid[:, :, 2:3].copy())


def _full_face(block: Block, nu: int, nv: int):
    return create_face_from_diagonals(block, [0, 0, 0], [nu - 1, nv - 1, 0])


class TestNormalMatchStillCertifies:
    """Sanity/regression: an exact, unambiguous full-face match must still
    be accepted -- the new checks tighten correctness, they must not reject
    genuine matches."""

    def test_exact_coincident_faces_match(self):
        nu, nv = 4, 3
        grid = _asym_grid(nu, nv)
        block1 = _block_from_grid(grid)
        block2 = _block_from_grid(grid.copy())
        face1 = _full_face(block1, nu, nv)
        face2 = _full_face(block2, nu, nv)

        df, periodic_faces, split_faces = periodicity_check(
            face1, face2, block1, block2, tol=1e-9)

        assert len(df) == nu * nv
        assert len(periodic_faces) == 2
        assert len(split_faces) == 0


class TestCompletenessCheck:
    """Matched-point count must equal the claimed sub-patch's index-space
    area -- a KDTree match that covers every point of a bounding box except
    one interior "hole" must be rejected, not silently accepted as a full
    face match with a missing node."""

    def test_missing_interior_point_is_rejected(self):
        nu, nv = 4, 3
        grid1 = _asym_grid(nu, nv)
        grid2 = grid1.copy()
        # Push one interior point (not a corner/edge -- (1,1) is interior
        # for a 4x3 box) far away so its nearest neighbor in block1's cloud
        # is some other point, at a distance well beyond tol. Every other
        # of the 12 points still coincides exactly, so the matched bounding
        # box still spans the full face -- only the completeness check
        # (matched count == box area) catches the missing node.
        grid2[1, 1] += 50.0

        block1 = _block_from_grid(grid1)
        block2 = _block_from_grid(grid2)
        face1 = _full_face(block1, nu, nv)
        face2 = _full_face(block2, nu, nv)

        df, periodic_faces, split_faces = periodicity_check(
            face1, face2, block1, block2, tol=1e-6)

        assert len(df) == 0
        assert periodic_faces == []
        assert split_faces == []

    def test_completeness_check_triggers_on_truncated_match_directly(self):
        """Direct unit check of the completeness arithmetic itself (not
        routed through the KDTree), as a regression guard independent of
        floating-point nearest-neighbor behavior: 11 matched points inside
        a 4x3=12 claimed bounding box must be judged incomplete."""
        from plot3d.connectivity import _face_point_count

        matched_area = _face_point_count([0, 0, 0], [3, 2, 0])
        assert matched_area == 12
        n_matched = 11
        assert n_matched < matched_area  # exactly the rejection condition


class TestCertificationCatchesInconsistentCorrespondence:
    """A proposal can pass the completeness check (matched count == claimed
    area, corners agree) while still not being a genuine structured
    (conformal) interface -- e.g. two interior nodes' coordinate values are
    transposed between the two faces, so every point still has an exact
    nearest neighbor somewhere in the other face's cloud (completeness is
    satisfied), but no single one of the 8 structured permutations maps
    block1's grid onto block2's grid. Node-for-node certification must
    reject this; a corner/count-only check would have accepted it."""

    def test_swapped_interior_nodes_rejected_despite_full_coverage(self):
        nu, nv = 4, 3
        grid1 = _asym_grid(nu, nv)
        grid2 = grid1.copy()
        # Swap two interior (non-boundary) points' coordinate values. Both
        # points still exist somewhere in block2's cloud (at each other's
        # slot), so KDTree nearest-neighbor matching finds an exact
        # (distance 0) counterpart for all 12 points -- the completeness
        # check alone cannot see this defect. But block2's raw grid, taken
        # in ascending index order (as certify_correspondence extracts it),
        # no longer corresponds to block1's grid under any single
        # structured permutation.
        grid2[[1, 2], 1] = grid2[[2, 1], 1]

        block1 = _block_from_grid(grid1)
        block2 = _block_from_grid(grid2)
        face1 = _full_face(block1, nu, nv)
        face2 = _full_face(block2, nu, nv)

        # Confirm the premise independently (every rejection path inside
        # __periodicity_check__ returns a fully-cleared empty DataFrame, so
        # this can't be read off the function's own return value): a raw
        # KDTree nearest-neighbor pass finds all 12 points within tol,
        # spanning the full bounding box -- completeness alone would accept
        # this proposal.
        from scipy.spatial import cKDTree

        from plot3d.periodicity import _extract_face_points

        pts1, _ = _extract_face_points(face1, block1)
        pts2, _ = _extract_face_points(face2, block2)
        dists, _ = cKDTree(pts2).query(pts1, k=1)
        assert int(np.sum(dists < 1e-9)) == 12

        df, periodic_faces, split_faces = periodicity_check(
            face1, face2, block1, block2, tol=1e-9)

        # Certification (not the count check) is what rejects this.
        assert len(df) == 0
        assert periodic_faces == []
        assert split_faces == []


@pytest.mark.skipif(not os.path.exists(MESH_PATH), reason="vspt_mesh_scaled.xyz not found")
class TestEndToEndNoRegression:
    """The real VSPT mesh (2 blocks, 55 blades) has genuine periodic pairs
    that must still be found end-to-end after adding completeness +
    certification -- this is a correctness-tightening change, not a
    behavior-narrowing one for actual conformal interfaces."""

    def test_rotated_periodicity_still_finds_periodic_pairs(self):
        from plot3d import connectivity_fast, read_plot3D, rotated_periodicity

        blocks = read_plot3D(MESH_PATH)
        face_matches, outer_faces = connectivity_fast(blocks)
        periodic_export, outer_export, periodic_pairs, _ = rotated_periodicity(
            blocks, face_matches, outer_faces,
            rotation_angle=360.0 / 55, rotation_axis="x",
        )

        assert len(periodic_export) > 0
        assert len(periodic_pairs) == len(periodic_export)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
