"""Distance-based point-coincidence primitive.

Mirrors plot3d-rs's geometry::PointGrid3/PointGrid2 fix (commit 0e1b1a1):
older coincidence checks bucketed points via `round(x/tol)` and compared
bucket-integer equality, which can miss two points that are genuinely
closer than `tol` but happen to straddle an integer bin boundary. This
file pins `plot3d.geometry.coincidence_count`'s distance-based behavior,
including the specific bin-boundary regression the bucket scheme got
wrong.
"""

import numpy as np

from plot3d.geometry import coincidence_count


class TestBasicSanity:
    def test_identical_point_sets_match_all(self):
        P = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-4.0, 5.0, 6.0]])
        assert coincidence_count(P, P.copy(), tol=1e-6) == 3

    def test_disjoint_far_apart_sets_match_none(self):
        P1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        P2 = np.array([[1000.0, 1000.0, 1000.0], [2000.0, 2000.0, 2000.0]])
        assert coincidence_count(P1, P2, tol=1e-6) == 0


class TestBinBoundaryRegression:
    """The actual bug this module exists to fix.

    With tol = 1e-4: point A sits at x = 4.4*tol, point B at x = 4.6*tol.
    Their separation is 0.2*tol (well within tol), but under the OLD
    `round(x/tol)` bucketing scheme:
        round(4.4) == 4
        round(4.6) == 5
    -- different integer bins, so the old bucket-equality test would
    report ZERO coincident points even though the points are genuinely
    within tol of each other. coincidence_count uses real Euclidean
    distance via cKDTree and must find the match.
    """

    def test_points_straddling_a_bin_boundary_are_found(self):
        tol = 1e-4
        a_x = 4.4 * tol
        b_x = 4.6 * tol
        assert round(a_x / tol) == 4
        assert round(b_x / tol) == 5  # confirms the old scheme would split these
        assert abs(b_x - a_x) < tol  # confirms they are genuinely coincident

        P1 = np.array([[a_x, 0.0, 0.0]])
        P2 = np.array([[b_x, 0.0, 0.0]])

        # Old round(x/tol)-bucket-equality scheme would return 0 here since
        # the two points quantize to different bins (4 vs 5).
        assert coincidence_count(P1, P2, tol=tol) == 1

    def test_points_just_outside_tol_are_not_matched(self):
        tol = 1e-4
        P1 = np.array([[0.0, 0.0, 0.0]])
        P2 = np.array([[1.5 * tol, 0.0, 0.0]])
        assert coincidence_count(P1, P2, tol=tol) == 0


class TestEmptyInputs:
    def test_empty_p1_returns_zero(self):
        P1 = np.zeros((0, 3))
        P2 = np.array([[0.0, 0.0, 0.0]])
        assert coincidence_count(P1, P2, tol=1e-6) == 0

    def test_empty_p2_returns_zero(self):
        P1 = np.array([[0.0, 0.0, 0.0]])
        P2 = np.zeros((0, 3))
        assert coincidence_count(P1, P2, tol=1e-6) == 0

    def test_both_empty_returns_zero(self):
        P1 = np.zeros((0, 3))
        P2 = np.zeros((0, 3))
        assert coincidence_count(P1, P2, tol=1e-6) == 0


class TestDimensionality:
    """face.py/periodicity.py call sites pass both (N,3) points (structured
    face grids) and (N,2) points (orthogonal-plane projections in
    periodicity.py's _orthogonal_precheck), so both must work."""

    def test_3d_points(self):
        P1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        P2 = np.array([[0.0, 0.0, 0.0], [5.0, 5.0, 5.0]])
        assert coincidence_count(P1, P2, tol=1e-6) == 1

    def test_2d_points(self):
        P1 = np.array([[0.0, 0.0], [1.0, 1.0]])
        P2 = np.array([[0.0, 0.0], [5.0, 5.0]])
        assert coincidence_count(P1, P2, tol=1e-6) == 1

    def test_partial_overlap_count(self):
        P1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        P2 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [99.0, 99.0, 99.0]])
        assert coincidence_count(P1, P2, tol=1e-6) == 2
