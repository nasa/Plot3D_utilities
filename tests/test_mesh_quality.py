"""Structured-grid mesh-quality battery tests.

Mirrors plot3d-rs ``src/mesh_quality.rs``'s ``#[cfg(test)] mod tests``
(commits ``9ebe7d3`` — DEGENERATE-vs-INVERTED cell split; ``ed01fd0`` —
element-type census; ``a6c9d9c`` — ``cell_aspect_ratio`` returns ``inf``
for a zero-length edge, not a huge finite artefact). Fixtures
(``_unit_cube``, ``_collapsed_line_block``) are fresh Python constructions
mirroring the Rust test module's ``unit_cube``/``collapsed_line_block``, not
reused from elsewhere in this codebase (no existing Python equivalent).
"""

import math

import numpy as np
import pytest

from plot3d.block import Block
from plot3d.mesh_quality import (
    BlockElementSummary,
    CellLocation,
    ElementInventory,
    ElementType,
    Handedness,
    MeshQualityReport,
    Severity,
    Thresholds,
    Violation,
    _rust_order_argmax,
    _rust_order_argmin,
    block_handedness,
    cell_aspect_ratio,
    cell_centroid,
    cell_dims,
    cell_distinct_node_count,
    cell_has_collapsed_edge,
    cell_signed_volume,
    cell_skewness,
    cell_volume_divergence,
    element_type_inventory,
    make_block_right_handed,
    run_all,
)


def _unit_cube(n: int) -> Block:
    """A right-handed n x n x n unit-cube block (n-1 cells/axis)."""
    h = 1.0 / (n - 1)
    idx = np.arange(n) * h
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    return Block(X.copy(), Y.copy(), Z.copy())


def _collapsed_line_block(n: int) -> Block:
    """A unit cube with node (i=1,j=0,k) snapped bit-identical onto
    (i=0,j=0,k) for every k — a deliberately collapsed grid line, mirroring
    an O-grid folded onto a blade-tip camber line (the Rust fixture of the
    same name).
    """
    b = _unit_cube(n)
    X, Y, Z = b.X.copy(), b.Y.copy(), b.Z.copy()
    X[1, 0, :] = X[0, 0, :]
    Y[1, 0, :] = Y[0, 0, :]
    Z[1, 0, :] = Z[0, 0, :]
    return Block(X, Y, Z)


# =============================================================================
# 1. Unit cube is clean
# =============================================================================


class TestUnitCubeIsClean:
    def test_right_handed(self):
        b = _unit_cube(5)
        assert block_handedness(b) == Handedness.RIGHT_HANDED

    def test_aspect_ratio_and_skewness_interior(self):
        b = _unit_cube(5)
        ar = cell_aspect_ratio(b)
        skew = cell_skewness(b)
        assert np.allclose(ar, 1.0, atol=1e-4)
        assert np.allclose(skew, 0.0, atol=1e-3)

    def test_run_all_passes_with_zero_violations(self):
        b = _unit_cube(5)
        report = run_all([b])
        assert report.handedness == [Handedness.RIGHT_HANDED]
        assert report.passes()
        assert report.n_error() == 0
        assert report.n_warn() == 0


# =============================================================================
# 2. The a6c9d9c regression: no aspect ratio for a collapsed cell
# =============================================================================


class TestCollapsedCellHasNoAspectRatio:
    def test_pinched_cell_aspect_ratio_is_inf(self):
        b = _collapsed_line_block(5)
        ar = cell_aspect_ratio(b)
        assert math.isinf(ar[0, 0, 0])
        assert cell_has_collapsed_edge(b)[0, 0, 0]

    def test_healthy_cell_elsewhere_is_finite(self):
        b = _collapsed_line_block(5)
        ar = cell_aspect_ratio(b)
        healthy = ar[2, 2, 2]
        assert np.isfinite(healthy)
        assert healthy >= 1.0

    def test_run_all_never_reports_a_non_finite_aspect_ratio(self):
        b = _collapsed_line_block(5)
        report = run_all([b])
        for v in report.violations:
            if v.check == "aspect_ratio":
                assert math.isfinite(v.actual), (
                    f"an aspect_ratio violation was raised with a non-finite "
                    f"value ({v.actual!r}) — the collapsed cell set the maximum"
                )


# =============================================================================
# 3. The 9ebe7d3 regression: collapsed line is degenerate, not inverted
# =============================================================================


class TestCollapsedLineIsDegenerateNotInverted:
    def test_no_negative_volume_violation(self):
        b = _collapsed_line_block(5)
        assert cell_signed_volume(b)[0, 0, 0] == 0.0
        assert cell_volume_divergence(b)[0, 0, 0] > 0.0
        assert cell_has_collapsed_edge(b)[0, 0, 0]

        report = run_all([b])
        assert not any(v.check == "negative_volume" for v in report.violations)

    def test_exactly_one_degenerate_cell_warning(self):
        b = _collapsed_line_block(5)
        report = run_all([b])
        degen = [v for v in report.violations if v.check == "degenerate_cell"]
        assert len(degen) == 1
        assert degen[0].severity == Severity.WARN

    def test_n_error_is_zero(self):
        b = _collapsed_line_block(5)
        report = run_all([b])
        assert report.n_error() == 0


# =============================================================================
# 4. A genuinely inverted cell is still fatal
# =============================================================================


class TestGenuinelyInvertedCellIsFatal:
    def test_negative_volume_error(self):
        b = _unit_cube(5)
        X, Y, Z = b.X.copy(), b.Y.copy(), b.Z.copy()
        X[0, 0, 0] = 10.0  # drag (0,0,0) far past the opposite face
        b = Block(X, Y, Z)

        assert cell_signed_volume(b)[0, 0, 0] < 0.0
        assert not cell_has_collapsed_edge(b)[0, 0, 0]

        report = run_all([b])
        assert any(
            v.check == "negative_volume" and v.severity == Severity.ERROR
            for v in report.violations
        )


# =============================================================================
# 5. Left-handed block detection + make_block_right_handed
# =============================================================================


class TestLeftHandedBlockFix:
    def test_flip_detected_and_fixed(self):
        cube = _unit_cube(5)
        lh = Block(cube.X[::-1, :, :].copy(), cube.Y[::-1, :, :].copy(), cube.Z[::-1, :, :].copy())
        assert block_handedness(lh) == Handedness.LEFT_HANDED

        fixed, axis = make_block_right_handed(lh)
        assert axis == 0
        assert block_handedness(fixed) == Handedness.RIGHT_HANDED

        report = run_all([fixed])
        assert report.passes()

    def test_noop_on_already_right_handed(self):
        cube = _unit_cube(5)
        same, axis = make_block_right_handed(cube)
        assert axis is None
        assert block_handedness(same) == Handedness.RIGHT_HANDED

    def test_original_block_not_mutated(self):
        cube = _unit_cube(5)
        lh_X = cube.X[::-1, :, :].copy()
        lh_Y = cube.Y[::-1, :, :].copy()
        lh_Z = cube.Z[::-1, :, :].copy()
        lh = Block(lh_X.copy(), lh_Y.copy(), lh_Z.copy())

        make_block_right_handed(lh)

        assert np.array_equal(lh.X, lh_X)
        assert np.array_equal(lh.Y, lh_Y)
        assert np.array_equal(lh.Z, lh_Z)
        assert block_handedness(lh) == Handedness.LEFT_HANDED


# =============================================================================
# 6. Thresholds.from_preset_name + exact numeric pinning
# =============================================================================


class TestThresholdsPreset:
    @pytest.mark.parametrize("name", ["strict", "STRICT", "Strict"])
    def test_case_insensitive_strict(self, name):
        assert Thresholds.from_preset_name(name) == Thresholds.STRICT

    @pytest.mark.parametrize("name", ["relaxed", "RELAXED", "Relaxed"])
    def test_case_insensitive_relaxed(self, name):
        assert Thresholds.from_preset_name(name) == Thresholds.RELAXED

    @pytest.mark.parametrize("name", ["standard", "STANDARD", "bogus", "", "nonsense"])
    def test_unknown_falls_back_to_standard(self, name):
        assert Thresholds.from_preset_name(name) == Thresholds.STANDARD

    def test_exact_current_values(self):
        # Regression guard against silent future threshold drift — pinned
        # from plot3d-rs src/mesh_quality.rs (post b6036d2/94c4def tightening).
        assert Thresholds.STRICT == Thresholds(
            skew_p99_deg=50.0,
            skew_max_deg=60.0,
            min_orthogonality_deg=30.0,
            max_ar_interior=100_000.0,
            max_ar_wall=200_000.0,
            min_cell_volume_ratio=1e-5,
            boundary_drop=2,
        )
        assert Thresholds.STANDARD == Thresholds(
            skew_p99_deg=80.0,
            skew_max_deg=75.0,
            min_orthogonality_deg=15.0,
            max_ar_interior=5_000.0,
            max_ar_wall=100_000.0,
            min_cell_volume_ratio=1e-6,
            boundary_drop=2,
        )
        relaxed = Thresholds.RELAXED
        assert relaxed.skew_p99_deg == 87.0
        assert relaxed.skew_max_deg == 90.0
        assert relaxed.min_orthogonality_deg == 1.0
        assert relaxed.max_ar_interior == 30_000.0
        assert math.isinf(relaxed.max_ar_wall)
        assert relaxed.min_cell_volume_ratio == 1e-9
        assert relaxed.boundary_drop == 2


# =============================================================================
# 7. element_type_inventory classifies pinched cells as PRISM
# =============================================================================


class TestElementTypeInventory:
    def test_collapsed_line_is_prism_not_hex(self):
        b = _collapsed_line_block(5)
        assert cell_distinct_node_count(b)[0, 0, 0] == 6
        assert ElementType.from_distinct_nodes(6) == ElementType.PRISM

        inv = element_type_inventory([b])
        assert inv.n_prism > 0
        # Every cell at j=0 along the collapsed line is pinched into a wedge.
        NI, NJ, NK = cell_dims(b)
        assert inv.per_block[0].n_prism == NK  # one pinched cell per k-layer
        assert inv.n_hex == NI * NJ * NK - NK

    def test_format_ads_style_smoke(self):
        b = _unit_cube(4)
        inv = element_type_inventory([b])
        text = inv.format_ads_style()
        assert "HEX" in text
        assert "ADVOLAREA" in text


# =============================================================================
# 8. np.percentile matches Rust's hand-rolled percentile
# =============================================================================


class TestPercentileEquivalence:
    def test_pinned_values_from_rust_test(self):
        # Rust's `percentile_linear_interpolation` test pins:
        #   percentile([0,1,2,3,4], 0.0) == 0.0
        #   percentile([0,1,2,3,4], 1.0) == 4.0
        #   percentile([0,1,2,3,4], 0.5) == 2.0
        # np.percentile uses a 0-100 scale for the same linear-interpolation
        # method, justifying delegation instead of reimplementing.
        v = [0.0, 1.0, 2.0, 3.0, 4.0]
        assert abs(np.percentile(v, 0) - 0.0) < 1e-6
        assert abs(np.percentile(v, 100) - 4.0) < 1e-6
        assert abs(np.percentile(v, 50) - 2.0) < 1e-6


# =============================================================================
# 9. _rust_order_argmax/argmin tie-break
# =============================================================================


class TestRustOrderArgExtrema:
    def test_argmax_breaks_ties_i_fastest(self):
        # Shape (NI,NJ,NK) = (2,2,2), tie the maximum at two cells:
        # (0,0,0) and (1,1,1). Rust's i-fastest scan (i,j,k all increasing,
        # i fastest) visits (0,0,0) before (1,1,1), so with a strict ">"
        # comparison (0,0,0) wins.
        field = np.zeros((2, 2, 2))
        field[0, 0, 0] = 5.0
        field[1, 1, 1] = 5.0
        value, (i, j, k) = _rust_order_argmax(field)
        assert value == 5.0
        assert (i, j, k) == (0, 0, 0)

    def test_argmin_breaks_ties_i_fastest(self):
        field = np.ones((2, 2, 2)) * 9.0
        field[1, 0, 0] = -3.0
        field[0, 1, 1] = -3.0
        value, (i, j, k) = _rust_order_argmin(field)
        assert value == -3.0
        # i-fastest scan order: (0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),... —
        # (1,0,0) is visited before (0,1,1).
        assert (i, j, k) == (1, 0, 0)

    def test_mask_excludes_cells(self):
        field = np.array([[[5.0, 1.0], [1.0, 1.0]], [[1.0, 1.0], [1.0, 1.0]]])
        mask = np.ones((2, 2, 2), dtype=bool)
        mask[0, 0, 0] = False  # exclude the true max
        value, _ = _rust_order_argmax(field, mask=mask)
        assert value == 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
