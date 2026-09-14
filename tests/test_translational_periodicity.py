"""Regression tests for the fabricated-1.0 spacing bug in
`translational_periodicity` (mirrors plot3d-rs's fix to its equivalent
in-plane-spacing helper: return None instead of a made-up 1.0 when a face
has too few points to measure a spacing from, and propagate "no tolerance
derivable" as "skip this pair" everywhere downstream).

Covers:
  1. `_median_inplane_spacing` returns None for a degenerate (<=1 point)
     face instead of a fabricated 1.0.
  2. `_combine_pair_spacings` (the `_pair_tol` combiner) returns None only
     when BOTH sides are None, and falls back to the available side's
     value when only one side is None.
  3. End to end: `translational_periodicity` never crashes and never
     manufactures a match out of a pair with no derivable tolerance.
"""

import numpy as np
import pytest

from plot3d.block import Block
from plot3d.facefunctions import create_face_from_diagonals
from plot3d.geometry import coincidence_count
from plot3d.periodicity import (
    _median_inplane_spacing,
    _combine_pair_spacings,
    translational_periodicity,
)


def _box_block(ni: int, nj: int, nk: int) -> Block:
    """A simple axis-aligned unit-spaced box block, (ni, nj, nk) points."""
    i, j, k = np.meshgrid(
        np.arange(ni, dtype=float),
        np.arange(nj, dtype=float),
        np.arange(nk, dtype=float),
        indexing="ij",
    )
    return Block(i, j, k)


# ---------------------------------------------------------------------
# 1) _median_inplane_spacing: degenerate faces return None, not 1.0
# ---------------------------------------------------------------------

def test_median_inplane_spacing_single_point_face_returns_none():
    block = _box_block(5, 5, 5)
    # K-constant face collapsed to a single point (I and J ranges both
    # degenerate) -- exactly 0 in-plane edges to measure.
    face = create_face_from_diagonals(block, [0, 0, 0], [0, 0, 0])
    assert _median_inplane_spacing(face, block) is None


def test_median_inplane_spacing_normal_face_returns_float():
    block = _box_block(5, 5, 5)
    # K-constant face spanning the full I,J extent -- plenty of in-plane
    # edges (unit spacing everywhere).
    face = create_face_from_diagonals(block, [0, 0, 0], [4, 4, 0])
    s = _median_inplane_spacing(face, block)
    assert s is not None
    assert s == pytest.approx(1.0)


# ---------------------------------------------------------------------
# 2) _combine_pair_spacings: the `_pair_tol` combiner logic
# ---------------------------------------------------------------------

def test_combine_pair_spacings_both_none_is_none():
    assert _combine_pair_spacings(None, None) is None


def test_combine_pair_spacings_one_none_falls_back_to_other_side():
    # Only B has data -> combiner must use B's value, not treat the pair
    # as unresolvable and not average None with a number.
    only_b = _combine_pair_spacings(None, 2.0)
    only_a = _combine_pair_spacings(2.0, None)
    assert only_b == only_a == max(0.03 * 2.0, 1e-4)


def test_combine_pair_spacings_both_present_uses_larger_side():
    assert _combine_pair_spacings(1.0, 3.0) == max(0.03 * 3.0, 1e-4)


# ---------------------------------------------------------------------
# 3) End-to-end: no crash, no fabricated-tolerance match
# ---------------------------------------------------------------------

def test_translational_periodicity_finds_real_pair():
    """Sanity baseline: a box periodic along z finds its K=0/K=max pair."""
    block = _box_block(4, 4, 5)  # z spans 0..4
    periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
        [block], outer_faces=[], translational_direction="z",
    )
    assert len(periodic_pairs) == 1


def test_translational_periodicity_no_spacing_data_skips_without_crash(monkeypatch):
    """If in-plane spacing can never be derived (simulating every face
    being degenerate), the pipeline must not crash with a TypeError from
    comparing against a None tolerance, and must simply find no periodic
    pairs rather than matching under a fabricated 1.0 tolerance."""
    # `plot3d/__init__.py` does `from .periodicity import periodicity`,
    # which shadows the `periodicity` submodule attribute on the `plot3d`
    # package with that function -- so `import plot3d.periodicity as x`
    # would resolve to the wrong object. Go through sys.modules instead.
    import sys
    periodicity_mod = sys.modules["plot3d.periodicity"]

    monkeypatch.setattr(periodicity_mod, "_median_inplane_spacing", lambda face, block: None)

    block = _box_block(4, 4, 5)
    periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
        [block], outer_faces=[], translational_direction="z",
    )
    assert periodic_pairs == []
    assert periodic_export == []


# ---------------------------------------------------------------------
# 4) Sub-patch location + node-for-node certification
#
# `translational_periodicity`'s coverage-fraction matcher (`faces_match`'s
# orthogonal precheck / `touches_by_nodes`) only ever answered "do these
# faces touch enough", never "does every node of an explicit claimed
# sub-patch actually correspond". These tests mirror
# `test_periodicity_certification.py`'s technique (applied there to
# rotational periodicity's `__periodicity_check__`) for the translational
# matcher: a clean multi-block pair must still be found (no regression),
# and a pair whose interior nodes are inconsistently ordered -- but whose
# points still individually coincide somewhere in the other face's cloud,
# fooling the coverage-fraction test -- must now be rejected.
# ---------------------------------------------------------------------

def _asym_layer(ni: int, nj: int) -> np.ndarray:
    """(ni, nj, 2) array of (x, y) with no reversal/swap symmetry.

    Mirrors `test_periodicity_certification.py`'s `_asym_grid`: every
    (i, j) maps to a geometrically distinct (x, y), so a deliberate
    interior-node swap (used below) is unambiguous -- it cannot be
    mistaken for one of the 8 structured permutations.
    """
    i = np.arange(ni, dtype=float)
    j = np.arange(nj, dtype=float)
    I, J = np.meshgrid(i, j, indexing="ij")
    X = I * 1.7 + 0.05 * J
    Y = J * 0.9 + 0.02 * I * I
    return np.stack([X, Y], axis=-1)


def _slab_block(layer: np.ndarray, z0: float, nk: int = 3) -> Block:
    """A `layer`-shaped (x, y) footprint extruded over `nk` unit-spaced
    z-layers starting at `z0` -- a simple axis-aligned slab block."""
    ni, nj = layer.shape[0], layer.shape[1]
    X = np.repeat(layer[:, :, 0:1], nk, axis=2)
    Y = np.repeat(layer[:, :, 1:2], nk, axis=2)
    zvals = z0 + np.arange(nk, dtype=float)
    Z = np.broadcast_to(zvals, (ni, nj, nk)).copy()
    return Block(X, Y, Z)


class TestNoRegressionMultiBlockPair:
    """A genuine 2-block periodic pair (distinct blocks, not one block
    self-periodic on its own two faces) must still be found after the
    sub-patch location + certification restructure."""

    def test_two_separate_blocks_periodic_pair_found(self):
        ni, nj = 4, 4
        layer = _asym_layer(ni, nj)
        d = 10.0

        block_lo = _slab_block(layer, 0.0)       # z = 0, 1, 2
        block_hi = _slab_block(layer, d - 2.0)   # z = 8, 9, 10

        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block_lo, block_hi], outer_faces=[], translational_direction="z",
        )

        assert len(periodic_pairs) == 1
        assert len(periodic_export) == 1
        rec = periodic_export[0]
        block_indices = {rec["block1"]["block_index"], rec["block2"]["block_index"]}
        assert block_indices == {0, 1}
        # A clean, fully-conformal pair must certify unambiguously.
        assert rec["orientation"]["permutation_index"] in range(-1, 8)


class TestCertificationCatchesInvalidCorrespondence:
    """A pair can pass the coverage-fraction test -- every point
    individually coincides somewhere in the other face's cloud -- without
    being a genuine structured (conformal) interface, e.g. two interior
    nodes' (x, y) values are swapped between (i, j) slots. The orthogonal
    precheck tests point-CLOUD coincidence (order independent), so it
    cannot see this; node-for-node certification of the located sub-patch
    must reject it."""

    def test_swapped_interior_nodes_rejected_despite_full_coverage(self):
        ni, nj, nk = 4, 4, 3
        swap = ((1, 1), (2, 1))  # both strictly interior for a 4x4 face

        layer = _asym_layer(ni, nj)
        (i1, j1), (i2, j2) = swap
        top_swapped = layer.copy()
        top_swapped[i1, j1], top_swapped[i2, j2] = (
            top_swapped[i2, j2].copy(), top_swapped[i1, j1].copy())

        # Premise check: the swap only reassigns which (i, j) owns which
        # (x, y) value -- the top layer's point SET is unchanged, so a
        # point-cloud coincidence test (exactly what the orthogonal
        # precheck runs) sees 100% overlap and would accept this pair.
        assert coincidence_count(
            layer.reshape(-1, 2), top_swapped.reshape(-1, 2), tol=1e-9,
        ) == ni * nj

        X = np.repeat(layer[:, :, 0:1], nk, axis=2)
        Y = np.repeat(layer[:, :, 1:2], nk, axis=2)
        Z = np.broadcast_to(np.arange(nk, dtype=float), (ni, nj, nk)).copy()
        X[i1, j1, nk - 1], X[i2, j2, nk - 1] = X[i2, j2, nk - 1], X[i1, j1, nk - 1]
        Y[i1, j1, nk - 1], Y[i2, j2, nk - 1] = Y[i2, j2, nk - 1], Y[i1, j1, nk - 1]
        block = Block(X, Y, Z)

        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block], outer_faces=[], translational_direction="z",
        )

        # Certification (not the coverage-fraction test) is what rejects
        # this -- the premise check above already confirmed coverage alone
        # would have accepted it.
        assert periodic_pairs == []
        assert periodic_export == []


class TestObliqueFallbackCertification:
    """The bladed-cascade oblique-fallback path (Phase 3: footprint filter,
    per-pair median offset, quantized-3D-intersection verification) builds
    its candidate matches differently from the main flat/global-extent path
    above, so it needs its own coverage of the shared
    `_locate_and_certify_periodic_patch` wiring.

    The oblique path only activates when the CALLER supplies a real
    `outer_faces` list -- `outer_face_dict_to_list([], ...)` returns `[]`
    unconditionally (unlike `find_bounding_faces`, which auto-detects outer
    faces from the blocks when given an empty list), so every other test in
    this file (which passes `outer_faces=[]`) always has an empty
    candidate pool for this path and never reaches it. No bundled fixture
    in `tests/data/` exercises it either (the only real-mesh caller,
    `test_weld.py`, needs a 1900-block proprietary mesh not present here),
    so this test forces the oblique path with the same simple, already-
    trusted synthetic geometry used above by monkeypatching
    `find_bounding_faces` to report empty lower/upper pools -- exactly as
    if the periodic faces were not at the block's global axis extreme --
    while supplying `outer_faces` so the oblique candidate pool is
    populated.
    """

    def test_oblique_path_finds_and_certifies_clean_pair(self, monkeypatch):
        import sys
        from plot3d.facefunctions import get_outer_faces

        periodicity_mod = sys.modules["plot3d.periodicity"]
        monkeypatch.setattr(
            periodicity_mod, "find_bounding_faces",
            lambda *a, **k: ([], [], [], []))

        ni, nj = 4, 4
        layer = _asym_layer(ni, nj)
        block = _slab_block(layer, 0.0, nk=3)  # z = 0, 1, 2

        outer_faces_list, _ = get_outer_faces(block)
        for f in outer_faces_list:
            f.set_block_index(0)
        outer_faces = [f.to_dict() for f in outer_faces_list]

        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block], outer_faces=outer_faces, translational_direction="z",
        )

        assert len(periodic_pairs) == 1
        assert len(periodic_export) == 1
        assert periodic_export[0]["mode"] == "z_oblique_pair"

    def test_oblique_path_rejects_swapped_interior_nodes(self, monkeypatch):
        """Same swap-defeats-coverage-but-not-certification technique as
        `TestCertificationCatchesInvalidCorrespondence`, routed through the
        oblique-fallback path instead of the main path."""
        import sys
        from plot3d.facefunctions import get_outer_faces

        periodicity_mod = sys.modules["plot3d.periodicity"]
        monkeypatch.setattr(
            periodicity_mod, "find_bounding_faces",
            lambda *a, **k: ([], [], [], []))

        ni, nj, nk = 4, 4, 3
        swap = ((1, 1), (2, 1))
        layer = _asym_layer(ni, nj)
        (i1, j1), (i2, j2) = swap

        X = np.repeat(layer[:, :, 0:1], nk, axis=2)
        Y = np.repeat(layer[:, :, 1:2], nk, axis=2)
        Z = np.broadcast_to(np.arange(nk, dtype=float), (ni, nj, nk)).copy()
        X[i1, j1, nk - 1], X[i2, j2, nk - 1] = X[i2, j2, nk - 1], X[i1, j1, nk - 1]
        Y[i1, j1, nk - 1], Y[i2, j2, nk - 1] = Y[i2, j2, nk - 1], Y[i1, j1, nk - 1]
        block = Block(X, Y, Z)

        outer_faces_list, _ = get_outer_faces(block)
        for f in outer_faces_list:
            f.set_block_index(0)
        outer_faces = [f.to_dict() for f in outer_faces_list]

        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block], outer_faces=outer_faces, translational_direction="z",
        )

        assert periodic_pairs == []
        assert periodic_export == []


# ---------------------------------------------------------------------
# A5: full-resolution re-validation after GCD reduction
#
# `translational_periodicity` GCD-reduces unconditionally (no flag). Unlike
# `rotated_periodicity`'s single global forward/backward rotation, each pair
# here carries its own per-pair axis shift, so re-validation is a direct
# inline certification loop (see `periodicity.py`, step "9b") rather than a
# shared `transforms` list. These tests cover:
#   1. A no-op regression on a synthetic box mesh sized so GCD reduction
#      actually applies (gcd > 1) -- nothing should be demoted.
#   2. A constructed case where a full-resolution-only interior node is
#      perturbed beyond tolerance -- the pair must be demoted, with a
#      RuntimeWarning, using this function's own per-pair-shift transform.
# ---------------------------------------------------------------------

import warnings as _warnings

from plot3d.blockfunctions import compute_min_gcd


def test_translational_periodicity_a5_no_op_when_gcd_greater_than_one():
    """GCD reduction must actually be in play (gcd > 1) for this check to
    be meaningful; every reduced-grid proposal found on a clean box mesh
    must survive full-resolution re-validation with no RuntimeWarning.
    """
    block = _box_block(5, 5, 9)  # (ni-1, nj-1, nk-1) = (4, 4, 8) -> gcd 4
    assert compute_min_gcd([block]) > 1

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block], outer_faces=[], translational_direction="z",
        )

    demotions = [w for w in caught
                 if issubclass(w.category, RuntimeWarning)
                 and "translational_periodicity" in str(w.message)]
    assert not demotions, f"unexpected demotions on clean mesh: {[str(w.message) for w in demotions]}"
    assert len(periodic_pairs) == 1
    assert len(periodic_export) == 1


def test_translational_periodicity_a5_demotes_full_resolution_only_perturbation():
    """A full-resolution-only interior node (I=1, J=1 on the z=0 face,
    skipped by GCD=4 reduction which only samples I,J in {0, 4}) is
    perturbed beyond tolerance. The reduced grid never sees the
    perturbation and finds the pair periodic; full-resolution
    re-validation must catch it via this function's own per-pair shift
    transform, demote the pair, and raise a RuntimeWarning naming
    `translational_periodicity`.
    """
    block = _box_block(5, 5, 9)
    assert compute_min_gcd([block]) == 4

    # Perturb an in-plane coordinate of a node that only exists at full
    # resolution (I=1 is skipped by the GCD=4 reduced sampling {0, 4}).
    block.X[1, 1, 0] += 0.05

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        periodic_export, periodic_pairs, remaining_outer = translational_periodicity(
            [block], outer_faces=[], translational_direction="z",
            node_tol_xyz=1e-3,
        )

    demotions = [w for w in caught
                 if issubclass(w.category, RuntimeWarning)
                 and "translational_periodicity" in str(w.message)]
    assert demotions, "expected a RuntimeWarning demoting the perturbed pair"
    assert periodic_pairs == []
    assert periodic_export == []
