"""Regression tests for the node-for-node certification sites added to
``connectivity.py`` (plot3d-rs port, mirrors Rust commit ``0e1b1a1``):

1. ``_try_permutations_with_transpose`` now scans all 8 permutations and
   accepts iff *exactly one* passes -- an ambiguous match (more than one
   permutation within tolerance) is rejected instead of silently accepting
   whichever permutation was tried first.
2. ``get_face_intersection``'s Step 3 (per-point geometric fallback) now
   certifies the claimed sub-patch via ``correspondence.certify_correspondence``
   after its existing head-count guard passes -- the head-count check alone
   cannot detect a matched point set that satisfies the count but does not
   form a valid, contiguous structured correspondence.
3. ``revalidate_full_resolution``/``demote_to_outer`` -- ``connectivity_fast``
   re-certifies every GCD-reduced-grid match proposal against the ORIGINAL
   full-resolution mesh before returning it, since GCD reduction can hide
   an interior-node perturbation that never survives into the reduced grid.
   A proposal that fails full-resolution certification is demoted: dropped
   from ``face_matches`` and both its faces added back to ``outer_faces``.
"""

import warnings

import numpy as np
import pytest

from plot3d import Block
from plot3d.connectivity import (
    _try_permutations_with_transpose,
    connectivity_fast,
    demote_to_outer,
    get_face_intersection,
    revalidate_full_resolution,
)
from plot3d.facefunctions import create_face_from_diagonals


def _asym_block(nu: int, nv: int) -> Block:
    """A single-k-layer block with no reversal/swap symmetry, so "exactly
    one permutation passes" isn't an accident of a too-symmetric fixture.
    Same construction as ``test_correspondence.py``'s ``_asym_block``.
    """
    i = np.arange(nu, dtype=float)
    j = np.arange(nv, dtype=float)
    I, J = np.meshgrid(i, j, indexing="ij")
    X = I * 1.7
    Y = J * 0.9
    Z = 0.05 * I * I + 0.03 * J + 0.01 * I * J
    return Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())


class TestTryPermutationsAmbiguity:
    def test_symmetric_face_is_rejected_not_first_match_accepted(self):
        # Coordinates depend only on j: every point is identical across the
        # i axis, so identity (perm 0) and i-reversed (perm 1) both certify
        # exactly. Before the fix, the first-found passing permutation
        # (perm 0) would have been accepted; now this must be rejected.
        nu, nv = 3, 2
        i = np.arange(nu, dtype=float)
        j = np.arange(nv, dtype=float)
        I, J = np.meshgrid(i, j, indexing="ij")
        X = J * 1.0
        Y = J * 2.0
        Z = np.zeros_like(J)
        block1 = Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())
        block2 = Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())

        lb = [0, 0, 0]
        ub = [nu - 1, nv - 1, 0]

        matched, df = _try_permutations_with_transpose(block1, lb, ub, block2, lb, ub, tol=1e-9)

        assert matched is False
        assert len(df) == 0
        assert list(df.columns) == ['i1', 'j1', 'k1', 'i2', 'j2', 'k2']

    def test_unambiguous_same_size_match_still_succeeds(self):
        # No regression on the ordinary case: an asymmetric grid matched to
        # an identical copy of itself has exactly one certifying
        # permutation (identity, perm 0) and must still be accepted.
        nu, nv = 4, 3
        block1 = _asym_block(nu, nv)
        block2 = _asym_block(nu, nv)

        lb = [0, 0, 0]
        ub = [nu - 1, nv - 1, 0]

        matched, df = _try_permutations_with_transpose(block1, lb, ub, block2, lb, ub, tol=1e-9)

        assert matched is True
        assert len(df) == nu * nv
        assert (df['i1'] == df['i2']).all()
        assert (df['j1'] == df['j2']).all()
        assert (df['k1'] == df['k2']).all()


class TestGetFaceIntersectionStep3Certification:
    def _build_blocks(self):
        # block1: a genuine 2x2 physical square, K-constant face.
        i = np.arange(2, dtype=float)
        j = np.arange(2, dtype=float)
        I, J = np.meshgrid(i, j, indexing="ij")
        X1 = I * 1.0
        Y1 = J * 1.0
        Z1 = np.zeros_like(I)
        block1 = Block(X1[:, :, None].copy(), Y1[:, :, None].copy(), Z1[:, :, None].copy())

        # block2: a 3x2 face whose I=0 and I=2 rows physically coincide with
        # block1's i=0 and i=1 rows (a genuine corner-to-corner match), but
        # whose middle row (I=1) is placed far away and corresponds to
        # nothing in block1 -- so the true target sub-region on block2 is
        # the *non-contiguous* index set {0, 2} x {0, 1}, not a real
        # rectangular sub-block.
        i2 = np.arange(3, dtype=float)
        j2 = np.arange(2, dtype=float)
        I2, J2 = np.meshgrid(i2, j2, indexing="ij")
        X2 = np.where(I2 == 1, 100.0 + J2, np.where(I2 == 2, 1.0, 0.0))
        Y2 = np.where(I2 == 1, 100.0 + J2, J2 * 1.0)
        Z2 = np.where(I2 == 1, 100.0, 0.0)
        block2 = Block(X2[:, :, None].copy(), Y2[:, :, None].copy(), Z2[:, :, None].copy())

        return block1, block2

    def test_noncontiguous_target_with_matching_headcount_is_rejected(self):
        block1, block2 = self._build_blocks()

        face1 = create_face_from_diagonals(block1, [0, 0, 0], [1, 1, 0])
        face2 = create_face_from_diagonals(block2, [0, 0, 0], [2, 1, 0])
        face1.set_block_index(0)
        face2.set_block_index(1)

        tol = 1e-6

        # Sanity: face sizes differ (n1=4, n2=6), so Step 1's same-size fast
        # path does not apply, and Step 2's subregion search must fail to
        # locate a same-size rectangular subregion (the only candidate
        # sub-block spans the full non-matching I=0..2 range), forcing
        # Step 3's per-point fallback to run. All 4 of face1's points do
        # have a true physical coincidence on face2 (at I in {0, 2}), so
        # the pre-existing head-count check (matched points == claimed
        # sub-patch bounding-box area) passes -- this is exactly the
        # necessary-but-not-sufficient case the certification pass exists
        # to catch.
        df, split1, split2 = get_face_intersection(face1, face2, block1, block2, tol)

        assert len(df) == 0, (
            "expected rejection: the matched points' bounding box on "
            "block2 spans a non-contiguous index range ({0,2} x {0,1}), "
            "not a valid structured correspondence, even though the "
            "head-count guard alone would have accepted it"
        )

    def test_certification_success_path_still_accepts_the_match(self, monkeypatch):
        # Same fixture as the rejection test above (a case that, by
        # construction, only reaches Step 3 when the target sub-region
        # genuinely is *not* a valid same-shape structured correspondence
        # -- see the analysis in that test's docstring: whenever Step 2's
        # own corner-to-corner bounding-box check would have found a
        # same-size subregion, Step 2 already handles the match via
        # ``_try_permutations_with_transpose``, so a "naturally occurring,
        # should-be-accepted" case that *also* reaches Step 3 with a
        # full-size headcount is not constructible here). Instead, this
        # confirms the "certification succeeds -> keep the match" branch
        # is wired correctly by stubbing ``certify_correspondence`` to
        # succeed and checking the match survives -- the rejection test
        # above already proves the "certification fails -> reject" branch
        # fires on a real, unmocked case; this proves the converse wiring
        # (success does not get accidentally discarded) without which the
        # try/except could be trivially (and wrongly) written to always
        # reject.
        block1, block2 = self._build_blocks()
        face1 = create_face_from_diagonals(block1, [0, 0, 0], [1, 1, 0])
        face2 = create_face_from_diagonals(block2, [0, 0, 0], [2, 1, 0])
        face1.set_block_index(0)
        face2.set_block_index(1)

        import plot3d.connectivity as connectivity_mod

        def _fake_certify_correspondence(*args, **kwargs):
            return "certified"  # any non-raising return value; result is unused

        monkeypatch.setattr(
            connectivity_mod.correspondence, "certify_correspondence",
            _fake_certify_correspondence)

        df, split1, split2 = get_face_intersection(face1, face2, block1, block2, 1e-6)

        assert len(df) == 4


def _cube_pair_blocks(perturb: bool = False):
    """Two 3x3x3 unit-spaced cube blocks glued face-to-face at x=2.

    Block 0 spans I in [0, 2] (x in [0, 2]); block 1 spans I in [0, 2] with
    an x offset of 2 (x in [2, 4]).  Block 0's I=2 face physically coincides
    with block 1's I=0 face -- a genuine 3x3-node interface.

    ``IMAX == JMAX == KMAX == 3`` on every block gives
    ``compute_min_gcd == gcd(2, gcd(2, 2)) == 2``, so ``connectivity_fast``
    reduces to indices {0, 2} in every axis -- the interior index 1 does not
    survive into the reduced grid.

    When ``perturb`` is True, block 0's interface-face interior node
    (I=2, J=1, K=1) is displaced by 1e-2 in Y -- far beyond the adaptive
    tolerance (~1e-5 at this coordinate magnitude), but invisible to the
    GCD-reduced grid since index 1 is dropped by the reduction.
    """
    i = np.arange(3, dtype=float)
    j = np.arange(3, dtype=float)
    k = np.arange(3, dtype=float)
    I, J, K = np.meshgrid(i, j, k, indexing='ij')

    X0 = I.copy()
    Y0 = J.copy()
    Z0 = K.copy()
    if perturb:
        Y0[2, 1, 1] += 1e-2
    block0 = Block(X0, Y0, Z0)

    X1 = I.copy() + 2.0
    Y1 = J.copy()
    Z1 = K.copy()
    block1 = Block(X1, Y1, Z1)

    return [block0, block1]


class TestConnectivityFastFullResolutionRevalidation:
    def test_clean_mesh_gcd_reduction_is_a_no_op(self):
        # GCD reduction *is* triggered (gcd_to_use == 2), and every
        # reduced-grid proposal survives full-resolution certification
        # since the mesh is clean -- confirms no spurious demotions.
        import warnings

        blocks = _cube_pair_blocks(perturb=False)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            face_matches, outer_faces = connectivity_fast(blocks)

        runtime_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
        assert runtime_warnings == []

        interface_matches = [
            m for m in face_matches
            if {m['block1']['block_index'], m['block2']['block_index']} == {0, 1}
        ]
        assert len(interface_matches) == 1

        m = interface_matches[0]
        b0_side = m['block1'] if m['block1']['block_index'] == 0 else m['block2']
        lb, ub = b0_side['lb'], b0_side['ub']
        assert min(lb[0], ub[0]) == max(lb[0], ub[0]) == 2
        assert sorted([min(lb[1], ub[1]), max(lb[1], ub[1])]) == [0, 2]
        assert sorted([min(lb[2], ub[2]), max(lb[2], ub[2])]) == [0, 2]

    def test_perturbed_interior_node_is_demoted_to_outer_faces(self):
        blocks = _cube_pair_blocks(perturb=True)

        with pytest.warns(RuntimeWarning, match="connectivity_fast"):
            face_matches, outer_faces = connectivity_fast(blocks)

        # (a) The interface match must NOT survive into face_matches.
        interface_matches = [
            m for m in face_matches
            if {m['block1']['block_index'], m['block2']['block_index']} == {0, 1}
        ]
        assert interface_matches == []

        # (b) Both faces of the demoted proposal must appear as outer faces.
        def _is_full_interface_face(o, block_index, i_plane):
            if o['block_index'] != block_index:
                return False
            lb, ub = o['lb'], o['ub']
            same_i = min(lb[0], ub[0]) == max(lb[0], ub[0]) == i_plane
            spans_j = sorted([min(lb[1], ub[1]), max(lb[1], ub[1])]) == [0, 2]
            spans_k = sorted([min(lb[2], ub[2]), max(lb[2], ub[2])]) == [0, 2]
            return same_i and spans_j and spans_k

        outer_block0 = [o for o in outer_faces if _is_full_interface_face(o, 0, 2)]
        outer_block1 = [o for o in outer_faces if _is_full_interface_face(o, 1, 0)]
        assert len(outer_block0) == 1, outer_faces
        assert len(outer_block1) == 1, outer_faces


class TestRevalidateFullResolutionUnit:
    """Direct unit coverage of revalidate_full_resolution/demote_to_outer
    against hand-built proposal dicts, independent of the full
    connectivity_fast pipeline."""

    def _blocks(self):
        return _cube_pair_blocks(perturb=False)

    def _good_proposal(self):
        return {
            'block1': {'block_index': 0, 'lb': [2, 0, 0], 'ub': [2, 2, 2], 'id': 1},
            'block2': {'block_index': 1, 'lb': [0, 0, 0], 'ub': [0, 2, 2], 'id': 2},
            'orientation': {
                'permutation_index': -1,
                'plane': 'in-plane',
                'permutation_matrix': [[1, 0], [0, 1]],
            },
        }

    def test_declared_orientation_that_certifies_is_kept(self):
        blocks = self._blocks()
        proposal = self._good_proposal()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kept, rejected = revalidate_full_resolution(
                blocks, [proposal], transforms=[lambda p: p], tol=1e-6,
                stage="unit-test",
            )

        assert kept == [proposal]
        assert rejected == []
        assert not any(issubclass(w.category, RuntimeWarning) for w in caught)

    def test_declared_orientation_that_fails_is_rejected_not_replaced_by_search(self):
        # An asymmetric grid whose only true certifying permutation is
        # identity (perm 0): block0 and block1 are identical copies of an
        # asymmetric face. The proposal DECLARES a reversed I direction on
        # block2's bounds (lb2/ub2 swapped along I) -- a wrong declared
        # correspondence, pinning down permutation 1 (u-reversed), which
        # does not hold for this asymmetric grid. certify_permutation must
        # fail as declared -- it must never silently fall back to a search
        # that WOULD find the correct permutation 0.
        nu, nv = 4, 3
        i = np.arange(nu, dtype=float)
        j = np.arange(nv, dtype=float)
        I, J = np.meshgrid(i, j, indexing='ij')
        X = I * 1.7
        Y = J * 0.9
        Z = 0.05 * I * I + 0.03 * J + 0.01 * I * J
        block0 = Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())
        block1 = Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())
        blocks = [block0, block1]

        proposal = {
            'block1': {'block_index': 0, 'lb': [0, 0, 0], 'ub': [nu - 1, nv - 1, 0], 'id': 1},
            'block2': {'block_index': 1, 'lb': [nu - 1, 0, 0], 'ub': [0, nv - 1, 0], 'id': 2},
            'orientation': {
                'permutation_index': 1,
                'plane': 'in-plane',
                'permutation_matrix': [[-1, 0], [0, 1]],
            },
        }

        with pytest.warns(RuntimeWarning, match="unit-test"):
            kept, rejected = revalidate_full_resolution(
                blocks, [proposal], transforms=[lambda p: p], tol=1e-9,
                stage="unit-test",
            )

        assert kept == []
        assert rejected == [proposal]

        # Sanity: an undeclared (searching) certification of the SAME
        # patches WOULD have succeeded, via permutation 0 -- proving the
        # declared path rejected on its own merits, not because no
        # permutation exists for this geometry at all.
        from plot3d import correspondence
        from plot3d.permutation import patch_from_bounds

        patch1 = patch_from_bounds(0, proposal['block1']['lb'], proposal['block1']['ub'])
        patch2 = patch_from_bounds(1, proposal['block2']['lb'], proposal['block2']['ub'])
        result = correspondence.certify_correspondence(block0, patch1, block1, patch2, 1e-9)
        assert result.permutation_index == 0

    def test_real_cross_plane_match_is_not_wrongly_demoted(self):
        # Regression guard: connectivity.py's exported
        # orientation['permutation_matrix'] uses _orient_vec_to_permutation's
        # own bit convention, which is keyed to face1's own axis order and
        # does NOT, in general, coincide with correspondence.py's canonical
        # (ascending, face2-keyed) bit convention for a swapped/cross-plane
        # match. Trusting that matrix's literal value as a
        # correspondence.py permutation index would wrongly reject this
        # genuine match (declared permutation 6 fails; the actual
        # certifying permutation is 5) -- revalidate_full_resolution must
        # derive the declared index from the proposal's own bounds instead
        # (see _perm_idx_from_declared_bounds), not from the matrix.
        import os

        from plot3d import read_plot3D

        mesh_path = os.path.join(
            os.path.dirname(__file__), "data", "cross_plane_pair.p3d")
        if not os.path.exists(mesh_path):
            pytest.skip("cross_plane_pair.p3d not found")

        blocks = read_plot3D(mesh_path)
        assert len(blocks) == 2

        # Exercise revalidate_full_resolution directly against the actual
        # declared match connectivity() finds on this mesh -- mirroring
        # exactly what connectivity_fast does when gcd_to_use > 1.
        from plot3d.connectivity import adaptive_tolerance, connectivity as _connectivity

        tol = adaptive_tolerance(blocks)
        proposed, _outer = _connectivity(blocks, tol=tol)
        assert len(proposed) == 1

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kept, rejected = revalidate_full_resolution(
                blocks, proposed, transforms=[lambda p: p], tol=tol,
                stage="regression-test",
            )

        assert rejected == [], (
            "the genuine cross-plane match must survive full-resolution "
            "re-validation -- it must not be wrongly demoted because the "
            "exported orientation.permutation_matrix uses a different bit "
            "convention than correspondence.py"
        )
        assert kept == proposed
        assert not any(issubclass(w.category, RuntimeWarning) for w in caught)

    def test_undeclared_proposal_is_searched_and_certifies(self):
        # No 'orientation' key at all (e.g. a self-match record) -> must be
        # searched via certify_correspondence rather than rejected outright.
        blocks = self._blocks()
        proposal = self._good_proposal()
        del proposal['orientation']

        kept, rejected = revalidate_full_resolution(
            blocks, [proposal], transforms=[lambda p: p], tol=1e-6,
            stage="unit-test",
        )

        assert kept == [proposal]
        assert rejected == []

    def test_perturbed_full_res_node_fails_certification(self):
        blocks = _cube_pair_blocks(perturb=True)
        proposal = self._good_proposal()

        with pytest.warns(RuntimeWarning, match=r"demoted 1 proposal"):
            kept, rejected = revalidate_full_resolution(
                blocks, [proposal], transforms=[lambda p: p], tol=1e-6,
                stage="unit-test",
            )

        assert kept == []
        assert rejected == [proposal]

    def test_demote_to_outer_appends_both_faces_with_fresh_ids(self):
        outer_faces = [{'block_index': 5, 'lb': [0, 0, 0], 'ub': [0, 1, 1], 'id': 3}]
        rejected = [self._good_proposal()]

        demote_to_outer(outer_faces, rejected)

        assert len(outer_faces) == 3
        new_ids = {o['id'] for o in outer_faces[1:]}
        assert new_ids == {4, 5}
        block_indices = {o['block_index'] for o in outer_faces[1:]}
        assert block_indices == {0, 1}

    def test_demote_to_outer_deduplicates_against_existing_and_within_batch(self):
        existing = {
            'block_index': 0, 'lb': [2, 0, 0], 'ub': [2, 2, 2], 'id': 7,
        }
        outer_faces = [existing]
        proposal_a = self._good_proposal()
        proposal_b = self._good_proposal()  # identical faces -> duplicate
        rejected = [proposal_a, proposal_b]

        demote_to_outer(outer_faces, rejected)

        # block1 side duplicates `existing` (already present) and is not
        # re-added; block2 side is new and appears exactly once despite
        # being rejected twice in this batch.
        assert len(outer_faces) == 2
        block1_side = [o for o in outer_faces if o['block_index'] == 1]
        assert len(block1_side) == 1
        assert block1_side[0]['id'] == 8


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
