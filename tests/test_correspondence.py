"""Node-for-node face-correspondence certification.

Mirrors plot3d-rs's ``correspondence.rs`` (commit ``0e1b1a1``): a conformal
structured interface is not "four corners agree" -- every node, interior
nodes included, must correspond under exactly one of the 8 structured
permutation mappings, within tolerance. Zero permutations passing is
``ExceedsTolerance``; more than one passing is ``Ambiguous`` (the patch
geometry is too symmetric/degenerate to certify a unique mapping); no
dimension-compatible permutation at all is ``IncompatibleDimensions``.
``certify_permutation`` certifies one declared permutation only, with no
fallback search -- a wrong declared orientation must fail even when a
different orientation would have certified.
"""

import numpy as np
import pytest

from plot3d import Block
from plot3d.permutation import patch_from_bounds, apply_permutation
from plot3d.correspondence import (
    Ambiguous,
    CertifiedMapping,
    ExceedsTolerance,
    IncompatibleDimensions,
    NodeDiscrepancy,
    certify_correspondence,
    certify_permutation,
)


def _asym_block(nu: int, nv: int) -> Block:
    """A single-k-layer block whose (i, j) face has no reversal/swap
    symmetry -- every permutation produces geometrically distinct points,
    so a certification test's "exactly one permutation passes" outcome
    isn't an accident of a too-symmetric fixture.
    """
    i = np.arange(nu, dtype=float)
    j = np.arange(nv, dtype=float)
    I, J = np.meshgrid(i, j, indexing="ij")
    X = I * 1.7
    Y = J * 0.9
    Z = 0.05 * I * I + 0.03 * J + 0.01 * I * J
    return Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())


def _full_face_patch(block_index: int, nu: int, nv: int):
    return patch_from_bounds(block_index, (0, 0, 0), (nu - 1, nv - 1, 0))


def _block_from_grid(grid: np.ndarray) -> Block:
    """Build a single-k-layer Block from an (nu, nv, 3) point grid."""
    return Block(grid[:, :, 0:1].copy(), grid[:, :, 1:2].copy(), grid[:, :, 2:3].copy())


def _face_grid(block: Block) -> np.ndarray:
    return np.stack([block.X[:, :, 0], block.Y[:, :, 0], block.Z[:, :, 0]], axis=-1)


class TestKnownUnambiguousTransform:
    """A patch related to itself by a known, unambiguous permutation."""

    def test_identity_match_certifies_permutation_zero(self):
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        mapping = certify_correspondence(block_a, patch_a, block_a, patch_a, tol=1e-9)

        assert isinstance(mapping, CertifiedMapping)
        assert mapping.permutation_index == 0
        assert mapping.plane == "in-plane"
        assert mapping.nodes_checked == nu * nv
        assert mapping.worst.distance == pytest.approx(0.0, abs=1e-12)

    def test_reversal_transform_certifies_correct_permutation(self):
        # perm 1 = u reversed; an involution, so applying it once to build
        # block B's raw grid means perm 1 is exactly what certifies it.
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        grid_b_raw = apply_permutation(_face_grid(block_a), 1)
        block_b = _block_from_grid(grid_b_raw)
        patch_b = _full_face_patch(1, nu, nv)

        mapping = certify_correspondence(block_a, patch_a, block_b, patch_b, tol=1e-9)

        assert mapping.permutation_index == 1
        assert mapping.plane == "in-plane"

    def test_swap_transform_certifies_correct_permutation(self):
        # perm 4 = axes swapped; also an involution. Rectangular (nu != nv)
        # patches so only the swap-family permutations (4-7) are even
        # dimension-compatible.
        nu, nv = 3, 4
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        grid_b_raw = apply_permutation(_face_grid(block_a), 4)  # shape (nv, nu, 3)
        block_b = _block_from_grid(grid_b_raw)
        patch_b = _full_face_patch(1, nv, nu)

        mapping = certify_correspondence(block_a, patch_a, block_b, patch_b, tol=1e-9)

        assert mapping.permutation_index == 4
        assert mapping.plane == "in-plane"


class TestExceedsTolerance:
    def test_single_node_nudged_outside_tolerance_is_rejected(self):
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        block_b = _asym_block(nu, nv)
        perturbed_ijk = (2, 1, 0)
        block_b.X[perturbed_ijk] += 10.0
        patch_b = _full_face_patch(1, nu, nv)

        with pytest.raises(ExceedsTolerance) as excinfo:
            certify_correspondence(block_a, patch_a, block_b, patch_b, tol=1e-6)

        err = excinfo.value
        # Every other permutation is grossly wrong across the whole patch
        # (the grid is asymmetric under every reversal/swap), so identity
        # (perm 0) is still the closest candidate, and its worst node is
        # exactly the one that was perturbed.
        assert err.best_permutation == 0
        assert isinstance(err.worst, NodeDiscrepancy)
        assert err.worst.node_a == perturbed_ijk
        assert err.worst.node_b == perturbed_ijk
        assert err.worst.distance == pytest.approx(10.0, rel=1e-9)


class TestAmbiguous:
    def test_symmetric_patch_certifies_under_two_permutations(self):
        # Coordinates depend only on j: every point is identical across the
        # u axis, so identity (perm 0) and u-reversed (perm 1) both certify
        # exactly -- the geometry cannot distinguish them.
        nu, nv = 3, 2
        i = np.arange(nu, dtype=float)
        j = np.arange(nv, dtype=float)
        I, J = np.meshgrid(i, j, indexing="ij")
        X = J * 1.0
        Y = J * 2.0
        Z = np.zeros_like(J)
        block = Block(X[:, :, None].copy(), Y[:, :, None].copy(), Z[:, :, None].copy())
        patch = _full_face_patch(0, nu, nv)

        with pytest.raises(Ambiguous) as excinfo:
            certify_correspondence(block, patch, block, patch, tol=1e-9)

        assert excinfo.value.permutations == [0, 1]


class TestIncompatibleDimensions:
    def test_no_permutation_reconciles_the_shapes(self):
        block_a = _asym_block(3, 4)
        patch_a = _full_face_patch(0, 3, 4)

        block_c = _asym_block(5, 2)
        patch_c = _full_face_patch(1, 5, 2)

        with pytest.raises(IncompatibleDimensions) as excinfo:
            certify_correspondence(block_a, patch_a, block_c, patch_c, tol=1e-9)

        assert excinfo.value.dims_a == (3, 4)
        assert excinfo.value.dims_b == (5, 2)


class TestCertifyPermutationNoFallback:
    def test_correct_declared_permutation_succeeds(self):
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        grid_b_raw = apply_permutation(_face_grid(block_a), 1)
        block_b = _block_from_grid(grid_b_raw)
        patch_b = _full_face_patch(1, nu, nv)

        mapping = certify_permutation(block_a, patch_a, block_b, patch_b, perm_idx=1, tol=1e-9)

        assert mapping.permutation_index == 1

    def test_wrong_declared_permutation_fails_with_no_silent_fallback(self):
        # Patch B truly corresponds under perm 1, but a discovery stage
        # wrongly declares perm 0 (identity). certify_permutation must
        # reject it outright -- never silently substitute perm 1.
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        grid_b_raw = apply_permutation(_face_grid(block_a), 1)
        block_b = _block_from_grid(grid_b_raw)
        patch_b = _full_face_patch(1, nu, nv)

        with pytest.raises(ExceedsTolerance) as excinfo:
            certify_permutation(block_a, patch_a, block_b, patch_b, perm_idx=0, tol=1e-9)

        assert excinfo.value.best_permutation == 0

        # Confirm perm 1 (the un-tried alternative) really would have
        # certified, proving the rejection above wasn't just a bad fixture.
        mapping = certify_permutation(block_a, patch_a, block_b, patch_b, perm_idx=1, tol=1e-9)
        assert mapping.permutation_index == 1

    def test_incompatible_dimensions_raised_for_declared_permutation(self):
        block_a = _asym_block(3, 4)
        patch_a = _full_face_patch(0, 3, 4)
        block_c = _asym_block(5, 2)
        patch_c = _full_face_patch(1, 5, 2)

        with pytest.raises(IncompatibleDimensions):
            certify_permutation(block_a, patch_a, block_c, patch_c, perm_idx=0, tol=1e-9)


class TestTransform:
    def test_transform_corrects_a_constant_offset(self):
        nu, nv = 4, 3
        block_a = _asym_block(nu, nv)
        patch_a = _full_face_patch(0, nu, nv)

        shift = np.array([100.0, 0.0, 0.0])
        block_b = _asym_block(nu, nv)
        block_b.X += shift[0]
        patch_b = _full_face_patch(1, nu, nv)

        # Without correcting for the offset, no permutation is within
        # tolerance.
        with pytest.raises(ExceedsTolerance):
            certify_correspondence(block_a, patch_a, block_b, patch_b, tol=1e-6)

        # transform is applied to patch B's grid before comparison, so
        # subtracting the shift makes it certify exactly.
        mapping = certify_correspondence(
            block_a, patch_a, block_b, patch_b, tol=1e-9,
            transform=lambda g: g - shift)

        assert mapping.permutation_index == 0
        assert mapping.worst.distance == pytest.approx(0.0, abs=1e-9)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
