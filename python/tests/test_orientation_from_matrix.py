"""Matrix→canonical-index round-trip for face-match orientation.

Mirrors plot3d-rs ``tests/test_orientation_from_matrix.rs`` (commit 920149f).
Asserts that the new ``_index_from_permutation_matrix`` helper recognises
all 8 canonical 2×2 signed-permutation matrices, and rejects everything
else with the ``-1`` sentinel.

The DECLARED-vs-UNDECLARED rule from
``plot3d.verify._resolve_match_perm`` depends on this lookup being
correct; an off-by-one here silently corrupts every cascade verifier.
"""

import numpy as np
import pytest

from plot3d.verify import _index_from_permutation_matrix, _resolve_declared_perm
from plot3d.connectivity import PERMUTATION_MATRICES


class TestIndexFromPermutationMatrix:
    """Direct round-trip on the 8 canonical matrices."""

    def test_all_eight_matrices_round_trip(self):
        for expected_idx in range(8):
            matrix = PERMUTATION_MATRICES[expected_idx]
            got_idx = _index_from_permutation_matrix(matrix.tolist())
            assert got_idx == expected_idx, (
                f"PERMUTATION_MATRICES[{expected_idx}].tolist() round-tripped "
                f"to index {got_idx}; expected {expected_idx}."
            )

    def test_ndarray_input_works(self):
        for expected_idx in range(8):
            matrix_arr = PERMUTATION_MATRICES[expected_idx]
            got_idx = _index_from_permutation_matrix(matrix_arr)
            assert got_idx == expected_idx

    def test_none_returns_minus_one(self):
        assert _index_from_permutation_matrix(None) == -1

    def test_wrong_shape_returns_minus_one(self):
        assert _index_from_permutation_matrix([[1, 0]]) == -1
        assert _index_from_permutation_matrix([[1, 0, 0], [0, 1, 0]]) == -1
        assert _index_from_permutation_matrix([1, 2, 3, 4]) == -1

    def test_non_canonical_matrix_returns_minus_one(self):
        # The 8 canonical matrices are signed permutations; anything
        # outside that set must be rejected.
        bogus_examples = [
            [[2, 0], [0, 1]],   # not a permutation (entry magnitude > 1)
            [[1, 1], [0, 1]],   # not a permutation (two 1s in a row)
            [[0, 0], [0, 0]],   # zero matrix
            [[1, 0], [1, 0]],   # collinear rows
        ]
        for m in bogus_examples:
            assert _index_from_permutation_matrix(m) == -1, f"matrix {m} should be rejected"


class TestResolveDeclaredPerm:
    """The orientation-dict resolver consumes ``permutation_matrix``
    preferentially, then falls back to ``permutation_index``, then
    returns ``-1`` for UNDECLARED.
    """

    def test_matrix_is_preferred(self):
        # Stored ``permutation_index = -1`` (in-plane sentinel) but
        # ``permutation_matrix`` is canonical index 1 (u-flip).
        orientation = {
            'permutation_index': -1,
            'permutation_matrix': PERMUTATION_MATRICES[1].tolist(),
            'plane': 'in-plane',
        }
        assert _resolve_declared_perm(orientation) == 1

    def test_index_used_when_matrix_absent(self):
        orientation = {'permutation_index': 5, 'plane': 'cross-plane'}
        assert _resolve_declared_perm(orientation) == 5

    def test_undeclared_returns_minus_one(self):
        # Empty dict, missing keys, or only the in-plane sentinel.
        assert _resolve_declared_perm({}) == -1
        assert _resolve_declared_perm(None) == -1
        assert _resolve_declared_perm({'plane': 'in-plane'}) == -1
        assert _resolve_declared_perm({'permutation_index': -1}) == -1

    def test_invalid_matrix_falls_back_to_index(self):
        # Bogus matrix but valid stored index — fall through to index.
        orientation = {
            'permutation_index': 3,
            'permutation_matrix': [[2, 0], [0, 1]],  # rejected by matrix lookup
            'plane': 'in-plane',
        }
        assert _resolve_declared_perm(orientation) == 3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
