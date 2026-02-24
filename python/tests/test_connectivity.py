"""Tests for connectivity detection between blocks."""
import numpy as np
import pytest
from plot3d import Block
from plot3d.connectivity import connectivity, connectivity_fast, face_matches_to_dict


class TestConnectivityTwoBlocks:
    def test_adjacent_blocks_share_one_face(self, two_adjacent_blocks):
        face_matches, outer_faces = connectivity(two_adjacent_blocks)
        assert len(face_matches) == 1

    def test_match_structure(self, two_adjacent_blocks):
        face_matches, outer_faces = connectivity(two_adjacent_blocks)
        m = face_matches[0]
        assert 'block1' in m
        assert 'block2' in m
        for side in ['block1', 'block2']:
            for key in ['block_index', 'IMIN', 'JMIN', 'KMIN', 'IMAX', 'JMAX', 'KMAX']:
                assert key in m[side]

    def test_outer_faces_count(self, two_adjacent_blocks):
        """Two cubes sharing one face -> 10 outer faces total (6+6-2)."""
        face_matches, outer_faces = connectivity(two_adjacent_blocks)
        assert len(outer_faces) == 10

    def test_matched_face_at_x_equals_1(self, two_adjacent_blocks):
        face_matches, _ = connectivity(two_adjacent_blocks)
        m = face_matches[0]
        # Block 0's imax face at i=4
        b1 = m['block1']
        b2 = m['block2']
        # One block has IMIN==IMAX (constant-i face)
        b1_const_i = b1['IMIN'] == b1['IMAX']
        b2_const_i = b2['IMIN'] == b2['IMAX']
        assert b1_const_i or b2_const_i


class TestConnectivityFourBlocks:
    def test_four_block_grid_connectivity(self, four_block_grid):
        face_matches, outer_faces = connectivity(four_block_grid)
        # 2x2 grid -> 4 internal shared faces
        assert len(face_matches) == 4

    def test_outer_faces_count(self, four_block_grid):
        """2x2 grid of cubes: 4*6=24 total - 4*2=8 internal = 16 outer."""
        face_matches, outer_faces = connectivity(four_block_grid)
        assert len(outer_faces) == 16


class TestConnectivityFast:
    def test_fast_matches_standard(self, two_adjacent_blocks):
        """connectivity_fast should produce same result count as connectivity."""
        fm_fast, of_fast = connectivity_fast(two_adjacent_blocks)
        fm_std, of_std = connectivity(two_adjacent_blocks)
        assert len(fm_fast) == len(fm_std)
        assert len(of_fast) == len(of_std)


class TestConnectivityResultStructure:
    def test_face_match_indices_consistent(self, two_adjacent_blocks):
        """Face match block indices should refer to valid blocks."""
        fm, of = connectivity(two_adjacent_blocks)
        for m in fm:
            assert 0 <= m['block1']['block_index'] < len(two_adjacent_blocks)
            assert 0 <= m['block2']['block_index'] < len(two_adjacent_blocks)

    def test_outer_face_has_required_keys(self, two_adjacent_blocks):
        """Each outer face dict should have index and dimension keys."""
        _, of = connectivity(two_adjacent_blocks)
        for o in of:
            for key in ['block_index', 'IMIN', 'JMIN', 'KMIN', 'IMAX', 'JMAX', 'KMAX']:
                assert key in o
