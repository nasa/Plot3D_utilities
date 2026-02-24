"""Tests for point matching."""
import numpy as np
import pytest
from plot3d.point_match import point_match


class TestPointMatch:
    def test_exact_match(self):
        """Point on a 2D face should return correct (i, j) indices."""
        # 3x3 face grid
        X = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]], dtype=float)
        Y = np.array([[0, 1, 2], [0, 1, 2], [0, 1, 2]], dtype=float)
        Z = np.zeros((3, 3))

        result = point_match(2.0, 2.0, 0.0, X, Y, Z)
        assert result[0] == 2  # i=2
        assert result[1] == 2  # j=2

    def test_within_tolerance(self):
        X = np.array([[0, 0], [1, 1]], dtype=float)
        Y = np.array([[0, 1], [0, 1]], dtype=float)
        Z = np.zeros((2, 2))
        result = point_match(1.0 + 1e-7, 1.0, 0.0, X, Y, Z, tol=1e-6)
        assert result[0] == 1
        assert result[1] == 1

    def test_no_match(self):
        X = np.array([[0, 0], [1, 1]], dtype=float)
        Y = np.array([[0, 1], [0, 1]], dtype=float)
        Z = np.zeros((2, 2))
        result = point_match(5.0, 5.0, 5.0, X, Y, Z)
        assert result[0] == -1
        assert result[1] == -1
