"""Point matching for structured grid face intersection.

Provides a vectorized point-to-grid matching function used by the connectivity
algorithm to find which grid node on face 2 corresponds to a given (x, y, z)
coordinate from face 1.
"""
import numpy as np


def point_match(x: float, y: float, z: float,
                X2: np.ndarray, Y2: np.ndarray, Z2: np.ndarray,
                tol: float = 1E-6):
    """Find the grid index on a face whose coordinates match a query point.

    Computes the Euclidean distance from ``(x, y, z)`` to every node in the
    2D grid defined by ``X2, Y2, Z2`` and returns the ``(p, q)`` indices of
    the closest node if it is within tolerance.

    Args:
        x: X coordinate of the query point.
        y: Y coordinate of the query point.
        z: Z coordinate of the query point.
        X2: 2D array of X coordinates on the target face, shape ``(Nu, Nv)``.
        Y2: 2D array of Y coordinates on the target face, shape ``(Nu, Nv)``.
        Z2: 2D array of Z coordinates on the target face, shape ``(Nu, Nv)``.
        tol: Maximum Euclidean distance for a valid match. Defaults to ``1e-6``.

    Returns:
        ``[p, q]`` indices into ``X2/Y2/Z2`` where the match was found, or
        ``[-1, -1]`` if no node is within tolerance.
    """
    dx = x - X2
    dy = y - Y2
    dz = z - Z2
    d = np.sqrt(dx * dx + dy * dy + dz * dz)
    val = np.amin(d)
    location = np.where(d == val)

    if val < tol:
        return [location[0][0], location[1][0]]

    return [-1, -1]
