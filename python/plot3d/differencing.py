"""Edge detection via finite differencing on structured grids.

Computes forward and backward differences at every grid node along each
computational direction. Useful for identifying parallel edges whose
vertices might intersect between two blocks.
"""
import numpy as np
import pandas as pd


def find_face_edges(X: np.ndarray, Y: np.ndarray, Z: np.ndarray):
    """Compute forward/backward difference vectors on a 2-D face grid.

    For each node ``(p, q)`` on a face, computes the backward and forward
    difference vectors along both the P and Q directions. This can be used
    to check if edges of two faces are parallel (a necessary condition for
    vertex intersection).

    Args:
        X: 2-D array of X coordinates with shape ``(PMAX, QMAX)``.
        Y: 2-D array of Y coordinates with shape ``(PMAX, QMAX)``.
        Z: 2-D array of Z coordinates with shape ``(PMAX, QMAX)``.

    Returns:
        pandas.DataFrame: One row per node with columns:

        - ``p``, ``q`` -- node indices.
        - ``dp`` -- tuple of backward and forward difference vectors
          ``((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))`` along the P direction.
        - ``dq`` -- same for the Q direction.

        Boundary differences are zero where no neighbor exists.
    """
    (PMAX, QMAX) = X.shape
    diffArray = list()

    for p in range(0, PMAX):
        for q in range(0, QMAX):
                dx_b = 0; dy_b = 0; dz_b = 0
                if p != 0:
                    dx_b = X[p-1, q] - X[p, q]
                    dy_b = Y[p-1, q] - Y[p, q]
                    dz_b = Z[p-1, q] - Z[p, q]

                dx_f = 0; dy_f = 0; dz_f = 0
                if p != PMAX-1:
                    dx_f = X[p+1, q] - X[p, q]
                    dy_f = Y[p+1, q] - Y[p, q]
                    dz_f = Z[p+1, q] - Z[p, q]

                dp = ((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))
                if q != 0:
                    dx_b = X[p, q-1] - X[p, q]
                    dy_b = Y[p, q-1] - Y[p, q]
                    dz_b = Z[p, q-1] - Z[p, q]

                if q != QMAX-1:
                    dx_f = X[p, q+1] - X[p, q]
                    dy_f = Y[p, q+1] - Y[p, q]
                    dz_f = Z[p, q+1] - Z[p, q]
                dq = ((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))

                diffArray.append({"p": p, "q": q, 'dp': dp, 'dq': dq})

    df = pd.DataFrame(data=diffArray)
    return df


def find_edges(X: np.ndarray, Y: np.ndarray, Z: np.ndarray):
    """Compute forward/backward difference vectors on a 3-D block grid.

    For each node ``(i, j, k)`` in a block, computes the backward and
    forward difference vectors along all three computational directions
    (I, J, K). This is the volumetric analog of :func:`find_face_edges`.

    Args:
        X: 3-D array of X coordinates with shape ``(IMAX, JMAX, KMAX)``.
        Y: 3-D array of Y coordinates with shape ``(IMAX, JMAX, KMAX)``.
        Z: 3-D array of Z coordinates with shape ``(IMAX, JMAX, KMAX)``.

    Returns:
        pandas.DataFrame: One row per node with columns:

        - ``i``, ``j``, ``k`` -- node indices.
        - ``di`` -- tuple of backward and forward difference vectors
          ``((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))`` along I.
        - ``dj`` -- same for J direction.
        - ``dk`` -- same for K direction.

        Boundary differences are zero where no neighbor exists.
    """
    (IMAX, JMAX, KMAX) = X.shape
    diffArray = list()

    for i in range(0, IMAX):
        for j in range(0, JMAX):
            for k in range(0, KMAX):
                dx_b = 0; dy_b = 0; dz_b = 0
                if i != 0:
                    dx_b = X[i-1, j, k] - X[i, j, k]
                    dy_b = Y[i-1, j, k] - Y[i, j, k]
                    dz_b = Z[i-1, j, k] - Z[i, j, k]

                dx_f = 0; dy_f = 0; dz_f = 0
                if i != IMAX-1:
                    dx_f = X[i, j, k] - X[i+1, j, k]
                    dy_f = Y[i, j, k] - Y[i+1, j, k]
                    dz_f = Z[i, j, k] - Z[i+1, j, k]

                di = ((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))
                if j != 0:
                    dx_b = X[i, j-1, k] - X[i, j, k]
                    dy_b = Y[i, j-1, k] - Y[i, j, k]
                    dz_b = Z[i, j-1, k] - Z[i, j, k]

                if j != JMAX-1:
                    dx_f = X[i, j, k] - X[i, j+1, k]
                    dy_f = Y[i, j, k] - Y[i, j+1, k]
                    dz_f = Z[i, j, k] - Z[i, j+1, k]
                dj = ((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))

                dx_f = 0; dy_f = 0; dz_f = 0
                if k != 0:
                    dx_b = X[i, j, k-1] - X[i, j, k]
                    dy_b = Y[i, j, k-1] - Y[i, j, k]
                    dz_b = Z[i, j, k-1] - Z[i, j, k]

                if k != KMAX-1:
                    dx_f = X[i, j, k] - X[i, j, k+1]
                    dy_f = Y[i, j, k] - Y[i, j, k+1]
                    dz_f = Z[i, j, k] - Z[i, j, k+1]
                dk = ((dx_b, dy_b, dz_b), (dx_f, dy_f, dz_f))

                diffArray.append({"i": i, "j": j, "k": k, 'di': di, 'dj': dj, 'dk': dk})
    df = pd.DataFrame(data=diffArray)
    return df
