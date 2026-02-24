"""Block class and related utilities for Plot3D structured grids.

A :class:`Block` represents a single structured grid block in a multi-block
Plot3D mesh.  Each block stores 3-D coordinate arrays ``X``, ``Y``, ``Z`` with
shape ``(IMAX, JMAX, KMAX)`` and provides basic geometric operations (scale,
shift, cylindrical conversion, cell volume computation).

The module also provides :func:`compute_gcd` and :func:`reduce_blocks` for
GCD-based mesh reduction used during connectivity and periodicity matching.
"""
import numpy as np
import math
from tqdm import trange
from typing import List
import numpy.typing as npt


class Block:
    """A single structured grid block in a Plot3D multi-block mesh.

    Attributes:
        X: 3-D array of X coordinates with shape ``(IMAX, JMAX, KMAX)``.
        Y: 3-D array of Y coordinates with shape ``(IMAX, JMAX, KMAX)``.
        Z: 3-D array of Z coordinates with shape ``(IMAX, JMAX, KMAX)``.
        IMAX: Number of nodes in the I direction.
        JMAX: Number of nodes in the J direction.
        KMAX: Number of nodes in the K direction.
        cx: X coordinate of the block centroid.
        cy: Y coordinate of the block centroid.
        cz: Z coordinate of the block centroid.
    """
    X: npt.NDArray
    Y: npt.NDArray
    Z: npt.NDArray
    IMAX: int
    JMAX: int
    KMAX: int

    def __init__(self, X: npt.NDArray, Y: npt.NDArray, Z: npt.NDArray):
        """Initialize a block from its coordinate arrays.

        When the module-level flag ``plot3d.use_single_precision`` is set,
        coordinates are cast to ``float32``; otherwise they remain as provided
        (typically ``float64``).

        Args:
            X: 3-D array of X coordinates, shape ``(IMAX, JMAX, KMAX)``.
            Y: 3-D array of Y coordinates, shape ``(IMAX, JMAX, KMAX)``.
            Z: 3-D array of Z coordinates, shape ``(IMAX, JMAX, KMAX)``.
        """
        import plot3d as _p3d
        if getattr(_p3d, 'use_single_precision', False):
            X = np.asarray(X, dtype=np.float32)
            Y = np.asarray(Y, dtype=np.float32)
            Z = np.asarray(Z, dtype=np.float32)
        self.IMAX, self.JMAX, self.KMAX = X.shape
        self.X = X
        self.Y = Y
        self.Z = Z
        self.cx = np.mean(X)
        self.cy = np.mean(Y)
        self.cz = np.mean(Z)

    def __repr__(self):
        return f"({self.IMAX},{self.JMAX},{self.KMAX})"

    def copy(self) -> 'Block':
        """Return an independent deep copy of this block."""
        return Block(self.X.copy(), self.Y.copy(), self.Z.copy())

    def scale(self, factor: float):
        """Uniformly scale all coordinates by a multiplicative factor.

        Args:
            factor: Scale factor applied to X, Y, and Z.
        """
        self.X *= factor
        self.Y *= factor
        self.Z *= factor

    def shift(self, shift_amount: float, direction: str = "z"):
        """Translate the block along a single axis.

        Args:
            shift_amount: Distance to shift.
            direction: Axis to shift along (``'x'``, ``'y'``, or ``'z'``).
                Defaults to ``'z'``.
        """
        if direction.lower() == 'z':
            self.Z += shift_amount
        elif direction.lower() == 'y':
            self.Y += shift_amount
        elif direction.lower() == 'x':
            self.X += shift_amount

    def cylindrical(self):
        """Convert to cylindrical coordinates assuming X is the rotation axis.

        Sets ``self.r`` (radius) and ``self.theta`` (azimuthal angle in radians)
        computed from Y and Z.
        """
        self.r = np.sqrt(self.Z * self.Z + self.Y * self.Y)
        self.theta = np.arctan2(self.Y, self.Z)

    def cell_volumes(self):
        """Compute the volume of every hexahedral cell in the block.

        Uses the method of Davies and Salmond (AIAA J., vol. 23, No. 6,
        pp. 954-956, 1985) which is exact for hexahedra whose faces are
        bilinear surfaces.

        Returns:
            numpy.ndarray: Array of cell volumes with shape
            ``(IMAX, JMAX, KMAX)``.  Volumes at boundary indices ``[0,:,:]``,
            ``[:,0,:]``, ``[:,:,0]`` are zero.
        """
        X = self.X
        Y = self.Y
        Z = self.Z
        a = [np.zeros(shape=(self.IMAX, self.JMAX, self.KMAX)) for _ in range(9)]
        # face csi=const
        for k in range(1, self.KMAX):
            for j in range(1, self.JMAX):
                for i in range(self.IMAX):
                    dx1 = X[i, j, k-1] - X[i, j-1, k]
                    dy1 = Y[i, j, k-1] - Y[i, j-1, k]
                    dz1 = Z[i, j, k-1] - Z[i, j-1, k]

                    dx2 = X[i, j, k] - X[i, j-1, k-1]
                    dy2 = Y[i, j, k] - Y[i, j-1, k-1]
                    dz2 = Z[i, j, k] - Z[i, j-1, k-1]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[0][i, j, k] = ax*0.5
                    a[1][i, j, k] = ay*0.5
                    a[2][i, j, k] = az*0.5

        # face eta=const
        for k in range(1, self.KMAX):
            for j in range(self.JMAX):
                for i in range(1, self.IMAX):
                    dx1 = X[i, j, k] - X[i-1, j, k-1]
                    dy1 = Y[i, j, k] - Y[i-1, j, k-1]
                    dz1 = Z[i, j, k] - Z[i-1, j, k-1]

                    dx2 = X[i, j, k-1] - X[i-1, j, k]
                    dy2 = Y[i, j, k-1] - Y[i-1, j, k]
                    dz2 = Z[i, j, k-1] - Z[i-1, j, k]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[3][i, j, k] = ax*0.5
                    a[4][i, j, k] = ay*0.5
                    a[5][i, j, k] = az*0.5

        # face zit=const
        for k in range(self.KMAX):
            for j in range(1, self.JMAX):
                for i in range(1, self.IMAX):
                    dx1 = X[i, j, k] - X[i-1, j-1, k]
                    dy1 = Y[i, j, k] - Y[i-1, j-1, k]
                    dz1 = Z[i, j, k] - Z[i-1, j-1, k]

                    dx2 = X[i-1, j, k] - X[i, j-1, k]
                    dy2 = Y[i-1, j, k] - Y[i, j-1, k]
                    dz2 = Z[i-1, j, k] - Z[i, j-1, k]

                    ax = dy1*dz2 - dz1*dy2
                    ay = dz1*dx2 - dx1*dz2
                    az = dx1*dy2 - dy1*dx2

                    a[6][i, j, k] = ax*0.5
                    a[7][i, j, k] = ay*0.5
                    a[8][i, j, k] = az*0.5

        cf = np.zeros(shape=(6, 3))
        v = np.zeros(shape=(self.IMAX, self.JMAX, self.KMAX))

        for k in trange(1, self.KMAX, desc='Calculating the volumes'):
            for j in range(1, self.JMAX):
                for i in range(1, self.IMAX):
                    cf[0, 0] = X[i-1, j-1, k-1] + X[i-1, j-1, k] + X[i-1, j, k-1] + X[i-1, j, k]
                    cf[0, 1] = Y[i-1, j-1, k-1] + Y[i-1, j-1, k] + Y[i-1, j, k-1] + Y[i-1, j, k]
                    cf[0, 2] = Z[i-1, j-1, k-1] + Z[i-1, j-1, k] + Z[i-1, j, k-1] + Z[i-1, j, k]
                    cf[1, 0] = X[i, j-1, k-1] + X[i, j-1, k] + X[i, j, k-1] + X[i, j, k]
                    cf[1, 1] = Y[i, j-1, k-1] + Y[i, j-1, k] + Y[i, j, k-1] + Y[i, j, k]
                    cf[1, 2] = Z[i, j-1, k-1] + Z[i, j-1, k] + Z[i, j, k-1] + Z[i, j, k]
                    cf[2, 0] = X[i-1, j-1, k-1] + X[i-1, j-1, k] + X[i, j-1, k-1] + X[i, j-1, k]
                    cf[2, 1] = Y[i-1, j-1, k-1] + Y[i-1, j-1, k] + Y[i, j-1, k-1] + Y[i, j-1, k]
                    cf[2, 2] = Z[i-1, j-1, k-1] + Z[i-1, j-1, k] + Z[i, j-1, k-1] + Z[i, j-1, k]
                    cf[3, 0] = X[i-1, j, k-1] + X[i-1, j, k] + X[i, j, k-1] + X[i, j, k]
                    cf[3, 1] = Y[i-1, j, k-1] + Y[i-1, j, k] + Y[i, j, k-1] + Y[i, j, k]
                    cf[3, 2] = Z[i-1, j, k-1] + Z[i-1, j, k] + Z[i, j, k-1] + Z[i, j, k]
                    cf[4, 0] = X[i-1, j-1, k-1] + X[i-1, j, k-1] + X[i, j-1, k-1] + X[i, j, k-1]
                    cf[4, 1] = Y[i-1, j-1, k-1] + Y[i-1, j, k-1] + Y[i, j-1, k-1] + Y[i, j, k-1]
                    cf[4, 2] = Z[i-1, j-1, k-1] + Z[i-1, j, k-1] + Z[i, j-1, k-1] + Z[i, j, k-1]
                    cf[5, 0] = X[i-1, j-1, k] + X[i-1, j, k] + X[i, j-1, k] + X[i, j, k]
                    cf[5, 1] = Y[i-1, j-1, k] + Y[i-1, j, k] + Y[i, j-1, k] + Y[i, j, k]
                    cf[5, 2] = Z[i-1, j-1, k] + Z[i-1, j, k] + Z[i, j-1, k] + Z[i, j, k]

                    vol12 = 0
                    for n in range(0, 2):
                        for l in range(0, 3):
                            vol12 += math.pow(-1, n+1) * (
                                        +cf[n, l]*a[l][i-1+n, j, k]
                                        +cf[2+n, l]*a[3+l][i, j-1+n, k]
                                        +cf[4+n, l]*a[6+l][i, j, k-1+n])
                    v[i, j, k] = vol12/12
        return v

    def get_faces(self):
        """Return the six boundary face coordinate arrays.

        Returns:
            dict: Keys are ``'imin'``, ``'imax'``, ``'jmin'``, ``'jmax'``,
            ``'kmin'``, ``'kmax'``. Values are ``(X_face, Y_face, Z_face)``
            tuples where each array is a 2-D slice of the block.
        """
        return {
            'imin': (self.X[0, :, :], self.Y[0, :, :], self.Z[0, :, :]),
            'imax': (self.X[-1, :, :], self.Y[-1, :, :], self.Z[-1, :, :]),
            'jmin': (self.X[:, 0, :], self.Y[:, 0, :], self.Z[:, 0, :]),
            'jmax': (self.X[:, -1, :], self.Y[:, -1, :], self.Z[:, -1, :]),
            'kmin': (self.X[:, :, 0], self.Y[:, :, 0], self.Z[:, :, 0]),
            'kmax': (self.X[:, :, -1], self.Y[:, :, -1], self.Z[:, :, -1]),
        }

    @property
    def size(self) -> int:
        """Total number of nodes in the block (``IMAX * JMAX * KMAX``)."""
        return self.IMAX * self.JMAX * self.KMAX


def compute_gcd(blocks: List[Block]) -> int:
    """Compute the minimum GCD across all blocks' dimensions.

    For each block, computes ``gcd(IMAX-1, JMAX-1, KMAX-1)`` and returns
    the minimum. This shared factor is used by :func:`reduce_blocks` to
    down-sample every block uniformly before connectivity or periodicity
    matching, then scale results back up.

    Args:
        blocks: List of blocks to analyze.

    Returns:
        The minimum GCD factor (always >= 1).
    """
    gcd_array = [math.gcd(b.IMAX - 1, math.gcd(b.JMAX - 1, b.KMAX - 1)) for b in blocks]
    return min(gcd_array)


def reduce_blocks(blocks: List[Block], factor: int):
    """Down-sample blocks by subsampling every *factor* indices.

    Each block's X, Y, Z arrays are sliced with stride ``factor`` along
    all three axes, and ``IMAX``, ``JMAX``, ``KMAX`` are updated to reflect
    the new shape.

    .. warning::
        This modifies blocks **in-place**. Pass copies
        (``[b.copy() for b in blocks]``) if the originals must be preserved.

    Args:
        blocks: List of blocks to reduce. Modified in-place.
        factor: Subsampling stride along each axis.

    Returns:
        The same list of blocks, now at reduced resolution.
    """
    for i in range(len(blocks)):
        blocks[i].X = blocks[i].X[::factor, ::factor, ::factor]
        blocks[i].Y = blocks[i].Y[::factor, ::factor, ::factor]
        blocks[i].Z = blocks[i].Z[::factor, ::factor, ::factor]
        blocks[i].IMAX, blocks[i].JMAX, blocks[i].KMAX = blocks[i].X.shape
    return blocks
