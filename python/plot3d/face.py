"""Face representation for Plot3D structured grid blocks.

This module defines the :class:`Face` class, which represents a single face
(quadrilateral or triangular surface) of a structured multi-block grid in the
Plot3D format.  A face is described by four (or three) corner vertices with
spatial coordinates ``(x, y, z)`` and structured indices ``(I, J, K)``.

The class provides:

* vertex-level geometry operations (centroid, normals, diagonal length),
* index-range queries (``IMIN`` ... ``KMAX``),
* equality / hashing based on block index and index ranges,
* polygon-based overlap detection (Sutherland-Hodgman clipping), and
* node-sharing predicates for structured connectivity analysis.
"""
from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np
import numpy.typing as npt
import math
from math import acos, degrees


class Face:
    """A face of a Plot3D structured grid block.

    A face is bounded by four corner vertices (quad) or three (triangle)
    and carries both physical coordinates ``(x, y, z)`` and structured
    indices ``(I, J, K)`` for each vertex.  The face may be I-constant,
    J-constant, or K-constant depending on which index is the same across
    all vertices.

    Attributes:
        x: X-coordinates of the corner vertices.
        y: Y-coordinates of the corner vertices.
        z: Z-coordinates of the corner vertices.
        I: Structured I-indices of the corner vertices.
        J: Structured J-indices of the corner vertices.
        K: Structured K-indices of the corner vertices.
        cx: X-component of the face centroid.
        cy: Y-component of the face centroid.
        cz: Z-component of the face centroid.
        nvertex: Number of vertices added so far (up to 4).
        blockIndex: Zero-based index of the parent block.
        id: Integer identifier for this face.
    """

    x: npt.NDArray
    y: npt.NDArray
    z: npt.NDArray
    I: npt.NDArray
    J: npt.NDArray
    K: npt.NDArray
    cx: float = 0.0
    cy: float = 0.0
    cz: float = 0.0
    nvertex: int = 0
    blockIndex: int = 0  # not really needed except in periodicity
    id: int = 0

    def __init__(self, nvertex: int = 4):
        """Initialise a face with pre-allocated vertex storage.

        Args:
            nvertex: Maximum number of vertices (4 for a quad, 3 for a
                triangle).  Storage is always allocated for 4 vertices
                regardless of this value.
        """
        self.x = np.zeros(4, dtype=float)
        self.y = np.zeros(4, dtype=float)
        self.z = np.zeros(4, dtype=float)
        self.I = np.zeros(4, dtype=np.int64)
        self.J = np.zeros(4, dtype=np.int64)
        self.K = np.zeros(4, dtype=np.int64)
        self.nvertex = 0
        self.cx = 0.0
        self.cy = 0.0
        self.cz = 0.0
        self.blockIndex = 0
        self.id = 0

    # ------------------------------------------------------------------
    # Basic utilities
    # ------------------------------------------------------------------

    def copy(self) -> 'Face':
        """Create an independent copy of this face.

        Uses ``numpy.ndarray.copy`` for the coordinate and index arrays,
        which is significantly faster than ``copy.deepcopy``.

        Returns:
            A new :class:`Face` instance with identical attribute values.
        """
        f = Face.__new__(Face)
        f.x = self.x.copy()
        f.y = self.y.copy()
        f.z = self.z.copy()
        f.I = self.I.copy()
        f.J = self.J.copy()
        f.K = self.K.copy()
        f.cx = self.cx
        f.cy = self.cy
        f.cz = self.cz
        f.nvertex = self.nvertex
        f.blockIndex = self.blockIndex
        f.id = self.id
        return f

    def to_dict(self) -> Dict[str, int]:
        """Return a dictionary representation of this face.

        The dictionary keys are ``IMIN``, ``JMIN``, ``KMIN``, ``IMAX``,
        ``JMAX``, ``KMAX``, ``id``, and ``block_index``.

        Returns:
            A dictionary mapping string keys to integer values that
            fully describe the face's index range and identity.
        """
        return {
            "IMIN": int(self.I.min()),
            "JMIN": int(self.J.min()),
            "KMIN": int(self.K.min()),
            "IMAX": int(self.I.max()),
            "JMAX": int(self.J.max()),
            "KMAX": int(self.K.max()),
            "id": int(self.id),
            "block_index": int(self.blockIndex),
        }

    @property
    def centroid(self) -> npt.NDArray[np.float64]:
        """Centroid of the face as a length-3 array ``[cx, cy, cz]``.

        Returns:
            A NumPy array of shape ``(3,)`` with dtype ``float64``.
        """
        # BUGFIX: previously returned (cx, cy, cx)
        return np.array([self.cx, self.cy, self.cz], dtype=np.float64)

    @property
    def IMIN(self) -> int:
        """Minimum I-index across all vertices.

        Returns:
            The smallest ``I`` value stored in this face.
        """
        return int(self.I.min())

    @property
    def JMIN(self) -> int:
        """Minimum J-index across all vertices.

        Returns:
            The smallest ``J`` value stored in this face.
        """
        return int(self.J.min())

    @property
    def KMIN(self) -> int:
        """Minimum K-index across all vertices.

        Returns:
            The smallest ``K`` value stored in this face.
        """
        return int(self.K.min())

    @property
    def IMAX(self) -> int:
        """Maximum I-index across all vertices.

        Returns:
            The largest ``I`` value stored in this face.
        """
        return int(self.I.max())

    @property
    def JMAX(self) -> int:
        """Maximum J-index across all vertices.

        Returns:
            The largest ``J`` value stored in this face.
        """
        return int(self.J.max())

    @property
    def KMAX(self) -> int:
        """Maximum K-index across all vertices.

        Returns:
            The largest ``K`` value stored in this face.
        """
        return int(self.K.max())

    @property
    def BlockIndex(self) -> int:
        """Zero-based block index of the parent block.

        Returns:
            The block index as a Python ``int``.
        """
        return int(self.blockIndex)

    @property
    def const_type(self) -> int:
        """Identify which structured index is constant across the face.

        Returns:
            ``0`` if I-constant (``IMIN == IMAX``), ``1`` if J-constant
            (``JMIN == JMAX``), ``2`` if K-constant (``KMIN == KMAX``),
            or ``-1`` if no single index is constant.
        """
        if self.IMIN == self.IMAX:
            return 0
        elif self.JMIN == self.JMAX:
            return 1
        elif self.KMIN == self.KMAX:
            return 2
        return -1

    @property
    def isEdge(self) -> bool:
        """Check whether the face degenerates to an edge.

        A face is an edge when two or more of its three index ranges
        collapse (i.e., min equals max).

        Returns:
            ``True`` if the face degenerates to an edge, ``False``
            otherwise.
        """
        return (
            int(self.IMIN == self.IMAX)
            + int(self.JMIN == self.JMAX)
            + int(self.KMIN == self.KMAX)
        ) > 1

    @property
    def isPoint(self) -> bool:
        """Check whether the face degenerates to a single point.

        A face is a point when all three index ranges collapse.

        Returns:
            ``True`` if the face degenerates to a point, ``False``
            otherwise.
        """
        return (
            int(self.IMIN == self.IMAX)
            + int(self.JMIN == self.JMAX)
            + int(self.KMIN == self.KMAX)
        ) > 2

    def get_val(self, i_val: int, j_val: int, k_val: int) -> Tuple[float, float, float]:
        """Retrieve the spatial coordinates for a specific vertex.

        Searches the stored vertices for one whose ``(I, J, K)`` indices
        match the requested values exactly.

        Args:
            i_val: I-index of the desired vertex.
            j_val: J-index of the desired vertex.
            k_val: K-index of the desired vertex.

        Returns:
            A tuple ``(x, y, z)`` of the matching vertex coordinates.

        Raises:
            KeyError: If no vertex with the given ``(I, J, K)`` exists in
                this face.
        """
        mask = (self.I == i_val) & (self.J == j_val) & (self.K == k_val)
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            raise KeyError("Face does not contain the requested (I,J,K) vertex.")
        d = int(idx[0])
        return float(self.x[d]), float(self.y[d]), float(self.z[d])

    def add_vertex(self, x: float, y: float, z: float, i: int, j: int, k: int):
        """Append a vertex to this face.

        Vertices are stored in insertion order.  Once the fourth vertex
        has been added the centroid ``(cx, cy, cz)`` is computed
        automatically as the mean of the four vertex positions.

        Args:
            x: X-coordinate of the vertex.
            y: Y-coordinate of the vertex.
            z: Z-coordinate of the vertex.
            i: Structured I-index of the vertex.
            j: Structured J-index of the vertex.
            k: Structured K-index of the vertex.
        """
        self.x[self.nvertex] = x
        self.y[self.nvertex] = y
        self.z[self.nvertex] = z
        self.I[self.nvertex] = i
        self.J[self.nvertex] = j
        self.K[self.nvertex] = k
        self.nvertex += 1
        if self.nvertex == 4:
            self.cx = float(self.x.mean())
            self.cy = float(self.y.mean())
            self.cz = float(self.z.mean())

    @property
    def size(self) -> int:
        """Parametric face size in index-space cells.

        Computed as the product of the two non-constant index spans.
        For a fully 3-D range (no constant index) the product of all
        three spans is returned.

        Returns:
            The number of index-space cells covered by this face.
        """
        if self.IMIN == self.IMAX:
            return (self.JMAX - self.JMIN) * (self.KMAX - self.KMIN)
        elif self.JMIN == self.JMAX:
            return (self.IMAX - self.IMIN) * (self.KMAX - self.KMIN)
        elif self.KMIN == self.KMAX:
            return (self.IMAX - self.IMIN) * (self.JMAX - self.JMIN)
        else:
            return (
                (self.IMAX - self.IMIN)
                * (self.JMAX - self.JMIN)
                * (self.KMAX - self.KMIN)
            )

    def set_block_index(self, val: int):
        """Set the parent block index for this face.

        Args:
            val: Zero-based block index.
        """
        self.blockIndex = int(val)

    def set_face_id(self, val: int):
        """Set the integer identifier for this face.

        Args:
            val: Face identifier.
        """
        self.id = int(val)

    # ------------------------------------------------------------------
    # Geometry (legacy normal)
    # ------------------------------------------------------------------

    def normal(self, block) -> npt.NDArray[np.float64]:
        """Compute the unnormalised geometric normal of this face.

        Uses three corner points from the parent ``block`` to form two
        edge vectors and returns their cross product.  The resulting
        vector is **not** unit-length.

        Args:
            block: The parent :class:`Block` object whose ``X``, ``Y``,
                ``Z`` arrays are indexed by ``(I, J, K)``.

        Returns:
            A NumPy array of shape ``(3,)`` representing the face normal
            vector (unnormalised).
        """
        if self.const_type == 0:  # I-constant: IMIN == IMAX
            p1 = np.array(
                [
                    block.X[self.IMIN, self.JMIN, self.KMIN],
                    block.Y[self.IMIN, self.JMIN, self.KMIN],
                    block.Z[self.IMIN, self.JMIN, self.KMIN],
                ]
            )
            p2 = np.array(
                [
                    block.X[self.IMIN, self.JMAX, self.KMIN],
                    block.Y[self.IMIN, self.JMAX, self.KMIN],
                    block.Z[self.IMIN, self.JMAX, self.KMIN],
                ]
            )
            p3 = np.array(
                [
                    block.X[self.IMIN, self.JMIN, self.KMAX],
                    block.Y[self.IMIN, self.JMIN, self.KMAX],
                    block.Z[self.IMIN, self.JMIN, self.KMAX],
                ]
            )
        elif self.const_type == 1:  # J-constant: JMIN == JMAX
            p1 = np.array(
                [
                    block.X[self.IMIN, self.JMIN, self.KMIN],
                    block.Y[self.IMIN, self.JMIN, self.KMIN],
                    block.Z[self.IMIN, self.JMIN, self.KMIN],
                ]
            )
            p2 = np.array(
                [
                    block.X[self.IMIN, self.JMIN, self.KMAX],
                    block.Y[self.IMIN, self.JMIN, self.KMAX],
                    block.Z[self.IMIN, self.JMIN, self.KMAX],
                ]
            )
            p3 = np.array(
                [
                    block.X[self.IMAX, self.JMIN, self.KMIN],
                    block.Y[self.IMAX, self.JMIN, self.KMIN],
                    block.Z[self.IMAX, self.JMIN, self.KMIN],
                ]
            )
        else:  # K-constant: KMIN == KMAX
            p1 = np.array(
                [
                    block.X[self.IMIN, self.JMIN, self.KMIN],
                    block.Y[self.IMIN, self.JMIN, self.KMIN],
                    block.Z[self.IMIN, self.JMIN, self.KMIN],
                ]
            )
            p2 = np.array(
                [
                    block.X[self.IMAX, self.JMIN, self.KMIN],
                    block.Y[self.IMAX, self.JMIN, self.KMIN],
                    block.Z[self.IMAX, self.JMIN, self.KMIN],
                ]
            )
            p3 = np.array(
                [
                    block.X[self.IMIN, self.JMAX, self.KMIN],
                    block.Y[self.IMIN, self.JMAX, self.KMIN],
                    block.Z[self.IMIN, self.JMAX, self.KMIN],
                ]
            )
        u = p2 - p1
        v = p3 - p1
        return np.cross(u, v)

    # ------------------------------------------------------------------
    # Equality / hashing
    # ------------------------------------------------------------------

    def match_indices(self, f: "Face") -> List[List[int]]:
        """Find matching vertex pairs between this face and another.

        Vertices are matched by spatial proximity using a tolerance of
        ``1e-6``.  Each vertex in ``self`` is matched to at most one
        vertex in *f*, and vice versa.

        Args:
            f: The other :class:`Face` to compare against.

        Returns:
            A list of ``[i, j]`` pairs where ``i`` is the vertex index
            in ``self`` and ``j`` is the matching vertex index in *f*.
        """
        matched_vertices: List[int] = []
        tol = 1e-6
        matchedIndices: List[List[int]] = []
        for i in range(self.nvertex):
            for j in range(f.nvertex):
                dx = abs(self.x[i] - f.x[j])
                dy = abs(self.y[i] - f.y[j])
                dz = abs(self.z[i] - f.z[j])
                if dx < tol and dy < tol and dz < tol and (j not in matched_vertices):
                    matchedIndices.append([i, j])
                    matched_vertices.append(j)
                    break
        return matchedIndices

    def __eq__(self, f: object) -> bool:
        """Test equality by block index and all six index-range bounds.

        Args:
            f: Another object to compare.  Must be a :class:`Face`.

        Returns:
            ``True`` if *f* is a :class:`Face` with identical
            ``BlockIndex``, ``IMIN``...``KMAX``.
        """
        if not isinstance(f, Face):
            return False
        return (
            (self.BlockIndex == f.BlockIndex)
            and (self.IMIN == f.IMIN)
            and (self.IMAX == f.IMAX)
            and (self.JMIN == f.JMIN)
            and (self.JMAX == f.JMAX)
            and (self.KMIN == f.KMIN)
            and (self.KMAX == f.KMAX)
        )

    def vertices_equals(self, f: "Face") -> bool:
        """Test whether all vertices of this face match another by coordinates.

        Uses :meth:`match_indices` with the default tolerance to check
        that every vertex in ``self`` has a spatial match in *f*.

        Args:
            f: The other :class:`Face` to compare against.

        Returns:
            ``True`` if the number of matched vertex pairs equals
            ``nvertex``.
        """
        matchedIndices = self.match_indices(f)
        return len(matchedIndices) == self.nvertex

    def __ne__(self, f: object) -> bool:
        """Test inequality (inverse of :meth:`__eq__`).

        Args:
            f: Another object to compare.

        Returns:
            ``True`` if ``self != f``.
        """
        return not self.__eq__(f)  # type: ignore

    def index_equals(self, f2: "Face") -> bool:
        """Test whether the six index-range bounds match exactly.

        Unlike :meth:`__eq__`, this does **not** compare
        ``blockIndex``.

        Args:
            f2: The other :class:`Face` to compare against.

        Returns:
            ``True`` if ``IMIN``...``KMAX`` are identical.
        """
        return (
            self.IMIN == f2.IMIN
            and self.JMIN == f2.JMIN
            and self.KMIN == f2.KMIN
            and self.IMAX == f2.IMAX
            and self.JMAX == f2.JMAX
            and self.KMAX == f2.KMAX
        )

    def __hash__(self) -> int:
        """Return a hash based on the first and last vertex indices.

        Returns:
            An integer hash value derived from the ``(I, J, K)``
            indices of the first and last stored vertices.
        """
        if len(self.I) > 0:
            return hash(
                (
                    int(self.I[0]),
                    int(self.J[0]),
                    int(self.K[0]),
                    int(self.I[-1]),
                    int(self.J[-1]),
                    int(self.K[-1]),
                )
            )
        else:
            return hash((0, 0, 0, 0, 0, 0))

    def __str__(self) -> str:
        """Return a human-readable string of block index and ranges.

        Returns:
            A string in the form
            ``blk: <blockIndex> [IMIN,JMIN,KMIN,IMAX,JMAX,KMAX]``.
        """
        if len(self.I) > 0:
            return "blk: {:d} [{:d},{:d},{:d},{:d},{:d},{:d}]".format(
                self.blockIndex,
                self.IMIN,
                self.JMIN,
                self.KMIN,
                self.IMAX,
                self.JMAX,
                self.KMAX,
            )
        else:
            return "blk: {:d} [{:d},{:d},{:d},{:d},{:d},{:d}]".format(
                self.blockIndex, 0, 0, 0, 0, 0, 0
            )

    def __repr__(self) -> str:
        """Return the same string as :meth:`__str__`.

        Returns:
            A human-readable representation identical to ``str(self)``.
        """
        return str(self)

    @property
    def diagonal_length(self) -> float:
        """Euclidean distance between the two extreme corners.

        Computes the distance between the vertex at
        ``(IMIN, JMIN, KMIN)`` and the vertex at ``(IMAX, JMAX, KMAX)``
        in physical ``(x, y, z)`` space.

        Returns:
            The diagonal length as a ``float``.
        """
        minIndx = 0
        maxIndx = 0
        for indx in range(len(self.I)):
            if (
                self.I[indx] == self.IMIN
                and self.J[indx] == self.JMIN
                and self.K[indx] == self.KMIN
            ):
                minIndx = indx
            if (
                self.I[indx] == self.IMAX
                and self.J[indx] == self.JMAX
                and self.K[indx] == self.KMAX
            ):
                maxIndx = indx
        dx = float(self.x[minIndx] - self.x[maxIndx])
        dy = float(self.y[minIndx] - self.y[maxIndx])
        dz = float(self.z[minIndx] - self.z[maxIndx])
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def get_corners(self) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Get the two extreme-corner coordinates of this face.

        Returns:
            A 2-tuple of ``(x, y, z)`` tuples.  The first element
            corresponds to the vertex at ``(IMIN, JMIN, KMIN)`` and the
            second to ``(IMAX, JMAX, KMAX)``.
        """
        minIndx = 0
        maxIndx = 0
        for indx in range(len(self.I)):
            if (
                self.I[indx] == self.IMIN
                and self.J[indx] == self.JMIN
                and self.K[indx] == self.KMIN
            ):
                minIndx = indx
            if (
                self.I[indx] == self.IMAX
                and self.J[indx] == self.JMAX
                and self.K[indx] == self.KMAX
            ):
                maxIndx = indx
        return (
            (float(self.x[minIndx]), float(self.y[minIndx]), float(self.z[minIndx])),
            (float(self.x[maxIndx]), float(self.y[maxIndx]), float(self.z[maxIndx])),
        )

    def is_connected(self, f: "Face", tol: float = 1e-8) -> bool:
        """Test connectivity by centroid proximity (legacy method).

        This is a simple heuristic retained for backward compatibility.
        For robust contact detection, prefer :meth:`touches` or
        :meth:`touches_by_nodes`.

        Args:
            f: The other :class:`Face` to test against.
            tol: Maximum Euclidean distance between centroids for the
                faces to be considered connected.

        Returns:
            ``True`` if the centroid-to-centroid distance is less than
            *tol*.
        """
        val = math.sqrt(
            (self.cx - f.cx) ** 2 + (self.cy - f.cy) ** 2 + (self.cz - f.cz) ** 2
        )
        return val < tol

    def shift(self, dx: float, dy: float, dz: float):
        """Translate the face vertices in-place.

        Args:
            dx: Translation along the X-axis.
            dy: Translation along the Y-axis.
            dz: Translation along the Z-axis.
        """
        self.x += dx
        self.y += dy
        self.z += dz  # BUGFIX: previously wrote self.dz

    # ------------------------------------------------------------------
    # Robust contact / overlap predicates (polygon & node sharing)
    # ------------------------------------------------------------------

    # ---- internal helpers (polygon) ----
    def _corner_points(self) -> np.ndarray:
        """Return the four corner points in a consistent winding order.

        The ordering is ``[p00, p10, p11, p01]`` where ``p00`` is the
        ``(IMIN, JMIN, KMIN)`` corner and ``p11`` is the
        ``(IMAX, JMAX, KMAX)`` corner.  The two remaining corners are
        resolved by index lookup with fallback heuristics.

        Returns:
            A NumPy array of shape ``(4, 3)`` containing the corner
            coordinates.
        """
        idx_map = {
            (int(self.I[p]), int(self.J[p]), int(self.K[p])): p for p in range(self.nvertex)
        }

        p00 = idx_map.get((self.IMIN, self.JMIN, self.KMIN), 0)
        p11 = idx_map.get((self.IMAX, self.JMAX, self.KMAX), -1)

        p10 = idx_map.get((self.IMAX, self.JMIN, self.KMIN), None)
        if p10 is None:
            p10 = idx_map.get((self.IMAX, self.JMIN, self.KMAX), None)
        if p10 is None:
            p10 = (set(range(self.nvertex)) - {p00, p11}).pop()

        p01 = idx_map.get((self.IMIN, self.JMAX, self.KMIN), None)
        if p01 is None:
            p01 = idx_map.get((self.IMIN, self.JMAX, self.KMAX), None)
        if p01 is None:
            p01 = (set(range(self.nvertex)) - {p00, p10, p11}).pop()

        P = np.array(
            [
                [self.x[p00], self.y[p00], self.z[p00]],
                [self.x[p10], self.y[p10], self.z[p10]],
                [self.x[p11], self.y[p11], self.z[p11]],
                [self.x[p01], self.y[p01], self.z[p01]],
            ],
            dtype=float,
        )
        return P

    @staticmethod
    def _unit_normal(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
        """Compute the unit normal of the triangle ``(p0, p1, p2)``.

        Args:
            p0: First vertex, shape ``(3,)``.
            p1: Second vertex, shape ``(3,)``.
            p2: Third vertex, shape ``(3,)``.

        Returns:
            A unit-length normal vector of shape ``(3,)``.  Falls back
            to ``[0, 0, 1]`` if the triangle is degenerate.
        """
        n = np.cross(p1 - p0, p2 - p0)
        ln = np.linalg.norm(n)
        return n / ln if ln > 0 else np.array([0.0, 0.0, 1.0])

    def _quad_normal(self) -> np.ndarray:
        """Compute a robust average normal from the face corners.

        Splits the quad into two triangles and averages their normals,
        which handles skew (non-planar) quads correctly.

        Returns:
            A unit-length normal vector of shape ``(3,)``.
        """
        q = self._corner_points()
        n1 = self._unit_normal(q[0], q[1], q[2])
        n2 = self._unit_normal(q[0], q[2], q[3])
        n = n1 + n2
        ln = np.linalg.norm(n)
        return (n / ln) if ln > 1e-12 else n1

    @staticmethod
    def _plane_distance(pts: np.ndarray, p0: np.ndarray, n: np.ndarray) -> np.ndarray:
        """Compute signed distances from points to a plane.

        The plane is defined by a point *p0* and normal *n*.

        Args:
            pts: Points to test, shape ``(N, 3)``.
            p0: A point on the plane, shape ``(3,)``.
            n: The plane normal, shape ``(3,)``.

        Returns:
            A 1-D array of signed distances, shape ``(N,)``.
        """
        return (pts - p0) @ n

    @staticmethod
    def _dominant_projection_axis(n: np.ndarray) -> int:
        """Return the axis index with the largest absolute normal component.

        This axis is dropped when projecting 3-D polygons to 2-D for
        intersection calculations.

        Args:
            n: Normal vector, shape ``(3,)``.

        Returns:
            ``0`` (X), ``1`` (Y), or ``2`` (Z).
        """
        return int(np.argmax(np.abs(n)))

    @staticmethod
    def _project_drop_axis(pts3: np.ndarray, drop_axis: int) -> np.ndarray:
        """Project 3-D points to 2-D by dropping one coordinate axis.

        Args:
            pts3: Points in 3-D, shape ``(N, 3)``.
            drop_axis: The axis index to remove (``0``, ``1``, or ``2``).

        Returns:
            A 2-D array of shape ``(N, 2)``.
        """
        keep = [0, 1, 2]
        keep.remove(drop_axis)
        return pts3[:, keep]  # (N,2)

    @staticmethod
    def _poly_area_2d(poly2d: np.ndarray) -> float:
        """Compute the signed area of a 2-D polygon using the shoelace formula.

        Args:
            poly2d: Vertices of the polygon, shape ``(N, 2)``, in order.

        Returns:
            The signed area (positive for counter-clockwise winding).
        """
        x, y = poly2d[:, 0], poly2d[:, 1]
        return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    @staticmethod
    def _clip_suth_hodg(subject: np.ndarray, clipper: np.ndarray) -> np.ndarray:
        """Clip a 2-D polygon against a convex clipper using Sutherland-Hodgman.

        Args:
            subject: Vertices of the polygon to be clipped, shape
                ``(N, 2)``.
            clipper: Vertices of the convex clipping polygon, shape
                ``(M, 2)``.

        Returns:
            A NumPy array of shape ``(K, 2)`` containing the vertices of
            the clipped polygon.  Returns an empty ``(0, 2)`` array if
            the intersection is empty.
        """
        def inside(p, a, b):
            return (b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]) >= 0.0
        def intersect(a, b, c, d):
            x1,y1 = a; x2,y2 = b; x3,y3 = c; x4,y4 = d
            den = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
            if abs(den) < 1e-15:
                return b
            px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)) / den
            py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)) / den
            return np.array([px, py], dtype=float)

        out = np.asarray(subject, dtype=float).reshape(-1, 2).tolist()
        if len(out) == 0:
            return np.zeros((0, 2), dtype=float)

        C = np.asarray(clipper, dtype=float).reshape(-1, 2)
        for i in range(len(C)):
            input_list = out
            out = []
            A = C[i]; B = C[(i+1) % len(C)]
            if len(input_list) == 0:
                break
            S = input_list[-1]
            for E in input_list:
                if inside(E, A, B):
                    if not inside(S, A, B):
                        out.append(intersect(S, E, A, B))
                    out.append(E)
                elif inside(S, A, B):
                    out.append(intersect(S, E, A, B))
                S = E
        return np.asarray(out, dtype=float).reshape(-1, 2)

    # ---- polygon-overlap API (relaxed/adaptive) ----
    def overlap_fraction(
        self,
        f: "Face",
        tol_angle_deg: float = 10.0,   # slightly looser
        tol_plane_dist: float = 1e-6,  # absolute floor; also scaled by size
    ) -> float:
        """Compute the fractional area overlap between two faces.

        Projects both faces onto the dominant 2-D plane, clips the
        polygons using Sutherland-Hodgman, and returns the ratio of the
        intersection area to the smaller face area.

        Faces whose normals differ by more than *tol_angle_deg* (parallel
        **or** anti-parallel) are treated as non-overlapping.  A
        plane-distance check scaled by the characteristic face length is
        also applied.

        Args:
            f: The other :class:`Face` to test overlap against.
            tol_angle_deg: Maximum angle (in degrees) between the face
                normals for the faces to be considered coplanar.
            tol_plane_dist: Absolute minimum tolerance for the
                plane-distance check.  The effective tolerance is
                ``max(tol_plane_dist, 1e-6 * max(1.0, Lc))`` where
                ``Lc`` is the average characteristic edge length.

        Returns:
            A ``float`` in ``[0.0, 1.0]`` representing the overlap
            fraction.  Returns ``0.0`` when the faces are not coplanar or
            do not intersect.
        """
        Q1 = self._corner_points()
        Q2 = f._corner_points()

        # orientation gate: accept either same or opposite orientation
        n1 = self._quad_normal()
        n2 = f._quad_normal()
        cosang_abs = float(np.clip(abs(np.dot(n1, n2)), -1.0, 1.0))
        ang = degrees(math.acos(cosang_abs))
        if ang > tol_angle_deg:
            return 0.0

        # characteristic length for scaling tolerances
        def char_len(Q):
            e = np.linalg.norm(np.roll(Q, -1, axis=0) - Q, axis=1)
            return float(e.mean())
        Lc = 0.5 * (char_len(Q1) + char_len(Q2))
        tol_plane = max(tol_plane_dist, 1e-6 * max(1.0, Lc))

        c1 = Q1.mean(axis=0); c2 = Q2.mean(axis=0)
        d1 = abs(float(self._plane_distance(c1[None, :], Q2[0], n2)))
        d2 = abs(float(self._plane_distance(c2[None, :], Q1[0], n1)))
        if d1 > tol_plane or d2 > tol_plane:
            return 0.0

        # quick 3D AABB prefilter
        pad = 1e-6 * max(1.0, Lc)
        mn1, mx1 = Q1.min(axis=0), Q1.max(axis=0)
        mn2, mx2 = Q2.min(axis=0), Q2.max(axis=0)
        if (mx1[0] + pad < mn2[0] - pad or mx2[0] + pad < mn1[0] - pad or
            mx1[1] + pad < mn2[1] - pad or mx2[1] + pad < mn1[1] - pad or
            mx1[2] + pad < mn2[2] - pad or mx2[2] + pad < mn1[2] - pad):
            return 0.0

        # project to dominant plane using average normal for stability
        navg = n1 + n2
        if np.linalg.norm(navg) < 1e-12:
            navg = n1
        drop = self._dominant_projection_axis(navg)
        P1 = self._project_drop_axis(Q1, drop)
        P2 = self._project_drop_axis(Q2, drop)

        A1 = self._poly_area_2d(P1);  A2 = self._poly_area_2d(P2)
        if A1 < 0: P1 = P1[::-1]; A1 = -A1
        if A2 < 0: P2 = P2[::-1]; A2 = -A2
        if min(A1, A2) == 0.0:
            return 0.0

        Pint = self._clip_suth_hodg(P1, P2)
        if Pint.size == 0:
            return 0.0
        Aint = abs(self._poly_area_2d(Pint))
        return float(Aint / min(A1, A2))

    def touches(
        self,
        f: "Face",
        tol_angle_deg: float = 10.0,
        tol_plane_dist: float = 1e-6,
        min_overlap_frac: float = 0.01,  # permissive; tune up once it works
    ) -> bool:
        """Test whether two faces are in surface contact.

        Two faces touch if they are nearly coplanar (normals parallel or
        anti-parallel within *tol_angle_deg*) and their projected area
        overlap fraction meets the *min_overlap_frac* threshold.

        Args:
            f: The other :class:`Face` to test against.
            tol_angle_deg: Maximum angular deviation (degrees) between
                normals.
            tol_plane_dist: Plane-distance tolerance passed to
                :meth:`overlap_fraction`.
            min_overlap_frac: Minimum overlap fraction (in ``[0, 1]``)
                for the faces to be considered touching.

        Returns:
            ``True`` if the faces are in surface contact, ``False``
            otherwise.
        """
        frac = self.overlap_fraction(f, tol_angle_deg=tol_angle_deg, tol_plane_dist=tol_plane_dist)
        return frac >= min_overlap_frac

    # ---- node-sharing helpers (structured faces) ----
    def _index_ranges(self) -> Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]:
        """Return the inclusive ``[min, max]`` index ranges for I, J, and K.

        Returns:
            A 3-tuple of ``(min, max)`` pairs for I, J, and K
            respectively.
        """
        return (self.IMIN, self.IMAX), (self.JMIN, self.JMAX), (self.KMIN, self.KMAX)

    def grid_points(self, block, stride_u: int = 1, stride_v: int = 1) -> np.ndarray:
        """Extract all XYZ grid points on this face from a parent block.

        Samples the structured face grid at the given strides in the
        two parametric directions.  Only works for faces with a single
        constant index (``const_type`` in ``{0, 1, 2}``); for
        unstructured faces, falls back to the stored corner vertices.

        Args:
            block: The parent :class:`Block` whose ``X``, ``Y``, ``Z``
                arrays are indexed by ``(I, J, K)``.
            stride_u: Sampling stride in the first parametric direction.
            stride_v: Sampling stride in the second parametric direction.

        Returns:
            A NumPy array of shape ``(N, 3)`` containing the sampled
            ``(x, y, z)`` coordinates.
        """
        if self.const_type == -1:
            # not a structured face; fall back to stored vertices
            return np.stack(
                [self.x[:self.nvertex], self.y[:self.nvertex], self.z[:self.nvertex]],
                axis=1,
            )

        (i0, i1), (j0, j1), (k0, k1) = self._index_ranges()
        su = max(1, int(stride_u))
        sv = max(1, int(stride_v))

        if self.const_type == 0:
            # I-constant: vary (j,k)
            I = self.IMIN
            js = np.arange(j0, j1 + 1, su, dtype=int)
            ks = np.arange(k0, k1 + 1, sv, dtype=int)
            JJ, KK = np.meshgrid(js, ks, indexing="ij")
            X = block.X[I, JJ, KK]
            Y = block.Y[I, JJ, KK]
            Z = block.Z[I, JJ, KK]

        elif self.const_type == 1:
            # J-constant: vary (i,k)
            J = self.JMIN
            is_ = np.arange(i0, i1 + 1, su, dtype=int)
            ks = np.arange(k0, k1 + 1, sv, dtype=int)
            II, KK = np.meshgrid(is_, ks, indexing="ij")
            X = block.X[II, J, KK]
            Y = block.Y[II, J, KK]
            Z = block.Z[II, J, KK]

        else:
            # K-constant: vary (i,j)
            K = self.KMIN
            is_ = np.arange(i0, i1 + 1, su, dtype=int)
            js = np.arange(j0, j1 + 1, sv, dtype=int)
            II, JJ = np.meshgrid(is_, js, indexing="ij")
            X = block.X[II, JJ, K]
            Y = block.Y[II, JJ, K]
            Z = block.Z[II, JJ, K]

        P = np.stack([X, Y, Z], axis=-1).reshape(-1, 3).astype(float)
        return P

    @staticmethod
    def _quantize_points(P: np.ndarray, tol: float) -> np.ndarray:
        """Snap 3-D points to an integer lattice for robust equality testing.

        Each coordinate is divided by *tol* and rounded to the nearest
        integer, so two points whose true distance is less than *tol*
        will map to the same lattice point.

        Args:
            P: Points to quantize, shape ``(N, 3)``.
            tol: Grid spacing for quantization.

        Returns:
            An ``int64`` array of shape ``(N, 3)`` with quantized
            coordinates.
        """
        s = tol if tol > 0 else 1e-12
        return np.round(P / s).astype(np.int64)

    @staticmethod
    def _row_view(a: np.ndarray) -> np.ndarray:
        """Create a structured 1-D view of rows for use with ``np.intersect1d``.

        Each row of the input array is treated as a single structured
        element, enabling fast set operations on multi-column arrays.

        Args:
            a: A 2-D array of shape ``(N, M)``.

        Returns:
            A 1-D structured array of length ``N`` suitable for
            ``np.intersect1d``.
        """
        if not a.flags["C_CONTIGUOUS"]:
            a = np.ascontiguousarray(a)
        return a.view([("", a.dtype)] * a.shape[1])

    def shared_point_fraction(
        self,
        other: "Face",
        block_self,
        block_other,
        tol_xyz: float = 1e-8,
        stride_u: int = 1,
        stride_v: int = 1,
    ) -> float:
        """Compute the fraction of coincident grid nodes between two faces.

        Samples both faces from their respective parent blocks, quantizes
        the coordinates, and counts how many nodes coincide within
        *tol_xyz*.

        Args:
            other: The other :class:`Face` to compare with.
            block_self: Parent block for ``self``.
            block_other: Parent block for *other*.
            tol_xyz: Spatial tolerance for node coincidence.
            stride_u: Sampling stride in the first parametric direction.
            stride_v: Sampling stride in the second parametric direction.

        Returns:
            ``(# coincident nodes) / min(N_self, N_other)`` as a
            ``float`` in ``[0.0, 1.0]``.  Returns ``0.0`` if either face
            has no points.
        """
        P1 = self.grid_points(block_self, stride_u=stride_u, stride_v=stride_v)
        P2 = other.grid_points(block_other, stride_u=stride_u, stride_v=stride_v)
        if P1.size == 0 or P2.size == 0:
            return 0.0

        Q1 = self._quantize_points(P1, tol_xyz)
        Q2 = self._quantize_points(P2, tol_xyz)

        v1 = self._row_view(Q1)
        v2 = self._row_view(Q2)
        inter = np.intersect1d(v1, v2, assume_unique=False)
        shared = int(inter.size)
        denom = min(Q1.shape[0], Q2.shape[0])
        return (shared / denom) if denom > 0 else 0.0

    def touches_by_nodes(
        self,
        other: "Face",
        block_self,
        block_other,
        tol_xyz: float = 1e-8,
        min_shared_frac: float = 0.02,
        min_shared_abs: int = 4,
        stride_u: int = 1,
        stride_v: int = 1,
    ) -> bool:
        """Test surface contact by counting shared structured grid nodes.

        Two faces are considered in contact when they share at least
        *min_shared_abs* coincident nodes **and** at least
        *min_shared_frac* of the smaller face's sampled nodes.  This
        method is more robust than polygon overlap for detecting partial
        contact between blocks that truly share grid nodes.

        Args:
            other: The other :class:`Face` to test against.
            block_self: Parent block for ``self``.
            block_other: Parent block for *other*.
            tol_xyz: Spatial tolerance for node coincidence.
            min_shared_frac: Minimum fraction of the smaller face's nodes
                that must coincide for contact to be detected.
            min_shared_abs: Minimum absolute number of coincident nodes
                required.
            stride_u: Sampling stride in the first parametric direction.
            stride_v: Sampling stride in the second parametric direction.

        Returns:
            ``True`` if both the absolute and fractional thresholds are
            met, ``False`` otherwise.
        """
        P1 = self.grid_points(block_self, stride_u=stride_u, stride_v=stride_v)
        P2 = other.grid_points(block_other, stride_u=stride_u, stride_v=stride_v)
        if P1.size == 0 or P2.size == 0:
            return False

        Q1 = self._quantize_points(P1, tol_xyz)
        Q2 = self._quantize_points(P2, tol_xyz)

        v1 = self._row_view(Q1)
        v2 = self._row_view(Q2)
        inter = np.intersect1d(v1, v2, assume_unique=False)
        shared = int(inter.size)

        denom = min(Q1.shape[0], Q2.shape[0])
        frac = (shared / denom) if denom > 0 else 0.0
        return (shared >= min_shared_abs) and (frac >= min_shared_frac)
