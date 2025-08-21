# face.py
from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np
import numpy.typing as npt
import math
from math import acos, degrees


class Face:
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

    """Defines a Face of a block (e.g., IMIN,JMIN,KMIN to IMAX,JMIN,KMIN)."""

    def __init__(self, nvertex: int = 4):
        """Define a face using nvertex (4 = quad; 3 = triangle)."""
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

    def to_dict(self) -> Dict[str, int]:
        """Return a dictionary representation of this face."""
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
        # BUGFIX: previously returned (cx, cy, cx)
        return np.array([self.cx, self.cy, self.cz], dtype=np.float64)

    @property
    def IMIN(self) -> int:
        return int(self.I.min())

    @property
    def JMIN(self) -> int:
        return int(self.J.min())

    @property
    def KMIN(self) -> int:
        return int(self.K.min())

    @property
    def IMAX(self) -> int:
        return int(self.I.max())

    @property
    def JMAX(self) -> int:
        return int(self.J.max())

    @property
    def KMAX(self) -> int:
        return int(self.K.max())

    @property
    def BlockIndex(self) -> int:
        return int(self.blockIndex)

    @property
    def const_type(self) -> int:
        """0 => I-constant, 1 => J-constant, 2 => K-constant, -1 => none."""
        if self.IMIN == self.IMAX:
            return 0
        elif self.JMIN == self.JMAX:
            return 1
        elif self.KMIN == self.KMAX:
            return 2
        return -1

    @property
    def isEdge(self) -> bool:
        """True if the face degenerates to an edge (two indices constant)."""
        return (
            int(self.IMIN == self.IMAX)
            + int(self.JMIN == self.JMAX)
            + int(self.KMIN == self.KMAX)
        ) > 1

    @property
    def isPoint(self) -> bool:
        """True if the face degenerates to a point (all three indices constant)."""
        return (
            int(self.IMIN == self.IMAX)
            + int(self.JMIN == self.JMAX)
            + int(self.KMIN == self.KMAX)
        ) > 2

    def get_val(self, i_val: int, j_val: int, k_val: int) -> Tuple[float, float, float]:
        """Get (x,y,z) where (I,J,K) == (i_val, j_val, k_val)."""
        mask = (self.I == i_val) & (self.J == j_val) & (self.K == k_val)
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            raise KeyError("Face does not contain the requested (I,J,K) vertex.")
        d = int(idx[0])
        return float(self.x[d]), float(self.y[d]), float(self.z[d])

    def add_vertex(self, x: float, y: float, z: float, i: int, j: int, k: int):
        """Add a vertex to define the face."""
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
        """Return parametric face 'size' in index-space cells."""
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
        self.blockIndex = int(val)

    def set_face_id(self, val: int):
        self.id = int(val)

    # ------------------------------------------------------------------
    # Geometry (legacy normal)
    # ------------------------------------------------------------------

    def normal(self, block) -> npt.NDArray[np.float64]:
        """Compute the (unnormalized) geometric normal using three points."""
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
        """Find vertex index matches between faces by coordinates (with tol)."""
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
        matchedIndices = self.match_indices(f)
        return len(matchedIndices) == self.nvertex

    def __ne__(self, f: object) -> bool:
        return not self.__eq__(f)  # type: ignore

    def index_equals(self, f2: "Face") -> bool:
        return (
            self.IMIN == f2.IMIN
            and self.JMIN == f2.JMIN
            and self.KMIN == f2.KMIN
            and self.IMAX == f2.IMAX
            and self.JMAX == f2.JMAX
            and self.KMAX == f2.KMAX
        )

    def __hash__(self) -> int:
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
        return str(self)

    @property
    def diagonal_length(self) -> float:
        """Diagonal length using (IMIN,JMIN,KMIN) and (IMAX,JMAX,KMAX)."""
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
        """Get the extreme corners (IMIN,JMIN,KMIN) and (IMAX,JMAX,KMAX)."""
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
        """Legacy: centroid proximity only (kept for backward compatibility)."""
        val = math.sqrt(
            (self.cx - f.cx) ** 2 + (self.cy - f.cy) ** 2 + (self.cz - f.cz) ** 2
        )
        return val < tol

    def shift(self, dx: float, dy: float, dz: float):
        """Shift the face vertices (in-place)."""
        self.x += dx
        self.y += dy
        self.z += dz  # BUGFIX: previously wrote self.dz

    # ------------------------------------------------------------------
    # Robust contact / overlap predicates (polygon & node sharing)
    # ------------------------------------------------------------------

    # ---- internal helpers (polygon) ----
    def _corner_points(self) -> np.ndarray:
        """Return four corner points (4,3) in a consistent order."""
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
        n = np.cross(p1 - p0, p2 - p0)
        ln = np.linalg.norm(n)
        return n / ln if ln > 0 else np.array([0.0, 0.0, 1.0])

    def _quad_normal(self) -> np.ndarray:
        """Robust average normal from face corners (handles skew quads)."""
        q = self._corner_points()
        n1 = self._unit_normal(q[0], q[1], q[2])
        n2 = self._unit_normal(q[0], q[2], q[3])
        n = n1 + n2
        ln = np.linalg.norm(n)
        return (n / ln) if ln > 1e-12 else n1

    @staticmethod
    def _plane_distance(pts: np.ndarray, p0: np.ndarray, n: np.ndarray) -> np.ndarray:
        return (pts - p0) @ n

    @staticmethod
    def _dominant_projection_axis(n: np.ndarray) -> int:
        return int(np.argmax(np.abs(n)))

    @staticmethod
    def _project_drop_axis(pts3: np.ndarray, drop_axis: int) -> np.ndarray:
        keep = [0, 1, 2]
        keep.remove(drop_axis)
        return pts3[:, keep]  # (N,2)

    @staticmethod
    def _poly_area_2d(poly2d: np.ndarray) -> float:
        x, y = poly2d[:, 0], poly2d[:, 1]
        return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    @staticmethod
    def _clip_suth_hodg(subject: np.ndarray, clipper: np.ndarray) -> np.ndarray:
        """Sutherland–Hodgman polygon clip (convex clipper). Always returns (N,2)."""
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
        """
        Return (intersection area) / min(area_self, area_f) after projecting to the
        dominant plane. Accept faces whose normals are parallel OR anti-parallel.
        Plane-distance tolerance scales with face size to handle large meshes.
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
        """Surface contact if nearly coplanar (parallel/anti-parallel normals)
        and projected overlap area fraction ≥ min_overlap_frac."""
        frac = self.overlap_fraction(f, tol_angle_deg=tol_angle_deg, tol_plane_dist=tol_plane_dist)
        return frac >= min_overlap_frac

    # ---- node-sharing helpers (structured faces) ----
    def _index_ranges(self) -> Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]:
        """Inclusive [min,max] index ranges for I,J,K."""
        return (self.IMIN, self.IMAX), (self.JMIN, self.JMAX), (self.KMIN, self.KMAX)

    def grid_points(self, block, stride_u: int = 1, stride_v: int = 1) -> np.ndarray:
        """
        Return all XYZ points on this face from the parent Block, sampled on the
        structured face grid. Shape: (N, 3).
        Works only for structured faces (const_type in {0,1,2}).
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
        """Quantize 3D points to an integer grid at spacing 'tol' for robust equality."""
        s = tol if tol > 0 else 1e-12
        return np.round(P / s).astype(np.int64)

    @staticmethod
    def _row_view(a: np.ndarray) -> np.ndarray:
        """Create a 1D view of rows to use with np.intersect1d."""
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
        """
        Fraction of shared nodes = (# coincident XYZ points within tol) / min(N_self, N_other),
        where points are sampled from each face's structured grid with the given strides.
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
        """
        True if faces share at least 'min_shared_abs' coincident nodes AND at least
        'min_shared_frac' of the smaller face's sampled nodes.

        Detects partial surface contact when the blocks truly share grid nodes.
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
