"""FacePool: Theta-bucketed spatial index for rotational periodicity matching.

Ported from the Rust ``face_pool.rs`` implementation in plot3d-rs. Provides
O(log N) candidate lookup via binary search on theta-sorted face centroids,
with axial/radial overlap pre-checks to reject non-overlapping faces.

The key idea: for rotational periodicity, we only need to check faces whose
angular position (theta) is within ``rotation_angle`` of the query face.
By sorting faces by theta and using binary search, we reduce candidate
selection from O(N) linear scan to O(log N + k) where k is the number of
candidates in the angular window.
"""
from __future__ import annotations
from typing import List, Set, Tuple, Optional
import bisect
import math
import numpy as np
from .face import Face
from .facefunctions import _to_theta, _to_radius
from .utils import face_key as _utils_face_key


class CylindricalFaceInfo:
    """Cylindrical coordinate extents for a single face."""
    __slots__ = (
        'theta_centroid', 'axial_centroid', 'radial_centroid',
        'theta_min', 'theta_max',
        'axial_min', 'axial_max',
        'radial_min', 'radial_max',
    )

    def __init__(self, face: Face, block, rotation_axis: str):
        n = face.nvertex
        xs = face.x[:n]
        ys = face.y[:n]
        zs = face.z[:n]

        # Centroid in cylindrical coords
        cx, cy, cz = face.cx, face.cy, face.cz
        self.theta_centroid = float(_to_theta(cx, cy, cz, rotation_axis))
        self.axial_centroid = float(_axial_coord(cx, cy, cz, rotation_axis))
        self.radial_centroid = float(_to_radius(cx, cy, cz, rotation_axis))

        # Extent from all face grid points (not just corners)
        pts = face.grid_points(block)
        if pts.size > 0:
            thetas = _to_theta(pts[:, 0], pts[:, 1], pts[:, 2], rotation_axis)
            axials = _axial_array(pts, rotation_axis)
            radials = _to_radius(pts[:, 0], pts[:, 1], pts[:, 2], rotation_axis)

            self.theta_min = float(np.min(thetas))
            self.theta_max = float(np.max(thetas))
            self.axial_min = float(np.min(axials))
            self.axial_max = float(np.max(axials))
            self.radial_min = float(np.min(radials))
            self.radial_max = float(np.max(radials))
        else:
            # Fallback to vertex-only extents
            thetas = _to_theta(xs, ys, zs, rotation_axis)
            axials = _axial_array_verts(xs, ys, zs, rotation_axis)
            radials = _to_radius(xs, ys, zs, rotation_axis)

            self.theta_min = float(np.min(thetas))
            self.theta_max = float(np.max(thetas))
            self.axial_min = float(np.min(axials))
            self.axial_max = float(np.max(axials))
            self.radial_min = float(np.min(radials))
            self.radial_max = float(np.max(radials))


def _axial_coord(x, y, z, rotation_axis: str) -> float:
    """Extract the axial coordinate along the rotation axis."""
    if rotation_axis == 'x':
        return float(x)
    elif rotation_axis == 'y':
        return float(y)
    else:
        return float(z)


def _axial_array(pts: np.ndarray, rotation_axis: str) -> np.ndarray:
    """Extract axial coordinate array from (N,3) points."""
    if rotation_axis == 'x':
        return pts[:, 0]
    elif rotation_axis == 'y':
        return pts[:, 1]
    else:
        return pts[:, 2]


def _axial_array_verts(x, y, z, rotation_axis: str) -> np.ndarray:
    """Extract axial coordinate array from separate x,y,z arrays."""
    if rotation_axis == 'x':
        return np.asarray(x)
    elif rotation_axis == 'y':
        return np.asarray(y)
    else:
        return np.asarray(z)


class FacePool:
    """Theta-sorted face pool for efficient rotational candidate lookup.

    Faces are indexed by their cylindrical theta centroid. Finding candidates
    for a given target theta is done via binary search + linear scan within
    the angular tolerance window, with axial/radial overlap checks.
    """

    def __init__(self, faces: List[Face], blocks, rotation_axis: str):
        """Build the face pool from a list of faces.

        Args:
            faces: List of Face objects to index
            blocks: List of Block objects (indexed by face.blockIndex)
            rotation_axis: 'x', 'y', or 'z'
        """
        self.faces: List[Face] = list(faces)
        self.rotation_axis = rotation_axis
        self._blocks = blocks

        # Compute cylindrical info for each face
        self.cyl_info: List[CylindricalFaceInfo] = []
        for f in self.faces:
            self.cyl_info.append(CylindricalFaceInfo(f, blocks[f.blockIndex], rotation_axis))

        # Build theta-sorted index
        self._theta_sorted: List[int] = sorted(
            range(len(self.faces)),
            key=lambda i: self.cyl_info[i].theta_centroid
        )
        self._theta_keys: List[float] = [
            self.cyl_info[i].theta_centroid for i in self._theta_sorted
        ]

        # Track consumed faces
        self.consumed: Set[Tuple] = set()

    def _face_key(self, f: Face) -> Tuple:
        return _utils_face_key(f)

    def is_consumed(self, idx: int) -> bool:
        return self._face_key(self.faces[idx]) in self.consumed

    def consume(self, idx: int):
        self.consumed.add(self._face_key(self.faces[idx]))

    def consume_face(self, face: Face):
        self.consumed.add(self._face_key(face))

    def active_indices(self) -> List[int]:
        """Return indices of faces that haven't been consumed."""
        return [i for i in range(len(self.faces)) if not self.is_consumed(i)]

    def add_face(self, face: Face):
        """Add a new face to the pool (e.g., a split remnant)."""
        idx = len(self.faces)
        self.faces.append(face)
        info = CylindricalFaceInfo(face, self._blocks[face.blockIndex], self.rotation_axis)
        self.cyl_info.append(info)

        # Insert into theta-sorted index (maintain sorted order)
        pos = bisect.bisect_left(self._theta_keys, info.theta_centroid)
        self._theta_sorted.insert(pos, idx)
        self._theta_keys.insert(pos, info.theta_centroid)

    def find_candidates(
        self,
        target_theta: float,
        axial_range: Tuple[float, float],
        radial_range: Tuple[float, float],
        theta_tol: float,
    ) -> List[int]:
        """Find face indices within theta_tol of target_theta with axial/radial overlap.

        Uses binary search for O(log N) initial lookup, then linear scan
        within the angular window. Handles theta wrapping at +/- pi.

        Args:
            target_theta: Target angular position (radians)
            axial_range: (min, max) axial extent of the query face
            radial_range: (min, max) radial extent of the query face
            theta_tol: Angular tolerance (radians)

        Returns:
            List of face pool indices that are candidates for matching
        """
        results = []

        def _scan_range(theta_lo: float, theta_hi: float):
            """Scan the theta-sorted array for faces in [theta_lo, theta_hi]."""
            start = bisect.bisect_left(self._theta_keys, theta_lo)
            i = start
            n = len(self._theta_sorted)
            while i < n and self._theta_keys[i] <= theta_hi:
                pool_idx = self._theta_sorted[i]
                if not self.is_consumed(pool_idx):
                    info = self.cyl_info[pool_idx]
                    # Check axial overlap (10% margin)
                    axial_margin = 0.1 * max(abs(axial_range[1] - axial_range[0]), 1e-12)
                    if not (info.axial_max < axial_range[0] - axial_margin or
                            info.axial_min > axial_range[1] + axial_margin):
                        # Check radial overlap (10% margin)
                        radial_margin = 0.1 * max(abs(radial_range[1] - radial_range[0]), 1e-12)
                        if not (info.radial_max < radial_range[0] - radial_margin or
                                info.radial_min > radial_range[1] + radial_margin):
                            results.append(pool_idx)
                i += 1

        theta_lo = target_theta - theta_tol
        theta_hi = target_theta + theta_tol

        # Main scan
        _scan_range(theta_lo, theta_hi)

        # Handle theta wrapping at ±π
        if theta_lo < -math.pi:
            _scan_range(theta_lo + 2 * math.pi, math.pi)
        if theta_hi > math.pi:
            _scan_range(-math.pi, theta_hi - 2 * math.pi)

        return results

    def find_rotational_candidates(
        self,
        face_idx: int,
        rotation_angle: float,
        theta_tol: float,
    ) -> List[int]:
        """Find candidate faces that could match face_idx after rotation.

        Looks for faces at theta ± rotation_angle from the given face's
        theta centroid, with axial/radial overlap filtering.

        Args:
            face_idx: Index of the source face in this pool
            rotation_angle: Rotation angle in radians
            theta_tol: Angular search tolerance in radians

        Returns:
            Sorted, deduplicated list of candidate face indices
        """
        info = self.cyl_info[face_idx]
        axial_range = (info.axial_min, info.axial_max)
        radial_range = (info.radial_min, info.radial_max)

        # Search forward rotation
        target_fwd = info.theta_centroid + rotation_angle
        cands_fwd = self.find_candidates(target_fwd, axial_range, radial_range, theta_tol)

        # Search backward rotation
        target_rev = info.theta_centroid - rotation_angle
        cands_rev = self.find_candidates(target_rev, axial_range, radial_range, theta_tol)

        # Merge, deduplicate, sort
        all_cands = set(cands_fwd) | set(cands_rev)
        all_cands.discard(face_idx)  # Don't match against self
        return sorted(all_cands)
