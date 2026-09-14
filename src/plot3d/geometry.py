# geometry.py
"""Distance-based point-coincidence primitives.

Mirrors plot3d-rs ``geometry::PointGrid3``/``PointGrid2`` (commit
0e1b1a1): earlier coincidence checks in this package bucketed points via
``round(x / tol)`` and compared bucket-integer equality, which can miss
two points that are genuinely closer than ``tol`` but happen to straddle
an integer bin boundary. ``coincidence_count`` replaces that with a real
Euclidean-distance nearest-neighbor query.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.spatial import cKDTree


def coincidence_count(P1: npt.NDArray, P2: npt.NDArray, tol: float) -> int:
    """Count of points in P1 that have an actual-Euclidean-distance match in P2
    within `tol`. Mirrors plot3d-rs's geometry::PointGrid3/PointGrid2 (commit
    0e1b1a1): a real distance test, not `round(x/tol)` bin-equality — two points
    closer than `tol` that straddle a bin boundary are no longer missed.
    """
    P1 = np.asarray(P1, dtype=float)
    P2 = np.asarray(P2, dtype=float)
    if P1.size == 0 or P2.size == 0:
        return 0

    tree = cKDTree(P2)
    dist, _ = tree.query(P1, k=1, distance_upper_bound=tol)
    # query() returns inf for points with no neighbor within
    # distance_upper_bound, so a plain <= correctly excludes those.
    return int(np.count_nonzero(dist <= tol))
