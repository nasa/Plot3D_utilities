"""Shared utility functions used by connectivity and periodicity modules."""
import math
from typing import List, Tuple, Dict


def euclidean_distance(x1, y1, z1, x2, y2, z2) -> float:
    """Euclidean distance between two 3D points."""
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def face_key(f) -> Tuple:
    """Return a hashable key identifying a face by block index and IJK bounds.

    Works with any object having BlockIndex/blockIndex and IMIN/JMIN/KMIN/IMAX/JMAX/KMAX
    attributes (Face objects from both face.py and facefunctions.py).
    """
    bi = getattr(f, 'BlockIndex', getattr(f, 'blockIndex', 0))
    return (bi, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)


def face_grid_dims(imin, imax, jmin, jmax, kmin, kmax) -> List[int]:
    """Return sorted list of non-zero face extents along each axis.

    E.g. a K-constant face from (0,0,5) to (10,20,5) returns [10, 20].
    """
    dims = []
    if imin != imax: dims.append(imax - imin)
    if jmin != jmax: dims.append(jmax - jmin)
    if kmin != kmax: dims.append(kmax - kmin)
    return sorted(dims)


def divide_face_dict_indices(records: List[Dict], gcd: int,
                              keys=('IMIN', 'JMIN', 'KMIN', 'IMAX', 'JMAX', 'KMAX'),
                              nested_sides=None):
    """Divide face dictionary indices by a GCD factor (integer division).

    Args:
        records: List of dicts containing face indices to scale down.
        gcd: Factor to divide each index by.
        keys: Which keys to scale within each record.
        nested_sides: If provided (e.g. ['block1', 'block2']), scale keys
            inside each nested dict rather than at the top level.
    """
    for rec in records:
        if nested_sides:
            for side in nested_sides:
                for k in keys:
                    rec[side][k] = rec[side][k] // gcd
        else:
            for k in keys:
                rec[k] = rec[k] // gcd


def enumerate_unique_corners(I: list, J: list, K: list) -> List[Tuple[int, int, int]]:
    """Enumerate unique (i, j, k) corners from min/max index lists.

    Args:
        I: [IMIN, IMAX]
        J: [JMIN, JMAX]
        K: [KMIN, KMAX]

    Returns:
        List of unique (i, j, k) tuples (up to 8 corners, fewer when indices coincide).
    """
    unique_corners = []
    seen = set()
    for i in I:
        for j in J:
            for k in K:
                key = (i, j, k)
                if key not in seen:
                    seen.add(key)
                    unique_corners.append(key)
    return unique_corners


def scale_face_dict_indices(records: List[Dict], gcd: int,
                            keys=('IMIN', 'JMIN', 'KMIN', 'IMAX', 'JMAX', 'KMAX'),
                            nested_sides=None):
    """Scale face dictionary indices by a GCD factor.

    Args:
        records: List of dicts containing face indices to scale.
        gcd: Factor to multiply each index by.
        keys: Which keys to scale within each record.
        nested_sides: If provided (e.g. ['block1', 'block2']), scale keys
            inside each nested dict rather than at the top level.
    """
    for rec in records:
        if nested_sides:
            for side in nested_sides:
                for k in keys:
                    rec[side][k] = int(rec[side][k] * gcd)
        else:
            for k in keys:
                rec[k] = int(rec[k] * gcd)
