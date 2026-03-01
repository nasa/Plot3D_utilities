"""Shared utility functions for the plot3d library.

Provides common helpers used across connectivity, periodicity, and verification
modules, including distance calculations, face key generation, index scaling,
and corner enumeration.
"""
import math
from typing import List, Tuple, Dict


def euclidean_distance(x1: float, y1: float, z1: float,
                       x2: float, y2: float, z2: float) -> float:
    """Compute the Euclidean distance between two 3D points.

    Args:
        x1: X coordinate of point 1.
        y1: Y coordinate of point 1.
        z1: Z coordinate of point 1.
        x2: X coordinate of point 2.
        y2: Y coordinate of point 2.
        z2: Z coordinate of point 2.

    Returns:
        The straight-line distance between the two points.
    """
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def face_key(f) -> Tuple:
    """Return a hashable key identifying a face by block index and IJK bounds.

    Works with any object that has ``BlockIndex`` (or ``blockIndex``) and
    ``IMIN``, ``JMIN``, ``KMIN``, ``IMAX``, ``JMAX``, ``KMAX`` attributes.

    Args:
        f: A Face-like object (from :class:`~plot3d.face.Face` or
            :func:`~plot3d.facefunctions.create_face_from_diagonals`).

    Returns:
        A 7-element tuple ``(block_index, IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)``.
    """
    bi = getattr(f, 'BlockIndex', getattr(f, 'blockIndex', 0))
    return (bi, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)


def face_grid_dims(imin: int, imax: int, jmin: int, jmax: int,
                   kmin: int, kmax: int) -> List[int]:
    """Return the sorted list of non-zero face extents along each axis.

    For a face defined by its IJK bounds, computes the extent
    (max - min) along each axis where the face is not constant,
    and returns them sorted smallest-first. Useful for checking whether
    two faces have compatible grid dimensions (allowing transposition).

    Args:
        imin: Minimum I index.
        imax: Maximum I index.
        jmin: Minimum J index.
        jmax: Maximum J index.
        kmin: Minimum K index.
        kmax: Maximum K index.

    Returns:
        Sorted list of non-zero extents. For example, a K-constant face
        from ``(0, 0, 5)`` to ``(10, 20, 5)`` returns ``[10, 20]``.
    """
    dims = []
    if imin != imax: dims.append(imax - imin)
    if jmin != jmax: dims.append(jmax - jmin)
    if kmin != kmax: dims.append(kmax - kmin)
    return sorted(dims)


def divide_face_dict_indices(records: List[Dict], gcd: int,
                              nested_sides=None):
    """Divide face dictionary lb/ub tuple indices by a GCD factor using integer division.

    This is the inverse of :func:`scale_face_dict_indices`. Used when
    reducing a mesh to GCD resolution before matching, to scale face
    indices down accordingly.

    Args:
        records: List of dicts containing ``'lb'`` and ``'ub'`` tuple indices
            to scale down.
        gcd: Divisor applied to each index via ``//``.
        nested_sides: If provided (e.g. ``['block1', 'block2']``), divides
            indices inside each nested sub-dict rather than at the top level.
    """
    for rec in records:
        if nested_sides:
            for side in nested_sides:
                rec[side]['lb'] = tuple(v // gcd for v in rec[side]['lb'])
                rec[side]['ub'] = tuple(v // gcd for v in rec[side]['ub'])
        else:
            rec['lb'] = tuple(v // gcd for v in rec['lb'])
            rec['ub'] = tuple(v // gcd for v in rec['ub'])


def enumerate_unique_corners(I: list, J: list, K: list) -> List[Tuple[int, int, int]]:
    """Enumerate unique ``(i, j, k)`` corners from min/max index lists.

    Given the ``[min, max]`` ranges along each axis, generates all
    combinations and deduplicates them. When an axis has ``min == max``,
    fewer than 8 corners are produced.

    Args:
        I: ``[IMIN, IMAX]`` index range.
        J: ``[JMIN, JMAX]`` index range.
        K: ``[KMIN, KMAX]`` index range.

    Returns:
        List of unique ``(i, j, k)`` tuples (up to 8 corners).
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
                            nested_sides=None):
    """Multiply face dictionary lb/ub tuple indices by a GCD factor.

    Used to scale indices back to the original grid resolution after
    performing matching at reduced (GCD) resolution.

    Args:
        records: List of dicts containing ``'lb'`` and ``'ub'`` tuple indices
            to scale up.
        gcd: Factor to multiply each index by.
        nested_sides: If provided (e.g. ``['block1', 'block2']``), scales
            indices inside each nested sub-dict rather than at the top level.
    """
    for rec in records:
        if nested_sides:
            for side in nested_sides:
                rec[side]['lb'] = tuple(int(v * gcd) for v in rec[side]['lb'])
                rec[side]['ub'] = tuple(int(v * gcd) for v in rec[side]['ub'])
        else:
            rec['lb'] = tuple(int(v * gcd) for v in rec['lb'])
            rec['ub'] = tuple(int(v * gcd) for v in rec['ub'])
