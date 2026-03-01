"""Face manipulation utilities for Plot3D structured grids.

This module provides functions for comparing, matching, splitting, and
searching block faces in multi-block structured grids. It includes
routines for identifying outer (boundary) faces, locating bounding faces
along translational or angular directions, and converting between
dictionary and ``Face`` object representations.
"""
from typing import Dict, List, Optional, Tuple
from .listfunctions import unique_pairs
from .block import Block, compute_gcd, reduce_blocks
from .face import Face
import numpy.typing as npt
import numpy as np
import math

def faces_match(face1: Tuple[npt.NDArray, npt.NDArray, npt.NDArray],face2: Tuple[npt.NDArray, npt.NDArray, npt.NDArray],tol: float = 1e-12) -> Tuple[bool, Optional[Tuple[bool, bool]]]:
    """Compare two block faces and determine whether they match geometrically.

    Each face is represented as a tuple of ``(X, Y, Z)`` 2-D arrays. The
    function tests all four combinations of up/down and left/right flips
    on ``face2`` and checks whether the corner coordinates agree with
    ``face1`` within the specified tolerance.

    Args:
        face1: Tuple of ``(X, Y, Z)`` arrays for the first face.
        face2: Tuple of ``(X, Y, Z)`` arrays for the second face.
        tol: Absolute tolerance for comparing corner coordinates.
            Defaults to ``1e-12``.

    Returns:
        A tuple ``(matched, flip_flags)`` where *matched* is ``True`` if
        the faces match and *flip_flags* is a ``(flip_ud, flip_lr)`` tuple
        indicating the flips applied to ``face2``. Returns
        ``(False, None)`` when no orientation produces a match.
    """
    def get_corners(X, Y, Z):
        return np.array([
            [X[0, 0], Y[0, 0], Z[0, 0]],
            [X[0, -1], Y[0, -1], Z[0, -1]],
            [X[-1, 0], Y[-1, 0], Z[-1, 0]],
            [X[-1, -1], Y[-1, -1], Z[-1, -1]],
        ])

    X1, Y1, Z1 = face1
    X2, Y2, Z2 = face2

    if X1.shape != X2.shape:
        return False, None

    corners1 = get_corners(X1, Y1, Z1)
    for flip_ud in [False, True]:
        for flip_lr in [False, True]:
            X2f, Y2f, Z2f = X2.copy(), Y2.copy(), Z2.copy()
            if flip_ud:
                X2f, Y2f, Z2f = np.flip(X2f, axis=0), np.flip(Y2f, axis=0), np.flip(Z2f, axis=0)
            if flip_lr:
                X2f, Y2f, Z2f = np.flip(X2f, axis=1), np.flip(Y2f, axis=1), np.flip(Z2f, axis=1)

            corners2 = get_corners(X2f, Y2f, Z2f)
            diffs = np.linalg.norm(corners1 - corners2, axis=1)
            if np.all(diffs <= tol):
                return True, (flip_ud, flip_lr)

    return False, None

def find_matching_faces(block1, block2, tol=1e-8):
    """Find a shared face between two blocks.

    Iterates over all six faces of each block and returns the first pair
    that matches geometrically (see ``faces_match``).

    Args:
        block1: First ``Block`` to compare.
        block2: Second ``Block`` to compare.
        tol: Absolute tolerance for the face comparison. Defaults to
            ``1e-8``.

    Returns:
        A tuple ``(face1_name, face2_name, flip_flags)`` for the first
        matching pair found. Returns ``(None, None, None)`` if no faces
        match.
    """
    faces1 = block1.get_faces()
    faces2 = block2.get_faces()
    for face1_name, face1_data in faces1.items():
        for face2_name, face2_data in faces2.items():
            match, flip_flags = faces_match(face1_data, face2_data, tol=tol)
            if match:
                return face1_name, face2_name, flip_flags
    return None, None, None


def get_outer_faces(block1:Block):
    """Identify the outer (boundary) faces of a single block.

    Builds the six bounding faces from the block's corner vertices and
    determines which faces are unique (outer) versus which pairs share
    identical vertex positions (interior/collapsed).

    Args:
        block1: A Plot3D ``Block`` whose outer faces are to be found.

    Returns:
        A tuple ``(non_matching, matching)`` where *non_matching* is a
        list of ``Face`` objects that have no duplicate within the block
        (guaranteed exterior faces) and *matching* is a list of
        ``(Face, Face)`` tuples representing interior face pairs that
        share the same vertex positions.
    """
    I = [0,block1.IMAX-1]               # Python index starts at 0, need to subtract 1 for it to get the i,j,k
    J = [0,block1.JMAX-1]
    K = [0,block1.KMAX-1]
    # Create the outer faces
    faces = list()
    face = Face(4)
    i=I[0]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    i=I[1]
    for j in J:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    j=J[0]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    j=J[1]
    for i in I:
        for k in K:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    k=K[0]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)

    faces.append(face)
    face = Face(4)
    k=K[1]
    for i in I:
        for j in J:
            face.add_vertex(block1.X[i,j,k], block1.Y[i,j,k], block1.Z[i,j,k],i,j,k)
    faces.append(face)

    # Check if faces match each other
    matching = list()
    non_matching = list()
    for i in range(len(faces)):
        matchFound = False
        for j in range(len(faces)):
            if (i!=j and faces[i].vertices_equals(faces[j])):
                matching.append((i,j))
                matchFound = True
        if not matchFound:
            non_matching.append(faces[i]) # these are guaranteed to be exterior
    matching = list(unique_pairs(matching))
    matching = [(faces[i],faces[j]) for i,j in matching]

    # Make sure normals do not intersect
    # block_center_to_face_center =  block1.cx
    return non_matching, matching # these should be the outer faces

def create_face_from_diagonals(block:Block,imin:int,jmin:int,kmin:int,imax:int,jmax:int,kmax:int) -> Face:
    """Create a ``Face`` on a block from two diagonal corner indices.

    Exactly one of the index pairs ``(imin, imax)``, ``(jmin, jmax)``, or
    ``(kmin, kmax)`` must be equal, which determines the orientation of the
    resulting 2-D face. The four vertices are generated by iterating over
    the two free index dimensions.

    Args:
        block: The ``Block`` on which to create the face.
        imin: Lower corner I index.
        jmin: Lower corner J index.
        kmin: Lower corner K index.
        imax: Upper corner I index.
        jmax: Upper corner J index.
        kmax: Upper corner K index.

    Returns:
        A ``Face`` object with four vertices populated from the block's
        coordinate arrays.
    """
    newFace = Face(4)           # This is because two of the corners either imin or imax can be equal
    if imin==imax:
        i = imin
        for j in [jmin,jmax]:
            for k in [kmin,kmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    elif jmin==jmax:
        j = jmin
        for i in [imin,imax]:
            for k in [kmin,kmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    elif kmin==kmax:
        k = kmin
        for i in [imin,imax]:
            for j in [jmin,jmax]:
                x = block.X[i,j,k]
                y = block.Y[i,j,k]
                z = block.Z[i,j,k]
                newFace.add_vertex(x,y,z,i,j,k)
    return newFace


AxisMap = {"x": 0, "y": 1, "z": 2}

def _face_axis_extreme(face: Face, axis: str) -> Tuple[float, float]:
    """Return the min and max coordinate of a face along a given axis.

    Args:
        face: The ``Face`` whose vertex coordinates are inspected.
        axis: Axis label, one of ``'x'``, ``'y'``, or ``'z'``.

    Returns:
        A tuple ``(vmin, vmax)`` of the minimum and maximum coordinate
        values along the specified axis.
    """
    if axis == "x":
        arr = face.x[:face.nvertex]
    elif axis == "y":
        arr = face.y[:face.nvertex]
    else:
        arr = face.z[:face.nvertex]
    return float(np.min(arr)), float(np.max(arr))

def _global_axis_extreme(blocks: List[Block], axis: str) -> Tuple[float, float]:
    """Compute the global coordinate range along an axis across all blocks.

    Args:
        blocks: List of ``Block`` objects to scan.
        axis: Axis label, one of ``'x'``, ``'y'``, or ``'z'``.

    Returns:
        A tuple ``(global_min, global_max)`` of the extreme coordinate
        values across every block.
    """
    vals = []
    idx = AxisMap[axis]
    for b in blocks:
        if idx == 0: vals.append(b.X.reshape(-1))
        elif idx == 1: vals.append(b.Y.reshape(-1))
        else: vals.append(b.Z.reshape(-1))
    cat = np.concatenate(vals, axis=0)
    return float(cat.min()), float(cat.max())

def _select_seed_faces(outer_faces: List[Face], blocks: List[Block],
                       axis: str, side: str, tol_abs: float) -> List[Face]:
    """Select outer faces whose extreme coordinate matches the global extreme.

    These seed faces serve as starting points for the BFS boundary
    collection in ``_bfs_collect_boundary``.

    Args:
        outer_faces: Candidate outer ``Face`` objects.
        blocks: List of all ``Block`` objects (used to compute global
            extremes).
        axis: Axis label, one of ``'x'``, ``'y'``, or ``'z'``.
        side: Which extreme to target: ``'min'`` or ``'max'``.
        tol_abs: Absolute tolerance for comparing the face extreme to
            the global extreme.

    Returns:
        A list of ``Face`` objects whose relevant extreme is within
        ``tol_abs`` of the global extreme along the given axis and side.
    """
    gmin, gmax = _global_axis_extreme(blocks, axis)
    target = gmin if side == "min" else gmax
    seeds: List[Face] = []
    for f in outer_faces:
        fmin, fmax = _face_axis_extreme(f, axis)
        face_ext = fmin if side == "min" else fmax
        if abs(face_ext - target) <= tol_abs:
            seeds.append(f)
    return seeds

def _bfs_collect_boundary(seed_faces: List[Face],
                          all_outer_faces: List[Face],
                          blocks: List[Block],
                          axis: str,
                          side: str,
                          tol_abs: float,
                          node_tol_xyz: float,
                          min_shared_abs: int = 2,
                          min_shared_frac: float = 0.005) -> List[Face]:
    """Flood-fill over outer faces to collect a connected boundary patch.

    Starting from the seed faces, this performs a breadth-first search
    through outer faces that lie on the same extreme plane and share
    grid nodes with at least one already-visited face.

    Args:
        seed_faces: Initial ``Face`` objects to start the flood fill.
        all_outer_faces: Complete pool of outer ``Face`` objects to
            search through.
        blocks: List of all ``Block`` objects (used for node coordinate
            lookups and global extreme computation).
        axis: Axis label, one of ``'x'``, ``'y'``, or ``'z'``.
        side: Which extreme plane to collect: ``'min'`` or ``'max'``.
        tol_abs: Absolute tolerance for determining whether a face lies
            on the extreme plane.
        node_tol_xyz: Absolute tolerance for node-level coordinate
            matching when testing face adjacency.
        min_shared_abs: Minimum number of shared nodes required for two
            faces to be considered neighbors. Defaults to ``2``.
        min_shared_frac: Minimum fraction of shared nodes (relative to
            the smaller face) required for adjacency. Defaults to
            ``0.005``.

    Returns:
        A list of ``Face`` objects forming the connected boundary patch
        reachable from the seed faces.
    """
    # index faces by (block, IMIN..KMAX) for visited bookkeeping
    def key(f: Face) -> Tuple[int,int,int,int,int,int,int]:
        return (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)

    # all faces that lie on the extreme plane
    on_plane = []
    for f in all_outer_faces:
        fmin, fmax = _face_axis_extreme(f, axis)
        v = fmin if side == "min" else fmax
        # Require the *entire* face to be on/above(below) the plane within tol to avoid picking X side
        gmin, gmax = _global_axis_extreme(blocks, axis)
        plane = gmin if side == "min" else gmax
        # Face must touch the plane and not protrude past it more than tol_abs
        touch_plane = abs(v - plane) <= tol_abs
        # The opposite extreme should not "cross" the plane the wrong way
        opp = fmax if side == "min" else fmin
        not_past = (opp - plane <= tol_abs) if side == "min" else (plane - opp <= tol_abs)
        if touch_plane and not_past:
            on_plane.append(f)

    pool = {key(f): f for f in on_plane}
    q = [f for f in seed_faces if key(f) in pool]
    visited = set()
    result: List[Face] = []

    while q:
        cur = q.pop()
        kcur = key(cur)
        if kcur in visited:
            continue
        visited.add(kcur)
        result.append(cur)

        # find neighbors by node-sharing (fast and reliable for structured grids)
        for k2, cand in list(pool.items()):
            if k2 in visited or k2 == kcur:
                continue
            if cur.touches_by_nodes(cand,
                                    blocks[cur.BlockIndex], blocks[cand.BlockIndex],
                                    tol_xyz=node_tol_xyz,
                                    min_shared_abs=min_shared_abs,
                                    min_shared_frac=min_shared_frac,
                                    stride_u=1, stride_v=1):
                q.append(cand)

    return result

def find_bounding_faces(blocks: List[Block],
                        outer_faces: List[Dict[str,int]],
                        direction: str = "z",
                        side: str = "both",
                        tol_rel: float = 1e-8,
                        node_tol_xyz: float = 1e-6) -> Tuple[List[Dict[str,int]], List[Dict[str,int]],List[Face], List[Face]]:
    """Find outer bounding faces at the global min/max of a given direction.

    Reduces the mesh by its GCD for alignment, identifies seed faces at
    the global extremes, then uses node-sharing BFS
    (``_bfs_collect_boundary``) to gather continuous boundary patches
    across multiple blocks. The resulting face indices are scaled back to
    the original grid resolution before being returned.

    Args:
        blocks: List of ``Block`` objects comprising the mesh.
        outer_faces: Precomputed outer faces in dictionary form. Pass an
            empty list to have them computed automatically.
        direction: Axis along which to find bounding faces, one of
            ``'x'``, ``'y'``, or ``'z'``. Defaults to ``'z'``.
        side: Which boundary to find: ``'both'``, ``'min'``, or
            ``'max'``. Defaults to ``'both'``.
        tol_rel: Relative tolerance for comparing face coordinates to
            global extremes. Defaults to ``1e-8``.
        node_tol_xyz: Absolute tolerance for node-level coordinate
            matching during BFS adjacency checks. Defaults to ``1e-6``.

    Returns:
        A tuple ``(lower_export, upper_export, lower_faces, upper_faces)``
        where the ``*_export`` entries are lists of face dictionaries
        (suitable for JSON serialization) and the ``*_faces`` entries are
        lists of ``Face`` objects at the original grid scale.
    """
    # 1) Reduce by GCD so grids line up
    gcd_to_use = compute_gcd(blocks)
    blocks_r = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # 2) Build outer face list at reduced resolution
    if len(outer_faces) == 0:
        outer_faces_all: List[Face] = []
        for bi, b in enumerate(blocks_r):
            outs, _ = get_outer_faces(b)
            for o in outs:
                o.set_block_index(bi)
            outer_faces_all.extend(outs)
    else:
        outer_faces_all = outer_face_dict_to_list(blocks_r, outer_faces, gcd_to_use)

    # 3) Absolute tolerance for plane selection
    gmin, gmax = _global_axis_extreme(blocks_r, direction)
    # scale tolerance by magnitude so meshes near origin work too
    tol_abs = max(1.0, abs(gmin) + abs(gmax)) * tol_rel

    # 4) Seeds on min/max planes
    want_min = (side in ("min", "both"))
    want_max = (side in ("max", "both"))

    lower_connected_faces: List[Face] = []
    upper_connected_faces: List[Face] = []

    if want_min:
        seeds_min = _select_seed_faces(outer_faces_all, blocks_r, direction, "min", tol_abs)
        lower_connected_faces = _bfs_collect_boundary(
            seed_faces=seeds_min,
            all_outer_faces=outer_faces_all,
            blocks=blocks_r,
            axis=direction,
            side="min",
            tol_abs=tol_abs,
            node_tol_xyz=node_tol_xyz,
            min_shared_abs=2,          # edge-wise continuity OK
            min_shared_frac=0.005      # small overlap fraction enough to chain
        )

    if want_max:
        seeds_max = _select_seed_faces(outer_faces_all, blocks_r, direction, "max", tol_abs)
        upper_connected_faces = _bfs_collect_boundary(
            seed_faces=seeds_max,
            all_outer_faces=outer_faces_all,
            blocks=blocks_r,
            axis=direction,
            side="max",
            tol_abs=tol_abs,
            node_tol_xyz=node_tol_xyz,
            min_shared_abs=2,
            min_shared_frac=0.005
        )

    # 5) Scale indices back up to the original grid
    #    (Make copies so we don't mutate faces held elsewhere.)
    def _rescale_faces(faces: List[Face]) -> List[Face]:
        out: List[Face] = []
        for f in faces:
            g = Face()
            g.x = f.x.copy(); g.y = f.y.copy(); g.z = f.z.copy()
            g.I = (f.I * gcd_to_use).astype(f.I.dtype)
            g.J = (f.J * gcd_to_use).astype(f.J.dtype)
            g.K = (f.K * gcd_to_use).astype(f.K.dtype)
            g.nvertex = f.nvertex
            g.cx, g.cy, g.cz = f.cx, f.cy, f.cz
            g.blockIndex = f.blockIndex
            g.id = f.id
            out.append(g)
        # de-duplicate by index ranges
        uniq = {}
        for f in out:
            k = (f.BlockIndex, f.IMIN, f.JMIN, f.KMIN, f.IMAX, f.JMAX, f.KMAX)
            uniq[k] = f
        return list(uniq.values())

    lower_connected_faces = _rescale_faces(lower_connected_faces)
    upper_connected_faces = _rescale_faces(upper_connected_faces)

    lower_connected_faces_export = [f.to_dict() for f in lower_connected_faces]
    upper_connected_faces_export = [f.to_dict() for f in upper_connected_faces]

    return (lower_connected_faces_export,
            upper_connected_faces_export,
            lower_connected_faces,
            upper_connected_faces)


def _to_theta(x, y, z, rotation_axis: str):
    """Compute the angular position (theta) about a given rotation axis.

    The convention matches ``Block.cylindrical()`` for each axis:

    - x-axis: ``theta = atan2(Y, Z)``
    - y-axis: ``theta = atan2(Z, X)``
    - z-axis: ``theta = atan2(Y, X)``

    Args:
        x: X coordinate(s), scalar or array.
        y: Y coordinate(s), scalar or array.
        z: Z coordinate(s), scalar or array.
        rotation_axis: Axis of rotation, one of ``'x'``, ``'y'``, or
            ``'z'``.

    Returns:
        The angular position(s) in radians, same shape as the inputs.
    """
    if rotation_axis == "x":
        return np.arctan2(y, z)
    elif rotation_axis == "y":
        return np.arctan2(z, x)
    else:  # "z"
        return np.arctan2(y, x)


def _to_radius(x, y, z, rotation_axis: str):
    """Compute the radial distance from a given rotation axis.

    Args:
        x: X coordinate(s), scalar or array.
        y: Y coordinate(s), scalar or array.
        z: Z coordinate(s), scalar or array.
        rotation_axis: Axis of rotation, one of ``'x'``, ``'y'``, or
            ``'z'``.

    Returns:
        The radial distance(s), same shape as the inputs.
    """
    if rotation_axis == "x":
        return np.sqrt(y * y + z * z)
    elif rotation_axis == "y":
        return np.sqrt(z * z + x * x)
    else:  # "z"
        return np.sqrt(y * y + x * x)


def _global_theta_extreme(blocks: List[Block], rotation_axis: str) -> Tuple[float, float]:
    """Compute the global angular range across all blocks.

    Args:
        blocks: List of ``Block`` objects to scan.
        rotation_axis: Axis of rotation, one of ``'x'``, ``'y'``, or
            ``'z'``.

    Returns:
        A tuple ``(theta_min, theta_max)`` in radians across every
        block's grid points.
    """
    thetas = []
    for b in blocks:
        theta = _to_theta(b.X.ravel(), b.Y.ravel(), b.Z.ravel(), rotation_axis)
        thetas.append(theta)
    all_theta = np.concatenate(thetas)
    return float(all_theta.min()), float(all_theta.max())


def _face_theta_extreme(face: Face, rotation_axis: str) -> Tuple[float, float]:
    """Return the angular range of a face's stored vertices.

    Args:
        face: The ``Face`` whose vertices are inspected.
        rotation_axis: Axis of rotation, one of ``'x'``, ``'y'``, or
            ``'z'``.

    Returns:
        A tuple ``(theta_min, theta_max)`` in radians for the face's
        vertices.
    """
    n = face.nvertex
    theta = _to_theta(face.x[:n], face.y[:n], face.z[:n], rotation_axis)
    return float(theta.min()), float(theta.max())


def find_angular_bounding_faces(
    blocks: List[Block],
    outer_faces: List[Dict[str, int]],
    rotation_axis: str = "x",
    tol_rel: float = 1e-6,
) -> Tuple[List[Dict[str, int]], List[Dict[str, int]], List[Face], List[Face]]:
    """Find outer faces on the angular min/max boundaries of an annular domain.

    This is the rotational analog of ``find_bounding_faces`` for
    translational periodicity. It classifies outer faces whose vertices
    all lie at the global theta minimum or maximum. Works for any
    rotation axis.

    Args:
        blocks: List of ``Block`` objects comprising the mesh.
        outer_faces: Outer faces in dictionary form (at original scale).
        rotation_axis: Axis of rotation, one of ``'x'``, ``'y'``, or
            ``'z'``. Defaults to ``'x'``.
        tol_rel: Relative tolerance for angular comparison, scaled by the
            global theta range. Defaults to ``1e-6``.

    Returns:
        A tuple ``(lower_export, upper_export, lower_faces, upper_faces)``
        where the ``*_export`` entries are lists of face dictionaries and
        the ``*_faces`` entries are lists of ``Face`` objects. Returns
        four empty lists if the domain is non-annular (theta range
        exceeds pi or is negligible).
    """
    axis = rotation_axis.lower().strip()
    assert axis in ("x", "y", "z")

    # Convert outer faces to Face objects (at original scale)
    outer_faces_all = outer_face_dict_to_list(blocks, outer_faces)

    # Global theta range from all block grid points
    theta_min, theta_max = _global_theta_extreme(blocks, axis)
    theta_range = theta_max - theta_min

    # Non-annular check: if domain spans more than half the circle, give up
    if theta_range > np.pi or theta_range < 1e-10:
        return [], [], [], []

    tol_abs = max(1e-8, tol_rel * theta_range)

    # Classify faces: a face is on the lower/upper angular boundary if ALL its
    # corner vertices are at the global theta_min/theta_max within tolerance.
    lower_faces: List[Face] = []
    upper_faces: List[Face] = []

    for f in outer_faces_all:
        f_theta_min, f_theta_max = _face_theta_extreme(f, axis)
        # ALL corners at theta_min
        if (f_theta_max - theta_min) <= tol_abs:
            lower_faces.append(f)
        # ALL corners at theta_max
        elif (theta_max - f_theta_min) <= tol_abs:
            upper_faces.append(f)

    lower_export = [f.to_dict() for f in lower_faces]
    upper_export = [f.to_dict() for f in upper_faces]

    return lower_export, upper_export, lower_faces, upper_faces


def find_closest_block(blocks:List[Block],x:np.ndarray,y:np.ndarray,z:np.ndarray,centroid:np.ndarray,translational_direction:str="x",minvalue:bool=True):
    """Find the block closest to a directional extreme and return a target point.

    Computes a target point offset beyond the global min or max extent
    along the specified axis, then selects the block whose centroid is
    nearest to that target. This is useful for identifying boundary blocks
    for translational periodicity.

    Args:
        blocks: List of ``Block`` objects in the mesh.
        x: Array of X centroid coordinates for all blocks.
        y: Array of Y centroid coordinates for all blocks.
        z: Array of Z centroid coordinates for all blocks.
        centroid: Reference centroid as ``(cx, cy, cz)``.
        translational_direction: Axis along which to search, one of
            ``'x'``, ``'y'``, or ``'z'``. Defaults to ``'x'``.
        minvalue: If ``True``, find the block nearest the minimum extent;
            if ``False``, the maximum. Defaults to ``True``.

    Returns:
        A tuple ``(selected_block_index, target_x, target_y, target_z)``
        where *selected_block_index* is the index of the closest block
        and ``target_*`` are the coordinates of the offset target point.
    """
    cx = centroid[0]; cy = centroid[1]; cz = centroid[2]
    target_x = cx; target_y = cy; target_z = cz
    if translational_direction=="x":
        xmins = [b.X.min() for b in blocks]
        xmaxes = [b.X.max() for b in blocks]
        xmin = min(xmins)
        xmax = max(xmaxes)
        dx = xmax - xmin
        if minvalue:
            target_x = xmin - dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))
        else:
            target_x = xmax + dx*0.5
            selected_block_indx = np.argmin(np.sqrt((target_x-x)**2 + (cy-y)**2 + (cz-z)**2))
        target_y= cy
        target_z = cz
    elif translational_direction=="y":
        ymins = [b.Y.min() for b in blocks]
        ymaxes = [b.Y.max() for b in blocks]
        ymin = min(ymins)
        ymax = max(ymaxes)
        dy = ymax-ymin
        if minvalue:
            target_y = ymin - dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        else:
            target_y = ymax + dy*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (target_y-y)**2 + (cz-z)**2))
        target_x = cx
        target_z = cz
    else: #  translational_direction=="z":
        zmins = [b.Z.min() for b in blocks]
        zmaxes = [b.Z.max() for b in blocks]
        zmin = min(zmins)
        zmax = max(zmaxes)
        dz = zmax - zmin
        if minvalue:
            target_z = zmin - dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))
        else:
            target_z = zmax + dz*0.5
            selected_block_indx = np.argmin(np.sqrt((cx-x)**2 + (cy-y)**2 + (target_z-z)**2))
        target_x = cx
        target_y = cy

    return selected_block_indx,target_x,target_y,target_z


def split_face(face_to_split:Face, block:Block,imin:int,jmin:int,kmin:int,imax:int,jmax:int,kmax:int):
    """Split a face into up to four sub-faces around an inner rectangular region.

    Given a face and an inner rectangle defined by diagonal indices, this
    produces the surrounding sub-faces (top, bottom, left, right) while
    excluding the center region and any degenerate edges. The split
    orientation depends on which index pair is constant (I, J, or K).

    ::

                    left face    top face     right face
         ________       __          __            __
        |   __   |     |  |        |__|          |  |     __
        |  |__|  |     |  |         __           |  |    |__|  face_to_split/center face
        |________|     |__|        |__|          |__|
                                bottom face

    Args:
        face_to_split: The ``Face`` to be subdivided.
        block: The ``Block`` on which the face resides.
        imin: Lower I index of the inner rectangle.
        jmin: Lower J index of the inner rectangle.
        kmin: Lower K index of the inner rectangle.
        imax: Upper I index of the inner rectangle.
        jmax: Upper J index of the inner rectangle.
        kmax: Upper K index of the inner rectangle.

    Returns:
        A list of ``Face`` objects representing the non-degenerate
        sub-faces surrounding the center rectangle. The center face
        itself is excluded.
    """
    center_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmin,kmax=kmax)

    if kmin == kmax:
        # In the picture above Horizontal = i, vertical = j
        left_face = create_face_from_diagonals(block,
                imin=face_to_split.IMIN,imax=imin,
                jmin=face_to_split.JMIN,jmax=face_to_split.JMAX,
                kmin=kmin, kmax=kmax)


        right_face = create_face_from_diagonals(block,
                imin=imax, imax=face_to_split.IMAX,
                jmin=face_to_split.JMIN, jmax=face_to_split.JMAX,
                kmin=kmin, kmax=kmax)

        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmax,jmax=face_to_split.JMAX,
            kmin=kmin,kmax=kmax)

        bottom_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=face_to_split.JMIN,jmax=jmin,
            kmin=kmin,kmax=kmax)

    elif (imin==imax):
        # In the picture above Horizontal = j, vertical = k
        left_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=face_to_split.JMIN, jmax=jmin,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        right_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmax, jmax=face_to_split.JMAX,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmax,kmax=face_to_split.KMAX)

        bottom_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=kmin)

    elif (jmin==jmax):
        # In the picture above Horizontal = i, vertical = k
        left_face = create_face_from_diagonals(block,
            imin=face_to_split.IMIN,imax=imin,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        right_face = create_face_from_diagonals(block,
            imin=imax,imax=face_to_split.IMAX,
            jmin=jmin,jmax=jmax,
            kmin=face_to_split.KMIN,kmax=face_to_split.KMAX)

        top_face = create_face_from_diagonals(block,
            imin=imin,imax=imax,
            jmin=jmin,jmax=jmax,
            kmin=kmax,kmax=face_to_split.KMAX)

        bottom_face = create_face_from_diagonals(block,
            imin=imin, imax=imax,
            jmin=jmin, jmax=jmax,
            kmin=face_to_split.KMIN, kmax=kmin)

    faces = [top_face,bottom_face,left_face,right_face]
    faces = [f for f in faces if not f.isEdge and not f.index_equals(center_face)] # Remove edges
    [f.set_block_index(face_to_split.blockIndex) for f in faces]
    return faces

def find_face_nearest_point(faces:List[Face], x:float,y:float,z:float):
    """Find the face whose centroid is nearest to a given point.

    Args:
        faces: List of ``Face`` objects to search.
        x: X coordinate of the reference point.
        y: Y coordinate of the reference point.
        z: Z coordinate of the reference point.

    Returns:
        The index into *faces* of the face with the smallest Euclidean
        distance from its centroid to ``(x, y, z)``.
    """
    n = list(range(len(faces)))
    dv = list()
    for i in n:
        dx = x-faces[i].cx
        dy = y-faces[i].cy
        dz = z-faces[i].cz
        dv.append(math.sqrt(dx*dx + dy*dy + dz*dz))
    face_index = np.argmin(np.array(dv))
    return face_index

def outer_face_dict_to_list(blocks:List[Block],outer_faces:List[Dict[str,int]],gcd:int=1) -> List[Face]:
    """Convert a list of outer-face dictionaries to ``Face`` objects.

    Each dictionary is expected to contain ``'block_index'``, ``'lb'``
    (tuple of 3 ints), ``'ub'`` (tuple of 3 ints), and optionally
    ``'id'``. Index values are divided by the GCD before face
    construction so that the resulting faces align with a GCD-reduced
    block list.

    Args:
        blocks: List of ``Block`` objects (possibly GCD-reduced).
        outer_faces: List of dictionaries, each describing one outer face.
        gcd: Greatest common divisor used to scale down face indices.
            Defaults to ``1`` (no scaling).

    Returns:
        A list of ``Face`` objects corresponding to the input
        dictionaries.
    """
    outer_faces_all = list()
    for o in outer_faces:
        face = create_face_from_diagonals(blocks[o['block_index']],
            int(o['lb'][0]/gcd), int(o['lb'][1]/gcd), int(o['lb'][2]/gcd),
            int(o['ub'][0]/gcd), int(o['ub'][1]/gcd), int(o['ub'][2]/gcd))
        if 'id' in o.keys():
            face.id = o['id']
        face.set_block_index(o['block_index'])
        outer_faces_all.append(face)

    return outer_faces_all

def match_faces_dict_to_list(blocks:List[Block],matched_faces:List[Dict[str,int]],gcd:int=1):
    """Convert a list of matched-face dictionaries to ``Face`` objects.

    Each dictionary is expected to contain ``'block1'`` and ``'block2'``
    sub-dictionaries, each with ``'block_index'``, ``'lb'`` (tuple of 3
    ints), ``'ub'`` (tuple of 3 ints), and optionally ``'id'``. The
    function creates two ``Face`` objects per entry (one for each side of
    the match) and returns them all in a flat list.

    Args:
        blocks: List of ``Block`` objects (possibly GCD-reduced).
        matched_faces: List of dictionaries, each describing a pair of
            matched faces.
        gcd: GCD factor to scale down face indices. Defaults to ``1``
            (no scaling).

    Returns:
        A flat list of ``Face`` objects, with two entries per matched-face
        dictionary (one for ``'block1'``, one for ``'block2'``).
    """
    matched_faces_all = list()
    for _,m in enumerate(matched_faces):
        face1 = create_face_from_diagonals(blocks[m['block1']['block_index']],
                            int(m['block1']['lb'][0]/gcd), int(m['block1']['lb'][1]/gcd), int(m['block1']['lb'][2]/gcd),
                            int(m['block1']['ub'][0]/gcd), int(m['block1']['ub'][1]/gcd), int(m['block1']['ub'][2]/gcd))
        face2 = create_face_from_diagonals(blocks[m['block2']['block_index']],
                            int(m['block2']['lb'][0]/gcd), int(m['block2']['lb'][1]/gcd), int(m['block2']['lb'][2]/gcd),
                            int(m['block2']['ub'][0]/gcd), int(m['block2']['ub'][1]/gcd), int(m['block2']['ub'][2]/gcd))
        face1.set_block_index(m['block1']['block_index'])
        if 'id' in m['block1'].keys():
            face1.id = m['block1']['id']
        face2.set_block_index(m['block2']['block_index'])
        if 'id' in m['block2'].keys():
            face2.id = m['block2']['id']
        matched_faces_all.append(face1)
        matched_faces_all.append(face2)
    return matched_faces_all
