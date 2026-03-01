"""Rotational and translational periodicity matching for structured multi-block meshes.

This module provides algorithms to detect periodic face pairs between blocks
in a Plot3D multi-block structured grid. Two types of periodicity are supported:

**Rotational periodicity** (``rotated_periodicity``):
    Identifies outer faces that match after rotation by a given angle about
    a coordinate axis (x, y, or z). Uses a three-phase algorithm with
    ``FacePool`` theta bucketing for efficient candidate selection:

    - Phase 1: Full-face corner matching (fast, O(N log N)).
    - Phase 2: Split-face matching with corner pre-checks (iterative).
    - Phase 3: Relaxed-tolerance matching for wavy-surface faces.

**Translational periodicity** (``translational_periodicity``):
    Identifies outer faces that match after translation by a fixed distance
    along a coordinate axis. Uses bounding-face detection, adaptive node
    tolerances, and orthogonal-plane pre-checks to efficiently pair faces.

Both algorithms support GCD-based mesh reduction for faster processing and
automatically scale results back to the original grid resolution.
"""
from typing import List, Dict, Tuple, Optional
import numpy as np
import math
from math import cos, radians, sin, sqrt, acos
from .block import Block, compute_gcd, reduce_blocks
from .blockfunctions import rotate_block
from .face import Face
from .facefunctions import (outer_face_dict_to_list, match_faces_dict_to_list,
    create_face_from_diagonals, find_bounding_faces, get_outer_faces,
    find_angular_bounding_faces, _to_theta, _to_radius)
from .connectivity import get_face_intersection, face_matches_to_dict
from .utils import (scale_face_dict_indices, face_key, face_grid_dims)
from tqdm import tqdm

def _scale_dicts_down(dicts: List[Dict], gcd: int):
    """Scale dictionary face lb/ub tuple indices down by a GCD factor.

    Divides each index in the ``'lb'`` and ``'ub'`` tuples by ``gcd``
    using integer division. This is used to convert face indices from the
    original grid resolution to the GCD-reduced resolution.

    Args:
        dicts: List of face dictionaries containing ``'lb'`` and ``'ub'``
            tuple indices.
        gcd: Greatest common divisor factor to divide indices by.
    """
    for d in dicts:
        if 'lb' in d:
            d['lb'] = tuple(int(v / gcd) for v in d['lb'])
        if 'ub' in d:
            d['ub'] = tuple(int(v / gcd) for v in d['ub'])


def create_rotation_matrix(rotation_angle:float, rotation_axis:str="x"):
    """Create a 3x3 rotation matrix for a given angle and axis.

    Constructs a standard right-hand rotation matrix that rotates
    coordinates about one of the three principal axes.

    Args:
        rotation_angle: Rotation angle in radians.
        rotation_axis: Axis of rotation: ``'x'``, ``'y'``, or ``'z'``.
            Defaults to ``'x'``.

    Returns:
        numpy.ndarray: A 3x3 rotation matrix.

    Raises:
        ValueError: If ``rotation_axis`` is not ``'x'``, ``'y'``, or ``'z'``.
    """

    if rotation_axis=='x':
        rotation_matrix = np.array([[1,0,0],
                            [0,cos(rotation_angle),-sin(rotation_angle)],
                            [0,sin(rotation_angle),cos(rotation_angle)]])

    elif rotation_axis=='y':
        rotation_matrix = np.array([[cos(rotation_angle),0,sin(rotation_angle)],
                            [0,1,0],
                            [-sin(rotation_angle),0,cos(rotation_angle)]])
    elif rotation_axis=='z':
        rotation_matrix = np.array([[cos(rotation_angle),-sin(rotation_angle), 0],
                            [sin(rotation_angle),cos(rotation_angle), 0],
                            [0, 0, 1]])
    else:
        raise ValueError(f"rotation_axis must be 'x', 'y', or 'z', got '{rotation_axis}'")

    return rotation_matrix


def _faces_could_match_rotationally(
    face1: Face, face2: Face,
    rotation_matrix: np.ndarray,
    rotation_axis: str,
    tol_rel: float = 0.1,
) -> bool:
    """Perform cheap geometric pre-checks to reject non-matching face pairs.

    Applies a series of inexpensive tests in order of increasing cost to
    quickly eliminate face pairs that cannot be periodic matches. The checks
    performed are:

    1. Both faces must be planar (have a valid ``const_type``).
    2. Axial extents along the rotation axis must overlap.
    3. Radial extents (distance from rotation axis) must overlap.
    4. Rotated centroid of ``face1`` must be within 1.5x the diagonal
       length of the target face.

    Args:
        face1: Source face to be rotated.
        face2: Target face to compare against.
        rotation_matrix: 3x3 rotation matrix to apply to ``face1``.
        rotation_axis: Axis of rotation: ``'x'``, ``'y'``, or ``'z'``.
        tol_rel: Relative tolerance as a fraction of the span used for
            axial and radial overlap checks. Defaults to ``0.1``.

    Returns:
        bool: ``True`` if the face pair passes all pre-checks and could
        potentially be a rotational match; ``False`` otherwise.
    """
    # 1. Each face must be planar (constant on some axis), but not necessarily the same one
    if face1.const_type == -1 or face2.const_type == -1:
        return False

    # 2. Axial extent overlap (along rotation axis)
    axis = rotation_axis.lower()
    n1, n2 = face1.nvertex, face2.nvertex
    if axis == "x":
        f1_ax = (float(face1.x[:n1].min()), float(face1.x[:n1].max()))
        f2_ax = (float(face2.x[:n2].min()), float(face2.x[:n2].max()))
    elif axis == "y":
        f1_ax = (float(face1.y[:n1].min()), float(face1.y[:n1].max()))
        f2_ax = (float(face2.y[:n2].min()), float(face2.y[:n2].max()))
    else:
        f1_ax = (float(face1.z[:n1].min()), float(face1.z[:n1].max()))
        f2_ax = (float(face2.z[:n2].min()), float(face2.z[:n2].max()))

    axial_span = max(f1_ax[1] - f1_ax[0], f2_ax[1] - f2_ax[0], 1e-12)
    tol_axial = tol_rel * axial_span
    if f1_ax[1] + tol_axial < f2_ax[0] or f2_ax[1] + tol_axial < f1_ax[0]:
        return False

    # 3. Radial extent overlap
    r1 = _to_radius(face1.x[:n1], face1.y[:n1], face1.z[:n1], axis)
    r2 = _to_radius(face2.x[:n2], face2.y[:n2], face2.z[:n2], axis)
    r1_range = (float(r1.min()), float(r1.max()))
    r2_range = (float(r2.min()), float(r2.max()))
    radial_span = max(r1_range[1] - r1_range[0], r2_range[1] - r2_range[0], 1e-12)
    tol_radial = tol_rel * radial_span
    if r1_range[1] + tol_radial < r2_range[0] or r2_range[1] + tol_radial < r1_range[0]:
        return False

    # 4. Rotated centroid proximity
    c1 = np.array([[face1.cx], [face1.cy], [face1.cz]])
    c1_rotated = (rotation_matrix @ c1).flatten()
    c2 = np.array([face2.cx, face2.cy, face2.cz])
    dist = float(np.sqrt(np.sum((c1_rotated - c2) ** 2)))
    max_diag = max(face1.diagonal_length, face2.diagonal_length)
    if dist > max_diag * 1.5:
        return False

    return True


def _count_rotated_corners_on_face(face_a: Face, face_b: Face, block_b: Block,
                                    rotation_matrix: np.ndarray, tol: float = 1e-6) -> int:
    """Count how many corners of a face land on another face's grid after rotation.

    Rotates the corner vertices of ``face_a`` using the given rotation matrix,
    then checks how many of those rotated corners are within ``tol`` of any
    grid point on ``face_b``. This serves as a cheap pre-check: if fewer than
    2 rotated corners hit ``face_b``, the expensive full intersection is skipped.

    Note:
        Ported from the Rust ``count_rotated_corners_on_face`` implementation.

    Args:
        face_a: Source face whose corners will be rotated.
        face_b: Target face to check against.
        block_b: Block containing ``face_b``, used to extract grid points.
        rotation_matrix: 3x3 rotation matrix to apply to ``face_a`` corners.
        tol: Euclidean distance tolerance for considering a corner matched.
            Defaults to ``1e-6``.

    Returns:
        int: Number of ``face_a`` corners (0 to 4) that match ``face_b``
        grid points after rotation.
    """
    pts_b = face_b.grid_points(block_b)
    if pts_b.size == 0:
        return 0

    # Get face_a's 4 corner coordinates
    n = face_a.nvertex
    corners_a = np.column_stack([face_a.x[:n], face_a.y[:n], face_a.z[:n]])

    # Rotate corners
    rotated = (rotation_matrix @ corners_a.T).T  # (n, 3)

    count = 0
    for rc in rotated:
        # Check if this rotated corner is close to any point on face_b
        dists = np.sqrt(np.sum((pts_b - rc) ** 2, axis=1))
        if np.min(dists) < tol:
            count += 1
    return count


def _full_face_match_rotated(face_a: Face, face_b: Face, rotation_matrix: np.ndarray,
                              tol: float = 1e-6) -> bool:
    """Check if all corners of two faces match bijectively after rotation.

    Rotates all corners of ``face_a`` and checks whether each rotated corner
    maps to a unique corner of ``face_b`` within the given tolerance. Both
    faces must have the same number of vertices and matching grid dimensions
    (allowing transposition).

    This is used in Phase 1 of periodic matching to quickly identify full-face
    periodic matches without expensive node-by-node intersection.

    Note:
        Ported from the Rust ``full_face_match_transformed`` implementation.

    Args:
        face_a: Source face whose corners will be rotated.
        face_b: Target face to compare against.
        rotation_matrix: 3x3 rotation matrix to apply to ``face_a`` corners.
        tol: Euclidean distance tolerance for corner matching.
            Defaults to ``1e-6``.

    Returns:
        bool: ``True`` if all corners match bijectively after rotation;
        ``False`` otherwise.
    """
    n_a = face_a.nvertex
    n_b = face_b.nvertex
    if n_a != n_b or n_a < 4:
        return False

    # Face dimensions must match (allowing transpose)
    from .utils import face_grid_dims
    if (face_grid_dims(face_a.IMIN, face_a.IMAX, face_a.JMIN, face_a.JMAX, face_a.KMIN, face_a.KMAX) !=
        face_grid_dims(face_b.IMIN, face_b.IMAX, face_b.JMIN, face_b.JMAX, face_b.KMIN, face_b.KMAX)):
        return False

    corners_a = np.column_stack([face_a.x[:n_a], face_a.y[:n_a], face_a.z[:n_a]])
    rotated_a = (rotation_matrix @ corners_a.T).T

    corners_b = np.column_stack([face_b.x[:n_b], face_b.y[:n_b], face_b.z[:n_b]])

    # Bijective matching: each rotated corner must match a unique face_b corner
    used = set()
    for rc in rotated_a:
        dists = np.sqrt(np.sum((corners_b - rc) ** 2, axis=1))
        best_idx = int(np.argmin(dists))
        if dists[best_idx] > tol or best_idx in used:
            return False
        used.add(best_idx)
    return len(used) == n_b


def _match_periodic_faces(
    blocks: List[Block],
    outer_faces_all: List[Face],
    matched_faces_all: List[Face],
    rotation_matrices: List[np.ndarray],
    rotation_axis: str,
    periodic_direction: Optional[str] = None,
) -> Tuple[List[Dict], List[Dict], List[Tuple], List[Face]]:
    """Run the three-phase rotational periodicity matching engine.

    This is the core matching algorithm used by ``rotated_periodicity``. It
    operates on GCD-reduced blocks and faces, using ``FacePool`` theta
    bucketing for O(log N) candidate selection per face.

    The three phases are:

    - **Phase 1** -- Full-face matching via corner comparison. For each face,
      find candidates at theta +/- ``rotation_angle`` and attempt a bijective
      corner match. No expensive node-by-node intersection is needed.
    - **Phase 2** -- Split-face matching with corner pre-check. Uses
      ``FacePool`` to find candidates, checks if at least 2 rotated corners
      land on the candidate face, then performs full intersection. Iterates
      until convergence (no iteration limit).
    - **Phase 3** -- Relaxed matching without pre-checks. Uses 5x tolerance
      (``5e-6``) to catch wavy-surface faces. Runs until convergence.

    Args:
        blocks: List of GCD-reduced blocks.
        outer_faces_all: Outer faces as ``Face`` objects at reduced scale.
        matched_faces_all: Already-matched connectivity faces as ``Face``
            objects at reduced scale, used for exclusion.
        rotation_matrices: List of rotation matrices to try (typically
            ``[+angle, -angle]``).
        rotation_axis: Axis of rotation: ``'x'``, ``'y'``, or ``'z'``.
        periodic_direction: If provided, only check faces whose constant
            index matches the given direction (``'i'``, ``'j'``, or ``'k'``).
            ``None`` means check all directions. Defaults to ``None``.

    Returns:
        tuple: A 4-tuple containing:

            - **periodic_faces_export** (*List[Dict]*): Periodic face pairs
              as export-ready dictionaries.
            - **outer_faces_export** (*List[Dict]*): Remaining unmatched
              outer faces as dictionaries.
            - **periodic_faces** (*List[Tuple]*): Periodic face pairs as
              tuples of ``Face`` objects.
            - **outer_faces_final** (*List[Face]*): Remaining unmatched
              outer faces as ``Face`` objects.
    """
    from .face_pool import FacePool

    # Build rotation caches (one per matrix, lazy)
    caches = {}
    for idx, mat in enumerate(rotation_matrices):
        caches[idx] = {}

    def get_rotated(matrix_idx: int, block_index: int) -> Block:
        """Retrieve a rotated block from cache, computing it lazily if needed.

        Args:
            matrix_idx: Index into ``rotation_matrices`` selecting which
                rotation to apply.
            block_index: Index of the block to rotate.

        Returns:
            Block: The rotated block.
        """
        cache = caches[matrix_idx]
        if block_index not in cache:
            cache[block_index] = rotate_block(blocks[block_index], rotation_matrices[matrix_idx])
        return cache[block_index]

    # Try to identify angular boundary faces
    outer_dicts = [f.to_dict() for f in outer_faces_all]
    lower_export, upper_export, lower_faces, upper_faces = find_angular_bounding_faces(
        blocks, outer_dicts, rotation_axis
    )
    del outer_dicts

    use_angular = len(lower_faces) > 0 and len(upper_faces) > 0

    # Use shared face_key utility (imported from utils)
    _face_key = face_key

    # Filter by periodic_direction and angular boundaries
    dir_const_map = {"i": 0, "j": 1, "k": 2}

    if use_angular:
        lower_keys = {_face_key(f) for f in lower_faces}
        upper_keys = {_face_key(f) for f in upper_faces}
        pool_faces = [f for f in outer_faces_all
                      if _face_key(f) in lower_keys or _face_key(f) in upper_keys]
    else:
        lower_keys = set()
        upper_keys = set()
        pool_faces = list(outer_faces_all)

    if periodic_direction is not None:
        required_const = dir_const_map.get(periodic_direction.lower())
        if required_const is not None:
            pool_faces = [f for f in pool_faces if f.const_type == required_const]

    # Compute rotation angle from the first rotation matrix
    # Extract angle from rotation matrix (works for any axis)
    rot_mat_0 = rotation_matrices[0]
    # The rotation angle can be extracted from the trace: trace = 1 + 2*cos(theta)
    trace = rot_mat_0[0, 0] + rot_mat_0[1, 1] + rot_mat_0[2, 2]
    rotation_angle_rad = math.acos(max(-1.0, min(1.0, (trace - 1.0) / 2.0)))
    theta_tol = abs(rotation_angle_rad) * 0.15 + 0.05  # Same as Rust

    # Build FacePool with theta bucketing
    pool = FacePool(pool_faces, blocks, rotation_axis)

    # Track results
    periodic_faces: List[Tuple] = []
    periodic_faces_export: List[Dict] = []
    split_faces_all: List[Face] = []
    non_matching: set = set()
    seen_pairs: set = set()

    MATCH_TOL = 1e-6

    # ===== PHASE 1: Full-face matching via corner comparison =====
    # Fast O(N log N) pass: for each face, check if all 4 rotated corners
    # match a candidate's corners bijectively. No expensive intersection needed.
    active = pool.active_indices()
    pbar_p1 = tqdm(total=len(active), desc="Periodic matching (Phase 1)", unit="face")
    for face_a_idx in active:
        pbar_p1.update(1)
        if pool.is_consumed(face_a_idx):
            continue
        face_a = pool.faces[face_a_idx]
        if face_a.const_type == -1:
            continue

        candidates = pool.find_rotational_candidates(
            face_a_idx, rotation_angle_rad, theta_tol
        )

        for cand_idx in candidates:
            if pool.is_consumed(cand_idx):
                continue
            face_b = pool.faces[cand_idx]
            if face_b.const_type == -1:
                continue

            for mat_idx, rot_mat in enumerate(rotation_matrices):
                if _full_face_match_rotated(face_a, face_b, rot_mat, MATCH_TOL):
                    pair_key = (min(_face_key(face_a), _face_key(face_b)),
                                max(_face_key(face_a), _face_key(face_b)))
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)

                    block_a_rot = get_rotated(mat_idx, face_a.blockIndex)
                    block_b = blocks[face_b.blockIndex]

                    periodic_faces.append((face_a, face_b))
                    periodic_faces_export.append(
                        face_matches_to_dict(face_a, face_b, block_a_rot, block_b)
                    )

                    pool.consume(face_a_idx)
                    pool.consume(cand_idx)

                    pbar_p1.set_description(
                        f"Phase 1: block {face_a.blockIndex} <-> block {face_b.blockIndex}")
                    pbar_p1.set_postfix(matches=len(periodic_faces))
                    break
            if pool.is_consumed(face_a_idx):
                break
    pbar_p1.close()

    # ===== PHASE 2: Split-face matching with corner pre-check =====
    pbar = tqdm(total=len(pool_faces), desc="Periodic matching (Phase 2)", unit="pair")

    changed = True
    iteration = 0
    while changed:
        changed = False
        iteration += 1
        active = pool.active_indices()

        for face_a_idx in active:
            if pool.is_consumed(face_a_idx):
                continue
            face_a = pool.faces[face_a_idx]

            # Skip non-planar faces
            if face_a.const_type == -1:
                continue

            # Use FacePool to find candidates via theta bucketing
            candidates = pool.find_rotational_candidates(
                face_a_idx, rotation_angle_rad, theta_tol
            )

            match_found = False
            for cand_idx in candidates:
                if pool.is_consumed(cand_idx):
                    continue
                face_b = pool.faces[cand_idx]

                if face_b.const_type == -1:
                    continue

                # Skip known non-matching pairs
                pair_id = (face_a_idx, cand_idx)
                if pair_id in non_matching:
                    continue

                # Try each rotation matrix
                for mat_idx, rot_mat in enumerate(rotation_matrices):
                    # Corner pre-check: at least 2 rotated corners must hit face_b
                    corners_hit = _count_rotated_corners_on_face(
                        face_a, face_b, blocks[face_b.blockIndex], rot_mat, MATCH_TOL
                    )
                    if corners_hit < 2:
                        continue

                    block_a_rot = get_rotated(mat_idx, face_a.blockIndex)
                    block_b = blocks[face_b.blockIndex]

                    _, periodic_temp, split_temp = __periodicity_check__(
                        face_a, face_b, block_a_rot, block_b, tol=MATCH_TOL
                    )

                    if len(periodic_temp) > 0:
                        # Check for duplicate pair
                        pair_key = (min(_face_key(periodic_temp[0]), _face_key(periodic_temp[1])),
                                    max(_face_key(periodic_temp[0]), _face_key(periodic_temp[1])))
                        if pair_key in seen_pairs:
                            continue
                        seen_pairs.add(pair_key)

                        periodic_faces.append(periodic_temp)
                        periodic_faces_export.append(
                            face_matches_to_dict(periodic_temp[0], periodic_temp[1], block_a_rot, block_b)
                        )

                        # Consume matched faces
                        pool.consume(face_a_idx)
                        pool.consume(cand_idx)

                        # Add split remnants back to pool
                        for sf in split_temp:
                            pool.add_face(sf)
                        split_faces_all.extend(split_temp)

                        changed = True
                        match_found = True
                        pbar.update(1)
                        pbar.set_description(
                            f"Phase 2: block {face_a.blockIndex} <-> block {face_b.blockIndex}")
                        pbar.set_postfix(matches=len(periodic_faces), iteration=iteration)
                        break

                if match_found:
                    break

                non_matching.add(pair_id)

            if match_found:
                break  # Restart search from beginning after a match

    pbar.close()

    # ===== PHASE 3: Relaxed matching for remaining wavy-surface faces =====
    relaxed_tol = 5e-6  # 5x the default
    active = pool.active_indices()

    phase3_found = True
    while phase3_found:
        phase3_found = False
        active = pool.active_indices()
        for i_a, face_a_idx in enumerate(active):
            if pool.is_consumed(face_a_idx):
                continue
            face_a = pool.faces[face_a_idx]
            if face_a.const_type == -1:
                continue

            for i_b in range(i_a + 1, len(active)):
                face_b_idx = active[i_b]
                if pool.is_consumed(face_b_idx):
                    continue
                face_b = pool.faces[face_b_idx]
                if face_b.const_type == -1:
                    continue

                for mat_idx, rot_mat in enumerate(rotation_matrices):
                    block_a_rot = get_rotated(mat_idx, face_a.blockIndex)
                    block_b = blocks[face_b.blockIndex]

                    _, periodic_temp, split_temp = __periodicity_check__(
                        face_a, face_b, block_a_rot, block_b, tol=relaxed_tol
                    )

                    if len(periodic_temp) > 0:
                        pair_key = (min(_face_key(periodic_temp[0]), _face_key(periodic_temp[1])),
                                    max(_face_key(periodic_temp[0]), _face_key(periodic_temp[1])))
                        if pair_key in seen_pairs:
                            continue
                        seen_pairs.add(pair_key)

                        periodic_faces.append(periodic_temp)
                        periodic_faces_export.append(
                            face_matches_to_dict(periodic_temp[0], periodic_temp[1], block_a_rot, block_b)
                        )

                        pool.consume(face_a_idx)
                        pool.consume(face_b_idx)

                        for sf in split_temp:
                            pool.add_face(sf)
                        split_faces_all.extend(split_temp)

                        phase3_found = True
                        break

                if phase3_found:
                    break
            if phase3_found:
                break

    # Free rotation caches
    caches.clear()

    # Build removal set for outer faces using lightweight tuple keys
    remove_keys = set()
    for p in periodic_faces:
        remove_keys.add(_face_key(p[0]))
        remove_keys.add(_face_key(p[1]))
    for m in matched_faces_all:
        remove_keys.add(_face_key(m))

    outer_faces_final = [p for p in outer_faces_all if _face_key(p) not in remove_keys]
    # Add non-matched split faces back
    final_keys = {_face_key(f) for f in outer_faces_final}
    for sf in split_faces_all:
        k = _face_key(sf)
        if k not in remove_keys and k not in final_keys:
            outer_faces_final.append(sf)
            final_keys.add(k)
    del remove_keys, final_keys

    # Build outer faces export
    outer_faces_export = [o.to_dict() for o in outer_faces_final]

    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_final


def _scale_results_up(periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all, gcd_to_use):
    """Scale GCD-reduced results back to the original grid resolution.

    Multiplies all face index values (IMIN, JMIN, KMIN, IMAX, JMAX, KMAX)
    and Face I/J/K attributes by the GCD factor to restore original-scale
    indices. Operates in-place on all four result containers.

    Args:
        periodic_faces_export: List of periodic face-pair dictionaries with
            nested ``'block1'`` and ``'block2'`` sub-dicts.
        outer_faces_export: List of remaining outer face dictionaries.
        periodic_faces: List of periodic face-pair tuples containing
            ``Face`` objects.
        outer_faces_all: List of remaining outer ``Face`` objects.
        gcd_to_use: GCD factor to multiply indices by. If ``<= 1``, this
            function returns immediately with no changes.
    """
    if gcd_to_use <= 1:
        return

    scale_face_dict_indices(periodic_faces_export, gcd_to_use, nested_sides=['block1', 'block2'])

    for pair in periodic_faces:
        for face in pair:
            face.I *= gcd_to_use
            face.J *= gcd_to_use
            face.K *= gcd_to_use

    scale_face_dict_indices(outer_faces_export, gcd_to_use)

    for face in outer_faces_all:
        face.I *= gcd_to_use
        face.J *= gcd_to_use
        face.K *= gcd_to_use


def rotated_periodicity(blocks:List[Block], matched_faces:List[Dict[str,int]], outer_faces:List[Dict[str,int]],
                        rotation_angle:float|None=None, rotation_axis:str = "x",
                        periodic_direction:str|None=None, nblades:int | None=None,
                        ReduceMesh:bool=True):
    """Find rotationally periodic face pairs between blocks.

    Detects which outer faces match after rotation by the given angle about
    the specified axis. Automatically reduces the mesh by GCD for faster
    matching, uses angular boundary detection to narrow candidates, and
    applies geometric pre-checks to reject non-matching pairs cheaply.

    Both ``+angle`` and ``-angle`` rotations are tried for robustness.

    Args:
        blocks: List of blocks. Do not duplicate before passing in; the
            function copies internally when reducing.
        matched_faces: Matched faces from connectivity, as dictionaries
            with keys ``'block_index'``, ``'lb'``, ``'ub'``, etc.
        outer_faces: Outer (unmatched) faces in dictionary form.
        rotation_angle: Rotation angle in degrees. Either this or
            ``nblades`` must be provided. If both are given,
            ``rotation_angle`` takes precedence.
        rotation_axis: Axis of rotation: ``'x'``, ``'y'``, or ``'z'``.
            Defaults to ``'x'``.
        periodic_direction: Filter to only check faces with a constant
            ``'i'``, ``'j'``, or ``'k'`` index. ``None`` means check all
            directions. Defaults to ``None``.
        nblades: Number of blades. Used to compute
            ``rotation_angle = 360.0 / nblades`` when ``rotation_angle``
            is not provided.
        ReduceMesh: If ``True``, reduces the mesh by GCD for faster
            matching. Defaults to ``True``.

    Returns:
        tuple: A 4-tuple containing:

            - **periodic_faces_export** (*List[Dict[str, int]]*): Periodic
              face pairs as export-ready dictionaries.
            - **outer_faces_export** (*List[Dict[str, int]]*): Remaining
              outer faces as dictionaries.
            - **periodic_faces** (*List[Tuple[Face, Face]]*): Periodic face
              pairs as ``Face`` objects.
            - **outer_faces_all** (*List[Face]*): Remaining outer faces as
              ``Face`` objects.

    Raises:
        ValueError: If neither ``rotation_angle`` nor ``nblades`` is
            provided.

    Example::

        # Using rotation_angle directly
        periodic_faces, outer_faces_export, _, _ = rotated_periodicity(
            blocks, face_matches, outer_faces,
            rotation_angle=6.545, rotation_axis="x")

        # Using nblades (computes 360/nblades)
        periodic_faces, outer_faces_export, _, _ = rotated_periodicity(
            blocks, face_matches, outer_faces,
            nblades=55, rotation_axis="x", periodic_direction="k")
    """
    # Resolve rotation angle
    if rotation_angle is None and nblades is not None:
        rotation_angle = 360.0 / nblades
    elif rotation_angle is None:
        raise ValueError("Either rotation_angle (degrees) or nblades must be provided.")

    original_blocks = blocks  # preserve reference before GCD reduction
    gcd_to_use = 1
    if ReduceMesh:
        gcd_to_use = compute_gcd(blocks)
        blocks = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # Build both +angle and -angle rotation matrices for robustness
    rotation_matrix_pos = create_rotation_matrix(radians(rotation_angle), rotation_axis)
    rotation_matrix_neg = create_rotation_matrix(radians(-rotation_angle), rotation_axis)

    outer_faces_all = outer_face_dict_to_list(blocks, outer_faces, gcd_to_use)
    matched_faces_all = match_faces_dict_to_list(blocks, matched_faces, gcd_to_use)

    periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all = _match_periodic_faces(
        blocks, outer_faces_all, matched_faces_all,
        [rotation_matrix_pos, rotation_matrix_neg], rotation_axis, periodic_direction
    )
    # Free reduced blocks and intermediate face lists — no longer needed after matching
    del blocks, matched_faces_all

    _scale_results_up(periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all, gcd_to_use)

    # Verify and correct diagonal corner ordering
    if periodic_faces_export:
        verified, mismatched = verify_periodicity(original_blocks, periodic_faces_export, rotation_angle, rotation_axis, tol=1e-4)
        if mismatched:
            print(f"verify_periodicity: {len(mismatched)} mismatched periodic faces could not be corrected")
        periodic_faces_export = verified + mismatched

    return periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all

def translational_periodicity(
    blocks: List[Block],
    outer_faces: List[Dict[str,int]],
    delta: float = None,
    translational_direction: str = "z",
    node_tol_xyz: float = None,        # global override; if None we compute per-pair adaptively
    min_shared_frac: float = 0.02,
    min_shared_abs: int = 4,
    stride_u: int = 1,
    stride_v: int = 1,
    ) -> Tuple[
               List[Dict[str, Dict[str, int]]],
               List[Tuple[Face, Face, Dict[str, str]]],
               List[Dict[str,int]]
            ]:
    """Find translationally periodic face pairs between blocks along a given axis.

    Detects which outer block faces are periodic counterparts across the domain
    in the specified translational direction (``'x'``, ``'y'``, or ``'z'``).
    The algorithm proceeds as follows:

    1. **Bounding faces**: Uses ``find_bounding_faces`` to identify candidate
       lower/upper faces for the given axis.
    2. **Grid reduction**: Reduces blocks to their GCD resolution for
       consistent indexing across blocks.
    3. **Shifting**: Creates shifted copies of all blocks in both positive
       and negative directions along the periodic axis.
    4. **Orthogonal-plane precheck**: Uses a fast projection test orthogonal
       to the periodic axis to determine whether two faces could match.
    5. **Node-based match**: Calls ``Face.touches_by_nodes`` on candidate
       pairs with an adaptive tolerance based on in-plane spacing.
    6. **Pairing**: Records each valid periodic pair with IJK index mappings
       and removes matched faces from the outer-face list.
    7. **Scaling back**: Rescales reduced indices to original grid spacing.

    Args:
        blocks: List of blocks.
        outer_faces: Outer faces as dictionaries with keys ``'block_index'``,
            ``'lb'`` (tuple of 3 ints), and ``'ub'`` (tuple of 3 ints).
        delta: Periodicity spacing along the chosen axis. If ``None``, it is
            inferred from the global block min/max extent.
        translational_direction: Axis to check: ``'x'``, ``'y'``, or
            ``'z'``. Defaults to ``'z'``.
        node_tol_xyz: Absolute coordinate tolerance for node matching. If
            ``None``, tolerance is computed adaptively based on median
            in-plane spacing of candidate faces.
        min_shared_frac: Minimum fraction of nodes that must overlap for
            two faces to be considered periodic. Defaults to ``0.02``.
        min_shared_abs: Minimum absolute number of shared nodes.
            Defaults to ``4``.
        stride_u: Subsampling stride along the first face index direction.
            Defaults to ``1`` (no skipping).
        stride_v: Subsampling stride along the second face index direction.
            Defaults to ``1`` (no skipping).

    Returns:
        tuple: A 3-tuple containing:

            - **periodic_faces_export** (*List[Dict[str, Dict[str, int]]]*):
              Export-ready dictionaries describing each periodic pair,
              including block indices, face extents, index mapping, and
              match mode.
            - **periodic_pairs** (*List[Tuple[Face, Face, Dict[str, str]]]*):
              Matched periodic face pairs as ``Face`` objects with IJK
              mapping dictionaries.
            - **outer_faces_remaining** (*List[Dict[str, int]]*): Updated
              list of outer faces with periodic ones removed, preserving
              any existing ``'id'`` fields.

    Note:
        The adaptive tolerance makes the method robust to small spacing
        differences between blocks. The orthogonal-plane precheck avoids
        expensive node comparisons when faces clearly do not align.
    """
    # 0) lower/upper via your finder (dicts at original scale)
    lower_connected_faces, upper_connected_faces, _, _ = find_bounding_faces(
        blocks, outer_faces, translational_direction, "both"
    )

    axis = translational_direction.lower().strip()
    assert axis in ("x","y","z")

    # 1) GCD reduce
    gcd_to_use = compute_gcd(blocks)

    lower_faces_r = outer_face_dict_to_list(blocks, lower_connected_faces, gcd_to_use)
    upper_faces_r = outer_face_dict_to_list(blocks, upper_connected_faces, gcd_to_use)
    blocks_r = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # 2) Δ along axis (if not provided)
    if axis == "x":
        a_min = min(b.X.min() for b in blocks_r); a_max = max(b.X.max() for b in blocks_r)
    elif axis == "y":
        a_min = min(b.Y.min() for b in blocks_r); a_max = max(b.Y.max() for b in blocks_r)
    else:
        a_min = min(b.Z.min() for b in blocks_r); a_max = max(b.Z.max() for b in blocks_r)
    d_axis = (a_max - a_min) if (delta is None) else float(delta)

    # 3) Shifted copies
    def shift_blocks(bb: List[Block], amount: float) -> List[Block]:
        """Create shifted copies of all blocks along the periodic axis.

        Args:
            bb: List of blocks to shift.
            amount: Distance to shift along the periodic axis.

        Returns:
            List[Block]: New block copies shifted by ``amount``.
        """
        cp = [b.copy() for b in bb]
        for b in cp:
            b.shift(amount, axis)
        return cp

    blocks_up = shift_blocks(blocks_r, +d_axis)
    blocks_dn = shift_blocks(blocks_r, -d_axis)

    def B(which: str, idx: int) -> Block:
        """Retrieve a block from the original, shifted-up, or shifted-down set.

        Args:
            which: Block set identifier: ``'orig'``, ``'up'``, or ``'dn'``.
            idx: Block index within the set.

        Returns:
            Block: The requested block.
        """
        return {"orig": blocks_r, "up": blocks_up, "dn": blocks_dn}[which][idx]

    # 4) Helpers for adaptive tolerance
    def _median_inplane_spacing(face: Face, block: Block) -> float:
        """Compute the median edge length on a face in its plane.

        Calculates edge lengths along both in-plane directions of the face
        and returns the median value, which is used to set adaptive
        matching tolerances.

        Args:
            face: The face to measure.
            block: The block containing the face.

        Returns:
            float: Median in-plane edge length. Returns ``1.0`` if the
            face has no edges (degenerate).
        """
        I0,I1,J0,J1,K0,K1 = face.IMIN,face.IMAX,face.JMIN,face.JMAX,face.KMIN,face.KMAX
        X,Y,Z = block.X, block.Y, block.Z
        if face.const_type == 0:  # I const -> vary (J,K)
            i = I0
            x = X[i,J0:J1+1,K0:K1+1]; y = Y[i,J0:J1+1,K0:K1+1]; z = Z[i,J0:J1+1,K0:K1+1]
        elif face.const_type == 1:  # J const -> vary (I,K)
            j = J0
            x = X[I0:I1+1,j,K0:K1+1]; y = Y[I0:I1+1,j,K0:K1+1]; z = Z[I0:I1+1,j,K0:K1+1]
        else:  # K const -> vary (I,J)
            k = K0
            x = X[I0:I1+1,J0:J1+1,k]; y = Y[I0:I1+1,J0:J1+1,k]; z = Z[I0:I1+1,J0:J1+1,k]
        s = []
        if x.shape[0] > 1:
            dx = np.diff(x, axis=0); dy = np.diff(y, axis=0); dz = np.diff(z, axis=0)
            s.append(np.sqrt(dx*dx + dy*dy + dz*dz))
        if x.shape[1] > 1:
            dx = np.diff(x, axis=1); dy = np.diff(y, axis=1); dz = np.diff(z, axis=1)
            s.append(np.sqrt(dx*dx + dy*dy + dz*dz))
        if not s: return 1.0
        return float(np.median(np.concatenate([v.ravel() for v in s])))

    def _pair_tol(fA: Face, fB: Face) -> float:
        """Compute adaptive absolute tolerance for a face pair.

        Uses approximately 3% of the larger face's median in-plane spacing,
        with a floor of ``1e-4``. If ``node_tol_xyz`` was provided to the
        outer function, that value is returned directly.

        Args:
            fA: First face of the candidate pair.
            fB: Second face of the candidate pair.

        Returns:
            float: Absolute coordinate tolerance for node matching.
        """
        if node_tol_xyz is not None:
            return float(node_tol_xyz)
        sA = _median_inplane_spacing(fA, B("orig", fA.BlockIndex))
        sB = _median_inplane_spacing(fB, B("orig", fB.BlockIndex))
        # ~3% of local in-plane spacing; floor at 1e-4 (tune if needed)
        return max(0.03 * max(sA, sB), 1e-4)

    # 5) General orthogonal-plane precheck (works for x/y/z periodicity)
    def _orthogonal_precheck(fA: Face, fB: Face, bA: Block, bB: Block,
                             d_axis_local: float, tol: float, axis_local: str) -> bool:
        """Perform a fast orthogonal-plane overlap test between two faces.

        Shifts ``fA`` along the periodic axis by ``d_axis_local``, then
        projects both faces onto the plane orthogonal to that axis. Grid
        points are quantized by ``tol`` and set-intersected. The test
        passes if at least ``min_shared_abs`` (or ``min_shared_frac`` of
        the smaller face) quantized points overlap.

        Args:
            fA: First face (will be shifted).
            fB: Second face (unchanged).
            bA: Block containing ``fA``.
            bB: Block containing ``fB``.
            d_axis_local: Distance to shift ``fA`` along the periodic axis.
            tol: Quantization tolerance for the orthogonal projection.
            axis_local: Periodic axis: ``'x'``, ``'y'``, or ``'z'``.

        Returns:
            bool: ``True`` if sufficient orthogonal-plane overlap exists.
        """
        PA = fA.grid_points(bA, stride_u=1, stride_v=1)
        PB = fB.grid_points(bB, stride_u=1, stride_v=1)
        if PA.size == 0 or PB.size == 0:
            return False

        if axis_local == "x":
            PA[:,0] += d_axis_local
            projA, projB = PA[:,1:], PB[:,1:]          # (y,z)
        elif axis_local == "y":
            PA[:,1] += d_axis_local
            projA, projB = PA[:,[0,2]], PB[:,[0,2]]    # (x,z)
        else:  # "z"
            PA[:,2] += d_axis_local
            projA, projB = PA[:,:2], PB[:,:2]          # (x,y)

        QA = np.round(projA / tol).astype(np.int64)
        QB = np.round(projB / tol).astype(np.int64)
        if not QA.flags["C_CONTIGUOUS"]: QA = np.ascontiguousarray(QA)
        if not QB.flags["C_CONTIGUOUS"]: QB = np.ascontiguousarray(QB)
        vA = QA.view([('', QA.dtype)] * QA.shape[1]).reshape(-1)
        vB = QB.view([('', QB.dtype)] * QB.shape[1]).reshape(-1)
        inter = np.intersect1d(vA, vB, assume_unique=False)
        return inter.size >= max(min_shared_abs, int(min_shared_frac * min(len(vA), len(vB))))

    # 6) Node-sharing matcher using per-pair tol + precheck
    def faces_match(fL: Face, fU: Face) -> Tuple[bool, str]:
        """Determine whether two faces are a translational periodic match.

        Tries multiple shift directions (lower-up, upper-down) with both
        orthogonal-plane prechecks and full node-based matching. Returns
        on the first successful match strategy.

        Args:
            fL: Lower-boundary face candidate.
            fU: Upper-boundary face candidate.

        Returns:
            tuple: A 2-tuple of:

                - **matched** (*bool*): ``True`` if the faces are periodic.
                - **mode** (*str*): Description of which matching strategy
                  succeeded (e.g., ``'lower_up_vs_upper_orig'``), or an
                  empty string if no match.
        """
        bl, bu = fL.BlockIndex, fU.BlockIndex
        tol_pair = _pair_tol(fL, fU)

        # Fast precheck on orthogonal plane (lower up vs upper orig)
        if _orthogonal_precheck(fL, fU, B("orig", bl), B("orig", bu), d_axis, tol_pair, axis):
            return True, f"{axis}_precheck_lower_up"

        # lower moved up vs upper orig
        if fL.touches_by_nodes(fU, B("up", bl), B("orig", bu),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "lower_up_vs_upper_orig"

        # lower orig vs upper moved down
        if fL.touches_by_nodes(fU, B("orig", bl), B("dn", bu),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "lower_orig_vs_upper_dn"

        # Symmetric guards
        if _orthogonal_precheck(fU, fL, B("orig", bu), B("orig", bl), d_axis, tol_pair, axis):
            return True, f"{axis}_precheck_upper_up"
        if fU.touches_by_nodes(fL, B("up", bu), B("orig", bl),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "upper_up_vs_lower_orig"
        if fU.touches_by_nodes(fL, B("orig", bu), B("dn", bl),
                               tol_xyz=tol_pair, min_shared_frac=min_shared_frac,
                               min_shared_abs=min_shared_abs, stride_u=stride_u, stride_v=stride_v):
            return True, "upper_orig_vs_lower_dn"

        return False, ""

    # 7) Index mapping
    def mapping_minmax(fA: Face, fB: Face) -> Dict[str, str]:
        """Determine the IJK index mapping direction between two periodic faces.

        For each axis (I, J, K), determines whether the min index of ``fA``
        maps to the min or max index of ``fB`` based on the closer alignment.

        Args:
            fA: First face of the periodic pair.
            fB: Second face of the periodic pair.

        Returns:
            Dict[str, str]: Mapping for each axis, e.g.,
            ``{'I': 'min->min', 'J': 'min->max', 'K': 'min->min'}``.
        """
        out = {}
        for ax in ("I","J","K"):
            Amin, Amax = getattr(fA, ax+"MIN"), getattr(fA, ax+"MAX")
            Bmin, Bmax = getattr(fB, ax+"MIN"), getattr(fB, ax+"MAX")
            if (Amin == Bmin) and (Amax == Bmax): out[ax] = "min->min"
            elif (Amin == Bmax) and (Amax == Bmin): out[ax] = "min->max"
            else:
                d_mm = abs(Amin-Bmin)+abs(Amax-Bmax)
                d_mM = abs(Amin-Bmax)+abs(Amax-Bmin)
                out[ax] = "min->min" if d_mm <= d_mM else "min->max"
        return out

    # 8) Greedy pairing and export
    lower_pool = list(dict.fromkeys(lower_faces_r))
    upper_pool = list(dict.fromkeys(upper_faces_r))
    periodic_pairs_r: List[Tuple[Face, Face, Dict[str,str]]] = []
    periodic_export: List[Dict[str, Dict[str,int]]] = []

    pbar = tqdm(list(lower_pool), desc="Translational periodicity", unit="face")
    for fL in pbar:
        for j, fU in enumerate(upper_pool):
            ok, mode = faces_match(fL, fU)
            if ok:
                pbar.set_description(
                    f"Translational match: block {fL.BlockIndex} <-> block {fU.BlockIndex}")
                pbar.set_postfix(matches=len(periodic_pairs_r) + 1, remaining=len(upper_pool) - 1)
                m = mapping_minmax(fL, fU)
                periodic_pairs_r.append((fL, fU, m))
                periodic_export.append({
                    "block1": {"block_index": fL.BlockIndex,
                               "lb": (fL.IMIN, fL.JMIN, fL.KMIN),
                               "ub": (fL.IMAX, fL.JMAX, fL.KMAX)},
                    "block2": {"block_index": fU.BlockIndex,
                               "lb": (fU.IMIN, fU.JMIN, fU.KMIN),
                               "ub": (fU.IMAX, fU.JMAX, fU.KMAX)},
                    "mapping": m,
                    "mode": mode
                    }) # type: ignore
                upper_pool.pop(j)
                break  # move to next fL

    # Free shifted block copies — no longer needed after matching
    del blocks_up, blocks_dn

    # 9) scale back up
    scale_face_dict_indices(periodic_export, gcd_to_use, nested_sides=["block1", "block2"])

    periodic_pairs: List[Tuple[Face, Face, Dict[str,str]]] = []
    for (fL, fU, m) in periodic_pairs_r:
        gL = fL.copy(); gU = fU.copy()
        gL.I *= gcd_to_use; gL.J *= gcd_to_use; gL.K *= gcd_to_use
        gU.I *= gcd_to_use; gU.J *= gcd_to_use; gU.K *= gcd_to_use
        periodic_pairs.append((gL, gU, m))

    # 10) remove periodic from outer_faces (keep 'id' on remaining)
    periodic_keys = set()
    for rec in periodic_export:
        for side in ("block1","block2"):
            bi = rec[side]["block_index"]
            key = (bi,) + tuple(rec[side]["lb"]) + tuple(rec[side]["ub"])
            periodic_keys.add(key)

    outer_faces_remaining = []
    for o in outer_faces:
        key = (o["block_index"],) + tuple(o["lb"]) + tuple(o["ub"])
        if key not in periodic_keys:
            outer_faces_remaining.append(o)

    return periodic_export, periodic_pairs, outer_faces_remaining

def linear_real_transform(face1:Face,face2:Face) -> Tuple:
    """Compute the rotation angle and matrix from one face to another.

    Calculates the rotation that transforms the diagonal vector of ``face1``
    into the diagonal vector of ``face2``. This is useful for verifying
    whether two faces within the same block are periodic. Assumes the
    rotation axis is in the x-direction.

    Note:
        Based on the Linear Real Transform (LRT) computation from GlennHT.
        See ``M_ccMBMesh.F`` / ``computeLRT`` at
        https://gitlab.grc.nasa.gov/lte-turbo/GlennHT.

    Args:
        face1: Source face whose diagonal defines the "from" direction.
        face2: Target face whose diagonal defines the "to" direction.

    Returns:
        tuple: A 2-tuple containing:

            - **ang** (*float*): Rotation angle in radians. Positive
              follows the right-hand rule about the x-axis. Zero if the
              faces are already aligned.
            - **rotation_matrix** (*numpy.ndarray*): 3x3 rotation matrix.
              A zero matrix if no rotation is needed (``ang == 0``).
    """

    cTo3,cTo1 = face1.get_corners()
    cFrom3,cFrom1 = face2.get_corners()

    dTo  = np.array(cTo3).transpose() - np.array(cTo1).transpose()                      # difference in corner points = diagonal vector for Face 1
    ldTo=np.sqrt(np.sum(dTo*dTo))
    if ldTo > 0:
        dTo=dTo/ldTo

    dFrom = np.array(cFrom3).transpose() - np.array(cFrom1).transpose()                 # difference in corner points = diagonal vector for Face 2
    ldFrom = np.sqrt(np.sum(dFrom*dFrom))
    if( ldFrom > 0 ):
        dFrom=dFrom/ldFrom

    dotprod = np.sum(dTo * dFrom)

    if( abs(dotprod-1) < 1E-10 ): # Case of no rotation
        ang = 0

        rotation_matrix  = np.zeros(shape=(3,3))
    else:
        #Compute the angle of rotation
        cosAng=(dTo[1]*dFrom[1]+dTo[2]*dFrom[2])/sqrt(dTo[1]*dTo[1]+dTo[2]*dTo[2])/sqrt(dFrom[1]*dFrom[1]+dFrom[2]*dFrom[2])
        sinAng=(dTo[2]*dFrom[1]-dTo[1]*dFrom[2])/sqrt(dTo[1]*dTo[1]+dTo[2]*dTo[2])/sqrt(dFrom[1]*dFrom[1]+dFrom[2]*dFrom[2])
        ang=acos(cosAng)

        rotation_matrix = [ [1, 0, 0],
                            [0, cosAng, -sinAng],
                            [0, sinAng, cosAng] ]
        if( sinAng < 0 ):
            ang*=-1
    return ang, rotation_matrix

def __periodicity_check__(face1:Face, face2:Face,block1:Block,block2:Block,tol:float=1E-6):
    """Check two faces for periodic intersection and produce matched/split faces.

    Orders the faces so that the shorter-diagonal face is checked against the
    longer one, then calls ``get_face_intersection`` to find node-level matches.
    If at least 4 matching nodes are found, constructs periodic face pairs from
    the matched region and returns any split (leftover) faces.

    Args:
        face1: First face to check.
        face2: Second face to check.
        block1: Block corresponding to ``face1``. For rotational periodicity,
            this should be the rotated copy of the block.
        block2: Block corresponding to ``face2``.
        tol: Euclidean distance tolerance for node matching.
            Defaults to ``1e-6``.

    Returns:
        tuple: A 3-tuple containing:

            - **match_rows** (*List[Dict]*): List of point-match dictionaries
              from ``get_face_intersection``, each containing ``'i1'``,
              ``'j1'``, ``'k1'``, ``'i2'``, ``'j2'``, ``'k2'`` keys.
            - **periodic_faces** (*List[Face]*): The two faces forming the
              periodic pair (in the original argument order), or an empty
              list if no match was found.
            - **split_faces** (*List[Face]*): Leftover face fragments that
              were not part of the matched region. These should be treated
              as outer faces for subsequent matching passes.
    """

    periodic_faces = list()
    split_faces = list()
    swapped = False
    if (face2.diagonal_length < face1.diagonal_length): # switch so that face 2 is always longer
        face1, face2 = face2, face1
        block1, block2 = block2, block1
        swapped = True

    match_rows, split_face1, split_face2 = get_face_intersection(face1, face2, block1, block2, tol)

    if len(match_rows) >= 4:
        f1 = create_face_from_diagonals(block1,
            imin=min(m['i1'] for m in match_rows), jmin=min(m['j1'] for m in match_rows), kmin=min(m['k1'] for m in match_rows),
            imax=max(m['i1'] for m in match_rows), jmax=max(m['j1'] for m in match_rows), kmax=max(m['k1'] for m in match_rows))
        f1.set_block_index(face1.blockIndex)
        f1.set_face_id(face1.id)

        f2 = create_face_from_diagonals(block2,
            imin=min(m['i2'] for m in match_rows), jmin=min(m['j2'] for m in match_rows), kmin=min(m['k2'] for m in match_rows),
            imax=max(m['i2'] for m in match_rows), jmax=max(m['j2'] for m in match_rows), kmax=max(m['k2'] for m in match_rows))
        f2.set_block_index(face2.blockIndex)
        f2.set_face_id(face2.id)

        split_faces.extend(split_face1)
        split_faces.extend(split_face2)
        if swapped:
            periodic_faces.append(f2)
            periodic_faces.append(f1)
        else:
            periodic_faces.append(f1)
            periodic_faces.append(f2)

    return match_rows, periodic_faces, split_faces


def _extract_face_points_per(block: Block, rec: dict) -> np.ndarray:
    """Extract face points in diagonal traversal order (lb -> ub)."""
    il, jl, kl = rec['lb']
    ih, jh, kh = rec['ub']
    il = min(il, block.IMAX - 1); ih = min(ih, block.IMAX - 1)
    jl = min(jl, block.JMAX - 1); jh = min(jh, block.JMAX - 1)
    kl = min(kl, block.KMAX - 1); kh = min(kh, block.KMAX - 1)

    def _directed(start, end):
        if start <= end: return range(start, end + 1)
        else: return range(start, end - 1, -1)

    pts = []
    for i in _directed(il, ih):
        for j in _directed(jl, jh):
            for k in _directed(kl, kh):
                pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
    return np.array(pts)


def verify_periodicity(blocks: List[Block], face_matches: list, theta: float,
                       rotation_axis: str = 'x', tol: float = 1E-6):
    """Verify that every grid point on each periodic face pair matches after rotation.

    For each match, extracts all grid points from both faces using the
    diagonal traversal convention (lb -> ub) on the full-resolution blocks,
    rotates face1 points by +/-theta, and checks element-by-element that
    every point pairs up within ``tol`` Euclidean distance.

    The diagonal convention means that points are extracted in lb->ub
    traversal order for both faces. Since matching faces share the same
    diagonal convention, no lexicographic sorting is needed -- points are
    compared directly in traversal order.

    Args:
        blocks (List[Block]): Full-resolution blocks.
        face_matches (list): List of face match dicts from
            ``rotated_periodicity``.  Each dict has ``'block1'`` and
            ``'block2'`` sub-dicts with ``'block_index'``, ``'lb'``
            (tuple), and ``'ub'`` (tuple).
        theta (float): Rotation angle in radians.
        rotation_axis (str): Axis of rotation: ``'x'``, ``'y'``, or ``'z'``.
        tol (float): Euclidean distance tolerance.  Defaults to 1e-6.

    Returns:
        tuple: A 2-tuple containing:

            - **verified** (list): Matches with all interior points coincident
              after rotation.
            - **mismatched** (list): Matches with at least one non-matching point.
    """
    rot_pos = create_rotation_matrix(theta, rotation_axis)
    rot_neg = create_rotation_matrix(-theta, rotation_axis)
    tol2 = tol * tol

    verified = []
    mismatched = []

    for fm in face_matches:
        b1 = fm['block1']
        b2 = fm['block2']

        if b1['block_index'] >= len(blocks) or b2['block_index'] >= len(blocks):
            mismatched.append(fm)
            continue

        pts1 = _extract_face_points_per(blocks[b1['block_index']], b1)
        pts2 = _extract_face_points_per(blocks[b2['block_index']], b2)

        if pts1.shape[0] != pts2.shape[0]:
            mismatched.append(fm)
            continue

        # Compare directly element-by-element (diagonal convention ensures matching order)
        found = False
        for rot in [rot_pos, rot_neg]:
            pts1_rot = pts1 @ rot.T  # Apply rotation
            dists2 = np.sum((pts1_rot - pts2) ** 2, axis=1)
            if np.all(dists2 < tol2):
                verified.append(fm)
                found = True
                break

        if not found:
            mismatched.append(fm)

    return verified, mismatched
