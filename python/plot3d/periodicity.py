from typing import List, Dict, Tuple, Optional
import warnings
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
from .utils import (euclidean_distance, enumerate_unique_corners, scale_face_dict_indices,
    face_key, divide_face_dict_indices)
from tqdm import tqdm

def _scale_dicts_down(dicts: List[Dict], gcd: int, keys=('IMIN','JMIN','KMIN','IMAX','JMAX','KMAX')):
    """Scale dictionary face indices down by GCD factor."""
    for d in dicts:
        for k in keys:
            if k in d:
                d[k] = int(d[k] / gcd)


def create_rotation_matrix(rotation_angle:float, rotation_axis:str="x"):
    """Creates a rotation matrix given an angle and axis 

    Args:
        rotation_angle (float): Rotation angle in radians
        rotation_axis (str, optional): Axis of rotation "x", "y", or "z". Defaults to "x".

    Returns:
        np.ndarray: 3x3 rotation matrix 
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
    """Cheap geometric pre-checks to reject obviously non-matching face pairs.

    Checks (in order of cost): same const_type, axial overlap,
    radial overlap, rotated centroid proximity.
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
    """Count how many corners of face_a, after rotation, land on face_b's grid.

    Ported from Rust ``count_rotated_corners_on_face``. Used as a cheap
    pre-check: if fewer than 2 rotated corners hit face_b, the expensive
    full intersection is skipped.

    Args:
        face_a: source face whose corners will be rotated
        face_b: target face to check against
        block_b: block containing face_b
        rotation_matrix: 3x3 rotation matrix
        tol: Euclidean distance tolerance

    Returns:
        Number of face_a corners (0-4) that match face_b grid points after rotation
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
    """Check if all corners of face_a after rotation match face_b's corners bijectively.

    Ported from Rust ``full_face_match_transformed``. Used in Phase 1 to quickly
    identify full-face periodic matches without expensive node-by-node intersection.

    Args:
        face_a: source face whose corners will be rotated
        face_b: target face to check against
        rotation_matrix: 3x3 rotation matrix
        tol: Euclidean distance tolerance

    Returns:
        True if all corners match bijectively after rotation
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
    """Core matching engine for rotated_periodicity().

    Uses a three-phase algorithm with FacePool theta bucketing for O(log N)
    candidate selection:

    **Phase 1**: Full-face matching using FacePool theta bucketing. For each
    face, find candidates at theta +/- rotation_angle, try full face match.

    **Phase 2**: Split-face matching with corner pre-check. Uses FacePool
    to find candidates, then checks if >= 2 rotated corners land on the
    candidate face before attempting expensive intersection. Iterates until
    no new matches are found (no iteration limit).

    **Phase 3**: Relaxed matching without pre-checks. Uses 5x tolerance
    to catch wavy-surface faces. Runs until convergence.

    Args:
        blocks: reduced blocks
        outer_faces_all: outer faces as Face objects (at reduced scale)
        matched_faces_all: matched faces as Face objects (at reduced scale)
        rotation_matrices: list of rotation matrices to try (1 or 2)
        rotation_axis: 'x', 'y', or 'z'
        periodic_direction: 'i', 'j', or 'k' to filter faces (None = any)

    Returns:
        (periodic_faces_export, outer_faces_export, periodic_faces, outer_faces_all)
    """
    from .face_pool import FacePool

    # Build rotation caches (one per matrix, lazy)
    caches = {}
    for idx, mat in enumerate(rotation_matrices):
        caches[idx] = {}

    def get_rotated(matrix_idx: int, block_index: int) -> Block:
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
    """Scale reduced-resolution results back to original grid resolution."""
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
    """Find rotational periodicity between outer faces by rotating blocks.

    Detects which outer faces match after rotation by the given angle about the
    specified axis. Automatically reduces the mesh by GCD for faster matching,
    uses angular boundary detection to narrow candidates, and applies geometric
    pre-checks to reject non-matching pairs cheaply.

    Tries both +angle and -angle rotations for robustness.

    Args:
        blocks (List[Block]): List of blocks (do not duplicate and pass in!).
        matched_faces (List[Dict[str,int]]): Matched faces from connectivity.
        outer_faces (List[Dict[str,int]]): Outer faces in dictionary form.
        rotation_angle (float, optional): Rotation angle in degrees. Either this
            or ``nblades`` must be provided. If both are given, ``rotation_angle``
            takes precedence.
        rotation_axis (str, optional): Axis of rotation: 'x', 'y', or 'z'. Defaults to 'x'.
        periodic_direction (str, optional): Filter to only check faces with a constant
            'i', 'j', or 'k' index. None means check all directions. Defaults to None.
        nblades (int, optional): Number of blades. Used to compute
            ``rotation_angle = 360.0 / nblades`` when ``rotation_angle`` is not provided.
        ReduceMesh (bool, optional): If True, reduces the mesh by GCD for faster matching.
            Defaults to True.

    Returns:
        (Tuple): containing

            - **periodic_faces_export** (List[Dict[str,int]]): periodic face pairs as dicts
            - **outer_faces_export** (List[Dict[str,int]]): remaining outer faces as dicts
            - **periodic_faces** (List[Tuple[Face,Face]]): periodic face pairs as Face objects
            - **outer_faces_all** (List[Face]): remaining outer faces as Face objects

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
    """
    Detect translational periodicity between block faces along a given axis.

    This function takes a set of outer block faces and attempts to identify 
    periodic counterparts across the domain in the specified translational 
    direction ('x', 'y', or 'z'). It works by:

      1. **Bounding faces:** Uses `find_bounding_faces` to identify candidate 
         lower/upper faces for the given axis.
      2. **Grid reduction:** Reduces blocks to their greatest common divisor 
         (GCD) resolution to make indexing consistent across blocks.
      3. **Shifting:** Creates shifted copies of all blocks in both positive 
         and negative directions along the periodic axis.
      4. **Precheck in orthogonal plane:** Uses a fast projection test 
         (orthogonal to the periodic axis) to determine whether two faces 
         could possibly match. This greatly reduces false negatives when 
         spacing/tolerances differ slightly.
      5. **Node-based match:** Calls `Face.touches_by_nodes` on candidate 
         pairs to check shared node positions, with an adaptive tolerance 
         based on the in-plane spacing of each face.
      6. **Pairing:** Records each valid pair of periodic faces, their 
         IJK index mappings (min→min or min→max), and removes matched faces 
         from the outer-face list.
      7. **Scaling back:** Rescales reduced indices back to the original grid 
         spacing so results are consistent with input block resolution.

    Args:
        blocks (List[Block]): List of blocks.
        outer_faces (List[Dict[str,int]]): Outer faces represented as 
            dictionaries (with IMIN, JMIN, KMIN, IMAX, JMAX, KMAX).
        delta (float, optional): Periodicity spacing along the chosen axis. 
            If None, it is inferred from the global block min/max extent.
        translational_direction (str, optional): Axis to check ('x','y','z'). 
            Default is 'z'.
        node_tol_xyz (float, optional): Absolute coordinate tolerance for 
            node-matching. If None, tolerance is computed adaptively based 
            on median in-plane spacing of candidate faces.
        min_shared_frac (float, optional): Minimum fraction of nodes that must 
            overlap for two faces to be considered periodic. Default 0.02.
        min_shared_abs (int, optional): Minimum absolute number of shared nodes. 
            Default 4.
        stride_u (int, optional): Subsampling stride along the first face index 
            direction. Default 1 (no skipping).
        stride_v (int, optional): Subsampling stride along the second face index 
            direction. Default 1 (no skipping).

    Returns:
        Tuple[
            List[Dict[str, Dict[str,int]]], 
            List[Tuple[Face, Face, Dict[str,str]]], 
            List[Dict[str,int]]
        ]:
            - **periodic_faces_export**: Export-ready dictionaries describing 
              each periodic pair (block indices, face extents, index mapping, 
              and match mode).
            - **periodic_pairs**: Matched periodic face pairs as `Face` objects 
              with IJK mapping.
            - **outer_faces_remaining**: Updated list of outer faces with 
              periodic ones removed (preserving any existing `id` fields).

    Notes:
        - Works for periodicity in **x**, **y**, or **z** directions.
        - The adaptive tolerance makes the method robust to small spacing 
          differences between blocks.
        - The orthogonal-plane precheck avoids expensive node comparisons 
          when faces clearly do not align.
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
        cp = [b.copy() for b in bb]
        for b in cp:
            b.shift(amount, axis)
        return cp

    blocks_up = shift_blocks(blocks_r, +d_axis)
    blocks_dn = shift_blocks(blocks_r, -d_axis)

    def B(which: str, idx: int) -> Block:
        return {"orig": blocks_r, "up": blocks_up, "dn": blocks_dn}[which][idx]

    # 4) Helpers for adaptive tolerance
    def _median_inplane_spacing(face: Face, block: Block) -> float:
        """Median edge length on the face (in-plane)."""
        I0,I1,J0,J1,K0,K1 = face.IMIN,face.IMAX,face.JMIN,face.JMAX,face.KMIN,face.KMAX
        X,Y,Z = block.X, block.Y, block.Z
        if face.const_type == 0:  # I const → vary (J,K)
            i = I0
            x = X[i,J0:J1+1,K0:K1+1]; y = Y[i,J0:J1+1,K0:K1+1]; z = Z[i,J0:J1+1,K0:K1+1]
        elif face.const_type == 1:  # J const → vary (I,K)
            j = J0
            x = X[I0:I1+1,j,K0:K1+1]; y = Y[I0:I1+1,j,K0:K1+1]; z = Z[I0:I1+1,j,K0:K1+1]
        else:  # K const → vary (I,J)
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
        """Adaptive absolute tolerance per pair (use global override if provided)."""
        if node_tol_xyz is not None:
            return float(node_tol_xyz)
        sA = _median_inplane_spacing(fA, B("orig", fA.BlockIndex))
        sB = _median_inplane_spacing(fB, B("orig", fB.BlockIndex))
        # ~3% of local in-plane spacing; floor at 1e-4 (tune if needed)
        return max(0.03 * max(sA, sB), 1e-4)

    # 5) General orthogonal-plane precheck (works for x/y/z periodicity)
    def _orthogonal_precheck(fA: Face, fB: Face, bA: Block, bB: Block,
                             d_axis_local: float, tol: float, axis_local: str) -> bool:
        """
        Shift face A along 'axis_local' by d_axis_local, then compare projections onto the
        orthogonal plane within tolerance. Requires both absolute and fractional overlap.
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
                               "IMIN": fL.IMIN, "JMIN": fL.JMIN, "KMIN": fL.KMIN,
                               "IMAX": fL.IMAX, "JMAX": fL.JMAX, "KMAX": fL.KMAX},
                    "block2": {"block_index": fU.BlockIndex,
                               "IMIN": fU.IMIN, "JMIN": fU.JMIN, "KMIN": fU.KMIN,
                               "IMAX": fU.IMAX, "JMAX": fU.JMAX, "KMAX": fU.KMAX},
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
            key = (bi, rec[side]["IMIN"], rec[side]["JMIN"], rec[side]["KMIN"],
                        rec[side]["IMAX"], rec[side]["JMAX"], rec[side]["KMAX"])
            periodic_keys.add(key)

    outer_faces_remaining = []
    for o in outer_faces:
        key = (o["block_index"], o["IMIN"], o["JMIN"], o["KMIN"],
                               o["IMAX"], o["JMAX"], o["KMAX"])
        if key not in periodic_keys:
            outer_faces_remaining.append(o)

    return periodic_export, periodic_pairs, outer_faces_remaining

def linear_real_transform(face1:Face,face2:Face) -> Tuple:
    """Computes the rotation angle from Face1 to Face2. This can be used to check if the faces are periodic 
        This function assumes the rotation axis is in the "x" direction. This is good for faces within the same block 

    Reference:
        - Linear Real Transforms from GlennHT https://gitlab.grc.nasa.gov/lte-turbo/GlennHT/-/blob/master/src/M_ccMBMesh.F See computeLRT
        
    Args:
        Face1 (Face): Face to rotate
        Face2 (Face): Face to rotate to

    Returns:
        (tuple): tuple containing:

            - **ang** (float): rotation angle
            - **rotation_matrix** (numpy.ndarray): Rotation matrix 3x3
        
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
    """General function to find periodicity within a given block.

    Steps:
        - 1: Take the face with the shorter diagonal.
        - 2: Rotate the shorter face by angle 360/nblades.
        - 3: Check to see if faces intersect

    Args:
        face1 (Face): An arbitrary face
        face2 (Face): An arbitrary face
        block1 (Block): block 1 cooresponding to face 1
        block2 (Block): block 2 cooresponding to face 2

    Returns:
        (tuple): containing

            - **match_rows** (List[Dict]): List of point matches for periodic surfaces
            - **periodic_surface** (List[Face]):  These are faces that are periodic
            - **split_surfaces** (List[Face]): Some blocks may have periodic faces with other blocks. But the faces may need to be split so say you pair a small face with a larger face. The split surfaces should be treated as an outer face

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


def periodicity(blocks, outer_faces, matched_faces, periodic_direction='k',
                rotation_axis='x', nblades=55):
    """Deprecated. Use rotated_periodicity() instead.

    Note: argument order differs from rotated_periodicity —
    this function takes (blocks, outer_faces, matched_faces)
    while rotated_periodicity takes (blocks, matched_faces, outer_faces).
    """
    warnings.warn(
        "periodicity() is deprecated. Use rotated_periodicity(blocks, matched_faces, "
        "outer_faces, nblades=N, periodic_direction='k') instead.",
        DeprecationWarning, stacklevel=2)
    return rotated_periodicity(blocks, matched_faces, outer_faces,
                               nblades=nblades, rotation_axis=rotation_axis,
                               periodic_direction=periodic_direction)


def periodicity_fast(blocks, outer_faces, matched_faces, periodic_direction='k',
                     rotation_axis='x', nblades=55):
    """Deprecated. Use rotated_periodicity() instead.

    Note: argument order differs from rotated_periodicity —
    this function takes (blocks, outer_faces, matched_faces)
    while rotated_periodicity takes (blocks, matched_faces, outer_faces).
    """
    warnings.warn(
        "periodicity_fast() is deprecated. Use rotated_periodicity(blocks, matched_faces, "
        "outer_faces, nblades=N, periodic_direction='k') instead.",
        DeprecationWarning, stacklevel=2)
    return rotated_periodicity(blocks, matched_faces, outer_faces,
                               nblades=nblades, rotation_axis=rotation_axis,
                               periodic_direction=periodic_direction,
                               ReduceMesh=True)


def verify_periodicity(blocks: List[Block], face_matches: list, theta: float, rotation_axis: str = 'x', tol: float = 1E-6):
    """Verifies that the diagonal corners of periodic face_matches are spatially
    consistent after rotating block1 by ±theta about the given axis.

    For each face_match, rotates block1 by +theta (and if needed -theta) and checks
    that the rotated lower/upper corner coordinates match block2's lower/upper
    corners within tolerance. If the stored diagonal doesn't match, tries all
    permutations of block2's face corners. If a valid permutation is found,
    the face_match is corrected and added to the verified list.

    Uses GCD reduction for efficient coordinate lookups.

    Args:
        blocks (List[Block]): List of all blocks (original full-resolution)
        face_matches (list): List of face_match dicts from periodicity or rotated_periodicity
        theta (float): Rotation angle in degrees
        rotation_axis (str, optional): Axis of rotation 'x', 'y', or 'z'. Defaults to 'x'.
        tol (float, optional): Euclidean distance tolerance. Defaults to 1E-6.

    Returns:
        (list): verified face_matches whose diagonals are confirmed or corrected
        (list): mismatched face_matches where no corner permutation matched
    """
    # Compute GCD and reduce blocks
    gcd_to_use = compute_gcd(blocks)
    reduced_blocks = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # Build rotation matrices for +theta and -theta
    rotation_matrix_pos = create_rotation_matrix(radians(theta), rotation_axis)
    rotation_matrix_neg = create_rotation_matrix(radians(-theta), rotation_axis)

    # Pre-rotate all reduced blocks in both directions
    rotated_blocks_pos = [rotate_block(b, rotation_matrix_pos) for b in reduced_blocks]
    rotated_blocks_neg = [rotate_block(b, rotation_matrix_neg) for b in reduced_blocks]

    # Scale down face_matches indices by GCD
    scaled_matches = [{k: (dict(v) if isinstance(v, dict) else v) for k, v in fm.items()} for fm in face_matches]
    divide_face_dict_indices(scaled_matches, gcd_to_use, nested_sides=['block1', 'block2'])

    verified = list()
    mismatched = list()

    for idx in range(len(scaled_matches)):
        fm = scaled_matches[idx]
        b1 = fm['block1']
        b2 = fm['block2']
        b1_idx = b1['block_index']
        b2_idx = b2['block_index']
        block2 = reduced_blocks[b2_idx]

        # Enumerate unique corners of block2's face
        unique_corners = enumerate_unique_corners(
            [b2['IMIN'], b2['IMAX']], [b2['JMIN'], b2['JMAX']], [b2['KMIN'], b2['KMAX']])

        found = False
        best_d_lower = float('inf')
        best_d_upper = float('inf')

        # Try +theta rotation first, then -theta
        for rotated_blocks in [rotated_blocks_pos, rotated_blocks_neg]:
            if found:
                break

            block1_rotated = rotated_blocks[b1_idx]

            # Block1 rotated diagonal coordinates
            x1_l = block1_rotated.X[b1['IMIN'], b1['JMIN'], b1['KMIN']]
            y1_l = block1_rotated.Y[b1['IMIN'], b1['JMIN'], b1['KMIN']]
            z1_l = block1_rotated.Z[b1['IMIN'], b1['JMIN'], b1['KMIN']]

            x1_u = block1_rotated.X[b1['IMAX'], b1['JMAX'], b1['KMAX']]
            y1_u = block1_rotated.Y[b1['IMAX'], b1['JMAX'], b1['KMAX']]
            z1_u = block1_rotated.Z[b1['IMAX'], b1['JMAX'], b1['KMAX']]

            # Check stored diagonal first
            x2_l = block2.X[b2['IMIN'], b2['JMIN'], b2['KMIN']]
            y2_l = block2.Y[b2['IMIN'], b2['JMIN'], b2['KMIN']]
            z2_l = block2.Z[b2['IMIN'], b2['JMIN'], b2['KMIN']]

            x2_u = block2.X[b2['IMAX'], b2['JMAX'], b2['KMAX']]
            y2_u = block2.Y[b2['IMAX'], b2['JMAX'], b2['KMAX']]
            z2_u = block2.Z[b2['IMAX'], b2['JMAX'], b2['KMAX']]

            d_lower = euclidean_distance(x1_l, y1_l, z1_l, x2_l, y2_l, z2_l)
            d_upper = euclidean_distance(x1_u, y1_u, z1_u, x2_u, y2_u, z2_u)

            if d_lower < best_d_lower:
                best_d_lower = d_lower
            if d_upper < best_d_upper:
                best_d_upper = d_upper

            if d_lower < tol and d_upper < tol:
                verified.append(face_matches[idx])
                found = True
                break

            # Try all permutations of block2's corners
            for corner_lower in unique_corners:
                for corner_upper in unique_corners:
                    if corner_lower == corner_upper:
                        continue

                    il, jl, kl = corner_lower
                    iu, ju, ku = corner_upper

                    dl = euclidean_distance(x1_l, y1_l, z1_l,
                                           block2.X[il, jl, kl], block2.Y[il, jl, kl], block2.Z[il, jl, kl])
                    du = euclidean_distance(x1_u, y1_u, z1_u,
                                           block2.X[iu, ju, ku], block2.Y[iu, ju, ku], block2.Z[iu, ju, ku])

                    if dl < best_d_lower:
                        best_d_lower = dl
                    if du < best_d_upper:
                        best_d_upper = du

                    if dl < tol and du < tol:
                        corrected = {k: (dict(v) if isinstance(v, dict) else v) for k, v in face_matches[idx].items()}
                        corrected['block2']['IMIN'] = il * gcd_to_use
                        corrected['block2']['JMIN'] = jl * gcd_to_use
                        corrected['block2']['KMIN'] = kl * gcd_to_use
                        corrected['block2']['IMAX'] = iu * gcd_to_use
                        corrected['block2']['JMAX'] = ju * gcd_to_use
                        corrected['block2']['KMAX'] = ku * gcd_to_use
                        verified.append(corrected)
                        found = True
                        break
                if found:
                    break

        if not found:
            orig = face_matches[idx]
            b1_orig = orig['block1']
            b2_orig = orig['block2']
            print(f"verify_periodicity: MISMATCH at face_match index {idx}")
            print(f"  block1 (block_index={b1_orig['block_index']}): "
                  f"lower=({b1_orig['IMIN']},{b1_orig['JMIN']},{b1_orig['KMIN']}) "
                  f"upper=({b1_orig['IMAX']},{b1_orig['JMAX']},{b1_orig['KMAX']})")
            print(f"  block2 (block_index={b2_orig['block_index']}): "
                  f"lower=({b2_orig['IMIN']},{b2_orig['JMIN']},{b2_orig['KMIN']}) "
                  f"upper=({b2_orig['IMAX']},{b2_orig['JMAX']},{b2_orig['KMAX']})")
            print(f"  Closest rotated block1 corner dist to block2 lower: {best_d_lower:.6e}")
            print(f"  Closest rotated block1 corner dist to block2 upper: {best_d_upper:.6e}")
            mismatched.append(face_matches[idx])

    return verified, mismatched
