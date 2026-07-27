"""Test rotated periodicity with a real VSPT mesh (2 blocks, 55 blades).

Validates that:
1. rotated_periodicity finds periodic face pairs
2. Rotating face1 by the transformation matrix geometrically matches face2
3. Connectivity face pairs share matching xyz points
4. Permutation matrices correctly map index offsets between paired faces
"""

import os
import pytest
import numpy as np
from scipy.spatial import cKDTree

MESH_PATH = os.path.join(os.path.dirname(__file__), "data", "vspt_mesh_scaled.xyz")
NBLADES = 55
ROTATION_AXIS = "x"
TOL = 1e-4

skip_no_mesh = pytest.mark.skipif(
    not os.path.exists(MESH_PATH),
    reason="vspt_mesh_scaled.xyz not found",
)


def _load():
    """Read mesh, compute connectivity and rotated periodicity."""
    from math import radians
    from plot3d import (
        read_plot3D,
        connectivity_fast,
        rotated_periodicity,
        create_rotation_matrix,
    )

    blocks = read_plot3D(MESH_PATH)
    face_matches, outer_faces = connectivity_fast(blocks)
    rotation_angle_deg = 360.0 / NBLADES
    periodic_export, outer_export, periodic_pairs, _ = rotated_periodicity(
        blocks, face_matches, outer_faces,
        rotation_angle=rotation_angle_deg,
        rotation_axis=ROTATION_AXIS,
    )
    rotation_matrix = create_rotation_matrix(radians(rotation_angle_deg), ROTATION_AXIS)
    return blocks, face_matches, periodic_export, rotation_matrix


def _extract_face_points(block, lb, ub):
    """Extract all xyz points on a face as (N, 3). Handles reversed lb/ub."""
    imin, imax = min(lb[0], ub[0]), max(lb[0], ub[0])
    jmin, jmax = min(lb[1], ub[1]), max(lb[1], ub[1])
    kmin, kmax = min(lb[2], ub[2]), max(lb[2], ub[2])
    sl = (slice(imin, imax + 1), slice(jmin, jmax + 1), slice(kmin, kmax + 1))
    X = block.X[sl].ravel()
    Y = block.Y[sl].ravel()
    Z = block.Z[sl].ravel()
    return np.column_stack([X, Y, Z])


def _extract_face_points_with_ijk(block, lb, ub):
    """Extract xyz points and (i,j,k) indices respecting lb->ub traversal."""
    def axis_range(lo, hi):
        return list(range(lo, hi + 1)) if lo <= hi else list(range(lo, hi - 1, -1))

    pts, ijk = [], []
    for i in axis_range(lb[0], ub[0]):
        for j in axis_range(lb[1], ub[1]):
            for k in axis_range(lb[2], ub[2]):
                pts.append([block.X[i, j, k], block.Y[i, j, k], block.Z[i, j, k]])
                ijk.append([i, j, k])
    return np.array(pts), np.array(ijk)


# ── Fixtures ──

@pytest.fixture(scope="module")
def mesh_data():
    pytest.importorskip("plot3d")
    if not os.path.exists(MESH_PATH):
        pytest.skip("vspt_mesh_scaled.xyz not found")
    return _load()


# ── Tests ──

@skip_no_mesh
def test_rotated_periodicity_finds_pairs(mesh_data):
    """rotated_periodicity must find at least one periodic face pair."""
    _, _, periodic_export, _ = mesh_data
    assert len(periodic_export) > 0, "No periodic faces found"


@skip_no_mesh
def test_periodic_faces_geometry(mesh_data):
    """Rotating face1 by the transformation matrix must match face2 points."""
    blocks, _, periodic_export, rotation_matrix = mesh_data

    for idx, pf in enumerate(periodic_export):
        b1 = pf["block1"]
        b2 = pf["block2"]
        pts1 = _extract_face_points(blocks[b1["block_index"]], b1["lb"], b1["ub"])
        pts2 = _extract_face_points(blocks[b2["block_index"]], b2["lb"], b2["ub"])

        pts1_rotated = (rotation_matrix @ pts1.T).T
        tree = cKDTree(pts2)
        dists, _ = tree.query(pts1_rotated, k=1)

        n_matched = int(np.sum(dists < TOL))
        assert n_matched == len(pts1), (
            f"periodic[{idx}]: only {n_matched}/{len(pts1)} points matched "
            f"(max_dist={dists.max():.2e})"
        )


@skip_no_mesh
def test_connectivity_faces_geometry(mesh_data):
    """Connectivity face pairs must share matching xyz points."""
    blocks, face_matches, _, _ = mesh_data

    for idx, fm in enumerate(face_matches):
        b1 = fm["block1"]
        b2 = fm["block2"]
        pts1 = _extract_face_points(blocks[b1["block_index"]], b1["lb"], b1["ub"])
        pts2 = _extract_face_points(blocks[b2["block_index"]], b2["lb"], b2["ub"])

        tree = cKDTree(pts2)
        dists, _ = tree.query(pts1, k=1)

        n_matched = int(np.sum(dists < TOL))
        assert n_matched == len(pts1), (
            f"connectivity[{idx}]: only {n_matched}/{len(pts1)} points matched "
            f"(max_dist={dists.max():.2e})"
        )


@skip_no_mesh
def test_permutation_matrix_connectivity(mesh_data):
    """For connectivity: P maps face1 index offsets to face2, same xyz."""
    blocks, face_matches, _, _ = mesh_data

    pairs_with_orient = [
        fm for fm in face_matches
        if "orientation" in fm and "permutation_matrix" in fm["orientation"]
    ]
    assert len(pairs_with_orient) > 0, "No connectivity pairs with permutation matrix"

    for idx, pair in enumerate(pairs_with_orient):
        _check_permutation(blocks, pair, rotation_matrix=None, tol=TOL, label=f"conn[{idx}]")


@skip_no_mesh
def test_permutation_matrix_periodic(mesh_data):
    """For periodicity: P maps face1 index offsets to face2, with rotation applied."""
    blocks, _, periodic_export, rotation_matrix = mesh_data

    pairs_with_orient = [
        pf for pf in periodic_export
        if "orientation" in pf and "permutation_matrix" in pf["orientation"]
    ]
    assert len(pairs_with_orient) > 0, "No periodic pairs with permutation matrix"

    for idx, pair in enumerate(pairs_with_orient):
        _check_permutation(blocks, pair, rotation_matrix=rotation_matrix, tol=TOL, label=f"peri[{idx}]")


def _check_permutation(blocks, pair, rotation_matrix, tol, label):
    """Verify permutation matrix maps face1 index offsets to face2.

    For connectivity (rotation_matrix=None): face1 xyz == face2 xyz.
    For periodicity: R @ face1 xyz == face2 xyz.
    """
    b1 = pair["block1"]
    b2 = pair["block2"]
    lb1, ub1 = b1["lb"], b1["ub"]
    lb2, ub2 = b2["lb"], b2["ub"]
    P = np.array(pair["orientation"]["permutation_matrix"])
    block1 = blocks[b1["block_index"]]
    block2 = blocks[b2["block_index"]]

    pts1, ijk1 = _extract_face_points_with_ijk(block1, lb1, ub1)
    varying1 = [d for d in range(3) if lb1[d] != ub1[d]]
    varying2 = [d for d in range(3) if lb2[d] != ub2[d]]

    if len(varying1) != 2 or len(varying2) != 2:
        return  # skip degenerate faces

    n_matched = 0
    for i in range(len(ijk1)):
        u1 = ijk1[i, varying1[0]] - lb1[varying1[0]]
        v1 = ijk1[i, varying1[1]] - lb1[varying1[1]]
        uv2 = P @ np.array([u1, v1])

        ijk2 = list(lb2)
        ijk2[varying2[0]] = lb2[varying2[0]] + int(uv2[0])
        ijk2[varying2[1]] = lb2[varying2[1]] + int(uv2[1])

        shape = block2.X.shape
        if not all(0 <= ijk2[d] < shape[d] for d in range(3)):
            continue

        pt1 = pts1[i]
        if rotation_matrix is not None:
            pt1 = rotation_matrix @ pt1

        pt2 = np.array([
            block2.X[ijk2[0], ijk2[1], ijk2[2]],
            block2.Y[ijk2[0], ijk2[1], ijk2[2]],
            block2.Z[ijk2[0], ijk2[1], ijk2[2]],
        ])
        if np.linalg.norm(pt1 - pt2) < tol:
            n_matched += 1

    assert n_matched == len(pts1), (
        f"{label}: permutation matrix matched {n_matched}/{len(pts1)} points"
    )
