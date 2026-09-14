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


# ---------------------------------------------------------------------
# A5: full-resolution re-validation after GCD reduction
#
# `rotated_periodicity` GCD-reduces `blocks` for speed, then finds periodic
# pairs on the coarse grid. A proposal that looks valid there can fail once
# every full-resolution node is checked (an interior node perturbed
# independently of its neighbours, invisible to the reduced grid). These
# tests cover:
#   1. A no-op regression on the known-good VSPT mesh -- re-validation must
#      not change anything when everything genuinely certifies.
#   2. A constructed case where a full-resolution-only interior node is
#      perturbed beyond tolerance -- the pair must be demoted from BOTH the
#      dict-form (`periodic_faces_export`/`outer_faces_export`) AND the
#      parallel Face-object-form (`periodic_faces`/`outer_faces_all`)
#      return channels, with a RuntimeWarning.
# ---------------------------------------------------------------------

import warnings as _warnings

from plot3d.block import Block
from plot3d.blockfunctions import compute_min_gcd
from plot3d.facefunctions import get_outer_faces
from plot3d.periodicity import rotated_periodicity


@skip_no_mesh
def test_rotated_periodicity_a5_no_op_on_known_good_mesh(mesh_data):
    """A5 re-validation must be a no-op on the known-good VSPT mesh: GCD
    reduction must actually be in play (gcd > 1) for this to be meaningful,
    and every reduced-grid proposal must survive full-resolution
    certification -- no RuntimeWarning, and the dict-form/Face-form channels
    must agree on counts (nothing silently demoted from one but not both).
    """
    blocks, face_matches, periodic_export_baseline, _ = mesh_data
    gcd = compute_min_gcd(blocks)
    assert gcd > 1, "fixture must exercise GCD reduction for this check to be meaningful"

    from plot3d import connectivity_fast
    _, outer_faces = connectivity_fast(blocks)
    rotation_angle_deg = 360.0 / NBLADES

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        periodic_export, outer_export, periodic_faces, outer_faces_all = rotated_periodicity(
            blocks, face_matches, outer_faces,
            rotation_angle=rotation_angle_deg, rotation_axis=ROTATION_AXIS,
        )

    demotions = [w for w in caught
                 if issubclass(w.category, RuntimeWarning)
                 and "rotated_periodicity" in str(w.message)]
    assert not demotions, f"unexpected demotions on known-good mesh: {[str(w.message) for w in demotions]}"
    assert len(periodic_export) == len(periodic_faces)
    assert len(outer_export) == len(outer_faces_all)
    assert len(periodic_export) == len(periodic_export_baseline)
    assert len(periodic_export) > 0


def _wedge_block(nx: int = 3, ntheta: int = 5, nr: int = 3,
                  rotation_deg: float = 20.0) -> Block:
    """A single-block angular wedge sector spanning exactly `rotation_deg`
    about the x-axis, self-periodic between its J=0 and J=ntheta-1 faces
    (rotating J=0 forward by `rotation_deg` about x lands exactly on
    J=ntheta-1). `(nx-1, ntheta-1, nr-1) = (2, 4, 2)` so
    `compute_min_gcd` == 2 -- GCD reduction skips the i=1/k=1 interior
    layers, leaving them full-resolution-only.
    """
    from math import radians
    i_idx = np.arange(nx, dtype=float)
    j_idx = np.arange(ntheta, dtype=float)
    k_idx = np.arange(nr, dtype=float)
    I, J, K = np.meshgrid(i_idx, j_idx, k_idx, indexing="ij")
    theta_full = radians(rotation_deg)
    theta = J * (theta_full / (ntheta - 1))
    r = 1.0 + K
    X = I
    Y = r * np.cos(theta)
    Z = r * np.sin(theta)
    return Block(X, Y, Z)


def test_rotated_periodicity_a5_demotes_full_resolution_only_perturbation():
    """A full-resolution-only interior node (i=1, skipped by GCD=2
    reduction) is perturbed beyond tolerance on the J=ntheta-1 face. The
    coarse (reduced) grid only samples i in {0, 2}, so it never sees the
    perturbation and finds the pair periodic; full-resolution
    re-validation must catch it and demote the pair from both return
    channels, with a RuntimeWarning naming `rotated_periodicity`.
    """
    rotation_deg = 20.0
    block = _wedge_block(rotation_deg=rotation_deg)
    assert compute_min_gcd([block]) == 2

    # Perturb the i=1 (full-resolution-only) node on the J=ntheta-1 face,
    # well beyond the default tol=1e-4.
    block.Y[1, -1, 0] += 0.05

    faces, _ = get_outer_faces(block)
    for f in faces:
        f.set_block_index(0)
    outer_faces = [f.to_dict() for f in faces]

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        periodic_export, outer_export, periodic_faces, outer_faces_all = rotated_periodicity(
            [block], matched_faces=[], outer_faces=outer_faces,
            rotation_angle=rotation_deg, rotation_axis="x",
            ReduceMesh=True, tol=1e-4,
        )

    demotions = [w for w in caught
                 if issubclass(w.category, RuntimeWarning)
                 and "rotated_periodicity" in str(w.message)]
    assert demotions, "expected a RuntimeWarning demoting the perturbed pair"

    # Dict-form channel: the pair must not appear as periodic.
    assert periodic_export == []
    # Face-object-form channel must agree: also empty, and both faces of
    # the demoted pair must have been appended to outer_faces_all instead
    # -- this is the dual-channel bookkeeping under test.
    assert periodic_faces == []
    j_const_outer = [f for f in outer_faces_all if f.JMIN == f.JMAX]
    assert len(j_const_outer) == 2, (
        "expected both J-constant faces of the demoted pair to land in "
        f"the Face-object outer_faces_all channel, got {len(j_const_outer)}"
    )
    # Dict-form outer channel must agree in count with the Face-form one.
    assert len(outer_export) == len(outer_faces_all)
