"""Tests for the flatmesh module (``plot3d.flatmesh``: ``flatten_mesh`` +
``FlatMesh``).

Fixture: the tutorial's 2-block annular wedge (see the ``wedge_block`` cell
in ``colab/Plot3D_Flatten.ipynb``) -- two 15x9x13-node
blocks spanning x in [0,0.5]/[0.5,1], r in [0.2,0.3], theta in [0,20deg].
20 degrees is exactly 1/18 of a full annulus (``nblades=18``), which is why
each block's own K=1/K=KMAX faces self-match under a +/-20deg rotation
(``rotated_periodicity``) instead of needing a second copy of the mesh.
``connectivity()``/``rotated_periodicity()`` run in well under a second on
this fixture (~1700 nodes/block, confirmed empirically), so nothing here is
gated on mesh size the way the real-stator round-trip in
``test_plot3d_flatten_deck.py`` is.

Volume/area correctness is deliberately NOT checked against the smooth
annular-sector analytic formula -- that formula assumes a smooth circular
arc, which is unreachable on a 13-point-per-block faceted theta
discretization. Instead:

- total volume is cross-checked against an INDEPENDENT 6-tet decomposition
  of every hex cell, built straight from ``points``/``cell_vertices`` (a
  completely different code path than ``Block.cell_volumes()``'s
  Davies & Salmond formula);
- the "cells are watertight" property (the key check) is verified by
  accumulating each cell's *own* signed face-area vectors -- rotating a
  periodic neighbor's ``face_area`` back into the owner's frame before
  subtracting -- and requiring the closed-surface sum to vanish, which is
  the divergence-theorem identity that must hold for every valid
  hexahedron, periodic or not.
"""

import copy
import math
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from plot3d import Block
from plot3d.connectivity import connectivity
from plot3d.periodicity import create_rotation_matrix, rotated_periodicity
from plot3d.glennht import Plot3DFlattenInletBC, Plot3DFlattenWallBC, tag_surfaces_geometric
from plot3d.flatmesh import (
    BC_INLET,
    BC_INTERIOR,
    BC_PERIODIC,
    BC_ROTATING_WALL,
    BC_WALL,
    FlatMesh,
    flatten_mesh,
)


# ---------------------------------------------------------------------------
# Fixture: 2-block annular wedge (matches
# colab/Plot3D_Flatten.ipynb's `wedge_block` cell).
# ---------------------------------------------------------------------------

DTHETA_DEG = 20.0
DTHETA_RAD = math.radians(DTHETA_DEG)
NBLADES = 18  # round(360 / 20)


def wedge_block(x0, x1, nx, r_hub=0.2, r_shroud=0.3, nr=9, dtheta=DTHETA_RAD, nt=13):
    x = np.linspace(x0, x1, nx)
    r = np.linspace(r_hub, r_shroud, nr)
    th = np.linspace(0.0, dtheta, nt)
    X = np.zeros((nx, nr, nt))
    Y = np.zeros_like(X)
    Z = np.zeros_like(X)
    for i in range(nx):
        for j in range(nr):
            for k in range(nt):
                X[i, j, k] = x[i]
                Y[i, j, k] = r[j] * np.cos(th[k])
                Z[i, j, k] = r[j] * np.sin(th[k])
    return Block(X, Y, Z)


def _make_wedge_blocks():
    # Two blocks stacked along x, sharing the x=0.5 face.
    return [wedge_block(0.0, 0.5, 15), wedge_block(0.5, 1.0, 15)]


@pytest.fixture(scope="module")
def wedge_blocks():
    return _make_wedge_blocks()


@pytest.fixture(scope="module")
def wedge_connectivity(wedge_blocks):
    """``(face_matches, periodic_faces, outer_faces, surface_ids)`` for the
    wedge, computed ONCE for the whole module.

    ``flatten_mesh`` never mutates any of its list/dict inputs (verified by
    reading its source), so every test below reuses these objects directly
    and read-only -- except the multi-surface test, which needs to
    reassign some ``outer_faces[].id`` values and therefore
    ``copy.deepcopy``s the list first, exactly as ``tag_surfaces_*`` tests
    do in ``test_plot3d_flatten_deck.py``.
    """
    face_matches, outer_faces = connectivity(wedge_blocks)
    periodic_faces, outer_faces, _, _ = rotated_periodicity(
        wedge_blocks, face_matches, outer_faces,
        rotation_angle=DTHETA_DEG, rotation_axis="x",
    )
    outer_faces, surface_ids = tag_surfaces_geometric(wedge_blocks, outer_faces, axis="x")
    return face_matches, periodic_faces, outer_faces, surface_ids


@pytest.fixture(scope="module")
def periodicity_meta():
    return {"rotation_angle_rad": DTHETA_RAD, "rotation_axis": "x"}


@pytest.fixture(scope="module")
def wedge_flatmesh(wedge_blocks, wedge_connectivity, periodicity_meta):
    """The default flattened wedge (no custom ``bcs``, so hub/shroud/blade
    resolve to plain ``wall`` via the name-based fallback). Module-scoped:
    every test that uses it only reads ``fm.<field>`` arrays, never
    mutates them in place, so sharing one instance across tests is safe
    and keeps the (already-fast) pipeline from re-running per test.
    """
    face_matches, periodic_faces, outer_faces, surface_ids = wedge_connectivity
    return flatten_mesh(
        wedge_blocks, face_matches, outer_faces,
        periodic_faces=periodic_faces, periodicity=dict(periodicity_meta),
        surface_ids=dict(surface_ids),
    )


def test_fixture_sanity(wedge_blocks, wedge_connectivity):
    assert [(b.IMAX, b.JMAX, b.KMAX) for b in wedge_blocks] == [(15, 9, 13), (15, 9, 13)]
    face_matches, periodic_faces, outer_faces, surface_ids = wedge_connectivity
    assert len(face_matches) == 1        # the shared x=0.5 interface
    assert len(periodic_faces) == 2      # each block's own K=1<->K=KMAX self-match
    ids_present = {f["id"] for f in outer_faces}
    assert ids_present == {1, 2, 4, 5}   # inlet, outlet, hub, shroud (no "blade" here)


# ---------------------------------------------------------------------------
# 1: independent volume (6-tet decomposition), never the smooth-sector formula
# ---------------------------------------------------------------------------

# 6 tets sharing the (0,6) space diagonal of the VTK_HEXAHEDRON corner order
# that flatmesh.py's `_CORNERS` / `cell_vertices` already use:
#   0=(i,j,k) 1=(i+1,j,k) 2=(i+1,j+1,k) 3=(i,j+1,k)
#   4=(i,j,k+1) 5=(i+1,j,k+1) 6=(i+1,j+1,k+1) 7=(i,j+1,k+1)
_HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def _independent_hex_volume_sum(points: np.ndarray, cell_vertices: np.ndarray) -> float:
    """Total volume of every hex cell via an independent 6-tet decomposition,
    V = |(B-A) . ((C-A) x (D-A))| / 6 per tet, summed in absolute value."""
    P = points[cell_vertices]  # (Nc, 8, 3)
    total = 0.0
    for a, b, c, d in _HEX_TETS:
        A, B, C, D = P[:, a], P[:, b], P[:, c], P[:, d]
        vol = np.abs(np.einsum("ij,ij->i", B - A, np.cross(C - A, D - A))) / 6.0
        total += float(vol.sum())
    return total


def test_independent_volume(wedge_flatmesh):
    fm = wedge_flatmesh
    independent_total = _independent_hex_volume_sum(fm.points, fm.cell_vertices)
    assert independent_total == pytest.approx(float(fm.cell_volume.sum()), rel=1e-9)
    assert np.all(fm.cell_volume > 0)


# ---------------------------------------------------------------------------
# 2: closed-cell sum(face_area) ~ 0, rotation-aware -- THE key check
# ---------------------------------------------------------------------------

def test_closed_cell_area_sum_is_rotation_aware(wedge_flatmesh):
    fm = wedge_flatmesh
    Nc = fm.cell_volume.shape[0]
    owner, neighbor, area, rot = (
        fm.face_owner, fm.face_neighbor, fm.face_area, fm.face_periodic_rotation,
    )

    has_neighbor = neighbor >= 0
    nb = neighbor[has_neighbor]
    S = area[has_neighbor]
    rot_nb = rot[has_neighbor]

    # Rotation-aware: a periodic face's neighbor cell is geometrically a
    # ROTATED copy, not the mirror image, of the owner's face -- so its own
    # outward area vector for that face is -(R @ S_owner), not just -S_owner
    # (R = identity, i.e. this collapses to the ordinary -S_owner, for
    # every non-periodic face, since face_periodic_rotation == 0 there).
    contrib = -S.copy()
    periodic_idx = np.nonzero(rot_nb != 0.0)[0]
    for i in periodic_idx:
        R = create_rotation_matrix(float(rot_nb[i]), "x")
        contrib[i] = -(R @ S[i])

    accum = np.zeros((Nc, 3))
    np.add.at(accum, owner, area)
    np.add.at(accum, nb, contrib)

    mean_area = float(np.mean(np.linalg.norm(area, axis=1)))
    max_resid = float(np.max(np.linalg.norm(accum, axis=1)))
    assert max_resid / mean_area < 1e-9

    # A NAIVE sum (no rotation correction) is NOT ~0 -- this confirms the
    # rotation correction is actually doing something (guards against a
    # regression where it silently became a no-op and both computations
    # "accidentally" agreed).
    assert periodic_idx.size > 0
    naive_accum = np.zeros((Nc, 3))
    np.add.at(naive_accum, owner, area)
    np.add.at(naive_accum, nb, -S)
    naive_max_resid = float(np.max(np.linalg.norm(naive_accum, axis=1)))
    assert naive_max_resid / mean_area > 1e-3


# ---------------------------------------------------------------------------
# 3: face_normal / face_area_mag invariant
# ---------------------------------------------------------------------------

def test_face_normal_and_area_mag_invariant(wedge_flatmesh):
    fm = wedge_flatmesh

    reconstructed = fm.face_normal * fm.face_area_mag[:, None]
    assert np.allclose(reconstructed, fm.face_area, atol=1e-10, rtol=1e-8)

    eps = 1e-12
    nonzero = fm.face_area_mag > eps
    assert np.all(nonzero), "no degenerate (zero-area) faces expected on this mesh"

    norms = np.linalg.norm(fm.face_normal[nonzero], axis=1)
    assert np.allclose(norms, 1.0, atol=1e-8)

    assert np.all(fm.face_area_mag > 0)


# ---------------------------------------------------------------------------
# 4: topology (owner/neighbor/boundary consistency)
# ---------------------------------------------------------------------------

def test_topology(wedge_flatmesh):
    fm = wedge_flatmesh
    Nc = fm.cell_volume.shape[0]
    Nf = fm.face_owner.shape[0]

    assert np.all(fm.face_owner >= 0) and np.all(fm.face_owner < Nc)

    interior_like = (fm.face_bc_type == BC_INTERIOR) | (fm.face_bc_type == BC_PERIODIC)
    has_neighbor = fm.face_neighbor >= 0
    # Interior, cross-block, and periodic faces always resolve a real
    # neighbor cell; true physical boundaries never do -- these two
    # partitions of the face set must coincide exactly.
    assert np.array_equal(interior_like, has_neighbor)
    assert np.all(fm.face_neighbor[has_neighbor] < Nc)

    boundary_mask = ~has_neighbor
    assert np.all(fm.face_surface_id[boundary_mask] >= 0)
    for sid in np.unique(fm.face_surface_id[boundary_mask]):
        assert int(sid) in fm.surface_ids

    # Counts are consistent: every face is either "has a neighbor"
    # (interior/cross-block/periodic) or "true boundary" -- exhaustively.
    n_with_neighbor = int(np.count_nonzero(has_neighbor))
    n_boundary = int(np.count_nonzero(boundary_mask))
    assert n_with_neighbor + n_boundary == Nf


# ---------------------------------------------------------------------------
# 5: periodic neighbors always resolved
# ---------------------------------------------------------------------------

def test_periodic_neighbors_always_resolved(wedge_flatmesh):
    fm = wedge_flatmesh
    periodic_mask = fm.face_periodic_rotation != 0.0
    assert np.count_nonzero(periodic_mask) > 0

    assert np.all(fm.face_neighbor[periodic_mask] >= 0)

    mags = np.abs(fm.face_periodic_rotation[periodic_mask])
    assert np.allclose(mags, DTHETA_RAD, atol=1e-9)


# ---------------------------------------------------------------------------
# 6: welding at the cross-block interface
# ---------------------------------------------------------------------------

def test_cross_block_welding(wedge_flatmesh):
    fm = wedge_flatmesh
    has_neighbor = fm.face_neighbor >= 0

    owner_block = fm.cell_block_id[fm.face_owner]
    neighbor_block = np.full(fm.face_neighbor.shape[0], -1, dtype=np.int64)
    neighbor_block[has_neighbor] = fm.cell_block_id[fm.face_neighbor[has_neighbor]]

    cross_mask = has_neighbor & (owner_block != neighbor_block)
    assert np.count_nonzero(cross_mask) > 0, "expected at least one cross-block interface face"

    for f in np.nonzero(cross_mask)[0]:
        owner_verts = set(int(v) for v in fm.cell_vertices[fm.face_owner[f]])
        neighbor_verts = set(int(v) for v in fm.cell_vertices[fm.face_neighbor[f]])
        for v in fm.face_vertices[f]:
            v = int(v)
            assert v in owner_verts, f"face {f} vertex {v} missing from owner cell_vertices"
            assert v in neighbor_verts, f"face {f} vertex {v} missing from neighbor cell_vertices"


# ---------------------------------------------------------------------------
# 7: node + face BC -- multiple surfaces (two inlets) & mixed walls
# ---------------------------------------------------------------------------

def test_multi_surface_and_mixed_wall_bcs(wedge_blocks, wedge_connectivity, periodicity_meta):
    face_matches, periodic_faces, outer_faces, surface_ids = wedge_connectivity
    outer_faces_custom = copy.deepcopy(outer_faces)

    # Split the single id=1 inlet face (block0's I=1 plane) in half along
    # theta into two distinct inlet surfaces: id=1 (main) and id=6 (a new
    # second inlet, e.g. a cooling-air feed) -- hub (4) and shroud (5) stay
    # untouched.
    inlet_entries = [f for f in outer_faces_custom if f["id"] == 1]
    assert len(inlet_entries) == 1
    inlet = inlet_entries[0]
    lb, ub = list(inlet["lb"]), list(inlet["ub"])
    assert lb[0] == ub[0], "expected an I-collapsed (constant-x) inlet face"
    k_lo, k_hi = sorted((lb[2], ub[2]))
    k_mid = (k_lo + k_hi) // 2
    assert k_lo < k_mid < k_hi

    outer_faces_custom.remove(inlet)
    outer_faces_custom.append({
        "block_index": inlet["block_index"],
        "lb": [lb[0], lb[1], k_lo], "ub": [lb[0], ub[1], k_mid],
        "id": 1,
    })
    outer_faces_custom.append({
        "block_index": inlet["block_index"],
        "lb": [lb[0], lb[1], k_mid], "ub": [lb[0], ub[1], k_hi],
        "id": 6,
    })

    bcs = [
        Plot3DFlattenInletBC(name="main_inlet", surfaces=[1],
                   total_pressure=350000.0, total_temperature=450.0),
        Plot3DFlattenInletBC(name="cooling_inlet", surfaces=[6],
                   total_pressure=360000.0, total_temperature=320.0),
        Plot3DFlattenWallBC(name="hub_wall", surfaces=[4], rotating=True, wall_rotation_rate=-6.0),
        Plot3DFlattenWallBC(name="shroud_wall", surfaces=[5]),
    ]

    fm = flatten_mesh(
        wedge_blocks, face_matches, outer_faces_custom,
        periodic_faces=periodic_faces, periodicity=dict(periodicity_meta),
        surface_ids=dict(surface_ids), bcs=bcs,
    )

    boundary_mask = fm.face_neighbor == -1
    sid_boundary = fm.face_surface_id[boundary_mask]
    bt_boundary = fm.face_bc_type[boundary_mask]

    # Both inlets keep distinct surface ids and are tagged `inlet`.
    for sid in (1, 6):
        m = sid_boundary == sid
        assert np.any(m), f"surface id {sid} missing from boundary faces"
        assert np.all(bt_boundary[m] == BC_INLET)

    # Hub is a rotating wall, shroud a stationary wall.
    m_hub = sid_boundary == 4
    assert np.any(m_hub)
    assert np.all(bt_boundary[m_hub] == BC_ROTATING_WALL)

    m_shroud = sid_boundary == 5
    assert np.any(m_shroud)
    assert np.all(bt_boundary[m_shroud] == BC_WALL)

    # A hub boundary node is tagged rotating_wall (rotating_wall has the
    # highest node-BC precedence, so this holds regardless of what else
    # that node's corner touches).
    hub_face_idx = int(np.nonzero(boundary_mask & (fm.face_surface_id == 4))[0][0])
    hub_node = int(fm.face_vertices[hub_face_idx][0])
    assert fm.point_bc_type[hub_node] == BC_ROTATING_WALL

    # boundary_conditions has one row per surface, and the hub row carries
    # wall_rotation_rate.
    for sid in (1, 6, 4, 5):
        assert sid in fm.boundary_conditions
    assert fm.boundary_conditions[4]["wall_rotation_rate"] == -6.0
    assert "wall_rotation_rate" not in fm.boundary_conditions[5]

    # A node lying on two surfaces has both ids in its CSR
    # (point_surf_off/point_surf_ids) entry.
    counts = fm.point_surf_off[1:] - fm.point_surf_off[:-1]
    multi_idx = np.nonzero(counts >= 2)[0]
    assert multi_idx.size > 0
    n0 = int(multi_idx[0])
    ids_n0 = [int(v) for v in fm.point_surf_ids[fm.point_surf_off[n0]:fm.point_surf_off[n0 + 1]]]
    assert len(ids_n0) == len(set(ids_n0)) >= 2


# ---------------------------------------------------------------------------
# 8: round-trip (hdf5 gated on h5py; npz always, loadable with bare numpy)
# ---------------------------------------------------------------------------

_ARRAY_FIELDS = (
    "points", "cell_volume", "cell_centroid", "cell_vertices", "cell_block_id", "cell_ijk",
    "face_owner", "face_neighbor", "face_area", "face_normal", "face_area_mag",
    "face_centroid", "face_vertices", "face_surface_id", "face_bc_type",
    "face_periodic_rotation", "face_periodic_translation",
    "point_bc_type", "point_surf_off", "point_surf_ids",
)
_DICT_FIELDS = ("surface_ids", "bc_type_legend", "boundary_conditions")


def test_hdf5_roundtrip(tmp_path, wedge_flatmesh):
    pytest.importorskip("h5py")
    fm = wedge_flatmesh
    path = tmp_path / "wedge.flatmesh.h5"
    fm.to_hdf5(str(path))
    fm2 = FlatMesh.from_hdf5(str(path))

    for name in _ARRAY_FIELDS:
        assert np.allclose(getattr(fm, name), getattr(fm2, name)), f"{name} mismatch after hdf5 round-trip"
    for name in _DICT_FIELDS:
        assert getattr(fm, name) == getattr(fm2, name), f"{name} mismatch after hdf5 round-trip"


def test_npz_roundtrip(tmp_path, wedge_flatmesh):
    fm = wedge_flatmesh
    path = tmp_path / "wedge.flatmesh.npz"
    fm.to_npz(str(path))

    # Loads with bare numpy.load -- no h5py import anywhere in this test.
    raw = np.load(str(path))
    for name in _ARRAY_FIELDS:
        assert name in raw.files
    assert "attributes_json" in raw.files

    fm2 = FlatMesh.from_npz(str(path))
    for name in _ARRAY_FIELDS:
        assert np.allclose(getattr(fm, name), getattr(fm2, name)), f"{name} mismatch after npz round-trip"
    for name in _DICT_FIELDS:
        assert getattr(fm, name) == getattr(fm2, name), f"{name} mismatch after npz round-trip"


# ---------------------------------------------------------------------------
# 9: to_vtu writes well-formed XML with the right cell/point counts
# ---------------------------------------------------------------------------

def test_to_vtu_writes_well_formed_xml(tmp_path, wedge_flatmesh):
    fm = wedge_flatmesh
    path = tmp_path / "wedge.flatmesh.vtu"
    fm.to_vtu(str(path))
    assert path.exists()

    tree = ET.parse(str(path))
    piece = tree.getroot().find(".//Piece")
    assert piece is not None
    assert int(piece.get("NumberOfPoints")) == fm.points.shape[0]
    assert int(piece.get("NumberOfCells")) == fm.cell_vertices.shape[0]



# ---------------------------------------------------------------------------
# 10: float32 grids weld at the default tolerance
#
# Regression: the default weld tolerance used to be a flat 1e-8. A float32
# grid stores each shared interface node twice, and the two copies can differ
# by an ulp -- ~2.5e-06 where coordinates reach 21, i.e. 250x that default.
# Those nodes then failed to weld, the interface node map came back
# incomplete, and _neighbor_grid died on a bare KeyError. The real VSPT
# stator (float32, coords to ~21) could not be flattened at all.
# ---------------------------------------------------------------------------

def _two_touching_blocks(dtype, x_offset=0.0):
    """Two blocks sharing their I=IMAX / I=1 plane, at an optional offset so
    the coordinate magnitude (and hence the float32 ulp) can be varied."""
    def blk(x0, x1):
        X, Y, Z = np.meshgrid(np.linspace(x0, x1, 5),
                              np.linspace(0.0, 1.0, 4),
                              np.linspace(0.0, 1.0, 3), indexing="ij")
        return Block(np.ascontiguousarray(X, dtype=dtype),
                     np.ascontiguousarray(Y, dtype=dtype),
                     np.ascontiguousarray(Z, dtype=dtype))
    return [blk(x_offset, x_offset + 1.0), blk(x_offset + 1.0, x_offset + 2.0)]


def test_auto_weld_tol_scales_with_dtype_and_magnitude():
    from plot3d.flatmesh import _auto_weld_tol

    # float64 stays at the historical floor -- previous behaviour preserved.
    assert _auto_weld_tol(_two_touching_blocks(np.float64, x_offset=20.0)) == 1e-8

    # float32 opens up in proportion to coordinate magnitude, and must clear
    # the ulp at that magnitude (which is what has to be bridged).
    tol = _auto_weld_tol(_two_touching_blocks(np.float32, x_offset=20.0))
    assert tol > 1e-8
    assert tol >= np.spacing(np.float32(22.0))


def test_flatten_float32_interface_welds_by_default():
    """A float32 interface whose two copies differ by one ulp must still
    weld -- and must produce the same graph as the exactly-coincident case."""
    blocks = _two_touching_blocks(np.float32, x_offset=20.0)
    face_matches, outer_faces = connectivity(blocks)

    exact = flatten_mesh(blocks, face_matches, outer_faces)

    # Nudge block1's copy of the shared plane by one ulp, as a real float32
    # grid written by a mesher would.
    blocks[1].X[0, :, :] = np.nextafter(blocks[1].X[0, :, :], np.float32(1e30))
    nudged = flatten_mesh(blocks, face_matches, outer_faces)

    assert nudged.points.shape == exact.points.shape, "one-ulp offset split the interface"
    # The interface is interior: both sides have a real neighbour.
    assert int((nudged.face_neighbor == -1).sum()) == int((exact.face_neighbor == -1).sum())


def test_incomplete_weld_raises_a_useful_error():
    """Too tight a tolerance must say so, not die on a bare KeyError."""
    blocks = _two_touching_blocks(np.float32, x_offset=20.0)
    face_matches, outer_faces = connectivity(blocks)
    blocks[1].X[0, :, :] = np.nextafter(blocks[1].X[0, :, :], np.float32(1e30))

    with pytest.raises(ValueError, match="weld_tol"):
        flatten_mesh(blocks, face_matches, outer_faces, weld_tol=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
