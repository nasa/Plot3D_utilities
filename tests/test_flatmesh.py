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

    # A dtype's own precision is not trusted as a bound below float32 --
    # a float64 array can still carry only single-precision-quantized data
    # (e.g. round-tripped through an ASCII grid file), so float64 gets the
    # same float32-scaled tolerance as float32 itself, not the bare 1e-8
    # floor. Both should be well above 1e-8 and roughly agree (same
    # coordinate magnitude drives both).
    tol64 = _auto_weld_tol(_two_touching_blocks(np.float64, x_offset=20.0))
    tol32 = _auto_weld_tol(_two_touching_blocks(np.float32, x_offset=20.0))
    assert tol64 > 1e-8
    assert tol32 > 1e-8
    assert tol64 == pytest.approx(tol32, rel=1e-9)

    # Must clear the actual float32 ulp at this magnitude, since that's
    # exactly the gap it exists to bridge.
    assert tol32 >= np.spacing(np.float32(22.0))

    # A tiny coordinate magnitude falls back to the absolute floor -- there's
    # no meaningful ulp to bridge worth opening the tolerance for.
    assert _auto_weld_tol(_two_touching_blocks(np.float64, x_offset=0.0)) == pytest.approx(
        max(1e-8, np.finfo(np.float32).eps * 2.0)
    )


def test_flatten_ascii_quantized_float64_interface_welds_by_default():
    """Regression: a *float64* array can still carry only single-precision
    fidelity if it was parsed from an ASCII grid file, whose fixed-width text
    formatting quantizes each coordinate independently of the runtime dtype.
    colab/VSPT_ASCII.xyz reads as float64 yet had 7 interface nodes 2.4e-07
    apart -- 24x the historical 1e-8 default -- which the dtype-only fix
    (float32 gets a looser tolerance, float64 does not) missed entirely."""
    blocks = _two_touching_blocks(np.float64, x_offset=20.0)
    face_matches, outer_faces = connectivity(blocks)

    exact = flatten_mesh(blocks, face_matches, outer_faces)

    # A few e-07 of quantization noise on the shared plane, as an ASCII
    # round-trip would introduce, while staying a float64 array throughout.
    blocks[1].X[0, :, :] = blocks[1].X[0, :, :] + 2.4e-7
    quantized = flatten_mesh(blocks, face_matches, outer_faces)

    assert quantized.points.shape == exact.points.shape, "quantization noise split the interface"
    assert int((quantized.face_neighbor == -1).sum()) == int((exact.face_neighbor == -1).sum())


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


# ---------------------------------------------------------------------------
# 11: welding is scoped to DECLARED interfaces, not a global radius
#
# Regression for the ECEF over-merge bug. `_auto_weld_tol` scales a float32
# noise floor by absolute coordinate magnitude (correctly -- ASCII/float32
# quantization is a fixed number of significant digits, so its absolute size
# grows with distance from the origin), which at ECEF scale (~6.4e6) gives
# weld_tol ~= 0.76. `_weld_all_points` used to feed that straight into ONE
# global `cKDTree.query_pairs(r=weld_tol)` over every point in the mesh, so
# every node of `_two_touching_blocks` -- real spacing 0.25-0.5 -- collapsed
# into a single point (verified against the pre-fix code: points.shape was
# (1, 3), 108 real nodes gone).
#
# The fix stops deciding "same node?" by distance alone: candidate pairs now
# come from `matched_faces`, i.e. interfaces `connectivity()` already vetted
# by shape/corner/orientation matching. A loose tolerance can then only blur
# nodes that are already declared to correspond.
# ---------------------------------------------------------------------------

_ECEF_OFFSET = 6.4e6  # metres from Earth's centre -- the reported failing scale


def test_ecef_scale_declared_interface_still_welds():
    """The reported bug, fixed: two blocks translated to ECEF scale whose
    shared interface `connectivity()` declares must flatten to exactly the
    same graph as the identical pair sitting at the origin -- both with a
    bit-identical interface and with realistic float64 round-off on it."""
    from plot3d.flatmesh import _auto_weld_tol

    reference = _two_touching_blocks(np.float64, x_offset=0.0)
    ref_matches, ref_outer = connectivity(reference)
    ref_fm = flatten_mesh(reference, ref_matches, ref_outer)
    # 5x4x3 nodes per block, the shared 4x3 I-plane welded once: 60+60-12.
    assert ref_fm.points.shape == (108, 3)

    blocks = _two_touching_blocks(np.float64, x_offset=_ECEF_OFFSET)
    face_matches, outer_faces = connectivity(blocks)
    assert len(face_matches) == 1, "the shared interface must still be detected"

    # The tolerance really is enormous in absolute terms at this magnitude --
    # far larger than the grid's own 0.25 spacing. That is exactly the input
    # that used to collapse the mesh, and must now be harmless.
    weld_tol = _auto_weld_tol(blocks)
    assert weld_tol > 0.5
    assert weld_tol > 0.25  # > the mesh's own node spacing

    exact = flatten_mesh(blocks, face_matches, outer_faces)
    assert exact.points.shape == ref_fm.points.shape
    assert exact.face_owner.shape == ref_fm.face_owner.shape
    assert int((exact.face_neighbor == -1).sum()) == int((ref_fm.face_neighbor == -1).sum())

    # ...and with a few ULPs of float64 round-off on the shared plane, the
    # noise a real translated grid actually carries at this magnitude
    # (np.spacing(6.4e6) ~ 9.3e-10).
    nudged_blocks = _two_touching_blocks(np.float64, x_offset=_ECEF_OFFSET)
    for _ in range(4):
        nudged_blocks[1].X[0, :, :] = np.nextafter(nudged_blocks[1].X[0, :, :], 1e30)
    nudged = flatten_mesh(nudged_blocks, face_matches, outer_faces)
    assert nudged.points.shape == ref_fm.points.shape, "ULP noise split the interface"
    assert int((nudged.face_neighbor == -1).sum()) == int((ref_fm.face_neighbor == -1).sum())


def test_ecef_scale_tolerance_bridges_noise_that_origin_scale_rejects():
    """The magnitude scaling is doing real work, not just riding along: a
    1e-06 interface discrepancy is well inside the ECEF-scale tolerance and
    well outside the origin-scale one."""
    at_origin = _two_touching_blocks(np.float64, x_offset=0.0)
    origin_matches, origin_outer = connectivity(at_origin)
    at_origin[1].X[0, :, :] += 1e-6
    with pytest.raises(ValueError, match="weld_tol"):
        flatten_mesh(at_origin, origin_matches, origin_outer)

    at_ecef = _two_touching_blocks(np.float64, x_offset=_ECEF_OFFSET)
    ecef_matches, ecef_outer = connectivity(at_ecef)
    at_ecef[1].X[0, :, :] += 1e-6
    assert flatten_mesh(at_ecef, ecef_matches, ecef_outer).points.shape == (108, 3)


def test_undeclared_nearby_surfaces_are_never_merged():
    """The failure mode no tolerance-only design could close: two surfaces
    that are geometrically close but that `connectivity()` never declared as
    a matching pair must stay distinct no matter how large weld_tol is,
    because they were never welding candidates in the first place.

    A blind global radius merges them by construction -- the pre-fix code
    reduced this fixture to points.shape == (1, 3) at weld_tol=1.0, and a
    *capped* tolerance would not have helped either, since the cap was
    derived from within-block structural spacing and says nothing about two
    unrelated blocks sitting near each other.
    """
    gap = 1e-3  # the two facing planes are 1e-03 apart -- 1000x under weld_tol
    blocks = _two_touching_blocks(np.float64, x_offset=0.0)
    blocks[1].X += gap  # slide block 1 off its neighbour: no shared interface

    face_matches, outer_faces = connectivity(blocks)
    assert face_matches == [], "premise: connectivity() declares no match here"

    n_nodes = sum(b.IMAX * b.JMAX * b.KMAX for b in blocks)
    assert n_nodes == 120

    for weld_tol in (1e-8, gap * 10, 1.0):
        fm = flatten_mesh(blocks, face_matches, outer_faces, weld_tol=weld_tol)
        assert fm.points.shape == (n_nodes, 3), (
            f"weld_tol={weld_tol} merged undeclared nodes"
        )
        # Every node really is geometrically distinct, so nothing was lost.
        assert np.unique(fm.cell_vertices).size == n_nodes


def test_periodic_partners_are_never_welded_even_when_coincident():
    """Periodic partners must stay distinct nodes. They are excluded
    structurally -- `_weld_all_points` only ever sees `matched_faces` -- not
    by happening to fall outside weld_tol, so even an *exactly coincident*
    pair declared as periodic survives as two separate points."""
    blocks = _two_touching_blocks(np.float64, x_offset=0.0)
    face_matches, outer_faces = connectivity(blocks)
    assert len(face_matches) == 1

    # Re-label the one genuine, exactly-coincident interface as periodic.
    periodic_faces = list(face_matches)
    fm = flatten_mesh(
        blocks, [], outer_faces,
        periodic_faces=periodic_faces,
        periodicity={"rotation_axis": "x"},
        weld_tol=1.0,
    )

    n_nodes = sum(b.IMAX * b.JMAX * b.KMAX for b in blocks)
    assert fm.points.shape == (n_nodes, 3), "periodic partners were welded"
    # The pairing itself still resolves (via the transform path), so the
    # interface is a real interior edge of the graph on both sides.
    assert int((fm.face_bc_type == BC_PERIODIC).sum()) > 0
    assert np.all(fm.face_neighbor[fm.face_bc_type == BC_PERIODIC] >= 0)


def test_ambiguous_face_internal_match_raises_instead_of_mismatching():
    """Scoping welding to a declared face pair keeps nearest-neighbour
    matching honest only while weld_tol stays well under that face's OWN
    in-plane node spacing. When it does not, the nearest face-2 node can be
    an adjacent node rather than the true twin -- that must fail loudly, not
    silently wire up the wrong correspondence.

    Here block 1's copy of the shared plane is slid half a node spacing
    (1/6) along j, so every interior face-1 node is exactly equidistant from
    two different face-2 nodes, and weld_tol (0.2) is comparable to the
    face's own spacing (1/3).
    """
    blocks = _two_touching_blocks(np.float64, x_offset=0.0)
    face_matches, outer_faces = connectivity(blocks)
    assert len(face_matches) == 1

    j_spacing = 1.0 / 3.0
    weld_tol = 0.6 * j_spacing  # 0.2 -- same order as the face's own spacing

    # Sanity: undisturbed, this interface welds cleanly even at that tolerance
    # (the twins are bit-identical, so the nearest neighbour is unambiguous no
    # matter how loose the radius is). The error below is caused by the
    # ambiguity, not by the loose tolerance on its own.
    assert flatten_mesh(
        blocks, face_matches, outer_faces, weld_tol=weld_tol
    ).points.shape == (108, 3)

    blocks[1].Y[0, :, :] += j_spacing / 2.0

    with pytest.raises(ValueError, match="ambiguous interface weld") as excinfo:
        flatten_mesh(blocks, face_matches, outer_faces, weld_tol=weld_tol)
    assert "weld_tol" in str(excinfo.value)

    # A tolerance tight enough that neither candidate is reachable reports the
    # other, pre-existing failure instead -- still a clear ValueError.
    with pytest.raises(ValueError, match="weld_tol"):
        flatten_mesh(blocks, face_matches, outer_faces, weld_tol=1e-9)


# ---------------------------------------------------------------------------
# 12: the individual guards inside `_weld_face_pair`.
#
# Section 11 exercises the guards only through `flatten_mesh`, where several
# of them shadow each other -- deleting any single one still leaves another
# raising on the same fixture, so those tests cannot tell them apart. These
# drive `_weld_face_pair` directly, one guard per test, so that removing or
# loosening one is caught.
# ---------------------------------------------------------------------------

def _plane(nu, nv, h=1.0):
    """(nu*nv, 3) grid of nodes on x=0 with in-plane spacing `h`."""
    u, v = np.meshgrid(np.arange(nu) * h, np.arange(nv) * h, indexing="ij")
    return np.stack([np.zeros(u.size), u.ravel(), v.ravel()], axis=-1)


def test_weld_face_pair_accepts_an_exact_interface():
    from plot3d.flatmesh import _weld_face_pair

    P = _plane(4, 3)
    partner = _weld_face_pair(P, P.copy(), 1e-9, "A", "B")
    assert np.array_equal(partner, np.arange(P.shape[0]))

    # Node order on the two sides is unrelated -- the map is geometric.
    perm = np.random.default_rng(0).permutation(P.shape[0])
    partner = _weld_face_pair(P, P[perm], 1e-9, "A", "B")
    assert np.allclose(P[perm][partner], P)


def test_weld_face_pair_exempts_round_off_copies_but_not_nearby_distinct_nodes():
    """The `candidate_sep > dup_tol` qualifier on the ambiguity guard, where
    `dup_tol` is a *round-off* radius (`_coincidence_tol`), not weld_tol.

    Two stored copies of one node (a collapsed/degenerate node) are exempt
    from the tie check because `_weld_all_points` will union them anyway.
    Anything further apart than round-off is two nodes as far as this module
    can tell, so a tie between them must raise -- even when both sit well
    inside weld_tol. (An earlier version exempted anything within weld_tol;
    that let a pole whose copies differ by more than round-off weld to an
    arbitrary copy per side and silently mis-wire the pole-adjacent cells'
    neighbours, see `test_pole_on_declared_interface_...` below.)
    """
    from plot3d.flatmesh import _weld_face_pair

    # (a) exempt: face 2 stores one node twice, bit-identically. With an
    # exact duplicate the ratio test alone lets it through (d1 < 2*d0 is
    # 0 < 0, false) -- so (a') is the case that actually exercises the
    # `candidate_sep > dup_tol` qualifier.
    P1 = np.array([[0.0, 0.1, 0.0]])
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert _weld_face_pair(P1, P2, 0.5, "A", "B")[0] in (0, 1)

    # (a') exempt: the two copies differ by a couple of float64 ULPs, the
    # accumulated round-off of e.g. r*cos(theta) evaluated at r == 0.
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 2 * np.spacing(0.1), 0.0]])
    assert _weld_face_pair(P1, P2, 0.5, "A", "B")[0] in (0, 1)

    # (b) NOT exempt: 1e-09 apart. Far inside weld_tol, but far outside
    # round-off -- two such nodes cannot be welded deterministically.
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 1e-9, 0.0]])
    with pytest.raises(ValueError, match="ambiguous interface weld") as excinfo:
        _weld_face_pair(P1, P2, 0.5, "A", "B")
    assert "round-off" in str(excinfo.value)

    # (c) not exempt: the two candidates are 1.0 apart -- genuinely different
    # nodes -- and the face-1 node sits exactly between them.
    P1 = np.array([[0.0, 0.5, 0.0]])
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    with pytest.raises(ValueError, match="ambiguous interface weld"):
        _weld_face_pair(P1, P2, 0.9, "A", "B")


def test_coincidence_tol_is_round_off_scaled_and_capped_by_weld_tol():
    from plot3d.flatmesh import _coincidence_tol, _COINCIDENT_ULPS

    P = np.array([[6.4e6, 0.0, 0.0]])
    tol = _coincidence_tol(0.76, P)
    assert tol == pytest.approx(_COINCIDENT_ULPS * np.spacing(6.4e6))
    assert tol < 1e-7  # ~1.5e-08: seven orders of magnitude under weld_tol
    # never looser than weld_tol itself
    assert _coincidence_tol(1e-20, P) == 1e-20
    # at the origin there is still a positive radius, so -0.0 / +0.0 group
    assert _coincidence_tol(1.0, np.zeros((2, 3))) > 0.0


def test_coincident_copy_pairs_groups_exact_and_round_off_copies_only():
    from plot3d.flatmesh import _coincident_copy_pairs

    def groups(P, tol):
        a, b = _coincident_copy_pairs(P, tol)
        parent = list(range(P.shape[0]))

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for i, j in zip(a.tolist(), b.tolist()):
            parent[find(i)] = find(j)
        return len({find(i) for i in range(P.shape[0])})

    h = 1.0
    P = np.array([[0, 0, 0], [0, 0, 0], [0, -0.0, np.spacing(1.0)],  # 3 copies
                  [0, h, 0], [0, 2 * h, 0]], dtype=float)              # 2 real nodes
    assert groups(P, 4 * np.spacing(1.0)) == 3
    # at tol == 0 only bit-identical rows group (the 1-ULP copy is its own
    # group); the real nodes a spacing apart are never chained either way
    assert groups(P, 0.0) == 4
    # a fully collapsed face of many bit-identical rows yields a linear
    # number of pairs (grouped by sort), not n^2 radius-query pairs
    big = np.zeros((5000, 3))
    a, _ = _coincident_copy_pairs(big, 1e-12)
    assert a.size == 4999
    assert groups(big, 1e-12) == 1


def test_weld_face_pair_ambiguity_margin_is_the_documented_factor():
    """`_WELD_AMBIGUITY_MARGIN` is 2.0: a runner-up less than 2x the winner's
    distance is a tie, one at more than 2x is a clear win. Pins the constant
    from both sides, so raising or lowering it fails here."""
    from plot3d.flatmesh import _weld_face_pair, _WELD_AMBIGUITY_MARGIN

    assert _WELD_AMBIGUITY_MARGIN == 2.0

    def ratio_case(r):
        # Winner at d, runner-up at r*d. weld_tol is set above both distances
        # (so `d1 <= weld_tol` holds) but below the candidates' separation (so
        # they count as genuinely different nodes) -- leaving the ratio as the
        # only thing that decides.
        d = 0.2
        # The second face-1 node sits exactly on the runner-up, so both face-2
        # nodes are reached and the surjectivity guard stays out of the way.
        P1 = np.array([[0.0, d, 0.0], [0.0, d + r * d, 0.0]])
        P2 = np.array([[0.0, 0.0, 0.0], [0.0, d + r * d, 0.0]])
        return _weld_face_pair(P1, P2, 0.5, "A", "B")

    with pytest.raises(ValueError, match="ambiguous interface weld"):
        ratio_case(1.9)  # inside the margin -> a tie
    assert ratio_case(2.1)[0] == 0  # outside it -> accepted


def test_weld_face_pair_guard_boundaries_are_exact():
    """The comparison operators themselves, on values chosen to be exactly
    representable so the boundary is hit dead-on rather than approached."""
    from plot3d.flatmesh import _weld_face_pair

    # `d1 <= weld_tol`: a runner-up sitting exactly ON the tolerance still
    # counts as a candidate. d0 = 5/16, d1 = 1/2 = weld_tol, ratio 1.6 < 2.
    P1 = np.array([[0.0, 0.3125, 0.0], [0.0, 0.8125, 0.0]])
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 0.8125, 0.0]])
    with pytest.raises(ValueError, match="ambiguous interface weld"):
        _weld_face_pair(P1, P2, 0.5, "A", "B")

    # `d1 < margin * d0`: a runner-up at exactly 2x the winner is a clear
    # enough win to accept. d0 = 1/4, d1 = 1/2 = 2 * d0 exactly.
    P1 = np.array([[0.0, 0.25, 0.0], [0.0, 0.75, 0.0]])
    P2 = np.array([[0.0, 0.0, 0.0], [0.0, 0.75, 0.0]])
    assert _weld_face_pair(P1, P2, 0.5, "A", "B").tolist() == [0, 1]

    # `d0 > weld_tol`: a gap exactly at the tolerance welds, one past it does
    # not. Pins that the acceptance radius really is weld_tol and not a
    # multiple of it.
    at = np.array([[0.0, 0.0, 0.0]])
    assert _weld_face_pair(at, np.array([[0.0, 0.5, 0.0]]), 0.5, "A", "B").tolist() == [0]
    with pytest.raises(ValueError, match="did not weld"):
        _weld_face_pair(at, np.array([[0.0, 0.5625, 0.0]]), 0.5, "A", "B")


def test_weld_face_pair_catches_a_whole_face_shifted_past_half_a_spacing():
    """On a real interface the guards are stronger than the ratio test alone.

    `_weld_face_pair`'s docstring notes that an isolated node pair whose noise
    exceeds 2/3 of the spacing gets a wrong partner with a confident-looking
    ratio. On a *face*, that cannot pass quietly: a wrong pick either collides
    with another face-1 node's partner (injectivity) or orphans a face-2 node
    (surjectivity). This pins that the layered guards cover the whole range
    above h/2, so no plausible shift welds silently wrong.
    """
    from plot3d.flatmesh import _weld_face_pair

    h = 1.0
    P2 = _plane(5, 4, h)
    for e in (0.0, 0.1 * h, 0.3 * h):  # below h/2: correct, and accepted
        P1 = P2.copy()
        P1[:, 1] += e
        assert np.array_equal(_weld_face_pair(P1, P2, 0.9 * h, "A", "B"),
                              np.arange(P2.shape[0])), f"e={e} mis-welded"

    # Above h/2 every shift is refused. Which guard fires varies with the
    # shift (ratio, then injectivity, then the plain out-of-tolerance check
    # once the far edge runs off the end of face 2) -- what matters is that
    # none of them lets a wrong correspondence through.
    for e in (0.45 * h, 0.5 * h, 0.6 * h, 0.7 * h, 0.9 * h):
        P1 = P2.copy()
        P1[:, 1] += e
        with pytest.raises(ValueError, match="ambiguous interface weld|did not weld"):
            _weld_face_pair(P1, P2, 0.9 * h, "A", "B")

    # Non-uniform noise behaves the same way.
    rng = np.random.default_rng(0)
    P1 = P2 + np.concatenate(
        [np.zeros((P2.shape[0], 1)), rng.uniform(-0.5 * h, 0.5 * h, (P2.shape[0], 2))],
        axis=1,
    )
    with pytest.raises(ValueError, match="ambiguous interface weld"):
        _weld_face_pair(P1, P2, 0.9 * h, "A", "B")


def test_weld_face_pair_rejects_a_many_to_one_collapse():
    """The injectivity guard: distinct face-1 nodes must not all weld onto one
    face-2 node unless they are themselves coincident."""
    from plot3d.flatmesh import _weld_face_pair

    # Two face-1 nodes 0.4 apart, one face-2 node between them.
    P1 = np.array([[0.0, 0.0, 0.0], [0.0, 0.4, 0.0]])
    P2 = np.array([[0.0, 0.2, 0.0]])
    with pytest.raises(ValueError, match="ambiguous interface weld"):
        _weld_face_pair(P1, P2, 0.3, "A", "B")

    # ...but coincident face-1 copies of one node are fine (degenerate node).
    P1 = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert np.array_equal(_weld_face_pair(P1, np.array([[0.0, 0.0, 0.0]]),
                                          1e-9, "A", "B"), [0, 0])


def test_weld_face_pair_rejects_unreached_face2_nodes():
    """The surjectivity guard. face1 -> face2 can be perfectly one-to-one
    while face 2 carries extra nodes nothing welds to; those would otherwise
    land in `points` as silent duplicates with no error raised anywhere."""
    from plot3d.flatmesh import _weld_face_pair

    P1 = _plane(3, 3)
    P2 = np.vstack([P1, [[0.0, 10.0, 10.0]]])  # one surplus, far away
    with pytest.raises(ValueError, match="no counterpart"):
        _weld_face_pair(P1, P2, 1e-9, "A", "B")

    # A surplus node coincident with a reached one is the degenerate case and
    # stays legal.
    P2 = np.vstack([P1, P1[0]])
    assert _weld_face_pair(P1, P2, 1e-9, "A", "B").shape == (P1.shape[0],)


def test_weld_face_pair_single_node_face():
    """The k == 1 path: a face-2 side with exactly one node has no runner-up
    to compare against, so the ambiguity guard must be skipped, not crash."""
    from plot3d.flatmesh import _weld_face_pair

    one = np.array([[1.0, 2.0, 3.0]])
    assert _weld_face_pair(one, one.copy(), 1e-9, "A", "B").tolist() == [0]
    with pytest.raises(ValueError, match="did not weld"):
        _weld_face_pair(np.array([[0.0, 0.0, 0.0]]), one, 1e-9, "A", "B")


def test_within_block_self_match_seam_is_welded():
    """An O-grid branch cut arrives from `connectivity()` as a `matched_faces`
    entry whose two sides name the SAME block. Welding is keyed on global flat
    indices, so it must merge exactly like a cross-block interface -- this is
    the within-block case real turbomachinery meshes actually depend on (the
    shipped VSPT mesh has one).
    """
    from plot3d.flatmesh import _weld_all_points

    # One block wrapped a full 360 degrees: its k=0 and k=KMAX planes are the
    # same physical surface.
    nr, nt, nx = 3, 9, 2
    r = np.linspace(1.0, 2.0, nr)
    th = np.linspace(0.0, 2.0 * np.pi, nt)  # th[-1] == th[0] (mod 2pi)
    x = np.linspace(0.0, 1.0, nx)
    X3, R3, T3 = np.meshgrid(x, r, th, indexing="ij")
    blk = Block(np.ascontiguousarray(X3),
                np.ascontiguousarray(R3 * np.cos(T3)),
                np.ascontiguousarray(R3 * np.sin(T3)))
    blocks = [blk]
    seam = [{"block1": {"block_index": 0, "lb": [0, 0, 0], "ub": [nx - 1, nr - 1, 0]},
             "block2": {"block_index": 0, "lb": [0, 0, nt - 1], "ub": [nx - 1, nr - 1, nt - 1]}}]

    n_nodes = nx * nr * nt
    points, ng = _weld_all_points(blocks, seam, 1e-9)
    seam_nodes = nx * nr
    assert points.shape == (n_nodes - seam_nodes, 3), "the branch cut did not weld"
    assert np.array_equal(ng[0][:, :, 0], ng[0][:, :, nt - 1])

    # Without the declaration nothing merges -- the seam is welded because it
    # was declared, not because the two planes happen to coincide.
    points_nodecl, _ = _weld_all_points(blocks, [], 1e-9)
    assert points_nodecl.shape == (n_nodes, 3)


def test_within_block_pole_line_is_documented_not_welded():
    """Known, deliberate scope limit (see `_weld_all_points`): a collapsed
    edge inside a block is a line, never a declared face, so its coincident
    nodes each keep their own id. Asserted so the behaviour is pinned rather
    than discovered downstream -- if a within-block collapse pass is ever
    added, this test is the one to update.
    """
    from plot3d.flatmesh import _weld_all_points

    nx, nr, nt = 2, 3, 4
    r = np.linspace(0.0, 1.0, nr)  # r == 0 at j == 0 -> a pole line
    th = np.linspace(0.0, 0.5, nt)
    X3, R3, T3 = np.meshgrid(np.linspace(0.0, 1.0, nx), r, th, indexing="ij")
    blk = Block(np.ascontiguousarray(X3),
                np.ascontiguousarray(R3 * np.cos(T3)),
                np.ascontiguousarray(R3 * np.sin(T3)))

    points, ng = _weld_all_points([blk], [], 1e-9)
    assert points.shape == (nx * nr * nt, 3)
    # The pole nodes really are coincident, and really are distinct ids.
    pole = np.array([[blk.X[0, 0, k], blk.Y[0, 0, k], blk.Z[0, 0, k]] for k in range(nt)])
    assert np.allclose(pole, pole[0])
    assert len({int(ng[0][0, 0, k]) for k in range(nt)}) == nt


# ---------------------------------------------------------------------------
# 13: a pole line lying IN a declared interface plane.
#
# Two blocks with an axis singularity (r == 0 at j == 0) stacked along x
# share an i-plane that contains nt stored copies of the pole point on EACH
# side. `connectivity()` declares this interface like any other. Before the
# coincident-copy union in `_weld_all_points`, every side-1 copy welded to
# whichever single side-2 copy cKDTree happened to return (k == 1 here);
# `_node_map_via_weld` then mapped all of them to that copy and
# `_neighbor_grid`'s min-over-corners resolved two of the four pole-ring
# interface faces to the WRONG block-2 cell -- silently, with a different
# wrong answer depending on which side was `block1`. The old global weld
# never had the problem because all copies shared one id.
# ---------------------------------------------------------------------------

def _stacked_pole_blocks(nx=3, nr=3, nt=5, noise=0.0, noisy_blocks=(0, 1)):
    def cyl(x0, x1, noisy):
        x = np.linspace(x0, x1, nx)
        r = np.linspace(0.0, 1.0, nr)  # pole at j == 0
        th = np.linspace(0.0, 0.5 * np.pi, nt)
        X3, R3, T3 = np.meshgrid(x, r, th, indexing="ij")
        Y = R3 * np.cos(T3)
        Z = R3 * np.sin(T3)
        if noise and noisy:
            # spread the pole copies so they are no longer identical
            Y[:, 0, :] += noise * np.arange(nt)
        return Block(np.ascontiguousarray(X3), np.ascontiguousarray(Y), np.ascontiguousarray(Z))
    return [cyl(0.0, 0.5, 0 in noisy_blocks), cyl(0.5, 1.0, 1 in noisy_blocks)]


def _every_interior_face_belongs_to_both_its_cells(fm):
    for f in np.flatnonzero(fm.face_neighbor >= 0):
        fv = set(fm.face_vertices[f].tolist())
        assert fv <= set(fm.cell_vertices[fm.face_owner[f]].tolist()), int(f)
        assert fv <= set(fm.cell_vertices[fm.face_neighbor[f]].tolist()), int(f)


def _unordered_interior_edges(fm):
    interior = np.flatnonzero(fm.face_neighbor >= 0)
    return {frozenset((int(fm.face_owner[f]), int(fm.face_neighbor[f]))) for f in interior}


def test_pole_on_declared_interface_resolves_the_same_graph_from_either_side():
    from plot3d.flatmesh import _weld_all_points

    nx, nr, nt = 3, 3, 5
    blocks = _stacked_pole_blocks(nx, nr, nt)
    matches, outer = connectivity(blocks)
    assert len(matches) == 1, "premise: connectivity() declares the shared plane"
    swapped = [{"block1": matches[0]["block2"], "block2": matches[0]["block1"]}]

    graphs = []
    for mf in (matches, swapped):
        fm = flatten_mesh(blocks, mf, outer, weld_tol=1e-9)
        _every_interior_face_belongs_to_both_its_cells(fm)
        graphs.append(fm)

    # The pole point on the interface is ONE id, shared by both sides, in
    # both orientations; the pole copies one node-layer away stay distinct
    # (the documented within-block scope limit).
    expected_points = 2 * (nx - 1) * nr * nt + (nr - 1) * nt + 1
    for fm in graphs:
        assert fm.points.shape == (expected_points, 3)
    assert graphs[0].points.shape == graphs[1].points.shape
    assert _unordered_interior_edges(graphs[0]) == _unordered_interior_edges(graphs[1])

    for mf in (matches, swapped):
        _, ng = _weld_all_points(blocks, mf, 1e-9)
        ids = {int(ng[0][nx - 1, 0, k]) for k in range(nt)} | {int(ng[1][0, 0, k]) for k in range(nt)}
        assert len(ids) == 1, ids


def test_pole_copies_differing_by_round_off_are_still_one_id():
    """Copies that differ by accumulated float64 round-off (not bit-identical)
    are grouped too -- the coincidence radius is ULP-scaled, not zero."""
    from plot3d.flatmesh import _weld_all_points

    nx, nr, nt = 3, 3, 5
    blocks = _stacked_pole_blocks(nx, nr, nt, noise=np.spacing(1.0))
    matches, outer = connectivity(blocks)
    assert len(matches) == 1
    fm = flatten_mesh(blocks, matches, outer, weld_tol=1e-9)
    _every_interior_face_belongs_to_both_its_cells(fm)
    _, ng = _weld_all_points(blocks, matches, 1e-9)
    ids = {int(ng[0][nx - 1, 0, k]) for k in range(nt)} | {int(ng[1][0, 0, k]) for k in range(nt)}
    assert len(ids) == 1


def test_pole_copies_differing_by_more_than_round_off_fail_loudly():
    """One side stores the pole bit-identically, the other with its copies
    spread 1e-09 apart -- inside weld_tol, far outside round-off. There is no
    deterministic pairing, so this must be a ValueError (whichever side the
    search runs from), not a silently mis-wired neighbour."""
    nx, nr, nt = 3, 3, 5
    blocks = _stacked_pole_blocks(nx, nr, nt, noise=1e-9, noisy_blocks=(1,))
    matches, outer = connectivity(blocks)
    assert len(matches) == 1
    swapped = [{"block1": matches[0]["block2"], "block2": matches[0]["block1"]}]
    # exact side -> spread side: the spread copies nobody reached are orphans
    with pytest.raises(ValueError, match="no counterpart"):
        flatten_mesh(blocks, matches, outer, weld_tol=1e-6)
    # spread side -> exact side: distinct nodes collapse onto one copy
    with pytest.raises(ValueError, match="not one-to-one"):
        flatten_mesh(blocks, swapped, outer, weld_tol=1e-6)

    # Whereas copies spread IDENTICALLY on both sides pair up one-to-one with
    # zero residual: that is an unambiguous (if odd) correspondence, and the
    # graph it yields is consistent.
    blocks = _stacked_pole_blocks(nx, nr, nt, noise=1e-9, noisy_blocks=(0, 1))
    matches, outer = connectivity(blocks)
    fm = flatten_mesh(blocks, matches, outer, weld_tol=1e-6)
    _every_interior_face_belongs_to_both_its_cells(fm)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
