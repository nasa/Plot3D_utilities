"""Tests for the ``.graph_p3d`` HDF5 bundle format (``plot3d.glennht.graph_p3d``,
re-exported from ``plot3d.glennht``): ``write_graph_p3d``, ``read_graph_p3d``,
``graph_p3d_to_files``, and the ``bundle=`` flag on
``plot3d.glennht.export_to_glennht_gpu``.

``h5py`` is an OPTIONAL dependency (the ``hdf5`` extra in ``pyproject.toml``)
-- this whole module is skipped cleanly (via the ``pytest.importorskip``
below) when it isn't installed, mirroring how ``graph_p3d.py`` itself only
ever imports h5py lazily, inside its own functions.

Uses the same tiny 2-block annular "wedge" mesh as the
``colab/Plot3D_GlennHT_GPU_GraphExport.ipynb`` tutorial (two 15x9x13-node
wedge blocks sharing the x=0.5 face, spanning a 20-degree sector) so
:func:`plot3d.connectivity.connectivity` -- which is slow on large meshes --
runs in well under a second.
"""

import json
import math
import os
from collections import Counter

import numpy as np
import pytest
import yaml

h5py = pytest.importorskip("h5py")

from plot3d import Block, read_plot3D
from plot3d.connectivity import connectivity
from plot3d.periodicity import create_rotation_matrix, rotated_periodicity
from plot3d.glennht import (
    GpuInletBC,
    GpuOutletBC,
    GpuWallBC,
    export_to_glennht_gpu,
    graph_p3d_to_files,
    read_graph_p3d,
    tag_surfaces_geometric,
    write_graph_p3d,
)


# ---------------------------------------------------------------------------
# Small synthetic fixture: two 15x9x13 annular "wedge" blocks sharing the
# x=0.5 face, spanning a 20-degree sector (18 blades for a full annulus).
# Copied verbatim from colab/Plot3D_GlennHT_GPU_GraphExport.ipynb.
# ---------------------------------------------------------------------------

def wedge_block(x0, x1, nx, r_hub=0.2, r_shroud=0.3, nr=9, dtheta=np.radians(20.0), nt=13):
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


DTHETA = np.radians(20.0)
NBLADES = round(2 * np.pi / DTHETA)  # 18


def _make_wedge_blocks():
    # Two blocks stacked along x, sharing the x=0.5 face.
    return [wedge_block(0.0, 0.5, 15), wedge_block(0.5, 1.0, 15)]


@pytest.fixture
def wedge_blocks():
    return _make_wedge_blocks()


def _build_wedge_connectivity(blocks):
    """(face_matches, outer_faces, periodic_faces, surface_ids, periodicity_meta)
    for the wedge fixture, computed in the SAME order
    ``export_to_glennht_gpu(..., tagging="geometric")`` uses internally:
    connectivity -> rotated_periodicity -> tag_surfaces_geometric.
    """
    face_matches, outer_faces = connectivity(blocks)
    periodic_faces, outer_faces, _, _ = rotated_periodicity(
        blocks, face_matches, outer_faces,
        rotation_angle=math.degrees(DTHETA), rotation_axis="x",
    )
    outer_faces, surface_ids = tag_surfaces_geometric(blocks, outer_faces, axis="x")

    transformation_matrix = create_rotation_matrix(DTHETA, "x")
    periodicity_meta = {
        "nblades": NBLADES,
        "rotation_axis": "x",
        "rotation_angle_rad": float(DTHETA),
        "rotation_angle_deg": float(math.degrees(DTHETA)),
        "transformation_matrix": transformation_matrix.tolist(),
        "convention": "Face_B_points = (transformation_matrix @ Face_A_points.T).T",
    }
    return face_matches, outer_faces, periodic_faces, surface_ids, periodicity_meta


def test_wedge_fixture_has_expected_counts(wedge_blocks):
    """Documents the fixture (mirrors test_glennht_gpu_export.py's analogous
    check) so a future drift in connectivity()/rotated_periodicity() shows
    up here first, rather than as a confusing failure in the tests below."""
    face_matches, outer_faces = connectivity(wedge_blocks)
    assert len(wedge_blocks) == 2
    assert len(face_matches) == 1
    assert len(outer_faces) == 10

    periodic_faces, outer_faces, _, _ = rotated_periodicity(
        wedge_blocks, face_matches, outer_faces,
        rotation_angle=math.degrees(DTHETA), rotation_axis="x",
    )
    assert len(periodic_faces) == 2
    assert len(outer_faces) == 6


# ---------------------------------------------------------------------------
# 1: round-trip coordinates + connectivity
# ---------------------------------------------------------------------------

def test_round_trip_coordinates_and_connectivity(tmp_path, wedge_blocks):
    blocks = wedge_blocks
    face_matches, outer_faces, periodic_faces, surface_ids, periodicity_meta = (
        _build_wedge_connectivity(blocks)
    )

    bundle_path = tmp_path / "wedge.graph_p3d"
    write_graph_p3d(
        str(bundle_path),
        blocks,
        matched_faces=face_matches,
        outer_faces=outer_faces,
        periodic_faces=periodic_faces,
        periodicity=periodicity_meta,
        surface_ids=surface_ids,
        case_name="wedge",
    )

    data = read_graph_p3d(str(bundle_path))
    rblocks = data["blocks"]
    conn = data["connectivity"]

    # --- coordinates + dimensions ---
    assert len(rblocks) == len(blocks) == 2
    for b, rb in zip(blocks, rblocks):
        assert (rb.IMAX, rb.JMAX, rb.KMAX) == (b.IMAX, b.JMAX, b.KMAX)
        assert np.allclose(rb.X, b.X)
        assert np.allclose(rb.Y, b.Y)
        assert np.allclose(rb.Z, b.Z)

    # --- connectivity counts ---
    assert conn["nblocks"] == 2
    assert len(conn["face_matches"]) == len(face_matches) == 1
    assert len(conn["outer_faces"]) == len(outer_faces) == 6
    assert len(conn["periodic_faces"]) == len(periodic_faces) == 2

    # --- surface_ids identical ---
    assert data["surface_ids"] == surface_ids
    assert conn["surface_ids"] == surface_ids

    # --- a face_match preserves DIRECTED lb/ub (not sorted) + permutation_index ---
    orig_fm = face_matches[0]
    got_fm = conn["face_matches"][0]
    assert got_fm["block1"]["lb"] == list(orig_fm["block1"]["lb"])
    assert got_fm["block1"]["ub"] == list(orig_fm["block1"]["ub"])
    assert got_fm["block2"]["lb"] == list(orig_fm["block2"]["lb"])
    assert got_fm["block2"]["ub"] == list(orig_fm["block2"]["ub"])
    assert got_fm["permutation_index"] == orig_fm["orientation"]["permutation_index"]


def test_write_graph_p3d_preserves_directed_face_match_diagonals(tmp_path, wedge_blocks):
    """Dedicated check that a deliberately BACKWARDS diagonal (lb > ub on one
    axis) round-trips byte-for-byte, i.e. is NOT normalised/sorted into
    min/max order. The real fixture's own face match happens to already be
    axis-sorted, so this uses a hand-crafted match to conclusively prove it
    -- mirrors
    test_glennht_gpu_export.py::test_face_match_diagonals_are_directed_ints.
    """
    blocks = wedge_blocks  # only need valid block_index 0/1; geometry unused here
    fake_fm = [{
        "block1": {"block_index": 0, "lb": [2, 0, 0], "ub": [0, 2, 2]},
        "block2": {"block_index": 1, "lb": [0, 0, 0], "ub": [0, 2, 2]},
        "orientation": {"permutation_index": -1, "permutation_matrix": [[1, 0], [0, 1]]},
    }]

    bundle_path = tmp_path / "backwards.graph_p3d"
    write_graph_p3d(str(bundle_path), blocks, matched_faces=fake_fm, outer_faces=[])

    conn = read_graph_p3d(str(bundle_path))["connectivity"]
    rec = conn["face_matches"][0]
    assert rec["block1"]["lb"] == [2, 0, 0]
    assert rec["block1"]["ub"] == [0, 2, 2]
    assert rec["block2"]["lb"] == [0, 0, 0]
    assert rec["block2"]["ub"] == [0, 2, 2]
    assert rec["permutation_index"] == -1


# ---------------------------------------------------------------------------
# 2: explode fidelity -- bundle vs. direct export must agree (the important test)
# ---------------------------------------------------------------------------

def _fm_keys(face_matches):
    """Hashable (b1,lb,ub,b2,lb,ub,permutation_index) key set for a
    face_matches (or periodic_faces) list, order-independent."""
    keys = set()
    for fm in face_matches:
        b1, b2 = fm["block1"], fm["block2"]
        keys.add((
            b1["block_index"], tuple(b1["lb"]), tuple(b1["ub"]),
            b2["block_index"], tuple(b2["lb"]), tuple(b2["ub"]),
            fm["permutation_index"],
        ))
    return keys


def test_graph_p3d_explode_matches_direct_export(tmp_path, wedge_blocks):
    """A .graph_p3d bundle, once exploded back to the loose
    {xyz, connectivity.json} trio via graph_p3d_to_files, must be
    structurally identical to what export_to_glennht_gpu(..., bundle=False)
    writes directly for the same mesh."""
    blocks = wedge_blocks

    bundle_dir = tmp_path / "bundle_out"
    direct_dir = tmp_path / "direct_out"
    explode_dir = tmp_path / "explode_out"

    bundle_paths = export_to_glennht_gpu(
        blocks, str(bundle_dir), "wedge",
        rotation_angle=DTHETA, rotation_axis="x", nblades=NBLADES,
        tagging="geometric", bundle=True,
    )
    direct_paths = export_to_glennht_gpu(
        blocks, str(direct_dir), "wedge",
        rotation_angle=DTHETA, rotation_axis="x", nblades=NBLADES,
        tagging="geometric", bundle=False,
    )

    assert "graph_p3d" in bundle_paths
    assert os.path.exists(bundle_paths["graph_p3d"])
    assert "graph_p3d" not in direct_paths

    exploded_paths = graph_p3d_to_files(bundle_paths["graph_p3d"], str(explode_dir))

    with open(exploded_paths["connectivity"]) as fh:
        exploded_conn = json.load(fh)
    with open(direct_paths["connectivity"]) as fh:
        direct_conn = json.load(fh)

    # face_matches: identical set of (b1,lb,ub,b2,lb,ub,permutation_index)
    assert _fm_keys(exploded_conn["face_matches"]) == _fm_keys(direct_conn["face_matches"])
    assert len(exploded_conn["face_matches"]) == len(direct_conn["face_matches"]) == 1

    # outer_faces: identical id distribution
    exploded_ids = Counter(f["id"] for f in exploded_conn["outer_faces"])
    direct_ids = Counter(f["id"] for f in direct_conn["outer_faces"])
    assert exploded_ids == direct_ids
    assert sum(exploded_ids.values()) == len(exploded_conn["outer_faces"]) == 6

    # periodicity
    assert np.allclose(
        exploded_conn["periodicity"]["transformation_matrix"],
        direct_conn["periodicity"]["transformation_matrix"],
    )

    # periodic_faces: same count (structurally identical set too, as a bonus)
    assert len(exploded_conn["periodic_faces"]) == len(direct_conn["periodic_faces"]) == 2
    assert _fm_keys(exploded_conn["periodic_faces"]) == _fm_keys(direct_conn["periodic_faces"])

    # the exploded .xyz reads back to the SAME coordinates as the original blocks
    exploded_blocks = read_plot3D(exploded_paths["grid"], binary=False)
    assert len(exploded_blocks) == len(blocks)
    for b, eb in zip(blocks, exploded_blocks):
        assert (eb.IMAX, eb.JMAX, eb.KMAX) == (b.IMAX, b.JMAX, b.KMAX)
        assert np.allclose(eb.X, b.X)
        assert np.allclose(eb.Y, b.Y)
        assert np.allclose(eb.Z, b.Z)


# ---------------------------------------------------------------------------
# 3: boundary-condition round trip
# ---------------------------------------------------------------------------

def test_boundary_conditions_round_trip(tmp_path, wedge_blocks):
    blocks = wedge_blocks
    bcs = [
        GpuInletBC(
            name="wedge_inlet", surfaces=[1],
            total_pressure=101325.0, total_temperature=288.15,
        ),
        GpuOutletBC(
            name="wedge_outlet", surfaces=[2],
            back_pressure=100000.0,
        ),
        GpuWallBC(name="wedge_walls", surfaces=[3, 4, 5]),
    ]

    paths = export_to_glennht_gpu(
        blocks, str(tmp_path / "bc_out"), "wedge",
        rotation_angle=DTHETA, rotation_axis="x", nblades=NBLADES,
        tagging="geometric", bundle=True, bcs=bcs,
    )

    data = read_graph_p3d(paths["graph_p3d"])
    assert data["boundary_conditions"]
    assert isinstance(data["boundary_conditions"], str)
    assert len(data["boundary_conditions"]) > 0

    exploded_paths = graph_p3d_to_files(paths["graph_p3d"], str(tmp_path / "bc_explode"))

    assert "boundary_conditions" in exploded_paths
    assert os.path.exists(exploded_paths["boundary_conditions"])

    with open(exploded_paths["boundary_conditions"]) as fh:
        reloaded = yaml.safe_load(fh)

    assert isinstance(reloaded["boundary_conditions"], list)
    assert len(reloaded["boundary_conditions"]) == 3

    by_name = {d["name"]: d for d in reloaded["boundary_conditions"]}
    assert set(by_name) == {"wedge_inlet", "wedge_outlet", "wedge_walls"}
    assert by_name["wedge_inlet"]["type"] == "inlet"
    assert by_name["wedge_outlet"]["type"] == "outlet"
    assert by_name["wedge_walls"]["type"] == "wall"
    assert by_name["wedge_inlet"]["surfaces"] == [1]
    assert by_name["wedge_outlet"]["surfaces"] == [2]
    assert by_name["wedge_walls"]["surfaces"] == [3, 4, 5]


# ---------------------------------------------------------------------------
# 4: export_to_glennht_gpu's bundle= flag
# ---------------------------------------------------------------------------

def test_export_to_glennht_gpu_bundle_true_includes_graph_p3d(tmp_path, wedge_blocks):
    paths = export_to_glennht_gpu(
        wedge_blocks, str(tmp_path), "wedge",
        tagging="geometric", bundle=True,
    )
    assert "graph_p3d" in paths
    assert os.path.exists(paths["graph_p3d"])

    data = read_graph_p3d(paths["graph_p3d"])
    assert data["meta"]["format"] == "graph_p3d"
    assert len(data["blocks"]) == 2


def test_export_to_glennht_gpu_bundle_false_by_default(tmp_path, wedge_blocks):
    paths = export_to_glennht_gpu(
        wedge_blocks, str(tmp_path), "wedge",
        tagging="geometric",
    )
    assert "graph_p3d" not in paths
    # sanity: the driver still did the normal (non-bundle) work
    assert "connectivity" in paths
    assert os.path.exists(paths["connectivity"])


# ---------------------------------------------------------------------------
# 5: HDF5 layout smoke check
# ---------------------------------------------------------------------------

def test_graph_p3d_hdf5_layout(tmp_path, wedge_blocks):
    blocks = wedge_blocks
    paths = export_to_glennht_gpu(
        blocks, str(tmp_path), "wedge",
        rotation_angle=DTHETA, rotation_axis="x", nblades=NBLADES,
        tagging="geometric", bundle=True,
    )

    with h5py.File(paths["graph_p3d"], "r") as f:
        fmt = f.attrs["format"]
        if isinstance(fmt, bytes):
            fmt = fmt.decode("utf-8")
        assert fmt == "graph_p3d"

        assert "grid" in f
        assert "block_0" in f["grid"]
        b0_grp = f["grid"]["block_0"]
        b0 = blocks[0]
        assert int(b0_grp.attrs["imax"]) == b0.IMAX
        assert int(b0_grp.attrs["jmax"]) == b0.JMAX
        assert int(b0_grp.attrs["kmax"]) == b0.KMAX
        for name in ("x", "y", "z"):
            assert b0_grp[name].shape == (b0.IMAX, b0.JMAX, b0.KMAX)

        assert "connectivity" in f
        assert "face_matches" in f["connectivity"]
        assert f["connectivity"]["face_matches"].shape[1] == 15


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
