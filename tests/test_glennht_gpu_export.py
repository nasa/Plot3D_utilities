"""Tests for the glennht-gpu export module (``plot3d.glennht.gpu`` +
``plot3d.glennht.gpu_boundary_conditions``).

Uses a tiny synthetic 2-block mesh (two adjacent 3x3x3-node unit-spaced
cubes sharing one face) so :func:`plot3d.connectivity` -- which is slow on
large meshes -- runs in well under a second. The fixture always produces
exactly 1 face match and 10 outer faces (5 per block, one side per block
consumed by the match).

A gated round-trip test at the bottom exercises a real code-tagged stator
case from an external test-data directory (configured via an environment
variable, not part of the repo) when available; it is skipped everywhere
else via ``pytest.mark.skipif``, mirroring the pattern in
``test_cross_plane.py`` / ``test_cmc009_connectivity_fortran.py``.
"""

import copy
import json
import math
import os

import numpy as np
import pytest
import yaml

from plot3d import Block, read_plot3D
from plot3d.connectivity import connectivity
from plot3d.glennht import (
    GpuInletBC,
    GpuOutletBC,
    GpuWallBC,
    export_to_glennht_gpu,
    merge_connectivity_json,
    tag_surfaces_from_bc_codes,
    tag_surfaces_from_diagonals,
    tag_surfaces_geometric,
    write_bc_codes_json,
    write_boundary_conditions_yaml,
    write_connectivity_json,
)
from plot3d.glennht.gpu import UNMATCHED_INTERFACE_ID


# ---------------------------------------------------------------------------
# Small synthetic fixture: two adjacent 3x3x3 cubes sharing the
# block0 I=IMAX / block1 I=1 face.
# ---------------------------------------------------------------------------

def _make_cube_block(x_offset: float) -> Block:
    n = 3
    x, y, z = np.meshgrid(
        np.arange(n) + x_offset, np.arange(n), np.arange(n), indexing="ij"
    )
    return Block(x.astype(float), y.astype(float), z.astype(float))


def _make_fixture_blocks():
    # block0 spans x in [0, 2], block1 spans x in [2, 4]: they share the
    # x=2 plane (block0's I=IMAX face == block1's I=1 face).
    return [_make_cube_block(0.0), _make_cube_block(2.0)]


@pytest.fixture
def fixture_blocks():
    return _make_fixture_blocks()


@pytest.fixture
def fixture_connectivity(fixture_blocks):
    """Fresh (blocks, face_matches, outer_faces) for each test that needs
    them -- tag_surfaces_* mutate outer_faces in place, so each test gets
    its own copy to avoid cross-test contamination."""
    face_matches, outer_faces = connectivity(fixture_blocks)
    return fixture_blocks, face_matches, outer_faces


def test_fixture_has_one_match_and_ten_outer_faces(fixture_connectivity):
    blocks, face_matches, outer_faces = fixture_connectivity
    assert len(blocks) == 2
    assert len(face_matches) == 1
    assert len(outer_faces) == 10


# ---------------------------------------------------------------------------
# 1 & 2: connectivity.json schema + directed diagonals
# ---------------------------------------------------------------------------

def test_write_connectivity_json_schema(tmp_path, fixture_connectivity):
    blocks, face_matches, outer_faces = fixture_connectivity
    out_path = tmp_path / "fixture.connectivity.json"

    write_connectivity_json(blocks, face_matches, copy.deepcopy(outer_faces), str(out_path))

    assert out_path.exists()
    with open(out_path) as fh:
        payload = json.load(fh)

    for key in (
        "mesh_file",
        "nblocks",
        "face_matches",
        "outer_faces",
        "permutation_matrices",
    ):
        assert key in payload

    assert payload["nblocks"] == len(blocks)
    assert len(payload["face_matches"]) == 1
    assert len(payload["outer_faces"]) == 10

    # permutation_matrices is 8x2x2
    pm = payload["permutation_matrices"]
    assert len(pm) == 8
    for row in pm:
        assert len(row) == 2
        for r in row:
            assert len(r) == 2

    fm = payload["face_matches"][0]
    for side in ("block1", "block2"):
        rec = fm[side]
        assert "block_index" in rec
        assert isinstance(rec["lb"], list) and len(rec["lb"]) == 3
        assert isinstance(rec["ub"], list) and len(rec["ub"]) == 3
    assert "permutation_index" in fm
    assert "permutation_matrix" in fm


def test_face_match_diagonals_are_directed_ints(fixture_connectivity):
    """lb/ub are length-3 int lists; schema allows lb > ub (directed
    traversal order), so don't assume they're sorted min/max."""
    blocks, face_matches, outer_faces = fixture_connectivity
    fm = face_matches[0]
    for side in ("block1", "block2"):
        rec = fm[side]
        lb, ub = rec["lb"], rec["ub"]
        assert len(lb) == 3 and len(ub) == 3
        assert all(isinstance(v, (int, np.integer)) for v in lb)
        assert all(isinstance(v, (int, np.integer)) for v in ub)

    # This fixture's match happens to be in-plane / axis-aligned (lb<=ub
    # on both sides), but the schema itself permits il > ih (a reversed
    # traversal) -- verify write_connectivity_json doesn't reorder/clamp
    # a deliberately "backwards" diagonal.
    fake_fm = [{
        "block1": {"block_index": 0, "lb": [2, 0, 0], "ub": [0, 2, 2]},
        "block2": {"block_index": 1, "lb": [0, 0, 0], "ub": [0, 2, 2]},
        "orientation": {"permutation_index": -1, "permutation_matrix": [[1, 0], [0, 1]]},
    }]
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "backwards.connectivity.json")
        payload = write_connectivity_json(blocks, fake_fm, [], path)
    rec = payload["face_matches"][0]["block1"]
    assert rec["lb"] == [2, 0, 0]
    assert rec["ub"] == [0, 2, 2]


# ---------------------------------------------------------------------------
# 3: tag_surfaces_from_bc_codes
# ---------------------------------------------------------------------------

def test_tag_surfaces_from_bc_codes(fixture_connectivity):
    blocks, face_matches, outer_faces = fixture_connectivity
    outer_faces = copy.deepcopy(outer_faces)

    # FACE_ORDER = [I=1, I=IMAX, J=1, J=JMAX, K=1, K=KMAX]
    # block0: I=IMAX is matched (unused code); block1: I=1 is matched.
    bc_codes = [
        [5, 0, 10, 9, 109, 6],    # block0
        [0, 6, 10, 9, 109, 0],    # block1
    ]

    tagged, surface_ids = tag_surfaces_from_bc_codes(blocks, outer_faces, bc_codes)

    expected_by_diag = {
        (0, (0, 0, 0), (0, 2, 2)): 1,   # block0 I=1 -> inlet
        (0, (0, 0, 0), (2, 0, 2)): 3,   # block0 J=1 -> blade
        (0, (0, 2, 0), (2, 2, 2)): 4,   # block0 J=JMAX -> hub
        (0, (0, 0, 0), (2, 2, 0)): 5,   # block0 K=1 -> shroud
        (0, (0, 0, 2), (2, 2, 2)): 2,   # block0 K=KMAX -> outlet
        (1, (2, 0, 0), (2, 2, 2)): 2,   # block1 I=IMAX -> outlet
        (1, (0, 0, 0), (2, 0, 2)): 3,   # block1 J=1 -> blade
        (1, (0, 2, 0), (2, 2, 2)): 4,   # block1 J=JMAX -> hub
        (1, (0, 0, 0), (2, 2, 0)): 5,   # block1 K=1 -> shroud
        (1, (0, 0, 2), (2, 2, 2)): UNMATCHED_INTERFACE_ID,  # block1 K=KMAX -> interface
    }

    seen = 0
    for face in tagged:
        key = (face["block_index"], tuple(face["lb"]), tuple(face["ub"]))
        assert key in expected_by_diag, f"unexpected outer face {face}"
        assert face["id"] == expected_by_diag[key], f"wrong id for {face}"
        seen += 1
    assert seen == 10

    assert surface_ids["1"] == "inlet"
    assert surface_ids["2"] == "outlet"
    assert surface_ids["3"] == "blade"
    assert surface_ids["4"] == "hub"
    assert surface_ids["5"] == "shroud"
    assert surface_ids[str(UNMATCHED_INTERFACE_ID)] == "unmatched_interface"


# ---------------------------------------------------------------------------
# 4: tag_surfaces_geometric
# ---------------------------------------------------------------------------

def test_tag_surfaces_geometric(fixture_connectivity):
    blocks, face_matches, outer_faces = fixture_connectivity
    outer_faces = copy.deepcopy(outer_faces)

    tagged, surface_ids = tag_surfaces_geometric(blocks, outer_faces, axis="x", band=0.1)

    valid_ids = {1, 2, 3, 4, 5}
    assert len(tagged) == 10
    for face in tagged:
        assert face["id"] in valid_ids

    assert surface_ids == {
        "1": "inlet",
        "2": "outlet",
        "3": "blade",
        "4": "hub",
        "5": "shroud",
    }
    # At least one face should land on each extreme axial side given the
    # fixture spans x in [0, 4].
    ids_present = {f["id"] for f in tagged}
    assert 1 in ids_present  # inlet (min x)
    assert 2 in ids_present  # outlet (max x)


# ---------------------------------------------------------------------------
# 5: tag_surfaces_from_diagonals
# ---------------------------------------------------------------------------

def test_tag_surfaces_from_diagonals(fixture_connectivity):
    blocks, face_matches, outer_faces = fixture_connectivity
    outer_faces = copy.deepcopy(outer_faces)

    # block0's I=1 face, exact diagonal.
    specs = [{
        "block_index": 0,
        "lb": [0, 0, 0],
        "ub": [0, 2, 2],
        "id": 42,
        "name": "test_inlet",
    }]

    tagged, surface_ids = tag_surfaces_from_diagonals(blocks, outer_faces, specs)

    matched = [f for f in tagged if f.get("id") == 42]
    assert len(matched) == 1
    face = matched[0]
    assert face["block_index"] == 0
    assert sorted(face["lb"]) == sorted([0, 0, 0])
    assert set(zip(face["lb"], face["ub"])) == {(0, 0), (0, 2), (0, 2)}

    assert surface_ids == {"42": "test_inlet"}

    # No other face should have been given id 42.
    others = [f for f in tagged if f is not face]
    assert all(o.get("id") != 42 for o in others)


# ---------------------------------------------------------------------------
# 6: write_bc_codes_json round trip
# ---------------------------------------------------------------------------

def test_write_bc_codes_json_roundtrip(tmp_path, fixture_blocks):
    bc_codes = [
        [5, 0, 10, 9, 109, 6],
        [0, 6, 10, 9, 109, 0],
    ]
    out_path = tmp_path / "fixture.bc_codes.json"

    write_bc_codes_json(fixture_blocks, bc_codes, str(out_path))

    with open(out_path) as fh:
        payload = json.load(fh)

    assert payload["face_order"] == ["I=1", "I=IMAX", "J=1", "J=JMAX", "K=1", "K=KMAX"]
    assert payload["blocks"] == bc_codes
    assert "face_code_legend" in payload
    assert payload["block_order"] == ["block_0", "block_1"]


# ---------------------------------------------------------------------------
# 7: merge_connectivity_json
# ---------------------------------------------------------------------------

def _tiny_connectivity_payload(nblocks, outer_ids):
    return {
        "mesh_file": "row.xyz",
        "nblocks": nblocks,
        "face_matches": [],
        "outer_faces": [
            {"block_index": 0, "lb": [0, 0, 0], "ub": [0, 1, 1], "id": outer_ids[0]},
            {"block_index": 1, "lb": [1, 0, 0], "ub": [1, 1, 1], "id": outer_ids[1]},
        ],
        "periodic_faces": [],
        "surface_ids": {str(outer_ids[0]): "a", str(outer_ids[1]): "b"},
        "permutation_matrices": [],
    }


def test_merge_connectivity_json(tmp_path):
    row0 = _tiny_connectivity_payload(2, [1, 2])
    row1 = _tiny_connectivity_payload(2, [1, 2])

    out_path = tmp_path / "merged.connectivity.json"
    merged = merge_connectivity_json([row0, row1], str(out_path), id_stride=100)

    assert out_path.exists()
    with open(out_path) as fh:
        reloaded = json.load(fh)
    assert reloaded == merged

    assert merged["nblocks"] == 4  # 2 + 2

    row0_outer = [f for f in merged["outer_faces"] if f["id"] in (1, 2)]
    row1_outer = [f for f in merged["outer_faces"] if f["id"] in (101, 102)]
    assert len(row0_outer) == 2
    assert len(row1_outer) == 2

    # Row 1's block_index is offset by row 0's block count (2).
    for f in row1_outer:
        assert f["block_index"] in (2, 3)
    for f in row0_outer:
        assert f["block_index"] in (0, 1)

    # Row 1's outer_faces[].id offset by id_stride=100.
    ids = sorted(f["id"] for f in merged["outer_faces"])
    assert ids == [1, 2, 101, 102]


def test_merge_connectivity_json_accepts_file_paths(tmp_path):
    row0_path = tmp_path / "row0.json"
    row1_path = tmp_path / "row1.json"
    with open(row0_path, "w") as fh:
        json.dump(_tiny_connectivity_payload(2, [1, 2]), fh)
    with open(row1_path, "w") as fh:
        json.dump(_tiny_connectivity_payload(3, [1, 2]), fh)

    out_path = tmp_path / "merged.connectivity.json"
    merged = merge_connectivity_json([str(row0_path), str(row1_path)], str(out_path), id_stride=100)

    assert merged["nblocks"] == 5  # 2 + 3
    row1_outer = [f for f in merged["outer_faces"] if f["id"] in (101, 102)]
    assert all(f["block_index"] in (2, 3, 4) for f in row1_outer)


# ---------------------------------------------------------------------------
# 8: BC YAML dataclasses
# ---------------------------------------------------------------------------

def test_write_boundary_conditions_yaml_roundtrip(tmp_path):
    inlet = GpuInletBC(
        name="m_inlet",
        surfaces=[1],
        total_pressure=101325.0,
        total_temperature=288.15,
    )
    outlet = GpuOutletBC(
        name="m_outlet",
        surfaces=[2],
        back_pressure=90000.0,
        extrapolation_order=0,
    )
    wall = GpuWallBC(
        name="blade_wall",
        surfaces=[3, 4, 5],
    )

    out_path = tmp_path / "bcs.yaml"
    rotation = {"omega_x": 0.0, "per_block": [{"blocks": [0, 1], "omega_x": -1792.7}]}

    text = write_boundary_conditions_yaml([inlet, outlet, wall], str(out_path), rotation=rotation)
    assert isinstance(text, str) and len(text) > 0

    with open(out_path) as fh:
        reloaded = yaml.safe_load(fh)

    assert isinstance(reloaded["boundary_conditions"], list)
    assert len(reloaded["boundary_conditions"]) == 3

    inlet_d, outlet_d, wall_d = reloaded["boundary_conditions"]

    # Inlet: only set (non-None) fields appear.
    assert inlet_d == {
        "name": "m_inlet",
        "type": "inlet",
        "surfaces": [1],
        "total_pressure": 101325.0,
        "total_temperature": 288.15,
    }
    assert "flow_angle_deg" not in inlet_d
    assert "turbulence_intensity" not in inlet_d
    assert "blades_full_ring" not in inlet_d

    # Outlet: extrapolation_order set, others absent.
    assert outlet_d == {
        "name": "m_outlet",
        "type": "outlet",
        "surfaces": [2],
        "back_pressure": 90000.0,
        "extrapolation_order": 0,
    }
    assert "mixing_plane_partner" not in outlet_d

    # Wall: only defaults (name/type/surfaces/thermal).
    assert wall_d == {
        "name": "blade_wall",
        "type": "wall",
        "surfaces": [3, 4, 5],
        "thermal": "adiabatic",
    }
    assert "wall_temperature" not in wall_d
    assert "rotating" not in wall_d

    assert reloaded["rotation"] == rotation


def test_gpu_bc_extra_fields_merge(tmp_path):
    inlet = GpuInletBC(
        name="downstream_inlet",
        surfaces=[101],
        total_pressure=200000.0,
        total_temperature=400.0,
        extra={"mach": 0.45},
    )
    out_path = tmp_path / "bcs_extra.yaml"
    write_boundary_conditions_yaml([inlet], str(out_path))
    with open(out_path) as fh:
        reloaded = yaml.safe_load(fh)
    assert reloaded["boundary_conditions"][0]["mach"] == 0.45
    assert "rotation" not in reloaded


# ---------------------------------------------------------------------------
# 9: export_to_glennht_gpu driver
# ---------------------------------------------------------------------------

def test_export_to_glennht_gpu_geometric(tmp_path, fixture_blocks):
    paths = export_to_glennht_gpu(
        fixture_blocks,
        str(tmp_path),
        "fixture",
        tagging="geometric",
        write_grid=True,
    )

    assert "grid" in paths
    assert "connectivity" in paths
    assert os.path.exists(paths["grid"])
    assert os.path.exists(paths["connectivity"])

    with open(paths["connectivity"]) as fh:
        payload = json.load(fh)
    assert payload["nblocks"] == 2
    assert len(payload["face_matches"]) == 1
    assert len(payload["outer_faces"]) == 10
    for f in payload["outer_faces"]:
        assert f["id"] in (1, 2, 3, 4, 5)


def test_export_to_glennht_gpu_requires_nblades_with_rotation_angle(tmp_path, fixture_blocks):
    with pytest.raises(ValueError):
        export_to_glennht_gpu(
            fixture_blocks,
            str(tmp_path),
            "fixture",
            rotation_angle=math.radians(10.0),
            tagging="geometric",
        )


def test_export_to_glennht_gpu_requires_bc_codes_for_bc_codes_tagging(tmp_path, fixture_blocks):
    with pytest.raises(ValueError):
        export_to_glennht_gpu(
            fixture_blocks,
            str(tmp_path),
            "fixture",
            tagging="bc_codes",
        )


def test_export_to_glennht_gpu_unknown_tagging(tmp_path, fixture_blocks):
    with pytest.raises(ValueError):
        export_to_glennht_gpu(
            fixture_blocks,
            str(tmp_path),
            "fixture",
            tagging="bogus",
        )


# ---------------------------------------------------------------------------
# 10: (gated) stator round-trip against real code-tagged mesh
# ---------------------------------------------------------------------------

STATOR_DIR = os.environ.get("PLOT3D_GLENNHT_GPU_STATOR_DIR", "")
STATOR_XYZ = os.path.join(STATOR_DIR, "stator.xyz") if STATOR_DIR else ""
STATOR_BC_CODES = os.path.join(STATOR_DIR, "stator.bc_codes.json") if STATOR_DIR else ""

skip_no_stator_data = pytest.mark.skipif(
    not (STATOR_DIR and os.path.exists(STATOR_XYZ) and os.path.exists(STATOR_BC_CODES)),
    reason=(
        "external code-tagged stator test data not found (set "
        "PLOT3D_GLENNHT_GPU_STATOR_DIR to enable this test)"
    ),
)


@skip_no_stator_data
def test_export_to_glennht_gpu_stator_roundtrip(tmp_path):
    """End-to-end against the real face-code-tagged stator mesh (external data).

    This calls the FULL connectivity() on a real ~205x21x81-ish 5-block
    mesh, so it is comparatively slow (order ~1 minute) -- kept as the
    only such test in this file, and skipped entirely on machines without
    the external test-case checkout.
    """
    blocks = read_plot3D(STATOR_XYZ)

    with open(STATOR_BC_CODES) as fh:
        bc_payload = json.load(fh)
    bc_codes = bc_payload["blocks"]

    paths = export_to_glennht_gpu(
        blocks,
        str(tmp_path),
        "stator",
        rotation_angle=math.radians(7.826086956521739),
        rotation_axis="x",
        nblades=46,
        tagging="bc_codes",
        bc_codes=bc_codes,
        write_grid=False,
    )

    assert "connectivity" in paths
    with open(paths["connectivity"]) as fh:
        payload = json.load(fh)

    assert len(payload["face_matches"]) == 10

    id_counts = {}
    for f in payload["outer_faces"]:
        fid = f.get("id")
        if fid is not None:
            id_counts[fid] = id_counts.get(fid, 0) + 1

    assert id_counts == {1: 3, 2: 3, 3: 1, 4: 5, 5: 5}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
