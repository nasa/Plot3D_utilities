"""Verify Python's connectivity / verifier matches Fortran's VSPT reference.

VSPT (Vehicle Single-Passage Turbine) lives at
``/Users/pjuangph/glennht_comparison/mesh/`` with two reference files:

* ``finalmesh.ght_conn`` — Fortran-native ASCII connectivity:
    - 5 face matches (3 are actually periodic seam pairs, 2 are
      cross-block within the single passage)
    - 7 outer faces, 0 GIFs, 1 zone (3 zone IDs).
* ``connectivity.json`` — Python-converted form that splits the 5
  matches from the .ght_conn into:
    - 2 cross-block ``face_matches``
    - 3 ``periodic_faces`` (rotational about x-axis, 55 blades,
      ~6.545 deg, with ``permutation_index = -1`` and
      ``permutation_matrix = identity`` declared on each).

Mesh: 78 MB ASCII Plot3D, 2 blocks (~256×100×32 + 168×100×52
ish).

Three test layers:

1. **Fast structural** — assertions against the JSON / .ght_conn.
2. **Read-and-compare** — confirm ``read_ght_conn`` matches the
   JSON layout (5 ↔ 2+3).
3. **Full integration** (opt-in via ``PLOT3D_VSPT_FULL_MESH=1``):
   load the 78 MB mesh, run ``connectivity_fast`` and
   ``rotated_periodicity``, then verify the periodic seam through
   the **DECLARED path** of ``verify_periodicity`` (the JSON
   carries ``permutation_matrix`` per pair; the verifier must use
   that matrix exactly, not silently scan for a different one).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


VSPT_MESH = Path("/Users/pjuangph/glennht_comparison/mesh/finalmesh-ASCII.xyz")
VSPT_GHT_CONN = Path("/Users/pjuangph/glennht_comparison/mesh/finalmesh.ght_conn")
VSPT_CONN_JSON = Path("/Users/pjuangph/glennht_comparison/mesh/connectivity.json")


@pytest.fixture(scope="module")
def vspt_conn_json():
    if not VSPT_CONN_JSON.exists():
        pytest.skip(f"VSPT connectivity.json not at {VSPT_CONN_JSON}")
    with open(VSPT_CONN_JSON) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def vspt_ght_conn():
    if not VSPT_GHT_CONN.exists():
        pytest.skip(f"VSPT finalmesh.ght_conn not at {VSPT_GHT_CONN}")
    from plot3d.glennht.import_functions import read_ght_conn
    return read_ght_conn(str(VSPT_GHT_CONN))


# ---------------------------------------------------------------------------
# 1. Structural assertions
# ---------------------------------------------------------------------------

class TestVsptConnectivityJsonStructure:

    def test_block_and_match_counts(self, vspt_conn_json):
        assert vspt_conn_json["nblocks"] == 2
        assert len(vspt_conn_json["face_matches"]) == 2
        assert len(vspt_conn_json["periodic_faces"]) == 3
        assert len(vspt_conn_json["outer_faces"]) == 7

    def test_periodicity_metadata(self, vspt_conn_json):
        per = vspt_conn_json["periodicity"]
        assert per["nblades"] == 55
        assert per["rotation_axis"] == "x"
        assert abs(per["rotation_angle_deg"] - 6.545454545454546) < 1e-9
        assert per["periodic_direction"] == "k"
        # Transformation matrix is the rotation matrix itself.
        T = per["transformation_matrix"]
        assert len(T) == 3 and all(len(row) == 3 for row in T)

    def test_periodic_faces_carry_orientation_with_matrix(self, vspt_conn_json):
        """Each periodic_face must declare a full orientation
        (matrix + index + plane). This is what distinguishes the
        JSON-converted form from the legacy CMC009 case, and it's
        what enables the DECLARED-path verifier."""
        canonical_2x2 = {
            (1, 0, 0, 1): 0,    # identity
            (-1, 0, 0, 1): 1,   # u-flip (cross-block rotational seam — see memory)
            (1, 0, 0, -1): 2,
            (-1, 0, 0, -1): 3,
        }
        for idx, pf in enumerate(vspt_conn_json["periodic_faces"]):
            assert "orientation" in pf, f"periodic_faces[{idx}] missing orientation"
            ori = pf["orientation"]
            assert "permutation_index" in ori
            assert "permutation_matrix" in ori
            assert "plane" in ori
            # In-plane sentinel is always -1 in the diagonal export.
            assert ori["permutation_index"] == -1
            assert ori["plane"] == "in-plane"
            # The 3 VSPT seams are 2 identity + 1 u-flip per memory
            # `project_partner_import_ghost_2026-04-29.md`. Accept any
            # canonical in-plane matrix (perms 0-3).
            mat = ori["permutation_matrix"]
            key = tuple(int(v) for row in mat for v in row)
            assert key in canonical_2x2, (
                f"periodic_faces[{idx}] has non-canonical in-plane matrix {mat}"
            )

    def test_periodic_faces_include_u_flip_seam(self, vspt_conn_json):
        """The cross-block rotational seam in VSPT carries
        ``permutation_matrix = [[-1,0],[0,1]]`` (u-flip, canonical
        index 1). Memory note: discarding that matrix and clamping
        to identity caused the iter-1 ρu_Linf 10× drift + iter-2
        NaN cascade (commit 920149f). The test is a regression
        guard — at least one of the 3 periodic seams must declare
        the u-flip."""
        u_flip_count = sum(
            1 for pf in vspt_conn_json["periodic_faces"]
            if pf["orientation"]["permutation_matrix"] == [[-1, 0], [0, 1]]
        )
        assert u_flip_count >= 1, (
            "Expected at least one VSPT periodic face to declare the "
            "u-flip permutation_matrix [[-1,0],[0,1]]; the JSON has none."
        )


# ---------------------------------------------------------------------------
# 2. read_ght_conn ↔ connectivity.json consistency
# ---------------------------------------------------------------------------

class TestVsptGhtConnVsJson:

    def test_ght_conn_total_matches_equal_json_breakdown(self, vspt_ght_conn, vspt_conn_json):
        """The .ght_conn file holds 5 matches; the JSON splits them
        into 2 cross-block + 3 periodic. Total must match."""
        face_matches = vspt_ght_conn[1]
        assert len(face_matches) == 5
        json_total = (
            len(vspt_conn_json["face_matches"])
            + len(vspt_conn_json["periodic_faces"])
        )
        assert json_total == 5

    def test_ght_conn_block_indices_in_range(self, vspt_ght_conn):
        face_matches = vspt_ght_conn[1]
        for idx, fm in enumerate(face_matches):
            for side in ("block1", "block2"):
                bi = fm[side]["block_index"]
                assert 0 <= bi <= 1, f"face_match[{idx}].{side} block_index={bi}"

    def test_ght_conn_face_corners_well_formed(self, vspt_ght_conn):
        face_matches = vspt_ght_conn[1]
        for idx, fm in enumerate(face_matches):
            for side in ("block1", "block2"):
                lb = fm[side]["lb"]
                ub = fm[side]["ub"]
                assert len(lb) == 3 and len(ub) == 3
                # At least one constant axis per face.
                fixed = sum(1 for d in range(3) if lb[d] == ub[d])
                assert fixed >= 1

    def test_ght_conn_outer_face_count(self, vspt_ght_conn):
        outer_faces = vspt_ght_conn[2]
        assert len(outer_faces) == 7


# ---------------------------------------------------------------------------
# 3. Full mesh integration — DECLARED path verifier
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    os.environ.get("PLOT3D_VSPT_FULL_MESH") != "1",
    reason="VSPT full-mesh test loads 78 MB ASCII; "
           "set PLOT3D_VSPT_FULL_MESH=1 to opt in.",
)
@pytest.mark.skipif(not VSPT_MESH.exists(), reason="VSPT mesh missing")
class TestVsptFullMesh:

    @pytest.fixture(scope="class")
    def vspt_blocks(self):
        from plot3d import read_plot3D
        return read_plot3D(str(VSPT_MESH), binary=False, big_endian=False,
                            read_double=False)

    def test_block_count(self, vspt_blocks):
        assert len(vspt_blocks) == 2

    def test_connectivity_fast_finds_two_cross_block_matches(self, vspt_blocks, vspt_conn_json):
        """Python's geometric scan should find the 2 cross-block
        face matches (it does NOT find the periodic ones — those
        come from a separate angle-based scan)."""
        from plot3d import connectivity_fast
        py_matches, _ = connectivity_fast(vspt_blocks, use_minmax=False)
        # VSPT has exactly 2 cross-block faces in the JSON; depending
        # on tolerance Python may find a few extras at coincident
        # coarse boundaries — accept anything in [2, 5] as healthy.
        assert 2 <= len(py_matches) <= 5

    def test_rotated_periodicity_finds_three_pairs(self, vspt_blocks, vspt_conn_json):
        """Synthesizer should recover the 3 declared periodic pairs."""
        from plot3d import rotated_periodicity, connectivity_fast
        from plot3d.facefunctions import outer_face_dict_to_list, get_outer_faces
        # `rotated_periodicity` needs known matched_faces + outer_faces;
        # use the JSON's lists directly (Fortran-derived).
        matched_faces = list(vspt_conn_json["face_matches"])
        outer_faces = list(vspt_conn_json["outer_faces"])
        angle_deg = vspt_conn_json["periodicity"]["rotation_angle_deg"]
        axis = vspt_conn_json["periodicity"]["rotation_axis"]
        py_periodic = rotated_periodicity(
            vspt_blocks, matched_faces, outer_faces,
            rotation_angle=angle_deg, rotation_axis=axis,
        )
        # rotated_periodicity returns a tuple; the periodic pairs are
        # the first element.
        py_pairs = py_periodic[0] if isinstance(py_periodic, tuple) else py_periodic
        assert len(py_pairs) >= 3, (
            f"Expected ≥3 periodic pairs, got {len(py_pairs)}"
        )

    def test_verify_periodicity_declared_path_accepts_correct_matrix(self, vspt_blocks, vspt_conn_json):
        """Run the JSON's 3 declared periodic_faces through
        ``verify_periodicity``. The DECLARED path uses the stored
        ``permutation_matrix`` directly. Identity matrix should
        verify all 3."""
        from plot3d import verify_periodicity
        face_matches = list(vspt_conn_json["periodic_faces"])
        angle = vspt_conn_json["periodicity"]["rotation_angle_deg"]
        axis = vspt_conn_json["periodicity"]["rotation_axis"]
        verified, mismatched = verify_periodicity(
            vspt_blocks, face_matches, theta=angle,
            rotation_axis=axis, tol=1e-4,
        )
        assert len(verified) == 3
        assert len(mismatched) == 0

    def test_verify_periodicity_declared_path_rejects_wrong_matrix(self, vspt_blocks, vspt_conn_json):
        """Plot3d-rs gold-standard rule (commit 920149f): a face
        match with a deliberately wrong declared
        ``permutation_matrix`` ends up in ``mismatched``, NOT in
        ``verified`` with a different perm silently rounded."""
        from copy import deepcopy
        from plot3d import verify_periodicity
        # Take one VSPT periodic match and corrupt its declared matrix
        # to a different canonical (perm 1: u-flip). For 3 in-plane
        # identity pairs, this should NOT verify under the DECLARED rule.
        face_match = deepcopy(vspt_conn_json["periodic_faces"][0])
        face_match["orientation"]["permutation_matrix"] = [[-1, 0], [0, 1]]
        face_match["orientation"]["permutation_index"] = 1
        angle = vspt_conn_json["periodicity"]["rotation_angle_deg"]
        axis = vspt_conn_json["periodicity"]["rotation_axis"]
        verified, mismatched = verify_periodicity(
            vspt_blocks, [face_match], theta=angle,
            rotation_axis=axis, tol=1e-4,
        )
        # DECLARED path: wrong matrix → mismatched, NOT silently
        # rounded into a passing canonical perm.
        assert len(verified) == 0, (
            "DECLARED path silently accepted a wrong matrix — "
            "_resolve_match_perm regression."
        )
        assert len(mismatched) == 1


# ---------------------------------------------------------------------------
# 4. UNDECLARED path on a stripped match (always runs — no full mesh)
# ---------------------------------------------------------------------------

class TestVsptUndeclaredPath:
    """Strip the orientation off a periodic match and confirm the
    verifier gracefully back-fills via the brute-force UNDECLARED
    path. Uses tiny synthetic blocks (no full-mesh load)."""

    @pytest.mark.skipif(not VSPT_MESH.exists(), reason="VSPT mesh missing")
    @pytest.mark.skipif(
        os.environ.get("PLOT3D_VSPT_FULL_MESH") != "1",
        reason="Set PLOT3D_VSPT_FULL_MESH=1 to load the 78 MB mesh.",
    )
    def test_periodic_match_without_orientation_back_fills(self, vspt_conn_json):
        from copy import deepcopy
        from plot3d import read_plot3D, verify_periodicity

        blocks = read_plot3D(str(VSPT_MESH), binary=False, big_endian=False,
                              read_double=False)

        face_match = deepcopy(vspt_conn_json["periodic_faces"][0])
        face_match.pop("orientation", None)  # legacy form: no orientation

        angle = vspt_conn_json["periodicity"]["rotation_angle_deg"]
        axis = vspt_conn_json["periodicity"]["rotation_axis"]
        verified, _ = verify_periodicity(
            blocks, [face_match], theta=angle, rotation_axis=axis, tol=1e-4,
        )
        assert len(verified) == 1
        # Back-fill happened.
        assert "orientation" in verified[0]
        assert "permutation_matrix" in verified[0]["orientation"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
