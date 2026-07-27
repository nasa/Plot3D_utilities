"""Verify Python's connectivity / verifier matches Fortran's CMC009 reference.

Fortran's reference connectivity for CMC009 lives at
``/Users/pjuangph/CMC009/shared/connectivity.json`` (converted from
the original ``LES009_444_new_diced_final.p3d_conn``). Counts:

* 593 blocks
* 1742 face_matches (block-to-block, no periodicity declared)
* 74 outer_faces
* 0 periodic_faces (CMC009 is a linear cascade; trans/rot periodicity
  is synthesised by ``plot3d.translational_periodicity`` if needed)

Two test layers:

1. **Fast**: assert structural counts and lb/ub well-formedness from
   the JSON alone — no mesh load.
2. **Slow / opt-in**: load the 1.1 GB ``RANS_009_refined2.p3d`` mesh
   and run ``connectivity_fast`` end-to-end. Gated on the env var
   ``PLOT3D_CMC009_FULL_MESH=1`` so default test runs (and CI without
   the fixture) skip it.

Even though Fortran's connectivity.json may not declare cross-plane
``permutation_index`` for every match (legacy format uses ``-1``
sentinel for in-plane), the verifier should still accept those
matches via the UNDECLARED brute-force path.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


CMC009_CONN_JSON = Path("/Users/pjuangph/CMC009/shared/connectivity.json")
CMC009_MESH = Path("/Users/pjuangph/CMC009/mesh/RANS_009_refined2.p3d")


@pytest.fixture(scope="module")
def cmc009_conn():
    """Loaded connectivity.json (JSON-converted Fortran reference)."""
    if not CMC009_CONN_JSON.exists():
        pytest.skip(f"CMC009 connectivity.json not at {CMC009_CONN_JSON}")
    with open(CMC009_CONN_JSON) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Structural assertions — fast, JSON-only
# ---------------------------------------------------------------------------

class TestCmc009ConnectivityJsonStructure:
    """Check the JSON Fortran-reference matches what we know about CMC009."""

    def test_block_and_match_counts(self, cmc009_conn):
        assert cmc009_conn["nblocks"] == 593
        assert len(cmc009_conn["face_matches"]) == 1742
        assert len(cmc009_conn["outer_faces"]) == 74
        # CMC009 is a linear cascade — Fortran's connectivity.json
        # leaves periodic_faces empty (no ROTATIONAL periodicity).
        # The 593 Z-span + 31 Y-pitch translational matches are
        # mixed into face_matches as same-block self-loops; the Rust
        # cascade verifier classifies them at load time. (See memory
        # note `project_partner_import_ghost_2026-04-29.md`.)
        assert len(cmc009_conn.get("periodic_faces", [])) == 0
        assert cmc009_conn["periodicity"]["nblades"] == 0

    def test_zones(self, cmc009_conn):
        zones = cmc009_conn["zones"]
        # CMC009 is single-zone (all fluid). The dumper's exact
        # interpretation of `ids` length varies — some emit one entry
        # per block, some per face, some per cell — but every entry
        # should be the fluid zone id (1).
        assert zones["count"] == 1
        assert len(zones["ids"]) > 0
        assert all(z == 1 for z in zones["ids"])
        assert all(t == 1 for t in zones["types"])

    def test_translational_self_loops_present(self, cmc009_conn):
        """CMC009's Fortran-derived face_matches mix in 593 Z-span
        same-block self-loops (block1.block_index == block2.block_index).
        The 31 Y-pitch matches live among the cross-block entries.

        Cascade-verified breakdown from memory
        `project_partner_import_ghost_2026-04-29.md`:
            1742 = 1118 cross-block + 31 Y-trans (cross-block) + 593 Z-trans (self-loops).
        """
        face_matches = cmc009_conn["face_matches"]
        self_loops = [
            fm for fm in face_matches
            if fm["block1"]["block_index"] == fm["block2"]["block_index"]
        ]
        cross_block = len(face_matches) - len(self_loops)
        # 593 Z-span self-loops (top↔bottom face within the same block).
        assert len(self_loops) == 593, (
            f"Expected 593 self-loop matches (Z-span periodicity); "
            f"got {len(self_loops)}."
        )
        # 1118 true cross-block + 31 Y-pitch cross-block = 1149.
        assert cross_block == 1149

    def test_face_match_corner_well_formedness(self, cmc009_conn):
        """Every face_match must have valid block_index + 3-component lb/ub."""
        nblocks = cmc009_conn["nblocks"]
        for idx, fm in enumerate(cmc009_conn["face_matches"]):
            for side_name in ("block1", "block2"):
                side = fm[side_name]
                assert "block_index" in side
                assert 0 <= side["block_index"] < nblocks, (
                    f"face_match[{idx}].{side_name}.block_index "
                    f"= {side['block_index']} out of [0, {nblocks})"
                )
                assert len(side["lb"]) == 3 and len(side["ub"]) == 3
                # CMC009 blocks are 33×33×33 → lb/ub in [0, 32]; allow up
                # to 32 inclusive on the high side.
                for d in range(3):
                    assert 0 <= min(side["lb"][d], side["ub"][d]) <= 32
                    assert 0 <= max(side["lb"][d], side["ub"][d]) <= 32

    def test_face_match_is_a_face_not_volume(self, cmc009_conn):
        """Every match must hold at least one constant axis (lb==ub on
        exactly one of i, j, k for both sides)."""
        for idx, fm in enumerate(cmc009_conn["face_matches"]):
            for side_name in ("block1", "block2"):
                lb = fm[side_name]["lb"]
                ub = fm[side_name]["ub"]
                fixed = sum(1 for d in range(3) if lb[d] == ub[d])
                assert fixed >= 1, (
                    f"face_match[{idx}].{side_name} has no constant axis: "
                    f"lb={lb} ub={ub}"
                )

    def test_cross_block_matches_count(self, cmc009_conn):
        """The 1149 cross-block (non-self-loop) matches sum of:
        1118 true cross-block (no periodicity) + 31 Y-pitch
        translational pairs across distinct blocks. The 593 Z-trans
        self-loops are handled in
        :meth:`test_translational_self_loops_present`.
        """
        cross_block_count = 0
        for fm in cmc009_conn["face_matches"]:
            if fm["block1"]["block_index"] != fm["block2"]["block_index"]:
                cross_block_count += 1
        assert cross_block_count == 1149


# ---------------------------------------------------------------------------
# UNDECLARED-path verifier — even without orientation, accept the match
# ---------------------------------------------------------------------------

class TestCmc009UndeclaredPath:
    """CMC009's JSON face_matches do NOT carry an `orientation` field.
    The verifier must treat them as UNDECLARED and back-fill via
    geometry-driven brute force.
    """

    def test_no_orientation_field_in_fortran_face_matches(self, cmc009_conn):
        """Confirm the precondition: Fortran-derived JSON didn't store
        orientation. The verifier MUST handle this."""
        no_orient_count = sum(
            1 for fm in cmc009_conn["face_matches"]
            if "orientation" not in fm
        )
        # Everything from Fortran's .p3d_conn comes through without
        # orientation by default (legacy format).
        assert no_orient_count == len(cmc009_conn["face_matches"]), (
            "Expected ALL CMC009 face_matches to lack orientation "
            "(Fortran ght_conn legacy format)."
        )


# ---------------------------------------------------------------------------
# Slow / opt-in — full mesh load + connectivity_fast
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    os.environ.get("PLOT3D_CMC009_FULL_MESH") != "1",
    reason="Full CMC009 mesh load is slow (1.1 GB ASCII). "
           "Set PLOT3D_CMC009_FULL_MESH=1 to opt in.",
)
@pytest.mark.skipif(not CMC009_MESH.exists(), reason="CMC009 mesh fixture missing")
class TestCmc009FullMeshIntegration:
    """End-to-end: load mesh, run connectivity_fast, compare to Fortran."""

    def test_connectivity_fast_finds_fortran_face_matches(self, cmc009_conn):
        from plot3d import read_plot3D, connectivity_fast

        blocks = read_plot3D(str(CMC009_MESH), binary=False, big_endian=False,
                              read_double=False)
        assert len(blocks) == 593

        # Discover face matches from geometry alone.
        py_matches, _ = connectivity_fast(blocks, use_minmax=False)

        # The number Python finds need not equal Fortran's 1742 byte-for-byte
        # (different tolerance + different match-selection rules), but it
        # should be in the same ballpark.
        assert 0.95 * 1742 <= len(py_matches) <= 1.05 * 1742, (
            f"Python connectivity_fast found {len(py_matches)} matches; "
            f"Fortran reference has 1742. >5% difference suggests a bug."
        )

        # Every Fortran match (by canonicalised block_index pair + face
        # location) should appear in Python's output.
        fortran_keys = _face_match_keys(cmc009_conn["face_matches"])
        python_keys = _face_match_keys(py_matches)
        missing = fortran_keys - python_keys
        assert not missing, (
            f"{len(missing)} Fortran face_matches have no Python equivalent. "
            f"First 3: {sorted(missing)[:3]}"
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _face_match_keys(face_matches):
    """Canonical key for a face_match: sorted (block, lb, ub) tuples.

    Two matches that differ only in which side is called "block1" vs
    "block2" map to the same key. Used for set-difference comparison
    between Python's discovered matches and Fortran's declared ones.
    """
    keys = set()
    for fm in face_matches:
        b1 = fm["block1"]
        b2 = fm["block2"]
        side_a = (b1["block_index"], tuple(b1["lb"]), tuple(b1["ub"]))
        side_b = (b2["block_index"], tuple(b2["lb"]), tuple(b2["ub"]))
        # Canonical: sort the two sides.
        keys.add(tuple(sorted([side_a, side_b])))
    return keys


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
