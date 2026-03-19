"""Integration test for the WELD mesh: connectivity + translational periodicity.

Reads /Volumes/T7/WELD/weld_ascii.xyz (1900 blocks), runs connectivity
and z-direction translational periodicity, then compares against reference JSON.

Skipped when the mesh file is not present.
"""

import json
import os
import pickle
import pytest
from pathlib import Path

MESH_PATH = "/Volumes/T7/WELD/weld_ascii.xyz"
CONN_JSON = "/Volumes/T7/WELD/weld_connectivity.json"
CONN_PERIOD_JSON = "/Volumes/T7/WELD/weld_connectivity-periodicity.json"
CACHE_PATH = "/Volumes/T7/WELD/cache.pkl"

skip_no_mesh = pytest.mark.skipif(
    not os.path.exists(MESH_PATH),
    reason=f"WELD mesh not found at {MESH_PATH}"
)


def canon_face_key(entry):
    """Canonical (block_index, lo, hi) with ascending lo/hi."""
    if isinstance(entry.get('lb'), list):
        lb, ub = entry['lb'], entry['ub']
        lo = tuple(min(a, b) for a, b in zip(lb, ub))
        hi = tuple(max(a, b) for a, b in zip(lb, ub))
    else:
        lo, hi = tuple(entry['lo']), tuple(entry['hi'])
    return (entry['block_index'], lo, hi)


def canon_match_key(match):
    """Order-independent match key for set comparison."""
    a = canon_face_key(match['block1'])
    b = canon_face_key(match['block2'])
    return (a, b) if a <= b else (b, a)


def load_blocks():
    """Load blocks from pickle cache (fast) or ASCII mesh (slow)."""
    if os.path.exists(CACHE_PATH):
        print(f"Loading blocks from cache {CACHE_PATH}...")
        with open(CACHE_PATH, 'rb') as f:
            data = pickle.load(f)
        # Cache may be a dict with 'blocks' key or a plain list
        if isinstance(data, dict):
            return data['blocks']
        return data
    else:
        from plot3d import read_plot3D
        print(f"Reading WELD mesh from {MESH_PATH}...")
        blocks = read_plot3D(MESH_PATH, binary=False)
        # Cache for next run
        with open(CACHE_PATH, 'wb') as f:
            pickle.dump(blocks, f)
        return blocks


@skip_no_mesh
def test_weld_connectivity_and_periodicity():
    from plot3d import (connectivity_fast, translational_periodicity,
                        verify_connectivity, extract_canonical_grid,
                        apply_permutation, verify_match, try_all_permutations,
                        get_bounds)

    # ── 1. Load blocks ──
    blocks = load_blocks()
    assert len(blocks) == 1900, f"Expected 1900 blocks, got {len(blocks)}"
    print(f"Loaded {len(blocks)} blocks")

    # ── 2. Run connectivity ──
    print("Running connectivity_fast...")
    face_matches, outer_faces = connectivity_fast(blocks)
    print(f"  {len(face_matches)} face matches, {len(outer_faces)} outer faces")

    # ── 3. Verify connectivity ──
    print("Verifying connectivity (point-by-point)...")
    verified, mismatched = verify_connectivity(blocks, face_matches, tol=1e-6)
    print(f"  Verified: {len(verified)}, Mismatched: {len(mismatched)}")

    # ── 4. Compare against reference JSON ──
    if os.path.exists(CONN_JSON):
        with open(CONN_JSON) as f:
            ref = json.load(f)
        ref_matches = ref['face_matches']

        rust_keys = set(canon_match_key(m) for m in verified)
        ref_keys = set(canon_match_key(m) for m in ref_matches)
        in_both = rust_keys & ref_keys

        print(f"  Reference comparison: {len(ref_keys)} ref, {len(in_both)} in both, "
              f"{len(ref_keys - rust_keys)} only-ref, {len(rust_keys - ref_keys)} only-Python")

        assert len(ref_keys - rust_keys) == 0, (
            f"{len(ref_keys - rust_keys)} matches in reference but not in Python result")
        assert len(rust_keys - ref_keys) == 0, (
            f"{len(rust_keys - ref_keys)} matches in Python but not in reference")

    # ── 5. Run translational periodicity (z) ──
    print("Running translational periodicity (z)...")
    periodic_matches, _, remaining_outer = translational_periodicity(
        blocks, outer_faces,
        delta=None,
        translational_direction='z',
    )
    print(f"  {len(periodic_matches)} periodic pairs, {len(remaining_outer)} remaining outer")

    # ── 6. Compare periodic matches against reference ──
    if os.path.exists(CONN_PERIOD_JSON):
        with open(CONN_PERIOD_JSON) as f:
            ref_p = json.load(f)
        ref_periodic = ref_p.get('periodic_faces', [])

        py_periodic_keys = set(canon_match_key(m) for m in periodic_matches)
        ref_periodic_keys = set(canon_match_key(m) for m in ref_periodic)
        in_both_p = py_periodic_keys & ref_periodic_keys

        print(f"  Reference comparison: {len(ref_periodic_keys)} ref, {len(in_both_p)} in both, "
              f"{len(ref_periodic_keys - py_periodic_keys)} only-ref, "
              f"{len(py_periodic_keys - ref_periodic_keys)} only-Python")

        assert len(ref_periodic_keys - py_periodic_keys) == 0, (
            f"{len(ref_periodic_keys - py_periodic_keys)} periodic matches in reference but not in Python")
        assert len(py_periodic_keys - ref_periodic_keys) == 0, (
            f"{len(py_periodic_keys - ref_periodic_keys)} periodic matches in Python but not in reference")

    # ── 7. Verify new helper functions work on a sample ──
    print("Testing extract_canonical_grid + apply_permutation + verify_match...")
    if verified:
        sample = verified[0]
        b1, b2 = sample['block1'], sample['block2']
        lo1, hi1 = get_bounds(b1)
        lo2, hi2 = get_bounds(b2)

        grid_a, nu_a, nv_a = extract_canonical_grid(blocks[b1['block_index']],
                                                      list(lo1), list(hi1))
        grid_b, nu_b, nv_b = extract_canonical_grid(blocks[b2['block_index']],
                                                      list(lo2), list(hi2))

        perm_idx = try_all_permutations(grid_a, grid_b, 1e-6)
        assert perm_idx >= 0, "Should find a valid permutation for a verified match"

        g = apply_permutation(grid_b, perm_idx)
        assert verify_match(grid_a, g, 1e-6), "Permuted grid should match"
        print(f"  Sample match: block {b1['block_index']}<->{b2['block_index']}, perm={perm_idx} OK")

    # ── Summary ──
    print(f"\n=== WELD Python Test Summary ===")
    print(f"  Blocks: {len(blocks)}")
    print(f"  Connectivity: {len(face_matches)} raw -> {len(verified)} verified ({len(mismatched)} mismatched)")
    print(f"  Periodicity: {len(periodic_matches)} z-periodic pairs")
    print(f"  Remaining outer: {len(remaining_outer)}")
    print(f"  PASS")
