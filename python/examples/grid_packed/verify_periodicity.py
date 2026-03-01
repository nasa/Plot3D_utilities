"""Verify periodic face matches using full directed traversal.

Checks:
  1. Face dimensions match (both faces must have the same shape).
  2. Every point from face1 (traversing lb→ub), rotated by ±theta about x,
     matches the corresponding point from face2 (traversing lb→ub) — no sorting.

The diagonal convention (lb/ub) encodes the traversal orientation.
Point n from face A's lb→ub traversal must equal point n from face B's
lb→ub traversal (after rotation). If they don't match, the solver will
map data to the wrong grid points.
"""

import json
import pickle
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent  # grid_packed dir
TOLERANCE = 1e-4  # periodicity tolerance (matches Rust verify_periodicity)
NBLADES = 44
THETA = np.radians(360.0 / NBLADES)


def rotate_x(pts, angle):
    """Rotate points about x-axis. pts is (N,3) or (3,)."""
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[1, 0, 0],
                  [0, c, -s],
                  [0, s,  c]])
    return pts @ R.T


def _directed(start, end):
    """Inclusive range from start to end, stepping +1 or -1."""
    if start <= end:
        return range(start, end + 1)
    else:
        return range(start, end - 1, -1)


def face_dims(lb, ub):
    """Return (di, dj, dk) face dimensions from diagonal corners."""
    return (abs(ub[0] - lb[0]) + 1,
            abs(ub[1] - lb[1]) + 1,
            abs(ub[2] - lb[2]) + 1)


def extract_face_points(blk, lb, ub):
    """Extract all face points in directed traversal order (lb → ub).

    Returns an (N, 3) array where point n from face A must match
    point n from face B — the diagonal convention preserves the
    node-to-node mapping between blocks.
    """
    pts = []
    for i in _directed(lb[0], ub[0]):
        for j in _directed(lb[1], ub[1]):
            for k in _directed(lb[2], ub[2]):
                pts.append([blk.X[i, j, k], blk.Y[i, j, k], blk.Z[i, j, k]])
    return np.array(pts)


# ── 1. Load blocks ──────────────────────────────────────────────────
pickle_path = SCRIPT_DIR / 'grid_packed_binary.tmp.tmp.pkl'
binary_path = SCRIPT_DIR / 'grid_packed_binary.tmp.tmp.p3d'

if pickle_path.exists():
    print(f"Loading blocks from {pickle_path.name}...")
    with open(pickle_path, 'rb') as f:
        blocks = pickle.load(f)
else:
    from plot3d import read_plot3D
    print(f"Pickle not found, reading from {binary_path.name}...")
    blocks = read_plot3D(str(binary_path), binary=True)
print(f"Loaded {len(blocks)} blocks")

# ── 2. Load connectivity JSON ──────────────────────────────────────
json_path = SCRIPT_DIR / 'connectivity.json'
print(f"Loading connectivity from {json_path.name}...")
with open(json_path) as f:
    data = json.load(f)

# ── 3. Verify connectivity face matches (full point-by-point) ─────
face_matches = data['connectivity_face_matches']
print(f"\nVerifying {len(face_matches)} connectivity face matches...")

conn_dim_fails = 0
conn_point_fails = 0
conn_passed = 0

for i, match in enumerate(face_matches):
    b1, b2 = match['block1'], match['block2']
    blk1, blk2 = blocks[b1['block_index']], blocks[b2['block_index']]

    dims1 = face_dims(b1['lb'], b1['ub'])
    dims2 = face_dims(b2['lb'], b2['ub'])

    if dims1 != dims2:
        conn_dim_fails += 1
        print(f"  ERROR dim mismatch: match #{i} block {b1['block_index']}<->{b2['block_index']} "
              f"dims {dims1} vs {dims2}")
        continue

    pts1 = extract_face_points(blk1, b1['lb'], b1['ub'])
    pts2 = extract_face_points(blk2, b2['lb'], b2['ub'])
    diffs = np.linalg.norm(pts1 - pts2, axis=1)

    if diffs.max() > 1e-6:
        conn_point_fails += 1
        n_bad = int(np.sum(diffs > 1e-6))
        print(f"  ERROR point mismatch: match #{i} block {b1['block_index']}<->{b2['block_index']} "
              f"{n_bad}/{len(pts1)} points, max_dist={diffs.max():.6e}")
    else:
        conn_passed += 1

print(f"Connectivity: {conn_passed} passed, {conn_dim_fails} dim failures, {conn_point_fails} point failures")

# ── 4. Verify periodic face matches (full point-by-point) ─────────
periodic_matches = data.get('periodic_face_matches', [])
print(f"\nVerifying {len(periodic_matches)} periodic face matches "
      f"(rotation {360.0 / NBLADES:.4f} deg about x)...")

dim_failures = []
point_failures = []
passed = 0
max_dist = 0.0

for i, match in enumerate(periodic_matches):
    b1, b2 = match['block1'], match['block2']
    blk1, blk2 = blocks[b1['block_index']], blocks[b2['block_index']]

    dims1 = face_dims(b1['lb'], b1['ub'])
    dims2 = face_dims(b2['lb'], b2['ub'])

    # ── Check 1: Face dimensions must match ──
    if dims1 != dims2:
        dim_failures.append({
            'index': i,
            'block1_idx': b1['block_index'],
            'block2_idx': b2['block_index'],
            'b1_lb': b1['lb'], 'b1_ub': b1['ub'],
            'b2_lb': b2['lb'], 'b2_ub': b2['ub'],
            'dims1': dims1, 'dims2': dims2,
        })
        continue

    # ── Check 2: Every point must match after rotation ──
    pts1 = extract_face_points(blk1, b1['lb'], b1['ub'])
    pts2 = extract_face_points(blk2, b2['lb'], b2['ub'])

    # Try +theta and -theta, pick whichever is better
    pts1_pos = rotate_x(pts1, THETA)
    pts1_neg = rotate_x(pts1, -THETA)

    diffs_pos = np.linalg.norm(pts1_pos - pts2, axis=1)
    diffs_neg = np.linalg.norm(pts1_neg - pts2, axis=1)

    worst_pos = float(diffs_pos.max())
    worst_neg = float(diffs_neg.max())

    if worst_pos <= worst_neg:
        best_diffs, best_worst, best_sign = diffs_pos, worst_pos, '+'
    else:
        best_diffs, best_worst, best_sign = diffs_neg, worst_neg, '-'

    max_dist = max(max_dist, best_worst)

    if best_worst > TOLERANCE:
        n_bad = int(np.sum(best_diffs > TOLERANCE))
        point_failures.append({
            'index': i,
            'block1_idx': b1['block_index'],
            'block2_idx': b2['block_index'],
            'b1_lb': b1['lb'], 'b1_ub': b1['ub'],
            'b2_lb': b2['lb'], 'b2_ub': b2['ub'],
            'dims': dims1,
            'total_points': len(pts1),
            'mismatched_points': n_bad,
            'max_dist': best_worst,
            'mean_dist': float(best_diffs[best_diffs > TOLERANCE].mean()),
            'sign': best_sign,
        })
    else:
        passed += 1


# ── 5. Report ──────────────────────────────────────────────────────
total = len(periodic_matches)
print(f"\n{'='*70}")
print(f"  Periodic Verification Results (full point-by-point)")
print(f"  Rotation: {360.0 / NBLADES:.4f} deg about x-axis")
print(f"{'='*70}")
print(f"  Total periodic matches:  {total}")
print(f"  Passed:                  {passed}")
print(f"  Dimension mismatches:    {len(dim_failures)}")
print(f"  Point mismatches:        {len(point_failures)}")
print(f"  Tolerance:               {TOLERANCE:.1e}")
print(f"  Max distance found:      {max_dist:.2e}")

if not dim_failures and not point_failures:
    print("\n  All periodic face points match within tolerance (directed traversal).")
else:
    if dim_failures:
        print(f"\n  ── ERROR: Face Dimension Mismatches ({len(dim_failures)}) ──")
        print(f"  The two faces in a periodic match must have the same shape.")
        for f in dim_failures:
            print(f"\n  Match #{f['index']:>5d}")
            print(f"    block {f['block1_idx']}: lb={f['b1_lb']} ub={f['b1_ub']} dims={f['dims1']}")
            print(f"    block {f['block2_idx']}: lb={f['b2_lb']} ub={f['b2_ub']} dims={f['dims2']}")

    if point_failures:
        print(f"\n  ── ERROR: Point Mismatches ({len(point_failures)}) ──")
        print(f"  Points traversed in lb→ub order (with rotation) do not match.")
        print(f"  This means the diagonal orientation is wrong or the faces")
        print(f"  do not share the same geometry after rotation.")
        for f in point_failures:
            pct = 100.0 * f['mismatched_points'] / f['total_points']
            print(f"\n  Match #{f['index']:>5d}")
            print(f"    block {f['block1_idx']}: lb={f['b1_lb']} ub={f['b1_ub']}")
            print(f"    block {f['block2_idx']}: lb={f['b2_lb']} ub={f['b2_ub']}")
            print(f"    face dims: {f['dims']}, total points: {f['total_points']}")
            print(f"    mismatched: {f['mismatched_points']} ({pct:.1f}%)")
            print(f"    max dist: {f['max_dist']:.6e}, mean dist: {f['mean_dist']:.6e}")
            print(f"    best rotation: {f['sign']}theta")
