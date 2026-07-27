"""GlennHT pre-run validation checks.

Validates exported mesh and connectivity against the same logic
GlennHT uses internally — catching failures before a crashed run.

Functions can be used individually or via ``run_all_checks()``.
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 1. Handedness check
# ---------------------------------------------------------------------------
def check_handedness(blocks: list, block_names: list[str]) -> bool:
    """Check block handedness (informational — GlennHT handles both).

    Returns:
        Always True (informational check).
    """
    print("\n" + "=" * 60)
    print("1. HANDEDNESS (informational — GlennHT handles both)")
    print("=" * 60)
    for i, b in enumerate(blocks):
        if hasattr(b, 'check_handedness'):
            rh = b.check_handedness()
        else:
            vol = b.cell_volumes()
            rh = bool(np.all(vol[1:, 1:, 1:] >= 0))
        status = "right-handed" if rh else "left-handed (OK — GlennHT handles this)"
        print(f"  {block_names[i]} ({b.IMAX}x{b.JMAX}x{b.KMAX}): {status}")
    return True


# ---------------------------------------------------------------------------
# 2. Corner matching
# ---------------------------------------------------------------------------
def check_corners(blocks: list, face_matches: list[dict],
                  periodic_matches: list[dict],
                  block_names: list[str]) -> bool:
    """Verify that lower-bound corners of matched faces coincide spatially.

    Returns:
        True (informational — verify_connectivity is authoritative).
    """
    print("\n" + "=" * 60)
    print("2. CORNER MATCHING (lb/ub corners)")
    print("=" * 60)
    tol = 1e-6
    for label, matches in [("Interface", face_matches), ("Periodic", periodic_matches)]:
        for i, m in enumerate(matches):
            b1, b2 = m["block1"], m["block2"]
            bi1, lb1, ub1 = b1["block_index"], b1["lb"], b1["ub"]
            bi2, lb2, ub2 = b2["block_index"], b2["lb"], b2["ub"]

            def _xyz(bi, ijk):
                return np.array([blocks[bi].X[ijk[0], ijk[1], ijk[2]],
                                 blocks[bi].Y[ijk[0], ijk[1], ijk[2]],
                                 blocks[bi].Z[ijk[0], ijk[1], ijk[2]]])

            c1_lb, c1_ub = _xyz(bi1, lb1), _xyz(bi1, ub1)
            c2_lb, c2_ub = _xyz(bi2, lb2), _xyz(bi2, ub2)
            d_direct = max(np.linalg.norm(c1_lb - c2_lb), np.linalg.norm(c1_ub - c2_ub))
            d_cross = max(np.linalg.norm(c1_lb - c2_ub), np.linalg.norm(c1_ub - c2_lb))
            n1, n2 = block_names[bi1], block_names[bi2]
            best = min(d_direct, d_cross)
            if best < tol:
                extra = " (reversed)" if d_cross < d_direct else ""
                print(f"  {label} {i+1}: {n1}<->{n2} OK{extra}")
            else:
                print(f"  {label} {i+1}: {n1}<->{n2} cross-axis (best corner dist={best:.2e})")
    return True


# ---------------------------------------------------------------------------
# 3. verify_connectivity
# ---------------------------------------------------------------------------
def check_verify_connectivity(blocks: list, face_matches: list[dict],
                              block_names: list[str]) -> bool:
    """Run plot3d's point-by-point connectivity verification.

    Returns:
        True if all face matches are verified with no mismatches.
    """
    from plot3d import verify_connectivity
    print("\n" + "=" * 60)
    print("3. VERIFY_CONNECTIVITY (point-by-point)")
    print("=" * 60)
    verified, mismatched = verify_connectivity(blocks, face_matches)
    print(f"  Verified: {len(verified)}, Mismatched: {len(mismatched)}")
    for m in mismatched:
        n1 = block_names[m["block1"]["block_index"]]
        n2 = block_names[m["block2"]["block_index"]]
        print(f"    FAIL: {n1} <-> {n2}")
    return len(mismatched) == 0


# ---------------------------------------------------------------------------
# 4. verify_periodicity
# ---------------------------------------------------------------------------
def check_verify_periodicity(blocks: list, periodic_matches: list[dict],
                             pitch_deg: float,
                             block_names: list[str]) -> bool:
    """Run plot3d's point-by-point periodicity verification with rotation.

    Returns:
        True if all periodic matches are verified with no mismatches.
    """
    from plot3d import verify_periodicity
    print("\n" + "=" * 60)
    print(f"4. VERIFY_PERIODICITY (pitch={pitch_deg:.4f} deg)")
    print("=" * 60)
    verified, mismatched = verify_periodicity(
        blocks, periodic_matches, pitch_deg, rotation_axis='x',
    )
    print(f"  Verified: {len(verified)}, Mismatched: {len(mismatched)}")
    for m in mismatched:
        n1 = block_names[m["block1"]["block_index"]]
        n2 = block_names[m["block2"]["block_index"]]
        print(f"    FAIL: {n1} <-> {n2}")
    return len(mismatched) == 0


# ---------------------------------------------------------------------------
# 5. GlennHT PM integer division check
# ---------------------------------------------------------------------------
def compute_pm(aA: list[int], bA: list[int],
               aB: list[int], bB: list[int]):
    """Compute the PM permutation/orientation matrix for a face match.

    Reproduces ``new_PTtransf`` from GlennHT's ``M_IndexVect.F``.

    Returns:
        (pm, (rem1, rem2)) on success, or (None, error_string) on failure.
    """
    va, vb = np.array(aA), np.array(bA)
    vc, vd = np.array(aB), np.array(bB)
    n1, n2 = np.zeros(3, dtype=int), np.zeros(3, dtype=int)
    for i in range(3):
        if va[i] == vb[i]:
            n1[i] = -1 if va[i] == 1 else 1
        if vc[i] == vd[i]:
            n2[i] = -1 if vc[i] == 1 else 1

    fVect = np.array([1, 2, 3])
    face1 = abs(np.sum(fVect * n1))
    face2 = abs(np.sum(fVect * n2))
    if face1 == 0 or face2 == 0:
        return None, "no constant axis"

    pf1f2 = -n1[face1 - 1] * n2[face2 - 1]
    s = ((-1) ** (face1 + face2)) * pf1f2
    d1 = vb - va
    d2 = vd - vc
    d = int(np.sum(d1 * d1))
    if d == 0:
        return None, "zero diagonal"

    idx = {1: (1, 2), 2: (0, 2), 3: (0, 1)}
    i11, i12 = idx[face1]
    i21, i22 = idx[face2]

    m = np.zeros((3, 3), dtype=int)
    m[i21, i11] = (d1[i11] * d2[i21] + s * d1[i12] * d2[i22]) // d
    m[i21, i12] = (d1[i12] * d2[i21] - s * d1[i11] * d2[i22]) // d
    m[i22, i11] = -s * m[i21, i12]
    m[i22, i12] = s * m[i21, i11]
    m[face2 - 1, face1 - 1] = pf1f2

    rem1 = (d1[i11] * d2[i21] + s * d1[i12] * d2[i22]) % d
    rem2 = (d1[i12] * d2[i21] - s * d1[i11] * d2[i22]) % d

    return m, (rem1, rem2)


def check_pm(blocks: list, face_matches: list[dict],
             periodic_matches: list[dict],
             block_names: list[str]) -> bool:
    """Validate GlennHT's PM transformation matrices for all face matches.

    Returns:
        True if all PM matrices are exact and map within bounds.
    """
    print("\n" + "=" * 60)
    print("5. INDEX TRANSFORM (GlennHT new_PTtransf)")
    print("=" * 60)
    all_ok = True
    all_matches = (
        [(m, "Interface") for m in face_matches] +
        [(m, "Periodic") for m in periodic_matches]
    )
    for idx, (m, label) in enumerate(all_matches):
        b1, b2 = m["block1"], m["block2"]
        aA = [b1["lb"][d] + 1 for d in range(3)]
        bA = [b1["ub"][d] + 1 for d in range(3)]
        aB = [b2["lb"][d] + 1 for d in range(3)]
        bB = [b2["ub"][d] + 1 for d in range(3)]
        n1 = block_names[b1["block_index"]]
        n2 = block_names[b2["block_index"]]

        pm, info = compute_pm(aA, bA, aB, bB)
        if pm is None:
            print(f"  {label} {idx+1}: {n1}<->{n2} PM FAILED: {info}")
            all_ok = False
            continue

        rem1, rem2 = info
        if rem1 != 0 or rem2 != 0:
            print(f"  {label} {idx+1}: {n1}<->{n2} FAIL: PM truncated (remainders={rem1},{rem2})")
            print(f"    PM = {pm.tolist()}")
            print(f"    Corners: {aA}->{bA} ↔ {aB}->{bB}")
            all_ok = False
        else:
            ilo = np.minimum(aA, bA)
            ihi = np.maximum(aA, bA)
            Iref, Jref = np.array(aA), np.array(aB)
            corners = [ilo, ihi]
            oob = False
            bi2 = b2["block_index"]
            bshape = [blocks[bi2].IMAX, blocks[bi2].JMAX, blocks[bi2].KMAX]
            for c in corners:
                mapped = Jref + pm @ (np.array(c) - Iref)
                for dd in range(3):
                    if mapped[dd] < 1 or mapped[dd] > bshape[dd]:
                        oob = True
            if oob:
                print(f"  {label} {idx+1}: {n1}<->{n2} PM OUT-OF-BOUNDS")
                all_ok = False
            else:
                print(f"  {label} {idx+1}: {n1}<->{n2} PM OK {pm.tolist()}")
    return all_ok


# ---------------------------------------------------------------------------
# 6. Face coverage
# ---------------------------------------------------------------------------
def check_face_coverage(blocks: list, face_matches: list[dict],
                        periodic_matches: list[dict],
                        outer_faces: list[dict],
                        block_names: list[str]) -> bool:
    """Check that every block face is accounted for by matches or BCs.

    Returns:
        True if no block face is completely uncovered.
    """
    print("\n" + "=" * 60)
    print("6. FACE COVERAGE (no orphan faces)")
    print("=" * 60)
    all_ok = True
    for bi, b in enumerate(blocks):
        ni, nj, nk = b.IMAX, b.JMAX, b.KMAX
        face_nodes = {
            "imin": nj * nk, "imax": nj * nk,
            "jmin": ni * nk, "jmax": ni * nk,
            "kmin": ni * nj, "kmax": ni * nj,
        }
        covered = {f: 0 for f in face_nodes}

        all_m = face_matches + periodic_matches
        for m_item in all_m:
            for side in ("block1", "block2"):
                f = m_item[side]
                if f["block_index"] != bi:
                    continue
                lb, ub = f["lb"], f["ub"]
                for axis in range(3):
                    mn, mx = min(lb[axis], ub[axis]), max(lb[axis], ub[axis])
                    dims = [ni, nj, nk]
                    if mn == mx:
                        if mn == 0:
                            fname = ["imin", "jmin", "kmin"][axis]
                        elif mn == dims[axis] - 1:
                            fname = ["imax", "jmax", "kmax"][axis]
                        else:
                            continue
                        other_axes = [d for d in range(3) if d != axis]
                        n1 = abs(ub[other_axes[0]] - lb[other_axes[0]]) + 1
                        n2 = abs(ub[other_axes[1]] - lb[other_axes[1]]) + 1
                        covered[fname] += n1 * n2

        for of in outer_faces:
            if of["block_index"] != bi:
                continue
            lb, ub = of["lb"], of["ub"]
            for axis in range(3):
                mn, mx = min(lb[axis], ub[axis]), max(lb[axis], ub[axis])
                dims = [ni, nj, nk]
                if mn == mx:
                    if mn == 0:
                        fname = ["imin", "jmin", "kmin"][axis]
                    elif mn == dims[axis] - 1:
                        fname = ["imax", "jmax", "kmax"][axis]
                    else:
                        continue
                    other_axes = [d for d in range(3) if d != axis]
                    n1 = abs(ub[other_axes[0]] - lb[other_axes[0]]) + 1
                    n2 = abs(ub[other_axes[1]] - lb[other_axes[1]]) + 1
                    covered[fname] += n1 * n2

        name = block_names[bi]
        uncovered = [f for f in face_nodes if covered[f] == 0]
        if uncovered:
            print(f"  {name}: UNCOVERED faces: {uncovered}")
            all_ok = False
        else:
            print(f"  {name}: all faces covered")
    return all_ok


# ---------------------------------------------------------------------------
# 7. Negative cell-volume check
# ---------------------------------------------------------------------------
def check_negative_volumes(blocks: list, block_names: list[str]) -> bool:
    """Compute cell volumes and report any with negative volume.

    Returns:
        True if no negative-volume cells are found.
    """
    print("\n" + "=" * 60)
    print("7. NEGATIVE CELL VOLUMES (Davies & Salmond)")
    print("=" * 60)
    all_ok = True
    for bi, b in enumerate(blocks):
        name = block_names[bi]
        print(f"  {name} ({b.IMAX}x{b.JMAX}x{b.KMAX}): computing volumes...", end="")
        vol = b.cell_volumes()
        interior = vol[1:, 1:, 1:]
        neg_mask = interior < 0
        n_neg = int(np.sum(neg_mask))
        total = interior.size
        if n_neg == 0:
            print(f" OK (min={interior.min():.4e})")
        elif n_neg == total:
            print(f" left-handed ({total} cells, GlennHT handles this)")
        else:
            all_ok = False
            pct = 100.0 * n_neg / total
            print(f" FAIL — {n_neg}/{total} cells negative ({pct:.1f}%)")
            locs = np.argwhere(neg_mask)
            for loc in locs[:5]:
                ii, jj, kk = loc[0] + 1, loc[1] + 1, loc[2] + 1
                print(f"    neg vol at (i={ii}, j={jj}, k={kk}): {interior[loc[0], loc[1], loc[2]]:.4e}")
            if len(locs) > 5:
                print(f"    ... and {len(locs) - 5} more")
    return all_ok


# ---------------------------------------------------------------------------
# 8. Face normals
# ---------------------------------------------------------------------------
def check_face_normals(blocks: list, block_names: list[str]) -> bool:
    """Print the average outward face normal for each block face.

    Returns:
        Always True (informational check).
    """
    print("\n" + "=" * 60)
    print("8. FACE NORMALS (average outward normal per face)")
    print("=" * 60)
    for bi, b in enumerate(blocks):
        name = block_names[bi]
        ni, nj, nk = b.IMAX, b.JMAX, b.KMAX
        X, Y, Z = b.X, b.Y, b.Z
        print(f"  {name} ({ni}x{nj}x{nk}):")

        def _avg_normal(pts1, pts2, pts3):
            v1 = np.stack([pts2[0] - pts1[0], pts2[1] - pts1[1], pts2[2] - pts1[2]], axis=-1)
            v2 = np.stack([pts3[0] - pts1[0], pts3[1] - pts1[1], pts3[2] - pts1[2]], axis=-1)
            normals = np.cross(v1, v2)
            avg = normals.mean(axis=tuple(range(normals.ndim - 1)))
            nrm = np.linalg.norm(avg)
            return avg / nrm if nrm > 0 else avg

        for idx, label in [(0, "imin"), (ni - 1, "imax")]:
            n = _avg_normal(
                (X[idx, :-1, :-1], Y[idx, :-1, :-1], Z[idx, :-1, :-1]),
                (X[idx, 1:, :-1], Y[idx, 1:, :-1], Z[idx, 1:, :-1]),
                (X[idx, :-1, 1:], Y[idx, :-1, 1:], Z[idx, :-1, 1:]),
            )
            print(f"    {label}: [{n[0]:+.4f}, {n[1]:+.4f}, {n[2]:+.4f}]")

        for idx, label in [(0, "jmin"), (nj - 1, "jmax")]:
            n = _avg_normal(
                (X[:-1, idx, :-1], Y[:-1, idx, :-1], Z[:-1, idx, :-1]),
                (X[:-1, idx, 1:], Y[:-1, idx, 1:], Z[:-1, idx, 1:]),
                (X[1:, idx, :-1], Y[1:, idx, :-1], Z[1:, idx, :-1]),
            )
            print(f"    {label}: [{n[0]:+.4f}, {n[1]:+.4f}, {n[2]:+.4f}]")

        for idx, label in [(0, "kmin"), (nk - 1, "kmax")]:
            n = _avg_normal(
                (X[:-1, :-1, idx], Y[:-1, :-1, idx], Z[:-1, :-1, idx]),
                (X[1:, :-1, idx], Y[1:, :-1, idx], Z[1:, :-1, idx]),
                (X[:-1, 1:, idx], Y[:-1, 1:, idx], Z[:-1, 1:, idx]),
            )
            print(f"    {label}: [{n[0]:+.4f}, {n[1]:+.4f}, {n[2]:+.4f}]")

    return True


# ---------------------------------------------------------------------------
# run_all_checks — convenience wrapper
# ---------------------------------------------------------------------------
def run_all_checks(blocks: list, block_names: list[str],
                   face_matches: list[dict],
                   periodic_matches: list[dict],
                   outer_faces: list[dict],
                   pitch_deg: float) -> dict[str, bool]:
    """Run all GlennHT pre-flight checks and return a pass/fail dict.

    Args:
        blocks: List of Plot3D Block objects.
        block_names: Names of each block.
        face_matches: Interface connectivity records.
        periodic_matches: Periodic connectivity records.
        outer_faces: Boundary-condition face records.
        pitch_deg: Blade pitch angle in degrees (360 / nblades).

    Returns:
        Dict mapping check name to True/False.
    """
    results = {}
    results["handedness"] = check_handedness(blocks, block_names)
    results["corners"] = check_corners(blocks, face_matches, periodic_matches, block_names)
    results["verify_conn"] = check_verify_connectivity(blocks, face_matches, block_names)
    results["verify_per"] = check_verify_periodicity(blocks, periodic_matches, pitch_deg, block_names)
    results["pm"] = check_pm(blocks, face_matches, periodic_matches, block_names)
    results["coverage"] = check_face_coverage(blocks, face_matches, periodic_matches, outer_faces, block_names)
    results["neg_volumes"] = check_negative_volumes(blocks, block_names)
    results["face_normals"] = check_face_normals(blocks, block_names)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    all_pass = True
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL"
        print(f"  {name}: {status}")
        if not ok:
            all_pass = False

    if all_pass:
        print("\nAll checks passed — ready for GlennHT!")
    else:
        print("\nSome checks FAILED — fix before running GlennHT.")

    return results
