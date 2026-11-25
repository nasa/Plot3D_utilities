"""Folded-grid-line untangling: detection, repair, refusal.

Mirrors plot3d-rs ``tests/test_untangle.rs``. Fixtures are cartesian blocks
with a thin cell layer (a wall-resolved spacing) under the ``jmax`` face,
where moving one node beyond its ``j+1`` neighbour reverses interior faces
without inverting a cell.
"""

import numpy as np
import pytest

from plot3d.block import Block
from plot3d.mesh_quality import (
    UNTANGLE_REPORT_LIMIT,
    RepairConstruction,
    ReversedFace,
    UntangleError,
    UntangleReport,
    block_handedness,
    Handedness,
    find_inverted_geometry,
    untangle_folded_lines,
    untangle_mesh,
    untangle_plot3d_file,
)

YS = [0.0, 0.4, 0.7, 0.99, 1.0]


def _block(ni, ys, nk):
    i, j, k = np.meshgrid(np.arange(ni), np.arange(len(ys)), np.arange(nk), indexing="ij")
    return Block(0.25 * i.astype(float), np.asarray(ys, dtype=float)[j], 0.25 * k.astype(float))


def _with_y(b, node, y):
    Y = b.Y.copy()
    Y[node] = y
    return Block(b.X.copy(), Y, b.Z.copy())


def _p(b, n):
    return np.array([b.X[n], b.Y[n], b.Z[n]])


def test_a_folded_node_next_to_a_face_is_rebuilt_from_its_neighbouring_lines():
    clean = _block(5, YS, 5)
    folded = _with_y(clean, (2, 3, 2), 1.02)

    faces, cells = find_inverted_geometry(folded, 7)
    assert faces, "fixture must reverse interior faces, else the repair path is never entered"
    assert not cells
    assert all(f.block == 7 and f.cos <= 0.0 for f in faces)

    repaired, report = untangle_folded_lines(folded, 7)
    assert len(report.folded_edges) == 1
    assert (report.folded_edges[0].axis, report.folded_edges[0].start) == (1, (2, 3, 2))
    assert len(report.moved) == 1
    m = report.moved[0]
    assert (m.block, m.node) == (7, (2, 3, 2))
    assert m.construction == RepairConstruction.NEIGHBOUR_MEAN

    # new = X(anchor) - mean of the j-edge on the four neighbouring lines.
    def edge(n):
        return _p(folded, (n[0], n[1] + 1, n[2])) - _p(folded, n)

    mean = sum(0.25 * edge(n) for n in [(1, 3, 2), (3, 3, 2), (2, 3, 1), (2, 3, 3)])
    expected = _p(folded, (2, 4, 2)) - mean
    np.testing.assert_allclose(m.new, expected, rtol=0, atol=1e-15)
    np.testing.assert_allclose(_p(repaired, (2, 3, 2)), expected, rtol=0, atol=1e-15)
    np.testing.assert_allclose(expected, _p(clean, (2, 3, 2)), rtol=0, atol=1e-12)
    assert m.old == (0.5, 1.02, 0.5)
    assert abs(m.displacement - 0.03) < 1e-12
    assert report.max_displacement == m.displacement

    faces_after, cells_after = find_inverted_geometry(repaired, 7)
    assert not faces_after and not cells_after
    changed = (repaired.X != folded.X) | (repaired.Y != folded.Y) | (repaired.Z != folded.Z)
    assert list(zip(*np.nonzero(changed))) == [(2, 3, 2)]

    assert report.blocks_scanned == 1
    assert report.rechecked_blocks == [7]
    text = str(report)
    assert text.startswith("untangle folded grid lines: scanned 1 block(s);")
    assert "\nblock 7: " in text
    f = report.reversed_faces[0]
    nb = list(f.owner)
    nb[f.axis] += 1
    assert (
        f"cell ({f.owner[0]},{f.owner[1]},{f.owner[2]}) -> ({nb[0]},{nb[1]},{nb[2]}), "
        f"cos(a,d) {f.cos:.4e}; folded j-edge (2,3,2)-(2,4,2), node (2,3,2) moved"
    ) in text
    assert "  moved node (2,3,2): (" in text
    assert "  re-check: 0 reversed face(s), 0 non-positive cell(s)" in text
    assert text.count("  reversed ") == len(report.reversed_faces)


def test_a_clean_block_is_returned_bit_identical_with_an_empty_report():
    b = _block(5, YS, 4)
    X = b.X + 0.1 * b.Y + 0.03 * b.Z * b.Z
    b = Block(X, b.Y, b.Z)
    out, report = untangle_folded_lines(b, 0)
    assert report.is_empty()
    assert report.max_displacement == 0.0
    for u, v in ((out.X, b.X), (out.Y, b.Y), (out.Z, b.Z)):
        assert u is not v
        assert np.array_equal(u.view(np.uint64), v.view(np.uint64))


def test_untangle_mesh_repairs_only_the_folded_block_and_labels_it():
    clean = _block(5, YS, 5)
    folded = _with_y(clean, (2, 3, 2), 1.02)
    out, report = untangle_mesh([clean, folded])
    assert np.array_equal(out[0].Y, clean.Y)
    assert [m.block for m in report.moved] == [1]
    assert all(f.block == 1 for f in report.reversed_faces)
    assert "1 node(s) moved" in report.format_report()


def test_a_left_handed_block_is_judged_relative_to_its_handedness():
    folded = _with_y(_block(5, YS, 5), (2, 3, 2), 1.02)
    lh = Block(folded.X[::-1].copy(), folded.Y[::-1].copy(), folded.Z[::-1].copy())
    assert block_handedness(lh) == Handedness.LEFT_HANDED
    assert len(find_inverted_geometry(lh)[0]) == len(find_inverted_geometry(folded)[0])
    _, report = untangle_folded_lines(lh, 0)
    assert [m.node for m in report.moved] == [(2, 3, 2)]


def test_a_dented_face_under_a_thin_layer_is_repaired_on_the_chord():
    ys = [0.0, 0.5, 0.9, 0.98, 1.0]
    folded = _with_y(_block(5, ys, 5), (2, 4, 2), 0.91)
    assert find_inverted_geometry(folded)[0], "fixture must reverse an interior face"
    repaired, report = untangle_folded_lines(folded, 0)
    assert len(report.moved) == 1
    m = report.moved[0]
    assert m.node == (2, 3, 2)
    assert m.construction == RepairConstruction.CHORD_REDISTRIBUTION
    u, w = ys[4] - ys[3], ys[3] - ys[2]
    r = w / (w + u)
    f, a = _p(folded, (2, 2, 2)), _p(folded, (2, 4, 2))
    np.testing.assert_allclose(m.new, f + r * (a - f), rtol=0, atol=1e-14)
    assert 0.9 < m.new[1] < 0.91
    assert find_inverted_geometry(repaired) == ([], [])
    assert np.array_equal(_p(repaired, (2, 4, 2)), a)


def test_a_fold_on_a_block_face_is_refused():
    clean = _block(5, YS, 5)
    folded = _with_y(clean, (0, 3, 2), 1.02)
    faces, _ = find_inverted_geometry(folded, 3)
    assert faces, "fixture must reverse an interior face"
    with pytest.raises(UntangleError) as exc:
        untangle_folded_lines(folded, 3)
    err = exc.value
    assert any(r.block == 3 and r.node == (0, 3, 2) and "both endpoints" in r.reason
               for r in err.refusals)
    text = str(err)
    assert "REFUSED" in text
    assert "block 3 (0,3,2): folded j-edge has both endpoints" in text
    assert "re-check: not run (refused before repair)" in text
    assert len(err.report.reversed_faces) == len(faces)
    with pytest.raises(UntangleError):
        untangle_mesh([clean, folded])


def test_an_interior_fold_with_no_unambiguous_anchor_is_refused():
    ys = [0.0, 0.3, 0.6, 0.89, 0.9, 0.95, 1.0]
    folded = _with_y(_block(7, ys, 7), (3, 3, 3), 0.92)
    assert find_inverted_geometry(folded)[0], "fixture must reverse an interior face"
    with pytest.raises(UntangleError) as exc:
        untangle_folded_lines(folded, 0)
    assert any(r.node == (3, 3, 3) and "unambiguous anchor" in r.reason
               for r in exc.value.refusals)


def test_a_long_report_prints_the_first_entries_and_counts_the_rest():
    n = UNTANGLE_REPORT_LIMIT + 5
    report = UntangleReport(
        blocks_scanned=1,
        reversed_faces=[ReversedFace(0, 0, (i, 0, 0), -1.0) for i in range(n)])
    text = report.format_report()
    assert text.count("  reversed i-face") == UNTANGLE_REPORT_LIMIT
    assert "... 5 more reversed face(s) not shown" in text
    assert f"{n} reversed interior face(s)" in text


def test_the_file_command_writes_the_untangled_mesh_and_report(tmp_path):
    from plot3d import read_plot3D, write_plot3D

    folded = _with_y(_block(5, YS, 5), (2, 3, 2), 1.02)
    src = tmp_path / "mesh.xyz"
    write_plot3D(str(src), [folded], binary=True)
    report = untangle_plot3d_file(str(src))
    assert len(report.moved) == 1
    out = read_plot3D(str(tmp_path / "mesh_untangled.xyz"))
    assert find_inverted_geometry(out[0]) == ([], [])
    assert "node (2,3,2) moved" in (tmp_path / "mesh_untangled_report.txt").read_text()
    with pytest.raises(FileExistsError):
        untangle_plot3d_file(str(src))
