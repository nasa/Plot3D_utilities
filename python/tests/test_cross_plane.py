"""Test cross-plane connectivity and directed diagonal reconstruction.

Uses cross_plane_pair.p3d (2 blocks) where block1's k=0 face matches
block2's j=0 face.  Expected GHT_CONN directed diagonal (1-based):

    1 1 1 1 25 409 1
    2 409 1 1 1 1 25
"""

import os
import tempfile
import pytest

MESH_PATH = os.path.join(os.path.dirname(__file__), "data", "cross_plane_pair.p3d")

skip_no_mesh = pytest.mark.skipif(
    not os.path.exists(MESH_PATH),
    reason="cross_plane_pair.p3d not found"
)


def _mesh_path():
    return MESH_PATH


@skip_no_mesh
def test_python_connectivity_produces_directed_diagonal():
    """Python connectivity_fast already produces directed lb/ub."""
    from plot3d import read_plot3D, connectivity_fast

    blocks = read_plot3D(_mesh_path(), binary=False)
    assert len(blocks) == 2

    face_matches, _outer = connectivity_fast(blocks)
    assert len(face_matches) == 1

    fm = face_matches[0]
    b1 = fm['block1']
    b2 = fm['block2']

    # 0-based expected values
    assert b1['lb'] == [0, 0, 0]
    assert b1['ub'] == [24, 408, 0]
    assert b2['lb'] == [408, 0, 0]
    assert b2['ub'] == [0, 0, 24]


@skip_no_mesh
def test_ght_conn_export_directed_diagonal():
    """End-to-end: connectivity + export produces correct ght_conn."""
    from plot3d import read_plot3D, connectivity_fast
    from plot3d.glennht.export_functions import export_to_glennht_conn

    blocks = read_plot3D(_mesh_path(), binary=False)
    face_matches, outer_faces = connectivity_fast(blocks)

    tmpf = tempfile.NamedTemporaryFile(suffix='.ght_conn', delete=False, mode='w')
    tmpf.close()
    try:
        export_to_glennht_conn(face_matches, outer_faces, tmpf.name, [], [], [])
        with open(tmpf.name) as f:
            lines = [l.strip() for l in f.readlines()]

        # Line 0: number of matches
        assert lines[0] == "1"
        # Lines 1-2: the match pair (1-based, tab/space separated)
        vals1 = lines[1].split()
        vals2 = lines[2].split()
        assert vals1 == ["1", "1", "1", "1", "25", "409", "1"]
        assert vals2 == ["2", "409", "1", "1", "1", "1", "25"]
    finally:
        os.unlink(tmpf.name)


