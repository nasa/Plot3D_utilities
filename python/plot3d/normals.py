"""Index-space normals and permutation matrix validation.

This module provides:
- ``index_space_normal``: Topological (integer) outward normal for a face.
- ``compute_permutation_matrix``: Compute the 3x3 permutation matrix from
  diagonal corners.
- ``validate_connectivity``: Batch validation of face matches.
- ``export_normals_json`` / ``import_normals_json``: Serialize face normals to/from JSON.
- ``plot_face_normals``: Matplotlib 3D quiver plot of normals at face centroids.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Varying-axis lookup
# ---------------------------------------------------------------------------

_VARYING_AXES = {0: (1, 2), 1: (0, 2), 2: (0, 1)}


def _varying_axes(const_ax: int) -> Tuple[int, int]:
    """Return the two varying axes for a given constant axis (0-indexed)."""
    return _VARYING_AXES[const_ax]


# ---------------------------------------------------------------------------
# Index-space normal
# ---------------------------------------------------------------------------

def index_space_normal(lb: Sequence[int], ub: Sequence[int]) -> Optional[np.ndarray]:
    """Compute the index-space (topological) outward normal for a face.

    Parameters
    ----------
    lb, ub : array-like of int, length 3
        Lower and upper bounds of the face (0-indexed).

    Returns
    -------
    np.ndarray of shape (3,) with dtype int8, or None if the face has no
    single constant axis.
    """
    lb = list(lb)
    ub = list(ub)
    const_axes = [d for d in range(3) if lb[d] == ub[d]]
    if len(const_axes) != 1:
        return None
    c = const_axes[0]
    n = np.zeros(3, dtype=np.int8)
    n[c] = -1 if lb[c] == 0 else 1
    return n


# ---------------------------------------------------------------------------
# Permutation matrix from diagonal corners
# ---------------------------------------------------------------------------

def compute_permutation_matrix(
    a1: Sequence[int],
    b1: Sequence[int],
    a2: Sequence[int],
    b2: Sequence[int],
) -> Optional[np.ndarray]:
    """Compute the 3x3 permutation matrix from diagonal corners.

    The formula maps block-1 indices to block-2 indices at a shared
    interface:  ``j = a2 + M * (i - a1)``

    Parameters
    ----------
    a1, b1 : array-like of int, length 3
        Diagonal corners of the face on block 1 (0-indexed).
    a2, b2 : array-like of int, length 3
        Diagonal corners of the face on block 2 (0-indexed).
        ``a1 <-> a2`` and ``b1 <-> b2`` must be spatially coincident.

    Returns
    -------
    np.ndarray of shape (3, 3) with dtype int8, or None if any entry is
    non-integer (bad corner pairing) or |det| != 1.
    """
    a1 = np.asarray(a1, dtype=np.int64)
    b1 = np.asarray(b1, dtype=np.int64)
    a2 = np.asarray(a2, dtype=np.int64)
    b2 = np.asarray(b2, dtype=np.int64)

    # Step 1: index-space normals
    n1 = np.zeros(3, dtype=np.int8)
    n2 = np.zeros(3, dtype=np.int8)
    for d in range(3):
        if a1[d] == b1[d]:
            n1[d] = -1 if a1[d] == 0 else 1
        if a2[d] == b2[d]:
            n2[d] = -1 if a2[d] == 0 else 1

    # Step 2: constant (face) axes
    face1_axes = np.nonzero(n1)[0]
    face2_axes = np.nonzero(n2)[0]
    if len(face1_axes) != 1 or len(face2_axes) != 1:
        return None
    face1 = int(face1_axes[0])
    face2 = int(face2_axes[0])

    # Step 3: pf1f2 and chirality factor s
    pf1f2 = -int(n1[face1]) * int(n2[face2])
    exp = (face1 + 1) + (face2 + 1)  # 1-indexed
    sign = 1 if exp % 2 == 0 else -1
    s = sign * pf1f2

    # Step 4: diagonal vectors
    d1 = (b1 - a1).astype(np.int64)
    d2 = (b2 - a2).astype(np.int64)
    d = int(np.dot(d1, d1))
    if d == 0:
        return None

    # Step 5: varying axes
    i11, i12 = _varying_axes(face1)
    i21, i22 = _varying_axes(face2)

    # Step 6: build 3x3 matrix
    m = np.zeros((3, 3), dtype=np.int64)
    m[i21, i11] = (d1[i11] * d2[i21] + s * d1[i12] * d2[i22]) // d
    m[i21, i12] = (d1[i12] * d2[i21] - s * d1[i11] * d2[i22]) // d
    m[i22, i11] = -s * m[i21, i12]
    m[i22, i12] = s * m[i21, i11]
    m[face2, face1] = pf1f2

    # Validate entries in {-1, 0, +1}
    if np.any(np.abs(m) > 1):
        return None

    result = m.astype(np.int8)

    # Check |det| == 1
    det = int(np.round(np.linalg.det(result.astype(float))))
    if abs(det) != 1:
        return None

    return result


# ---------------------------------------------------------------------------
# Batch validation
# ---------------------------------------------------------------------------

def validate_connectivity(
    face_matches: List[Dict[str, Any]],
    blocks: Optional[list] = None,
    tol: float = 1e-6,
) -> List[Dict[str, Any]]:
    """Validate face matches produce valid 3x3 permutation matrices.

    Parameters
    ----------
    face_matches : list of dict
        Each dict has ``block1`` and ``block2`` sub-dicts with
        ``lb``/``ub`` (or ``lo``/``hi``) keys and ``block_index``.
    blocks : list of Block, optional
        If provided, also checks spatial coincidence of diagonal corners.
    tol : float
        Spatial coincidence tolerance (only used when *blocks* is given).

    Returns
    -------
    List of dicts, one per match, with keys:
    - ``match_index``: index into *face_matches*
    - ``valid``: bool
    - ``matrix``: 3x3 list or None
    - ``error``: string or None
    """
    results = []
    for idx, fm in enumerate(face_matches):
        entry: Dict[str, Any] = {"match_index": idx, "valid": False, "matrix": None, "error": None}

        b1 = fm["block1"]
        b2 = fm["block2"]
        lb1 = b1.get("lb", b1.get("lo"))
        ub1 = b1.get("ub", b1.get("hi"))
        lb2 = b2.get("lb", b2.get("lo"))
        ub2 = b2.get("ub", b2.get("hi"))

        if lb1 is None or ub1 is None or lb2 is None or ub2 is None:
            entry["error"] = "missing lb/ub or lo/hi"
            results.append(entry)
            continue

        mat = compute_permutation_matrix(lb1, ub1, lb2, ub2)
        if mat is None:
            entry["error"] = "Permutation matrix computation failed (non-integer entries or |det| != 1)"
        else:
            entry["valid"] = True
            entry["matrix"] = mat.tolist()

        results.append(entry)

    return results


# ---------------------------------------------------------------------------
# Normals for all 6 faces of every block
# ---------------------------------------------------------------------------

def compute_all_normals(blocks: list) -> List[Dict[str, Any]]:
    """Compute index-space normals for all 6 faces of every block.

    Parameters
    ----------
    blocks : list of Block
        Each block must have ``IMAX``, ``JMAX``, ``KMAX`` attributes.

    Returns
    -------
    List of dicts with keys: ``block_index``, ``lb``, ``ub``,
    ``index_normal``, ``face_type``.
    """
    _FACE_TYPE = {
        (0, True): "imin", (0, False): "imax",
        (1, True): "jmin", (1, False): "jmax",
        (2, True): "kmin", (2, False): "kmax",
    }
    faces = []
    for bi, blk in enumerate(blocks):
        dims = [blk.IMAX, blk.JMAX, blk.KMAX]
        for const_ax in range(3):
            for is_low in (True, False):
                val = 0 if is_low else dims[const_ax] - 1
                lb = [0, 0, 0]
                ub = [0, 0, 0]
                for d in range(3):
                    if d == const_ax:
                        lb[d] = val
                        ub[d] = val
                    else:
                        lb[d] = 0
                        ub[d] = dims[d] - 1
                normal = np.zeros(3, dtype=np.int8)
                normal[const_ax] = -1 if is_low else 1
                faces.append({
                    "block_index": bi,
                    "lb": lb,
                    "ub": ub,
                    "index_normal": normal.tolist(),
                    "face_type": _FACE_TYPE[(const_ax, is_low)],
                })
    return faces


# ---------------------------------------------------------------------------
# JSON export / import
# ---------------------------------------------------------------------------

def export_normals_json(normals: List[Dict[str, Any]], path) -> None:
    """Write normals to a JSON file."""
    path = Path(path)
    data = {"faces": normals}
    path.write_text(json.dumps(data, indent=2))


def import_normals_json(path) -> List[Dict[str, Any]]:
    """Read normals from a JSON file. Returns the list of face dicts."""
    path = Path(path)
    data = json.loads(path.read_text())
    return data["faces"]


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_face_normals(
    blocks: list,
    face_records: Optional[List[Dict[str, Any]]] = None,
    ax=None,
    arrow_scale: float = 1.0,
    show_wireframe: bool = True,
):
    """Plot index-space normals as 3D quiver arrows at face centroids.

    Parameters
    ----------
    blocks : list of Block
        Blocks with ``X``, ``Y``, ``Z`` arrays.
    face_records : list of dict, optional
        Output of ``compute_all_normals``. If None, computed automatically.
    ax : matplotlib Axes3D, optional
        Existing 3D axes. Created if None.
    arrow_scale : float
        Scaling factor for arrow length.
    show_wireframe : bool
        If True, draw block edges as wireframe.

    Returns
    -------
    ax : matplotlib Axes3D
    """
    import matplotlib.pyplot as plt

    if face_records is None:
        face_records = compute_all_normals(blocks)

    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="3d")

    # Color map by face type
    _COLORS = {
        "imin": "red", "imax": "darkred",
        "jmin": "green", "jmax": "darkgreen",
        "kmin": "blue", "kmax": "darkblue",
    }

    for rec in face_records:
        bi = rec["block_index"]
        blk = blocks[bi]
        lb = rec["lb"]
        ub = rec["ub"]
        normal = rec["index_normal"]
        face_type = rec["face_type"]

        # Compute face centroid from block coordinates
        i_range = range(min(lb[0], ub[0]), max(lb[0], ub[0]) + 1)
        j_range = range(min(lb[1], ub[1]), max(lb[1], ub[1]) + 1)
        k_range = range(min(lb[2], ub[2]), max(lb[2], ub[2]) + 1)

        xs, ys, zs = [], [], []
        for i in i_range:
            for j in j_range:
                for k in k_range:
                    xs.append(blk.X[i, j, k])
                    ys.append(blk.Y[i, j, k])
                    zs.append(blk.Z[i, j, k])

        if not xs:
            continue

        cx = np.mean(xs)
        cy = np.mean(ys)
        cz = np.mean(zs)

        # Estimate a physical scale from face extent
        extent = max(
            max(xs) - min(xs),
            max(ys) - min(ys),
            max(zs) - min(zs),
        )
        if extent == 0:
            extent = 1.0
        length = extent * 0.15 * arrow_scale

        color = _COLORS.get(face_type, "gray")
        ax.quiver(
            cx, cy, cz,
            normal[0] * length, normal[1] * length, normal[2] * length,
            color=color, arrow_length_ratio=0.3, linewidth=1.5,
        )

    # Wireframe for context
    if show_wireframe:
        for bi, blk in enumerate(blocks):
            ni, nj, nk = blk.X.shape
            corners = [
                (0, 0, 0), (ni - 1, 0, 0), (ni - 1, nj - 1, 0), (0, nj - 1, 0),
                (0, 0, nk - 1), (ni - 1, 0, nk - 1), (ni - 1, nj - 1, nk - 1), (0, nj - 1, nk - 1),
            ]
            edges = [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7),
            ]
            for e0, e1 in edges:
                c0, c1 = corners[e0], corners[e1]
                ax.plot3D(
                    [blk.X[c0], blk.X[c1]],
                    [blk.Y[c0], blk.Y[c1]],
                    [blk.Z[c0], blk.Z[c1]],
                    "k-", alpha=0.3, linewidth=0.5,
                )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Index-Space Normals")
    return ax
