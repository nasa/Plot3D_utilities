# glennht_gpu.py
"""Export Plot3D multi-block connectivity to the ``connectivity.json``
format consumed by the ``glennht-gpu`` solver.

This module produces the same ``connectivity.json`` schema that
glennht-gpu's own upstream mesh-preprocessing tooling emits, letting
``plot3d_utilities`` -- which already computes block connectivity, outer
faces and rotational periodicity -- produce that JSON directly, without a
separate preprocessing round-trip.

Permutation-index convention
-----------------------------
:func:`plot3d.connectivity_fast` (and :func:`plot3d.connectivity`) follow
the *diagonal* (lb/ub) convention documented in ``connectivity.py``:

- **in-plane** matches (permutation 0-3) export ``permutation_index = -1``
  because the traversal direction is already fully encoded in block2's
  lb/ub ordering;
- **cross-plane** matches (permutation 4-7) export the real ``0..7`` index
  because lb/ub alone cannot represent an axis swap.

The glennht-gpu connectivity reader prefers a ``permutation_matrix`` over
the bare ``permutation_index`` whenever both are present: it derives the
canonical ``0..7`` index directly from the matrix, and only falls back to
``permutation_index.clamp(0, 7)`` when no matrix is supplied. A bare ``-1``
with **no** matrix would silently clamp to ``0`` (identity) -- which is
*wrong* for in-plane permutations 1-3 (u and/or v reversed with no axis
swap). To avoid that silent corruption, :func:`write_connectivity_json`
always emits **both** ``permutation_index`` (as-is, including ``-1``) *and*
``permutation_matrix`` at the top level of every ``face_matches`` /
``periodic_faces`` entry, so the reader's matrix-preferring path -- the one
it always tries first -- sees the real orientation regardless of what
``permutation_index`` says.

This differs from connectivity JSON produced by glennht-gpu's own upstream
mesh-preprocessing tooling, which carries only a flat ``permutation_index``
and no ``permutation_matrix``. That form's indices happen to be small real
values (0, 1, ...) rather than -1 because its face-match representation
doesn't have a "diagonal implies direction" shortcut. Both forms are
accepted by the glennht-gpu reader; we keep the richer one (with the
matrix) since it is unambiguous.
"""
from __future__ import annotations

import json
import math
import os
from copy import deepcopy
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from ..block import Block
from ..connectivity import PERMUTATION_MATRICES, connectivity
from ..facefunctions import create_face_from_diagonals
from ..periodicity import create_rotation_matrix, rotated_periodicity
from ..write import write_plot3D
from .gpu_boundary_conditions import write_boundary_conditions_yaml

__all__ = [
    "write_connectivity_json",
    "tag_surfaces_from_diagonals",
    "tag_surfaces_from_bc_codes",
    "tag_surfaces_geometric",
    "write_bc_codes_json",
    "merge_connectivity_json",
    "export_to_glennht_gpu",
]

#: Provenance string stamped into the ``periodicity`` block's ``source``
#: field when the caller doesn't supply one.
_SOURCE = "plot3d_utilities.glennht.gpu.write_connectivity_json"

#: Fixed solver face order used by ``bc_codes`` tables: one code per block
#: side, in this order. Matches the face-side convention used by
#: glennht-gpu's connectivity format.
FACE_ORDER: List[str] = ["I=1", "I=IMAX", "J=1", "J=JMAX", "K=1", "K=KMAX"]

#: Face-code -> surface-id legend (matches the ``bc_codes.json``
#: schema's ``face_code_legend`` used by glennht-gpu).
BC_CODE_LEGEND: Dict[str, str] = {
    "5": "inlet",
    "6": "outlet",
    "10": "blade",
    "9": "hub",
    "109": "shroud",
    "13": "periodic_lo",
    "14": "periodic_hi",
    "0": "interface",
}

#: Surface id used for an interface face (code 0) that plot3d_utilities
#: failed to block-match. These are CONFORMAL interior interfaces, NOT
#: walls -- matches glennht-gpu's convention for unmatched interfaces.
UNMATCHED_INTERFACE_ID = 8

#: Default surface-id -> name map used by :func:`tag_surfaces_geometric`
#: and :func:`tag_surfaces_from_bc_codes` (ids 1-5 mirror glennht-gpu's
#: default surface-id convention; 6 and 8 cover the remaining face-code
#: outcomes).
_DEFAULT_SURFACE_IDS: Dict[str, str] = {
    "1": "inlet",
    "2": "outlet",
    "3": "blade",
    "4": "hub",
    "5": "shroud",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_face_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Reduce a face-record dict to the schema's ``{block_index, lb, ub}``.

    Drops any extra keys (e.g. ``id``) that ``connectivity_fast`` /
    ``rotated_periodicity`` attach to ``block1``/``block2`` sub-dicts but
    which glennht-gpu's ``connectivity.json`` format does not carry on
    match entries.
    """
    return {
        "block_index": int(rec["block_index"]),
        "lb": [int(v) for v in rec["lb"]],
        "ub": [int(v) for v in rec["ub"]],
    }


def _clean_face_match(fm: Dict[str, Any]) -> Dict[str, Any]:
    """Reduce one ``connectivity_fast``/``rotated_periodicity`` match dict
    to the glennht-gpu ``face_matches`` / ``periodic_faces`` entry schema.

    Handles both known shapes:

    - matches with a nested ``orientation`` dict (the normal case emitted
      by :func:`plot3d.connectivity` / :func:`plot3d.connectivity_fast` /
      :func:`plot3d.rotated_periodicity`);
    - self-matches (O-grid branch-cut) which carry no ``orientation`` key
      at all -- these fall back to the identity permutation.
    """
    orientation = fm.get("orientation") or {}
    perm_idx = orientation.get("permutation_index", fm.get("permutation_index", -1))
    perm_mat = orientation.get("permutation_matrix", fm.get("permutation_matrix"))
    if perm_mat is None:
        perm_mat = PERMUTATION_MATRICES[0].tolist()  # identity fallback

    return {
        "block1": _clean_face_record(fm["block1"]),
        "block2": _clean_face_record(fm["block2"]),
        "permutation_index": int(perm_idx),
        "permutation_matrix": [[int(c) for c in row] for row in perm_mat],
    }


def _clean_outer_face(face: Dict[str, Any]) -> Dict[str, Any]:
    """Reduce an outer-face dict to ``{block_index, lb, ub, id}``."""
    fid = face.get("id")
    return {
        "block_index": int(face["block_index"]),
        "lb": [int(v) for v in face["lb"]],
        "ub": [int(v) for v in face["ub"]],
        "id": int(fid) if fid is not None else None,
    }


def _collapsed_axis_value(lb: Sequence[int], ub: Sequence[int]) -> Optional[Tuple[int, int]]:
    """Return ``(axis, value)`` for the constant axis of a normalised
    (min/max) ``lb``/``ub`` diagonal, or ``None`` if none of the three axes
    is constant (should not happen for a structured outer/matched face)."""
    if lb[0] == ub[0]:
        return 0, lb[0]
    if lb[1] == ub[1]:
        return 1, lb[1]
    if lb[2] == ub[2]:
        return 2, lb[2]
    return None


def _block_side_of_face(face: Dict[str, Any]) -> Optional[int]:
    """Which of the 6 block sides an outer face lies on, as an index into
    :data:`FACE_ORDER` (``[I=1, I=IMAX, J=1, J=JMAX, K=1, K=KMAX]``).

    Matches glennht-gpu's face-side convention: the low/high side is
    decided purely by whether the constant index is 0 (low) or not (high).
    Returns ``None`` if the face isn't collapsed on a single axis.
    """
    lb, ub = face["lb"], face["ub"]
    il, jl, kl = lb
    ih, jh, kh = ub
    if il == ih:
        return 0 if il == 0 else 1
    if jl == jh:
        return 2 if jl == 0 else 3
    if kl == kh:
        return 4 if kl == 0 else 5
    return None


def _face_code_to_surface_id(code: int) -> Optional[int]:
    """Map an integer face code to a surface id.

    Implements glennht-gpu's face-code surface-code convention:

    - 5 -> 1 (inlet), 6 -> 2 (outlet), 10 -> 3 (blade), 9 -> 4 (hub),
      109 -> 5 (shroud);
    - 13 or 14 (leftover periodic) -> 6;
    - 0 (interface) -> :data:`UNMATCHED_INTERFACE_ID` (8) -- a
      CONFORMAL interior interface that plot3d_utilities failed to
      block-match, NOT a wall;
    - anything else -> ``None`` (fall back to geometry).
    """
    mapping = {
        5: 1,
        6: 2,
        10: 3,
        9: 4,
        109: 5,
        13: 6,
        14: 6,
        0: UNMATCHED_INTERFACE_ID,
    }
    return mapping.get(int(code))


def _axial_and_radius(block: Block, i: int, j: int, k: int, axis: str) -> Tuple[float, float]:
    """Axial coordinate and radius of node (i,j,k) for a given rotation axis."""
    x = float(block.X[i, j, k])
    y = float(block.Y[i, j, k])
    z = float(block.Z[i, j, k])
    if axis == "x":
        return x, math.sqrt(y * y + z * z)
    elif axis == "y":
        return y, math.sqrt(x * x + z * z)
    else:
        return z, math.sqrt(x * x + y * y)


def _face_axial_radius(block: Block, face: Dict[str, Any], axis: str) -> Tuple[float, float]:
    """Mean axial coordinate and mean radius over a face's node span."""
    lb, ub = face["lb"], face["ub"]
    i0, i1 = sorted((lb[0], ub[0]))
    j0, j1 = sorted((lb[1], ub[1]))
    k0, k1 = sorted((lb[2], ub[2]))
    sa = sr = n = 0.0
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            for k in range(k0, k1 + 1):
                a, r = _axial_and_radius(block, i, j, k, axis)
                sa += a
                sr += r
                n += 1.0
    if n == 0.0:
        return 0.0, 0.0
    return sa / n, sr / n


# ---------------------------------------------------------------------------
# B1: connectivity.json writer
# ---------------------------------------------------------------------------

def write_connectivity_json(
    blocks: List[Block],
    matched_faces: List[Dict[str, Any]],
    outer_faces: List[Dict[str, Any]],
    filename: str,
    *,
    periodic_faces: Optional[List[Dict[str, Any]]] = None,
    periodicity: Optional[Dict[str, Any]] = None,
    surface_ids: Optional[Dict[Any, str]] = None,
    mesh_file: Optional[str] = None,
) -> Dict[str, Any]:
    """Write a glennht-gpu ``connectivity.json`` file.

    Args:
        blocks (List[Block]): The multi-block mesh (only ``len(blocks)`` is
            used, for ``nblocks``).
        matched_faces (List[Dict]): Face matches as returned by
            :func:`plot3d.connectivity_fast` / :func:`plot3d.connectivity`.
        outer_faces (List[Dict]): Outer (unmatched) faces, as returned by
            :func:`plot3d.connectivity_fast` or by
            :func:`tag_surfaces_geometric` / :func:`tag_surfaces_from_bc_codes`
            / :func:`tag_surfaces_from_diagonals` (id-tagged).
        filename (str): Output path, e.g. ``"stator.connectivity.json"``.
        periodic_faces (List[Dict], optional): Rotationally-periodic face
            pairs, as returned by :func:`plot3d.rotated_periodicity`.
        periodicity (Dict, optional): Global periodicity metadata:
            ``{nblades, rotation_axis, rotation_angle_rad,
            rotation_angle_deg, transformation_matrix, convention, source}``.
            Required (by glennht-gpu) whenever ``periodic_faces`` is
            non-empty. If given without a ``source`` key, one is filled in.
        surface_ids (Dict, optional): ``{id: name}`` map, e.g.
            ``{1: "inlet", 2: "outlet", ...}``. Written with string keys.
        mesh_file (str, optional): Value for the ``mesh_file`` field
            (typically the basename of the ``.xyz`` grid this connectivity
            was computed from). Defaults to the basename of ``filename``
            with its ``.connectivity.json`` suffix swapped for ``.xyz``.

    Returns:
        (Dict): The payload that was written (useful for tests).
    """
    if mesh_file is None:
        base = os.path.basename(filename)
        for suffix in (".connectivity.json", ".json"):
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break
        mesh_file = base + ".xyz"

    if periodic_faces and periodicity is None:
        raise ValueError(
            "periodic_faces was given but periodicity is missing; the "
            "glennht-gpu parser requires a 'periodicity' block whenever "
            "'periodic_faces' is non-empty."
        )

    payload: Dict[str, Any] = {
        "mesh_file": mesh_file,
        "nblocks": len(blocks),
        "face_matches": [_clean_face_match(fm) for fm in matched_faces],
        "outer_faces": [_clean_outer_face(f) for f in outer_faces],
        "periodic_faces": [_clean_face_match(fm) for fm in (periodic_faces or [])],
    }

    if periodicity is not None:
        periodicity = dict(periodicity)
        periodicity.setdefault("source", _SOURCE)
        payload["periodicity"] = periodicity

    payload["surface_ids"] = {str(k): v for k, v in (surface_ids or _DEFAULT_SURFACE_IDS).items()}
    payload["permutation_matrices"] = PERMUTATION_MATRICES.tolist()

    with open(filename, "w") as fh:
        json.dump(payload, fh, indent=2)

    return payload


# ---------------------------------------------------------------------------
# B2: surface tagging
# ---------------------------------------------------------------------------

def tag_surfaces_from_diagonals(
    blocks: List[Block],
    outer_faces: List[Dict[str, Any]],
    specs: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    """Tag outer faces by explicit index-diagonal specification.

    Args:
        blocks (List[Block]): The multi-block mesh.
        outer_faces (List[Dict]): Outer faces to tag (mutated in place).
        specs (List[Dict]): Each spec is
            ``{"lb": [i,j,k], "ub": [i,j,k], "block_index": int,
            "id": int, "name": str}``. A :class:`~plot3d.face.Face` is
            built via :func:`create_face_from_diagonals` (mainly to
            validate the diagonal against the block) and matched to the
            outer face on the same block whose collapsed (constant) axis
            and index value coincide, and whose index range contains the
            spec's range.

    Returns:
        (Tuple): containing

        - **outer_faces** (List[Dict]): The same list, with ``id`` set on
          every outer face that matched a spec.
        - **surface_ids** (Dict[str, str]): ``{"<id>": name}`` map built
          from the specs.
    """
    surface_ids: Dict[str, str] = {}

    for spec in specs:
        block_index = int(spec["block_index"])
        lb = list(spec["lb"])
        ub = list(spec["ub"])
        surf_id = int(spec["id"])
        name = str(spec.get("name", surf_id))

        # Validate the diagonal actually describes a face on this block
        # (raises/degenerates naturally if it doesn't collapse to a plane).
        create_face_from_diagonals(blocks[block_index], lb, ub)

        n_lb = [min(lb[d], ub[d]) for d in range(3)]
        n_ub = [max(lb[d], ub[d]) for d in range(3)]
        spec_axis_val = _collapsed_axis_value(n_lb, n_ub)

        for face in outer_faces:
            if int(face.get("block_index", -1)) != block_index:
                continue
            f_lb = [min(face["lb"][d], face["ub"][d]) for d in range(3)]
            f_ub = [max(face["lb"][d], face["ub"][d]) for d in range(3)]
            face_axis_val = _collapsed_axis_value(f_lb, f_ub)

            if face_axis_val is None or spec_axis_val is None:
                continue
            if face_axis_val != spec_axis_val:
                continue
            # spec's index range must lie within the outer face's range
            if all(f_lb[d] <= n_lb[d] and n_ub[d] <= f_ub[d] for d in range(3)):
                face["id"] = surf_id
                break

        surface_ids[str(surf_id)] = name

    return outer_faces, surface_ids


def tag_surfaces_from_bc_codes(
    blocks: List[Block],
    outer_faces: List[Dict[str, Any]],
    bc_codes: List[List[int]],
) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    """Tag outer faces using per-block solver face codes (authoritative).

    Args:
        blocks (List[Block]): The multi-block mesh.
        outer_faces (List[Dict]): Outer faces to tag (mutated in place).
        bc_codes (List[List[int]]): One 6-int row per block, in
            :data:`FACE_ORDER` order (``[I=1, I=IMAX, J=1, J=JMAX, K=1,
            K=KMAX]``) -- the same shape as ``stator.bc_codes.json``'s
            ``"blocks"`` array.

    Returns:
        (Tuple): containing

        - **outer_faces** (List[Dict]): The same list, with ``id`` set on
          every outer face whose (block, side) has a recognised code.
        - **surface_ids** (Dict[str, str]): The default id -> name map
          (:data:`_DEFAULT_SURFACE_IDS` plus ``6``/``8``).
    """
    for face in outer_faces:
        block_index = int(face.get("block_index", -1))
        if block_index < 0 or block_index >= len(bc_codes):
            continue
        side = _block_side_of_face(face)
        if side is None:
            continue
        code = bc_codes[block_index][side]
        surf_id = _face_code_to_surface_id(code)
        if surf_id is not None:
            face["id"] = surf_id

    surface_ids = dict(_DEFAULT_SURFACE_IDS)
    surface_ids["6"] = "periodic"
    surface_ids[str(UNMATCHED_INTERFACE_ID)] = "unmatched_interface"
    return outer_faces, surface_ids


def tag_surfaces_geometric(
    blocks: List[Block],
    outer_faces: List[Dict[str, Any]],
    axis: str = "x",
    band: float = 0.1,
    overrides: Optional[Dict[Tuple[int, int], int]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    """Tag outer faces by geometric position (fallback when no face codes
    are available).

    Matches glennht-gpu's geometric fallback convention: the global axial
    and radial extent of all outer-face centroids is computed first; a
    face within ``band`` (a fraction of that extent) of the -axial extreme
    is the inlet, of the +axial extreme the outlet, of min-radius the hub,
    of max-radius the shroud; everything else is tagged blade.

    Args:
        blocks (List[Block]): The multi-block mesh.
        outer_faces (List[Dict]): Outer faces to tag (mutated in place).
        axis (str): Rotation/machine axis, ``"x"``, ``"y"``, or ``"z"``.
        band (float): Fraction of the axial/radial range within which a
            face is classed as touching that extreme.
        overrides (Dict[(block_index, side), id], optional): Forces
            specific ``(block_index, side)`` pairs (``side`` from
            :func:`_block_side_of_face` / :data:`FACE_ORDER`) to a given
            surface id, bypassing the geometric heuristic.

    Returns:
        (Tuple): containing

        - **outer_faces** (List[Dict]): The same list, with ``id`` set on
          every face.
        - **surface_ids** (Dict[str, str]): ``{"1": "inlet", ..., "5":
          "shroud"}``.
    """
    overrides = overrides or {}

    ax_min = math.inf
    ax_max = -math.inf
    r_min = math.inf
    r_max = -math.inf
    centroids: List[Tuple[float, float]] = []
    for face in outer_faces:
        a, r = _face_axial_radius(blocks[int(face["block_index"])], face, axis)
        centroids.append((a, r))
        ax_min = min(ax_min, a)
        ax_max = max(ax_max, a)
        r_min = min(r_min, r)
        r_max = max(r_max, r)

    ax_range = max(ax_max - ax_min, 1e-30)
    r_range = max(r_max - r_min, 1e-30)

    for face, (a, r) in zip(outer_faces, centroids):
        side = _block_side_of_face(face)
        override_key = (int(face["block_index"]), side)
        if override_key in overrides:
            face["id"] = int(overrides[override_key])
            continue

        dx_lo = (a - ax_min) / ax_range  # 0 at inlet (-axial)
        dx_hi = (ax_max - a) / ax_range  # 0 at outlet (+axial)
        dr_lo = (r - r_min) / r_range    # 0 at hub
        dr_hi = (r_max - r) / r_range    # 0 at shroud

        if dx_lo <= band:
            face["id"] = 1
        elif dx_hi <= band:
            face["id"] = 2
        elif dr_lo <= band:
            face["id"] = 4
        elif dr_hi <= band:
            face["id"] = 5
        else:
            face["id"] = 3

    return outer_faces, dict(_DEFAULT_SURFACE_IDS)


# ---------------------------------------------------------------------------
# B3: bc_codes.json sidecar writer
# ---------------------------------------------------------------------------

def write_bc_codes_json(
    blocks: List[Block],
    bc_codes: List[List[int]],
    filename: str,
    *,
    block_order: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Write a ``<mesh>.bc_codes.json`` sidecar (solver per-face codes).

    Matches glennht-gpu's ``bc_codes.json`` schema: ``face_order``,
    ``face_code_legend``, ``blocks`` (the per-block 6-int arrays),
    ``block_order``.

    Args:
        blocks (List[Block]): The multi-block mesh (used only for
            ``len(blocks)`` when ``block_order`` isn't given).
        bc_codes (List[List[int]]): One 6-int row per block, in
            :data:`FACE_ORDER` order.
        filename (str): Output path, e.g. ``"stator.bc_codes.json"``.
        block_order (List[str], optional): Human-readable name per block
            (e.g. ``["S35#0001", ...]``). Defaults to ``["block_0", ...]``.

    Returns:
        (Dict): The payload that was written (useful for tests).
    """
    if block_order is None:
        block_order = [f"block_{i}" for i in range(len(blocks))]

    payload = {
        "face_order": list(FACE_ORDER),
        "face_code_legend": dict(BC_CODE_LEGEND),
        "blocks": [[int(c) for c in row] for row in bc_codes],
        "block_order": list(block_order),
    }

    with open(filename, "w") as fh:
        json.dump(payload, fh, indent=2)

    return payload


# ---------------------------------------------------------------------------
# B5: multi-row merge
# ---------------------------------------------------------------------------

def _offset_match_block_indices(fm: Dict[str, Any], off: int) -> None:
    """Shift ``block1``/``block2``'s ``block_index`` by ``off`` (in place).

    Mirrors how glennht-gpu's own connectivity-merge step shifts block
    indices when concatenating rows. Note it only touches ``block_index``
    -- NOT any per-side ``id`` that some producers attach to ``block1``/
    ``block2`` in ``periodic_faces`` entries: those per-side ids are left
    unchanged by row offsetting even though ``block_index`` is correctly
    shifted to match the merged mesh's numbering. Only
    ``outer_faces[].id`` is ever shifted by ``id_stride``.
    """
    for side in ("block1", "block2"):
        b = fm.get(side)
        if b is not None and "block_index" in b:
            b["block_index"] = int(b["block_index"]) + off


def merge_connectivity_json(
    row_files: List[Union[str, "os.PathLike[str]", Dict[str, Any]]],
    out_filename: str,
    *,
    id_stride: int = 100,
    mesh_file: Optional[str] = None,
) -> Dict[str, Any]:
    """Merge several per-row ``connectivity.json`` payloads into one.

    Companion to a concatenated multi-row ``.xyz`` (e.g. rotor blocks then
    stator blocks): mirrors how glennht-gpu merges several single-row
    ``connectivity.json`` payloads (e.g. a rotor row and a stator row) into
    one multi-row mesh's connectivity file.

    For row ``r`` (0-indexed, in flow order matching the concatenated grid):

    - every ``block_index`` in ``face_matches``, ``outer_faces`` and
      ``periodic_faces`` is shifted by the cumulative block count of the
      prior rows;
    - every ``outer_faces[].id`` is shifted by ``id_stride * r`` (the
      ``id_stride=100`` convention: row 0 keeps ids 1-5, row 1 becomes
      101-105, etc. -- see the ``boundary_conditions:`` list of your
      glennht-gpu run configuration);
    - each row's own ``periodicity.rotation_angle_rad`` (its blade pitch)
      is stamped onto that row's ``periodic_faces`` entries, so a merged
      mesh with different blade counts per row still tags each seam's
      pitch correctly;
    - the merged top-level ``periodicity`` block is row 0's (a
      default/fallback -- the authoritative per-seam angle rides on each
      ``periodic_faces`` entry's own ``rotation_angle_rad``), matching
      glennht-gpu's own merge behaviour.

    NOTE on ``surface_ids``: glennht-gpu's own connectivity merge does not
    carry this field at all -- a merged file it produces has no
    ``surface_ids`` key, even though both per-row inputs do. This function
    instead merges and offsets ``surface_ids`` (by ``id_stride * r`` on the
    integer keys, matching ``outer_faces[].id``) because it is useful and
    harmless (nothing downstream reads it; the run configuration hardcodes
    its own surface-id map), but a byte-for-byte match against glennht-gpu's
    own merged output is therefore neither possible nor expected.

    Args:
        row_files (List[str | Dict]): Per-row connectivity payloads, in
            flow order. Each entry is either a path to a
            ``*.connectivity.json`` file, or an already-loaded dict with
            the same schema.
        out_filename (str): Output path for the merged JSON.
        id_stride (int): Per-row surface-id offset. Defaults to 100.
        mesh_file (str, optional): Value for the merged payload's
            ``mesh_file`` field. Defaults to the basename of
            ``out_filename`` with its ``.connectivity.json``/``.json``
            suffix swapped for ``.xyz``.

    Returns:
        (Dict): The merged payload that was written (useful for tests).
    """
    all_fm: List[Dict[str, Any]] = []
    all_outer: List[Dict[str, Any]] = []
    all_periodic: List[Dict[str, Any]] = []
    merged_surface_ids: Dict[str, str] = {}
    global_periodicity: Optional[Dict[str, Any]] = None
    block_off = 0

    for row, rf in enumerate(row_files):
        if isinstance(rf, dict):
            d = rf
        else:
            with open(rf) as fh:
                d = json.load(fh)

        nblocks = int(d.get("nblocks", 0))
        id_off = id_stride * row
        row_periodicity = d.get("periodicity")
        row_angle = (row_periodicity or {}).get("rotation_angle_rad")
        if global_periodicity is None and row_periodicity is not None:
            global_periodicity = deepcopy(row_periodicity)

        for fm in d.get("face_matches", []):
            fm = deepcopy(fm)
            _offset_match_block_indices(fm, block_off)
            all_fm.append(fm)

        for of in d.get("outer_faces", []):
            of = deepcopy(of)
            of["block_index"] = int(of["block_index"]) + block_off
            if of.get("id") is not None:
                of["id"] = int(of["id"]) + id_off
            all_outer.append(of)

        for pf in d.get("periodic_faces", []):
            pf = deepcopy(pf)
            _offset_match_block_indices(pf, block_off)
            if row_angle is not None:
                pf["rotation_angle_rad"] = row_angle
            all_periodic.append(pf)

        for k, v in (d.get("surface_ids") or {}).items():
            merged_surface_ids[str(int(k) + id_off)] = v

        block_off += nblocks

    if mesh_file is None:
        base = os.path.basename(out_filename)
        for suffix in (".connectivity.json", ".json"):
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break
        mesh_file = base + ".xyz"

    payload: Dict[str, Any] = {
        "mesh_file": mesh_file,
        "nblocks": block_off,
        "face_matches": all_fm,
        "outer_faces": all_outer,
        "periodic_faces": all_periodic,
    }
    if global_periodicity is not None:
        payload["periodicity"] = global_periodicity
    if merged_surface_ids:
        payload["surface_ids"] = merged_surface_ids
    payload["permutation_matrices"] = PERMUTATION_MATRICES.tolist()

    with open(out_filename, "w") as fh:
        json.dump(payload, fh, indent=2)

    return payload


# ---------------------------------------------------------------------------
# B6: top-level driver
# ---------------------------------------------------------------------------

def export_to_glennht_gpu(
    blocks: List[Block],
    out_dir: str,
    case_name: str,
    *,
    rotation_angle: Optional[float] = None,
    rotation_axis: str = "x",
    nblades: Optional[int] = None,
    tagging: str = "geometric",
    bc_codes: Optional[List[List[int]]] = None,
    surface_specs: Optional[List[Dict[str, Any]]] = None,
    bcs: Optional[List[Any]] = None,
    rotation: Optional[Dict[str, Any]] = None,
    write_grid: bool = True,
    tagging_kwargs: Optional[Dict[str, Any]] = None,
    bundle: bool = False,
) -> Dict[str, str]:
    """Single-call happy path: Plot3D blocks -> glennht-gpu case files.

    Runs the FULL :func:`plot3d.connectivity.connectivity` (never
    ``connectivity_fast`` -- see the module docstring / verification
    notes: on ``stator.xyz`` the fast GCD-reduced path drops 5 of 10 real
    face matches, mis-tagging them as unmatched interfaces), optionally
    layers rotational periodicity on top via
    :func:`plot3d.periodicity.rotated_periodicity`, tags the remaining
    outer faces with the requested method, then writes whichever of the
    grid / connectivity / bc_codes / boundary-conditions files the caller
    asked for.

    Args:
        blocks (List[Block]): The multi-block mesh (one row/case).
        out_dir (str): Output directory (created if missing).
        case_name (str): Base name for every output file, e.g.
            ``"stator"`` -> ``stator.xyz`` / ``stator.connectivity.json``
            / ``stator.bc_codes.json`` / ``stator_boundary_conditions.yaml``.
        rotation_angle (float, optional): Blade pitch angle, in
            **radians**, e.g. ``math.radians(7.826)`` for a 46-blade row.
            Must be given together with ``nblades`` to compute rotational
            periodicity; if omitted, no ``periodic_faces``/``periodicity``
            block is produced (all periodic faces stay tagged as ordinary
            outer faces).
        rotation_axis (str): ``"x"``, ``"y"``, or ``"z"``. Defaults to
            ``"x"``.
        nblades (int, optional): Blade count for a full annulus. Must be
            given together with ``rotation_angle``.
        tagging (str): One of ``"geometric"``, ``"bc_codes"``,
            ``"diagonals"`` -- which :func:`tag_surfaces_*` helper to run
            on the outer faces before writing ``connectivity.json``.
        bc_codes (List[List[int]], optional): Required when
            ``tagging="bc_codes"`` (passed to
            :func:`tag_surfaces_from_bc_codes`); also written to
            ``{case_name}.bc_codes.json`` whenever given, regardless of
            ``tagging``.
        surface_specs (List[Dict], optional): Required when
            ``tagging="diagonals"`` (passed to
            :func:`tag_surfaces_from_diagonals`).
        bcs (List[GpuInletBC | GpuOutletBC | GpuWallBC], optional): If
            given, written to ``{case_name}_boundary_conditions.yaml`` via
            :func:`write_boundary_conditions_yaml`.
        rotation (Dict, optional): Forwarded to
            :func:`write_boundary_conditions_yaml`'s ``rotation=`` (the
            run yaml's ``rotation:`` block); only used when ``bcs`` is
            given.
        write_grid (bool): If True (default), also write
            ``{case_name}.xyz`` (formatted/ASCII) via
            :func:`plot3d.write.write_plot3D`.
        tagging_kwargs (Dict, optional): Extra keyword arguments forwarded
            to the chosen ``tag_surfaces_*`` function (e.g. ``{"axis":
            "x", "band": 0.1}`` for ``"geometric"``, or ``{"overrides":
            {...}}``).
        bundle (bool): If True, also write a single-file HDF5
            ``{case_name}.graph_p3d`` bundle (grid + connectivity +
            surface tags + ``bcs``) via
            :func:`plot3d.glennht.graph_p3d.write_graph_p3d`. Requires the
            optional ``h5py`` dependency (the ``hdf5`` extra). glennht-gpu
            cannot read this bundle directly -- explode it back into the
            loose trio with :func:`plot3d.glennht.graph_p3d.graph_p3d_to_files`
            first. Defaults to False, so existing callers are unaffected.

    Returns:
        (Dict[str, str]): ``{"connectivity": ..., "grid": ...,
        "bc_codes": ..., "boundary_conditions": ..., "graph_p3d": ...}``
        -- only the keys for files actually written are present.
    """
    if (rotation_angle is None) != (nblades is None):
        raise ValueError(
            "rotation_angle and nblades must be given together (both or "
            "neither) -- got rotation_angle=%r, nblades=%r" % (rotation_angle, nblades)
        )

    os.makedirs(out_dir, exist_ok=True)
    tagging_kwargs = dict(tagging_kwargs or {})

    # FULL connectivity -- see module/function docstring: connectivity_fast's
    # GCD reduction silently drops non-aligned O-grid/H-split face matches.
    face_matches, outer_faces = connectivity(blocks)

    periodic_faces: List[Dict[str, Any]] = []
    periodicity_meta: Optional[Dict[str, Any]] = None
    if rotation_angle is not None:
        periodic_faces, outer_faces, _, _ = rotated_periodicity(
            blocks,
            face_matches,
            outer_faces,
            rotation_angle=math.degrees(rotation_angle),
            rotation_axis=rotation_axis,
        )
        transformation_matrix = create_rotation_matrix(rotation_angle, rotation_axis)
        periodicity_meta = {
            "nblades": int(nblades),
            "rotation_axis": rotation_axis,
            "rotation_angle_rad": float(rotation_angle),
            "rotation_angle_deg": float(math.degrees(rotation_angle)),
            "transformation_matrix": transformation_matrix.tolist(),
            "convention": "Face_B_points = (transformation_matrix @ Face_A_points.T).T",
            "source": _SOURCE,
        }

    if tagging == "bc_codes":
        if bc_codes is None:
            raise ValueError("tagging='bc_codes' requires bc_codes=")
        outer_faces, surface_ids = tag_surfaces_from_bc_codes(
            blocks, outer_faces, bc_codes, **tagging_kwargs
        )
    elif tagging == "diagonals":
        if surface_specs is None:
            raise ValueError("tagging='diagonals' requires surface_specs=")
        outer_faces, surface_ids = tag_surfaces_from_diagonals(
            blocks, outer_faces, surface_specs, **tagging_kwargs
        )
    elif tagging == "geometric":
        tagging_kwargs.setdefault("axis", rotation_axis)
        outer_faces, surface_ids = tag_surfaces_geometric(
            blocks, outer_faces, **tagging_kwargs
        )
    else:
        raise ValueError(
            f"unknown tagging={tagging!r}; expected 'geometric', 'bc_codes', or 'diagonals'"
        )

    paths: Dict[str, str] = {}

    if write_grid:
        grid_path = os.path.join(out_dir, f"{case_name}.xyz")
        write_plot3D(grid_path, blocks, binary=False)
        paths["grid"] = grid_path

    conn_path = os.path.join(out_dir, f"{case_name}.connectivity.json")
    write_connectivity_json(
        blocks,
        face_matches,
        outer_faces,
        conn_path,
        periodic_faces=periodic_faces,
        periodicity=periodicity_meta,
        surface_ids=surface_ids,
        mesh_file=f"{case_name}.xyz",
    )
    paths["connectivity"] = conn_path

    if bc_codes is not None:
        bc_codes_path = os.path.join(out_dir, f"{case_name}.bc_codes.json")
        write_bc_codes_json(blocks, bc_codes, bc_codes_path)
        paths["bc_codes"] = bc_codes_path

    if bcs is not None:
        bcs_path = os.path.join(out_dir, f"{case_name}_boundary_conditions.yaml")
        write_boundary_conditions_yaml(bcs, bcs_path, rotation=rotation)
        paths["boundary_conditions"] = bcs_path

    if bundle:
        # Local import: graph_p3d.py imports helpers FROM this module at
        # its own top level, so importing it back here at gpu.py's module
        # top would be circular. Deferring it into this branch also keeps
        # h5py (optional dependency, imported lazily inside graph_p3d.py)
        # out of the import chain entirely unless bundle=True is used.
        from .graph_p3d import write_graph_p3d

        graph_p3d_path = os.path.join(out_dir, f"{case_name}.graph_p3d")
        write_graph_p3d(
            graph_p3d_path,
            blocks,
            matched_faces=face_matches,
            outer_faces=outer_faces,
            periodic_faces=periodic_faces,
            periodicity=periodicity_meta,
            surface_ids=surface_ids,
            bcs=bcs,
            rotation=rotation,
            case_name=case_name,
            mesh_file=f"{case_name}.xyz",
        )
        paths["graph_p3d"] = graph_p3d_path

    return paths
