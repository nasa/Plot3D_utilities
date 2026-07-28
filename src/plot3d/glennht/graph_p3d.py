# graph_p3d.py
"""Single-file HDF5 bundle format (``.graph_p3d``) for a Plot3D multi-block
mesh + its glennht-gpu block-connectivity graph + boundary-condition
surface tags.

``connectivity.json`` (see :mod:`plot3d.glennht.gpu`) always ships as a
*trio* of loose files: the ``.xyz`` grid, the ``.connectivity.json``, and
(optionally) a ``_boundary_conditions.yaml`` fragment. ``.graph_p3d``
bundles all three into a single portable HDF5 file -- convenient for
storing/moving one case as one artifact.

glennht-gpu CANNOT read a ``.graph_p3d`` bundle directly -- it only
understands the loose trio. Use :func:`graph_p3d_to_files` to *explode* a
bundle back into ``{case}.xyz`` / ``{case}.connectivity.json`` /
``{case}_boundary_conditions.yaml`` right before handing files to
glennht-gpu.

``h5py`` is an OPTIONAL dependency (the ``hdf5`` extra in
``pyproject.toml``). It is imported lazily, inside :func:`_require_h5py`
only -- so ``import plot3d`` / ``import plot3d.glennht`` keep working with
h5py absent; only actually *calling* one of this module's functions
without h5py installed raises a clear ``ImportError``.

HDF5 layout
-----------
::

    <case>.graph_p3d  (HDF5)
    - root attrs: format="graph_p3d", format_version="1.0",
                  created_by="plot3d <version>",
                  source="plot3d_utilities.glennht.graph_p3d",
                  nblocks, case_name, mesh_file
    - /grid/block_<n>/            one subgroup per structured block
        attrs: imax, jmax, kmax
        x, y, z    (imax,jmax,kmax) float64, gzip-compressed             # (i,j,k)
    - /connectivity/
        face_matches                 (N,15) int64
        face_match_permutations      (N,2,2) int64
        outer_faces                  (M,8)  int64
        periodic_faces                (P,17) int64            [only if P>0]
        periodic_face_permutations    (P,2,2) int64           [only if P>0]
        periodic_rotation_angles      (P,) float64            [only if any
                                                                pair carries one]
        permutation_matrices         (8,2,2) int64
        surface_ids/  (group)        attrs {"<id>": name}
    - /periodicity/  (created ONLY if periodicity metadata was given)
        attrs: nblades, rotation_axis, rotation_angle_rad, rotation_angle_deg,
               convention, source
        transformation_matrix (3,3) float64
    - /boundary_conditions   (created ONLY if bcs given): scalar string
        dataset holding the BC YAML fragment text.

Column layout (see :func:`write_graph_p3d` for the full encoding notes):

- ``face_matches`` row: ``[b1, b1_lb(3), b1_ub(3), b2, b2_lb(3), b2_ub(3),
  perm_index]``.
- ``outer_faces`` row: ``[block, lb(3), ub(3), id]`` (``id=-1`` when
  null/untagged).
- ``periodic_faces`` row: ``[b1, b1_id, b1_lb(3), b1_ub(3), b2, b2_id,
  b2_lb(3), b2_ub(3), perm_index]`` (``*_id=-1`` when unset). Unlike
  ``face_matches``, each side's ``id`` (as attached by
  :func:`plot3d.rotated_periodicity`) is preserved here even though
  ``connectivity.json`` itself drops it on match entries.

Directed ``lb``/``ub`` diagonals are stored exactly as given (never
sorted), so orientation round-trips through :func:`read_graph_p3d`.
"""
from __future__ import annotations

import math
import os
from typing import Any, Dict, List, Optional

import numpy as np

from ..block import Block
from ..connectivity import PERMUTATION_MATRICES
from ..write import write_plot3D
from .gpu import (
    _DEFAULT_SURFACE_IDS,
    _clean_face_match,
    _clean_outer_face,
    write_connectivity_json,
)
from .gpu_boundary_conditions import _boundary_conditions_yaml_text

__all__ = [
    "write_graph_p3d",
    "read_graph_p3d",
    "graph_p3d_to_files",
]

#: Format tag / version stamped into the bundle's root attrs.
_FORMAT = "graph_p3d"
_FORMAT_VERSION = "1.0"

#: Provenance string stamped into the bundle's root ``source`` attr (and
#: used as a fallback for the ``/periodicity`` group's ``source`` attr).
_SOURCE = "plot3d_utilities.glennht.graph_p3d"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_h5py():
    """Import and return the ``h5py`` module, or raise a clear error.

    ``h5py`` is optional (the ``hdf5`` extra) -- every public function in
    this module calls this helper before touching HDF5, so the failure
    mode is a single, actionable ``ImportError`` instead of a bare one
    from deep inside ``h5py``'s own import machinery.
    """
    try:
        import h5py  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without h5py
        raise ImportError(
            "graph_p3d support requires h5py. Install with "
            "`pip install \"plot3d[hdf5]\"` (or `pip install h5py`)."
        ) from exc
    return h5py


def _package_version() -> str:
    """Best-effort ``plot3d`` package version for the ``created_by`` attr."""
    try:
        from importlib.metadata import PackageNotFoundError, version
        try:
            return version("plot3d")
        except PackageNotFoundError:
            return "unknown"
    except Exception:  # pragma: no cover - defensive
        return "unknown"


def _native(value: Any) -> Any:
    """Convert an h5py attr/scalar (numpy scalar or bytes) to native Python."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _side_id(rec: Dict[str, Any]) -> int:
    """Extract a face-record's ``id`` (used on periodic-face sides).

    Returns ``-1`` when unset (``None`` or absent), matching the
    ``outer_faces`` ``id=-1`` null convention.
    """
    fid = rec.get("id")
    return int(fid) if fid is not None else -1


# ---------------------------------------------------------------------------
# write_graph_p3d
# ---------------------------------------------------------------------------

def write_graph_p3d(
    filename: str,
    blocks: List[Block],
    *,
    matched_faces: List[Dict[str, Any]],
    outer_faces: List[Dict[str, Any]],
    periodic_faces: Optional[List[Dict[str, Any]]] = None,
    periodicity: Optional[Dict[str, Any]] = None,
    surface_ids: Optional[Dict[Any, str]] = None,
    bcs: Optional[List[Any]] = None,
    rotation: Optional[Dict[str, Any]] = None,
    case_name: Optional[str] = None,
    mesh_file: Optional[str] = None,
    compression: Optional[str] = "gzip",
) -> None:
    """Write a single-file HDF5 ``.graph_p3d`` bundle.

    Bundles a Plot3D multi-block mesh's node coordinates, its
    block-connectivity graph (as produced by :func:`plot3d.connectivity` /
    :func:`plot3d.rotated_periodicity`), and boundary-condition surface
    tags into one portable HDF5 file -- a single-file alternative to the
    loose ``.xyz`` + ``.connectivity.json`` + boundary-conditions-YAML trio
    that :func:`plot3d.glennht.export_to_glennht_gpu` writes.

    glennht-gpu CANNOT read this bundle directly -- explode it back into
    the three loose files with :func:`graph_p3d_to_files` first.

    Args:
        filename (str): Output path, e.g. ``"wedge.graph_p3d"``.
        blocks (List[Block]): The multi-block mesh.
        matched_faces (List[Dict]): Face matches, in the same shape
            accepted by :func:`plot3d.glennht.write_connectivity_json`
            (i.e. straight from :func:`plot3d.connectivity` /
            :func:`plot3d.connectivity_fast`).
        outer_faces (List[Dict]): Outer (unmatched, possibly id-tagged)
            faces, same shape as :func:`write_connectivity_json`'s.
        periodic_faces (List[Dict], optional): Rotationally-periodic face
            pairs, as returned by :func:`plot3d.rotated_periodicity`. Each
            side's ``block1``/``block2`` sub-dict may carry an ``id`` (as
            ``rotated_periodicity`` attaches from the source ``Face``);
            it is preserved in the bundle even though ``connectivity.json``
            itself drops it (see :func:`plot3d.glennht.gpu._clean_face_record`).
        periodicity (Dict, optional): Global periodicity metadata (see
            :func:`write_connectivity_json`). Required whenever
            ``periodic_faces`` is non-empty.
        surface_ids (Dict, optional): ``{id: name}`` map. Defaults to the
            same ``{1: "inlet", ...}`` map :func:`write_connectivity_json`
            defaults to.
        bcs (List[GpuInletBC | GpuOutletBC | GpuWallBC], optional): If
            given, serialized to YAML text (the same text
            :func:`plot3d.glennht.write_boundary_conditions_yaml` would
            write to a file) and stored under ``/boundary_conditions``.
        rotation (Dict, optional): The run-yaml ``rotation:`` block (rotor
            speeds), included at the top of the stored BC YAML when ``bcs``
            is given -- keeps the exploded ``_boundary_conditions.yaml`` in
            parity with the direct :func:`export_to_glennht_gpu` output.
        case_name (str, optional): Stored as the root ``case_name`` attr
            (informational).
        mesh_file (str, optional): Value for the root ``mesh_file`` attr
            and the reconstructed connectivity dict's ``mesh_file`` field.
            Defaults to ``f"{case_name}.xyz"`` when ``case_name`` is given,
            else ``""``.
        compression (str, optional): h5py compression filter used for the
            grid ``x``/``y``/``z`` datasets. Defaults to ``"gzip"``; pass
            ``None`` to disable.

    Raises:
        ImportError: If ``h5py`` is not installed.
        ValueError: If ``periodic_faces`` is non-empty but ``periodicity``
            is ``None`` (mirrors :func:`write_connectivity_json`).
    """
    h5py = _require_h5py()

    periodic_faces = periodic_faces or []
    if periodic_faces and periodicity is None:
        raise ValueError(
            "periodic_faces was given but periodicity is missing; a "
            "graph_p3d bundle requires 'periodicity' whenever "
            "'periodic_faces' is non-empty (matches write_connectivity_json)."
        )

    if mesh_file is None:
        mesh_file = f"{case_name}.xyz" if case_name else ""

    with h5py.File(filename, "w") as f:
        f.attrs["format"] = _FORMAT
        f.attrs["format_version"] = _FORMAT_VERSION
        f.attrs["created_by"] = f"plot3d {_package_version()}"
        f.attrs["source"] = _SOURCE
        f.attrs["nblocks"] = int(len(blocks))
        f.attrs["case_name"] = case_name or ""
        f.attrs["mesh_file"] = mesh_file

        # --- /grid ---
        grid_grp = f.create_group("grid")
        for i, blk in enumerate(blocks):
            bgrp = grid_grp.create_group(f"block_{i}")
            bgrp.attrs["imax"] = int(blk.IMAX)
            bgrp.attrs["jmax"] = int(blk.JMAX)
            bgrp.attrs["kmax"] = int(blk.KMAX)
            for name, arr in (("x", blk.X), ("y", blk.Y), ("z", blk.Z)):
                bgrp.create_dataset(
                    name,
                    data=np.asarray(arr, dtype=np.float64),
                    compression=compression,
                )

        # --- /connectivity ---
        conn_grp = f.create_group("connectivity")

        fm_rows: List[List[int]] = []
        fm_perms: List[List[List[int]]] = []
        for fm in matched_faces:
            cleaned = _clean_face_match(fm)
            b1, b2 = cleaned["block1"], cleaned["block2"]
            fm_rows.append([
                b1["block_index"], *b1["lb"], *b1["ub"],
                b2["block_index"], *b2["lb"], *b2["ub"],
                cleaned["permutation_index"],
            ])
            fm_perms.append(cleaned["permutation_matrix"])
        conn_grp.create_dataset(
            "face_matches",
            data=np.array(fm_rows, dtype=np.int64).reshape(len(fm_rows), 15),
        )
        conn_grp.create_dataset(
            "face_match_permutations",
            data=np.array(fm_perms, dtype=np.int64).reshape(len(fm_perms), 2, 2),
        )

        of_rows: List[List[int]] = []
        for face in outer_faces:
            cleaned = _clean_outer_face(face)
            fid = cleaned["id"] if cleaned["id"] is not None else -1
            of_rows.append([cleaned["block_index"], *cleaned["lb"], *cleaned["ub"], fid])
        conn_grp.create_dataset(
            "outer_faces",
            data=np.array(of_rows, dtype=np.int64).reshape(len(of_rows), 8),
        )

        if periodic_faces:
            pf_rows: List[List[int]] = []
            pf_perms: List[List[List[int]]] = []
            rates: List[float] = []
            any_rate = False
            for fm in periodic_faces:
                cleaned = _clean_face_match(fm)
                b1, b2 = cleaned["block1"], cleaned["block2"]
                b1_id = _side_id(fm.get("block1") or {})
                b2_id = _side_id(fm.get("block2") or {})
                pf_rows.append([
                    b1["block_index"], b1_id, *b1["lb"], *b1["ub"],
                    b2["block_index"], b2_id, *b2["lb"], *b2["ub"],
                    cleaned["permutation_index"],
                ])
                pf_perms.append(cleaned["permutation_matrix"])
                rate = fm.get("rotation_angle_rad")
                if rate is not None:
                    any_rate = True
                    rates.append(float(rate))
                else:
                    rates.append(float("nan"))

            conn_grp.create_dataset(
                "periodic_faces",
                data=np.array(pf_rows, dtype=np.int64).reshape(len(pf_rows), 17),
            )
            conn_grp.create_dataset(
                "periodic_face_permutations",
                data=np.array(pf_perms, dtype=np.int64).reshape(len(pf_perms), 2, 2),
            )
            if any_rate:
                conn_grp.create_dataset(
                    "periodic_rotation_angles",
                    data=np.array(rates, dtype=np.float64),
                )

        conn_grp.create_dataset(
            "permutation_matrices",
            data=PERMUTATION_MATRICES.astype(np.int64),
        )

        sids_grp = conn_grp.create_group("surface_ids")
        for k, v in (surface_ids or _DEFAULT_SURFACE_IDS).items():
            sids_grp.attrs[str(k)] = str(v)

        # --- /periodicity ---
        if periodicity is not None:
            per_grp = f.create_group("periodicity")
            per_grp.attrs["nblades"] = int(periodicity.get("nblades", 0))
            per_grp.attrs["rotation_axis"] = str(periodicity.get("rotation_axis", ""))
            per_grp.attrs["rotation_angle_rad"] = float(periodicity.get("rotation_angle_rad", 0.0))
            per_grp.attrs["rotation_angle_deg"] = float(periodicity.get("rotation_angle_deg", 0.0))
            per_grp.attrs["convention"] = str(periodicity.get("convention", ""))
            per_grp.attrs["source"] = str(periodicity.get("source", _SOURCE))
            tm = np.asarray(
                periodicity.get("transformation_matrix", np.eye(3).tolist()),
                dtype=np.float64,
            )
            per_grp.create_dataset("transformation_matrix", data=tm)

        # --- /boundary_conditions ---
        if bcs is not None:
            text = _boundary_conditions_yaml_text(bcs, rotation=rotation)
            dt = h5py.string_dtype(encoding="utf-8")
            f.create_dataset("boundary_conditions", data=text, dtype=dt)


# ---------------------------------------------------------------------------
# read_graph_p3d
# ---------------------------------------------------------------------------

def read_graph_p3d(filename: str) -> Dict[str, Any]:
    """Read a ``.graph_p3d`` HDF5 bundle back into Python objects.

    Args:
        filename (str): Path to the ``.graph_p3d`` bundle.

    Returns:
        (Dict[str, Any]): with keys:

        - **blocks** (List[Block]): The multi-block mesh.
        - **connectivity** (Dict): Same schema as ``connectivity.json``
          (``mesh_file``, ``nblocks``, ``face_matches``, ``outer_faces``,
          ``periodic_faces``, ``periodicity`` (only if present),
          ``surface_ids``, ``permutation_matrices``) -- pass straight to
          :func:`write_connectivity_json`.
        - **surface_ids** (Dict[str, str]): ``{"<id>": name}``.
        - **periodicity** (Dict or None): Periodicity metadata, if present.
        - **boundary_conditions** (str or None): The stored BC YAML text.
        - **meta** (Dict): The root HDF5 attrs (``format``,
          ``format_version``, ``created_by``, ``source``, ``nblocks``,
          ``case_name``, ``mesh_file``).

    Raises:
        ImportError: If ``h5py`` is not installed.
    """
    h5py = _require_h5py()

    with h5py.File(filename, "r") as f:
        meta = {k: _native(v) for k, v in f.attrs.items()}
        nblocks = int(meta.get("nblocks", 0))

        # --- /grid ---
        blocks: List[Block] = []
        grid_grp = f["grid"]
        for i in range(nblocks):
            bgrp = grid_grp[f"block_{i}"]
            X = np.asarray(bgrp["x"][()], dtype=np.float64)
            Y = np.asarray(bgrp["y"][()], dtype=np.float64)
            Z = np.asarray(bgrp["z"][()], dtype=np.float64)
            blocks.append(Block(X, Y, Z))

        # --- /connectivity ---
        conn_grp = f["connectivity"]

        fm_arr = np.asarray(conn_grp["face_matches"][()], dtype=np.int64)
        fm_perm_arr = np.asarray(conn_grp["face_match_permutations"][()], dtype=np.int64)
        face_matches: List[Dict[str, Any]] = []
        for row, perm in zip(fm_arr.tolist(), fm_perm_arr.tolist()):
            (b1, b1lb0, b1lb1, b1lb2, b1ub0, b1ub1, b1ub2,
             b2, b2lb0, b2lb1, b2lb2, b2ub0, b2ub1, b2ub2, perm_idx) = row
            face_matches.append({
                "block1": {"block_index": int(b1),
                           "lb": [int(b1lb0), int(b1lb1), int(b1lb2)],
                           "ub": [int(b1ub0), int(b1ub1), int(b1ub2)]},
                "block2": {"block_index": int(b2),
                           "lb": [int(b2lb0), int(b2lb1), int(b2lb2)],
                           "ub": [int(b2ub0), int(b2ub1), int(b2ub2)]},
                "permutation_index": int(perm_idx),
                "permutation_matrix": [[int(c) for c in r] for r in perm],
            })

        of_arr = np.asarray(conn_grp["outer_faces"][()], dtype=np.int64)
        outer_faces: List[Dict[str, Any]] = []
        for row in of_arr.tolist():
            block_idx, lb0, lb1, lb2, ub0, ub1, ub2, fid = row
            outer_faces.append({
                "block_index": int(block_idx),
                "lb": [int(lb0), int(lb1), int(lb2)],
                "ub": [int(ub0), int(ub1), int(ub2)],
                "id": None if fid == -1 else int(fid),
            })

        periodic_faces: List[Dict[str, Any]] = []
        if "periodic_faces" in conn_grp:
            pf_arr = np.asarray(conn_grp["periodic_faces"][()], dtype=np.int64)
            pf_perm_arr = np.asarray(conn_grp["periodic_face_permutations"][()], dtype=np.int64)
            pf_rates = (
                np.asarray(conn_grp["periodic_rotation_angles"][()], dtype=np.float64)
                if "periodic_rotation_angles" in conn_grp else None
            )
            for idx, (row, perm) in enumerate(zip(pf_arr.tolist(), pf_perm_arr.tolist())):
                (b1, b1id, b1lb0, b1lb1, b1lb2, b1ub0, b1ub1, b1ub2,
                 b2, b2id, b2lb0, b2lb1, b2lb2, b2ub0, b2ub1, b2ub2, perm_idx) = row
                entry: Dict[str, Any] = {
                    "block1": {"block_index": int(b1),
                               "lb": [int(b1lb0), int(b1lb1), int(b1lb2)],
                               "ub": [int(b1ub0), int(b1ub1), int(b1ub2)],
                               "id": None if b1id == -1 else int(b1id)},
                    "block2": {"block_index": int(b2),
                               "lb": [int(b2lb0), int(b2lb1), int(b2lb2)],
                               "ub": [int(b2ub0), int(b2ub1), int(b2ub2)],
                               "id": None if b2id == -1 else int(b2id)},
                    "permutation_index": int(perm_idx),
                    "permutation_matrix": [[int(c) for c in r] for r in perm],
                }
                if pf_rates is not None:
                    rate = float(pf_rates[idx])
                    if not math.isnan(rate):
                        entry["rotation_angle_rad"] = rate
                periodic_faces.append(entry)

        permutation_matrices = np.asarray(
            conn_grp["permutation_matrices"][()], dtype=np.int64
        ).tolist()

        surface_ids: Dict[str, str] = {}
        if "surface_ids" in conn_grp:
            sids_grp = conn_grp["surface_ids"]
            for k, v in sids_grp.attrs.items():
                surface_ids[str(k)] = _native(v)

        # --- /periodicity ---
        periodicity: Optional[Dict[str, Any]] = None
        if "periodicity" in f:
            per_grp = f["periodicity"]
            periodicity = {
                "nblades": int(_native(per_grp.attrs["nblades"])),
                "rotation_axis": _native(per_grp.attrs["rotation_axis"]),
                "rotation_angle_rad": float(_native(per_grp.attrs["rotation_angle_rad"])),
                "rotation_angle_deg": float(_native(per_grp.attrs["rotation_angle_deg"])),
                "transformation_matrix": np.asarray(
                    per_grp["transformation_matrix"][()], dtype=np.float64
                ).tolist(),
                "convention": _native(per_grp.attrs["convention"]),
                "source": _native(per_grp.attrs["source"]),
            }

        # --- /boundary_conditions ---
        boundary_conditions: Optional[str] = None
        if "boundary_conditions" in f:
            raw = f["boundary_conditions"][()]
            boundary_conditions = raw.decode("utf-8") if isinstance(raw, bytes) else str(raw)

        connectivity: Dict[str, Any] = {
            "mesh_file": meta.get("mesh_file", ""),
            "nblocks": nblocks,
            "face_matches": face_matches,
            "outer_faces": outer_faces,
            "periodic_faces": periodic_faces,
            "surface_ids": surface_ids,
            "permutation_matrices": permutation_matrices,
        }
        if periodicity is not None:
            connectivity["periodicity"] = periodicity

    return {
        "blocks": blocks,
        "connectivity": connectivity,
        "surface_ids": surface_ids,
        "periodicity": periodicity,
        "boundary_conditions": boundary_conditions,
        "meta": meta,
    }


# ---------------------------------------------------------------------------
# graph_p3d_to_files (explode)
# ---------------------------------------------------------------------------

def graph_p3d_to_files(
    filename: str,
    out_dir: str,
    case_name: Optional[str] = None,
) -> Dict[str, str]:
    """Explode a ``.graph_p3d`` bundle back into the loose files
    glennht-gpu reads: ``{case}.xyz`` + ``{case}.connectivity.json`` (+
    ``{case}_boundary_conditions.yaml`` if the bundle carries BCs).

    Args:
        filename (str): Path to the ``.graph_p3d`` bundle.
        out_dir (str): Output directory (created if missing).
        case_name (str, optional): Base name for the written files.
            Defaults to the bundle's stored ``case_name`` attr, falling
            back to ``filename``'s basename with its extension stripped.

    Returns:
        (Dict[str, str]): ``{"grid": ..., "connectivity": ...,
        "boundary_conditions": ...}`` -- only the keys for files actually
        written are present (``boundary_conditions`` only when the bundle
        has one).
    """
    data = read_graph_p3d(filename)
    blocks = data["blocks"]
    conn = data["connectivity"]

    if case_name is None:
        case_name = data["meta"].get("case_name") or None
    if not case_name:
        case_name = os.path.splitext(os.path.basename(filename))[0]

    os.makedirs(out_dir, exist_ok=True)
    paths: Dict[str, str] = {}

    grid_path = os.path.join(out_dir, f"{case_name}.xyz")
    write_plot3D(grid_path, blocks, binary=False)
    paths["grid"] = grid_path

    conn_path = os.path.join(out_dir, f"{case_name}.connectivity.json")
    write_connectivity_json(
        blocks,
        conn["face_matches"],
        conn["outer_faces"],
        conn_path,
        periodic_faces=conn.get("periodic_faces") or [],
        periodicity=conn.get("periodicity"),
        surface_ids=conn.get("surface_ids"),
        mesh_file=f"{case_name}.xyz",
    )
    paths["connectivity"] = conn_path

    if data.get("boundary_conditions"):
        bc_path = os.path.join(out_dir, f"{case_name}_boundary_conditions.yaml")
        with open(bc_path, "w") as fh:
            fh.write(data["boundary_conditions"])
        paths["boundary_conditions"] = bc_path

    return paths
