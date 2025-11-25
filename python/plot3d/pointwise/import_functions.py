from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from plot3d import Block, connectivity_fast
from plot3d.glennht.class_definitions import (
    BCGroup,
    BoundaryConditionType,
    GIF,
    InletBC,
    Job,
    OutletBC,
    SymmetricSlipBC,
    VolumeZone,
    WallBC,
)
from plot3d.glennht.export_functions import (
    export_to_boundary_condition,
    export_to_glennht_conn,
)


@dataclass
class Patch:
    """Patch entry from a Pointwise .inp file."""

    i1: int
    i2: int
    j1: int
    j2: int
    k1: int
    k2: int
    tag: int


@dataclass
class BlockDecl:
    """Block declaration from a Pointwise .inp file."""

    name: str
    ni: int
    nj: int
    nk: int
    npatch_declared: int
    patches: List[Patch]


@dataclass
class Plot3DInp:
    header_flag: int
    nblocks: int
    blocks: List[BlockDecl]


@dataclass
class FvbndPatch:
    """Single boundary entry from a Pointwise .fvbnd file."""

    bc_id: int
    block_index: int
    i1: int
    i2: int
    j1: int
    j2: int
    k1: int
    k2: int

    def to_face(
        self,
        *,
        bc_name: str,
        bc_type: str,
    ) -> Dict[str, Any]:
        return {
            "block_index": self.block_index,
            "IMIN": self.i1,
            "JMIN": self.j1,
            "KMIN": self.k1,
            "IMAX": self.i2,
            "JMAX": self.j2,
            "KMAX": self.k2,
            "id": self.bc_id,
            "bc_name": bc_name,
            "bc_type": bc_type,
        }


def read_pointwise_inp(
    path: str | Path,
    *,
    index_zero_based_in_file: bool = False,
) -> Plot3DInp:
    """Parse a Pointwise .inp connectivity file.

    Args:
        path: Location of the .inp file.
        index_zero_based_in_file: Set to True if the patch indices in the
            file are already zero-based. Pointwise typically writes them
            one-based, so the default converts to zero-based.

    Returns:
        Parsed :class:`Plot3DInp` structure.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]

    if len(lines) < 2:
        raise ValueError(f"{path} does not look like a Pointwise .inp file")

    header_flag = int(lines[0].split()[0])
    nblocks = int(lines[1].split()[0])
    idx_offset = 0 if index_zero_based_in_file else -1

    def _offset(val: int) -> int:
        if index_zero_based_in_file or val < 0:
            return val
        return val + idx_offset

    blocks: List[BlockDecl] = []
    i = 2
    while i < len(lines):
        dims = [int(x) for x in lines[i].split()]
        if len(dims) != 3:
            raise ValueError(f"Expected 3 integers for dims at line {i}: {lines[i]}")
        ni, nj, nk = dims

        name = lines[i + 1].strip()
        npatch_declared = int(lines[i + 2].split()[0])

        patches: List[Patch] = []
        j = i + 3
        while j < len(lines):
            tokens = lines[j].split()
            try:
                ints = [int(x) for x in tokens]
            except ValueError:
                break

            if len(ints) == 7:
                i1, i2, j1, j2, k1, k2, tag = ints
                patches.append(
                    Patch(
                        _offset(i1),
                        _offset(i2),
                        _offset(j1),
                        _offset(j2),
                        _offset(k1),
                        _offset(k2),
                        tag,
                    )
                )
                j += 1
            elif len(ints) == 3:
                break
            else:
                j += 1

        blocks.append(
            BlockDecl(
                name=name,
                ni=ni,
                nj=nj,
                nk=nk,
                npatch_declared=npatch_declared,
                patches=patches,
            )
        )
        i = j

    return Plot3DInp(header_flag=header_flag, nblocks=nblocks, blocks=blocks)


def read_pointwise_fvbnd(
    path: str | Path,
    *,
    block_zero_based_in_file: bool = False,
    index_zero_based_in_file: bool = False,
) -> Tuple[Dict[int, str], List[FvbndPatch]]:
    """Parse a Pointwise .fvbnd file into a list of boundary patches.

    Args:
        path: Location of the .fvbnd file.
        block_zero_based_in_file: Set True if block ids are already zero-based.
        index_zero_based_in_file: Set True if i/j/k bounds are already zero-based.

    Returns:
        Tuple of (bc_id -> bc_name mapping, list of :class:`FvbndPatch`).
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        raw_lines = [ln.strip() for ln in f if ln.strip()]

    if not raw_lines:
        raise ValueError(f"{path} is empty")

    names: Dict[int, str] = {}
    patches: List[FvbndPatch] = []

    # Skip the first header line (e.g., "FVBND 1 3")
    start = 1
    while start < len(raw_lines):
        line = raw_lines[start]
        if line.upper().startswith("BOUNDARIES"):
            start += 1
            break
        names[len(names) + 1] = line.strip()
        start += 1

    block_offset = 0 if block_zero_based_in_file else -1
    idx_offset = 0 if index_zero_based_in_file else -1

    for line in raw_lines[start:]:
        tokens = line.split()
        if len(tokens) != 8:
            raise ValueError(
                f"Expected 8 integer columns for a boundary row, got {len(tokens)}: {line}"
            )

        bc_id = int(tokens[0])
        block_id = int(tokens[1]) + block_offset
        i1, i2, j1, j2, k1, k2 = (int(tok) + idx_offset for tok in tokens[2:])

        patches.append(
            FvbndPatch(
                bc_id=bc_id,
                block_index=block_id,
                i1=i1,
                i2=i2,
                j1=j1,
                j2=j2,
                k1=k1,
                k2=k2,
            )
        )

    return names, patches


def _guess_bc_type(name: str, overrides: Optional[Dict[int | str, str]], bc_id: int) -> str:
    """Infer a GlennHT BC type keyword from a Pointwise BC name or override map."""
    if overrides:
        if bc_id in overrides:
            return overrides[bc_id]
        lowered = name.lower()
        if lowered in overrides:
            return overrides[lowered]

    lowered_name = name.lower()
    if "inlet" in lowered_name or "inflow" in lowered_name:
        return "inlet"
    if "outlet" in lowered_name or "outflow" in lowered_name:
        return "outlet"
    if "sym" in lowered_name or "slip" in lowered_name:
        return "symm_slip"
    if "far" in lowered_name:
        return "symm_slip"
    return "wall"


def build_pointwise_bc_group(
    boundary_patches: Iterable[FvbndPatch],
    bc_names: Dict[int, str],
    *,
    bc_type_overrides: Optional[Dict[int | str, str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Build a GridPro-style ``bc_group`` mapping from Pointwise boundary data."""
    bc_group: Dict[str, Dict[str, Any]] = {}

    for patch in boundary_patches:
        bc_name = bc_names.get(patch.bc_id, f"BC_{patch.bc_id}")
        bc_type = _guess_bc_type(bc_name, bc_type_overrides, patch.bc_id)
        face = patch.to_face(bc_name=bc_name, bc_type=bc_type)

        if bc_name not in bc_group:
            bc_group[bc_name] = {"type": bc_type, "pty_ids": [patch.bc_id], "faces": []}
        entry = bc_group[bc_name]
        faces = entry.setdefault("faces", [])
        faces.append(face)

        if "pty_ids" in entry:
            ids = set(entry["pty_ids"])
            ids.add(patch.bc_id)
            entry["pty_ids"] = sorted(ids)

    return bc_group


def bcgroup_from_pointwise(bc_group: Dict[str, Dict[str, Any]]) -> BCGroup:
    """Convert a ``bc_group`` mapping into GlennHT ``BCGroup`` dataclasses."""
    inlets: List[InletBC] = []
    outlets: List[OutletBC] = []
    slips: List[SymmetricSlipBC] = []
    walls: List[WallBC] = []

    for name, entry in bc_group.items():
        bc_type = entry.get("type", "wall")
        faces = entry.get("faces", [])
        ids: List[int] = sorted(
            {face.get("id") for face in faces if isinstance(face, dict) and face.get("id") is not None}
        )
        if not ids:
            ids = entry.get("pty_ids", []) or [0]

        for sid in ids:
            if bc_type == "inlet":
                inlets.append(
                    InletBC(
                        BCType=BoundaryConditionType.Inlet,
                        SurfaceID=int(sid),
                        Name=name,
                        inlet_subType=int(sid),
                    )
                )
            elif bc_type == "outlet":
                outlets.append(
                    OutletBC(
                        BCType=BoundaryConditionType.Outlet,
                        SurfaceID=int(sid),
                        Name=name,
                        outlet_subType=int(sid),
                    )
                )
            elif bc_type == "symm_slip":
                slips.append(
                    SymmetricSlipBC(
                        BCType=BoundaryConditionType.SymmetryOrSlip,
                        SurfaceID=int(sid),
                        Name=name,
                        slip_subType=int(sid),
                    )
                )
            else:
                walls.append(
                    WallBC(
                        BCType=BoundaryConditionType.Wall,
                        SurfaceID=int(sid),
                        Name=name,
                        wall_subType=int(sid),
                    )
                )

    return BCGroup(Inlets=inlets, Outlets=outlets, SymmetricSlips=slips, Walls=walls)


def build_pointwise_connectivity(
    blocks: List[Block],
    fvbnd_path: str | Path,
    *,
    inp_path: str | Path | None = None,
    bc_type_overrides: Optional[Dict[int | str, str]] = None,
) -> Dict[str, Any]:
    """Create a connectivity bundle from Pointwise inputs."""
    bc_names, boundary_patches = read_pointwise_fvbnd(fvbnd_path)
    bc_group = build_pointwise_bc_group(boundary_patches, bc_names, bc_type_overrides=bc_type_overrides)

    face_matches, _ = connectivity_fast(blocks)

    outer_faces: List[Dict[str, Any]] = []
    for entry in bc_group.values():
        faces = entry.get("faces", [])
        outer_faces.extend(faces)

    block_sizes: List[int] = []
    if inp_path:
        inp = read_pointwise_inp(inp_path)
        block_sizes = [(blk.ni - 1) * (blk.nj - 1) * (blk.nk - 1) for blk in inp.blocks]
    else:
        block_sizes = [
            (blk.IMAX - 1) * (blk.JMAX - 1) * (blk.KMAX - 1) for blk in blocks  # type: ignore[attr-defined]
        ]

    volume_zones = [{"block_index": i, "zone_type": "fluid", "contiguous_index": 1} for i in range(len(blocks))]

    return {
        "face_matches": face_matches,
        "bc_group": bc_group,
        "gif_faces": [],
        "gif_pairs": [],
        "volume_zones": volume_zones,
        "blocksizes": block_sizes,
        "bc_names": bc_names,
        "boundary_patches": boundary_patches,
    }


def export_pointwise_to_glennht(
    blocks: List[Block],
    fvbnd_path: str | Path,
    *,
    inp_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    conn_filename: str = "connectivity.ght_conn",
    bc_filename: str = "boundary_conditions.bcs",
    job: Job | None = None,
    bc_type_overrides: Optional[Dict[int | str, str]] = None,
    gif_pairs: Optional[List[GIF | Dict[str, Any]]] = None,
    volume_zones: Optional[List[VolumeZone | Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Convenience wrapper to write GlennHT connectivity + BC files for Pointwise meshes."""
    bundle = build_pointwise_connectivity(
        blocks,
        fvbnd_path,
        inp_path=inp_path,
        bc_type_overrides=bc_type_overrides,
    )

    out_dir = Path(output_dir) if output_dir else Path(fvbnd_path).resolve().parent
    conn_path = out_dir / conn_filename
    bc_path = out_dir / bc_filename

    gif_pairs = gif_pairs or []
    gif_faces = bundle.get("gif_faces", [])
    volume_zones = volume_zones or bundle.get("volume_zones", [])

    outer_faces: List[Dict[str, Any]] = []
    for entry in bundle["bc_group"].values():
        faces = entry.get("faces", [])
        outer_faces.extend(faces)
    outer_faces = sorted(
        outer_faces,
        key=lambda f: (
            f.get("id", 0),
            f.get("block_index", 0),
            f.get("IMIN", 0),
            f.get("JMIN", 0),
            f.get("KMIN", 0),
        ),
    )

    export_to_glennht_conn(
        bundle["face_matches"],
        outer_faces,
        str(conn_path),
        gif_pairs,
        gif_faces,
        volume_zones,
    )

    bcgroup_dataclass = bcgroup_from_pointwise(bundle["bc_group"])
    if job is None:
        job = Job()
    job.JobFiles.ConnFILE = conn_path.name
    job.JobFiles.BCSpecFILE = bc_path.name
    export_to_boundary_condition(
        str(bc_path),
        job_settings=job,
        bc_group=bcgroup_dataclass,
        gif_pairs=gif_pairs,
        volume_zones=volume_zones,
    )

    bundle["conn_path"] = conn_path
    bundle["bc_path"] = bc_path
    bundle["job"] = job
    return bundle
