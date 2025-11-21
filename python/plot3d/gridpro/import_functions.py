from typing import Dict, List, Tuple, Optional, Union, Set
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List, Tuple, Optional, Any
from plot3d import Block
from collections import defaultdict

def _parse_header(line: str) -> Optional[Tuple[int,int,int]]:
    parts = line.strip().split()
    if len(parts) < 3:
        return None
    try:
        i, j, k = int(float(parts[0])), int(float(parts[1])), int(float(parts[2]))
        return (i, j, k) if (i > 0 and j > 0 and k > 0) else None
    except ValueError:
        return None

def read_gridpro_to_blocks(
    filename: str,
    encoding: str = "utf-8",
    comment_prefixes: Tuple[str, ...] = ("#", "//"),
) -> List[Block]:
    """Read a structured GridPro text grid into Plot3D Block objects.

    Args:
        filename: Path to the GridPro grid file.
        encoding: Encoding to use while reading the text file.
        comment_prefixes: Line prefixes that should be treated as comments.

    Raises:
        ValueError: If a block does not contain the expected number of floats.

    Returns:
        List[Block]: A list of Block objects representing each GridPro block.
    """
    blocks = []

    with open(filename, "r", encoding=encoding, newline="") as f:
        # seek first header
        line = f.readline()
        while line:
            st = line.strip()
            if st and not any(st.startswith(pfx) for pfx in comment_prefixes):
                if _parse_header(st) is not None:
                    break
            line = f.readline()

        # main block loop
        while line:
            hdr = _parse_header(line)
            if hdr is None:
                line = f.readline()
                continue

            IMAX, JMAX, KMAX = hdr
            total_pts = IMAX * JMAX * KMAX
            expected = 3 * total_pts

            # use tqdm for this block
            with tqdm(total=total_pts, unit="pts", desc=f"Block ({IMAX},{JMAX},{KMAX})") as pbar:
                # np.fromfile reads ALL floats in one go
                arr = np.fromfile(f, sep=" ", count=expected, dtype=np.float64)
                if arr.size != expected:
                    raise ValueError(
                        f"Expected {expected} floats for block {hdr}, got {arr.size}."
                    )
                # update tqdm to 100%
                pbar.update(total_pts)

            coords = arr.reshape(total_pts, 3)
            X = coords[:, 0].reshape(IMAX, JMAX, KMAX)
            Y = coords[:, 1].reshape(IMAX, JMAX, KMAX)
            Z = coords[:, 2].reshape(IMAX, JMAX, KMAX)

            blocks.append(Block(X, Y, Z))  # type: ignore[name-defined]

            # advance to next header (skip blanks/comments)
            line = f.readline()
            while line and (not line.strip() or any(line.strip().startswith(p) for p in comment_prefixes)):
                line = f.readline()

    return blocks


def read_gridpro_connectivity(
    file_path: str,
    sb_zero_based_in_file: bool = True,     # False if sb1/sb2 are 1-based in file
    index_zero_based_in_file: bool = True,  # False if IMIN..KMAX are 1-based in file
    inlet_ids: Optional[Union[List[int], Dict[str, Any]]] = None,
    outlet_ids: Optional[Union[List[int], Dict[str, Any]]] = None,
    wall_ids: Optional[Union[List[int], Dict[str, Any]]] = None,
    symm_slip_ids: Optional[Union[List[int], Dict[str, Any]]] = None,
    custom_bc_ids: Optional[Dict[str, List[int]]] = None,
) -> Dict[str, object]:
    """Parse a GridPro connectivity file into plot3d.connectivity_fast-like structures.

    Expected P-line layout (as provided):
    ``P pid sb1 sf1 sb2 sf2 fmap L1i L1j L1k H1i H1j H1k L2i L2j L2k H2i H2j H2k pty lbid``

    Boundary-condition IDs default to GridPro's PTY values but can be overridden
    via the provided *_ids arguments (either a list of PTY ints or
    ``{"name": str, "ids": [...]}``) or extended via
    ``custom_bc_ids={"name": [ids]}``.

    Args:
        file_path: Path to the GridPro connectivity file.
        sb_zero_based_in_file: Set False if ``sb1/sb2`` are 1-based in the file.
        index_zero_based_in_file: Set False if IMIN..KMAX are 1-based in the file.
        inlet_ids: PTY ids to treat as inlet; defaults to ``[5]``.
            Pass a dict ``{"name": "...", "ids": [...]}`` to rename the group, or
            provide a mapping like ``{"Inlet": [...], "Pressure Inlet": [...]}``
            to emit multiple named inlet groups. Each named group remembers its
            base type in the returned ``bc_group`` mapping.
        outlet_ids: PTY ids to treat as outlet; defaults to ``[6]`` (same dict
            options as ``inlet_ids``).
        wall_ids: PTY ids to treat as wall; defaults to ``[2]`` (same dict options).
        symm_slip_ids: PTY ids to treat as symmetry/slip; defaults to ``[4]`` (same dict options).
        custom_bc_ids: Optional mapping of ``group_name -> [pty ids]`` for any
            additional boundary-condition groupings (e.g., ``{"cooling": [500]}``).

    Returns:
        Dict[str, object]: A dictionary with keys:
            - ``face_matches``: list of face-pair dictionaries for connected blocks.
            - ``all_surfaces``: list of faces on the exterior (no neighbor).
            - ``ungrouped_surfaces``: exterior faces whose PTY does not belong to any BC group.
            - ``bc_group``: mapping of bc name to a metadata dict
              ``{"type": base_type, "pty_ids": [...], "faces": [...]}``. Each
              face dict also includes ``bc_name`` (group label) and ``bc_type``.
            - ``gif_faces``: list of faces tagged as GIF (pty 12..21 or 1000).
            - ``periodic_faces``: list of paired periodic faces (pty 3).
            - ``volume_zones``: list describing volume zone type per superblock.
            - ``blocksizes``: list of per-superblock cell counts (I*J*K).
            - ``patches``: pandas DataFrame of raw patch records.

    Examples:
        Basic usage with defaults::

            data = read_gridpro_connectivity("connectivity.dat")
            inlet_faces = data["bc_group"]["inlet"]["faces"]

        Override built-in BC ids and add a custom group::

            data = read_gridpro_connectivity(
                "connectivity.dat",
                inlet_ids=[5, 105],
                custom_bc_ids={"cooling_hole1": [500]},
            )
            cooling_faces = data["bc_group"]["cooling_hole1"]["faces"]
    """
    # ---------------------------- parsing ----------------------------
    superblock_ptys: List[int] = []
    superblock_sizes: List[int] = []
    patch_rows: List[Dict[str, object]] = []

    sb_offset = 0 if sb_zero_based_in_file else -1
    idx_offset = 0 if index_zero_based_in_file else -1  # convert to 0-based if file is 1-based

    def parse_patch(tokens: List[str]) -> None:
        # Expect at least 21 tokens: P + 20 fields
        # P pid sb1 sf1 sb2 sf2 fmap  L1i L1j L1k  H1i H1j H1k  L2i L2j L2k  H2i H2j H2k  pty lbid
        if len(tokens) < 21:
            raise ValueError(f"Patch line too short: {' '.join(tokens)}")

        try:
            pid  = int(tokens[1])
            sb1_raw = int(tokens[2])
            sf1  = int(tokens[3])
            sb2_raw = int(tokens[4])
            sf2  = int(tokens[5])
            
            sb1 = sb1_raw - 1            # 1..76 -> 0..75
            sb2 = (sb2_raw - 1) if sb2_raw > 0 else -1   # 0 -> -1 (no neighbor), 1..76 -> 0..75

            fmap = tokens[6]  # keep as string like '012'

            # Triplets
            L1i, L1j, L1k = (int(tokens[7])  + idx_offset,
                              int(tokens[8])  + idx_offset,
                              int(tokens[9])  + idx_offset)
            H1i, H1j, H1k = (int(tokens[10]) + idx_offset,
                              int(tokens[11]) + idx_offset,
                              int(tokens[12]) + idx_offset)
            L2i, L2j, L2k = (int(tokens[13]) + idx_offset,
                              int(tokens[14]) + idx_offset,
                              int(tokens[15]) + idx_offset)
            H2i, H2j, H2k = (int(tokens[16]) + idx_offset,
                              int(tokens[17]) + idx_offset,
                              int(tokens[18]) + idx_offset)

            pty  = int(tokens[19])
            lbid = int(tokens[20])
        except Exception as e:
            raise ValueError(f"Failed to parse patch line: {' '.join(tokens)}") from e

        patch_rows.append({
            "pid": pid, "sb1": sb1, "sf1": sf1, "sb2": sb2, "sf2": sf2,
            "fmap": fmap,
            "L1i": L1i, "L1j": L1j, "L1k": L1k, "H1i": H1i, "H1j": H1j, "H1k": H1k,
            "L2i": L2i, "L2j": L2j, "L2k": L2k, "H2i": H2i, "H2j": H2j, "H2k": H2k,
            "pty": pty, "lbid": lbid
        })

    with open(file_path, "r", encoding="utf-8", newline="") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            toks = line.split()
            tag = toks[0]
            if tag == "SB":
                superblock_ptys.append(int(toks[-2]))
                try:
                    I_dim = int(toks[2])
                    J_dim = int(toks[3])
                    K_dim = int(toks[4])
                    superblock_sizes.append(I_dim * J_dim * K_dim)
                except (IndexError, ValueError):
                    superblock_sizes.append(0)
            elif tag == "P":
                parse_patch(toks)

    patches_df = pd.DataFrame(patch_rows)

    # ---------------------------- helpers ----------------------------
    def face_dict(sb: int, imin: int, jmin: int, kmin: int, imax: int, jmax: int, kmax: int, pty:int=-1) -> Dict[str, int]:
        # Matches plot3d.connectivity_fast single-face dict shape
        return {
            "block_index": sb,
            "IMIN": imin, "JMIN": jmin, "KMIN": kmin,
            "IMAX": imax, "JMAX": jmax, "KMAX": kmax,
            "id": pty
        }

    def pair_dict(sb1: int, r1: Tuple[int,int,int,int,int,int],
                  sb2: int, r2: Tuple[int,int,int,int,int,int]) -> Dict[str, Dict[str, int]]:
        L1i,L1j,L1k,H1i,H1j,H1k = r1
        L2i,L2j,L2k,H2i,H2j,H2k = r2
        return {
            "block1": face_dict(sb1, L1i, L1j, L1k, H1i, H1j, H1k),
            "block2": face_dict(sb2, L2i, L2j, L2k, H2i, H2j, H2k),
        }

    # ------------------------- build outputs -------------------------
    # Connections: sb2 != -1 (i.e., not "none") and pty in {1,3}
    connections: List[Dict[str, Dict[str, int]]] = []
    for p in patches_df.itertuples(index=False):
        if p.sb2 != -1 and (p.pty in (1, 3)):
            connections.append(
                pair_dict(
                    p.sb1, (p.L1i, p.L1j, p.L1k, p.H1i, p.H1j, p.H1k),
                    p.sb2, (p.L2i, p.L2j, p.L2k, p.H2i, p.H2j, p.H2k),
                )
            )

    # Boundary-condition groups (by pty)
    def bc_faces_for(name: str, pty_vals: List[int], bc_type: str) -> List[Dict[str, int]]:
        targets = set(pty_vals)
        faces: List[Dict[str, int]] = []
        for p in patches_df.itertuples(index=False):
            if p.pty in targets:
                face = face_dict(p.sb1, p.L1i, p.L1j, p.L1k, p.H1i, p.H1j, p.H1k, p.pty)
                face["bc_name"] = name  # type: ignore
                face["bc_type"] = bc_type  # type: ignore
                faces.append(face)
        return faces

    def normalize_bc_arg(
        arg: Optional[Union[List[int], Dict[str, Any]]],
        bc_type: str,
        default_ids: List[int],
    ) -> List[Tuple[str, List[int]]]:
        def coerce_ids(values: Any) -> List[int]:
            if values is None:
                return []
            # Allow nested {"ids": [...]} dictionaries to be passed directly
            if isinstance(values, dict):
                return coerce_ids(values.get("ids"))
            if isinstance(values, (list, tuple, set)):
                return [int(v) for v in values if v is not None]
            return [int(values)]

        if arg is None:
            return [(bc_type, coerce_ids(default_ids))]

        if isinstance(arg, dict):
            # Support legacy {"name": str, "ids": [...]} as well as
            # {"CustomName": [...], "AnotherName": [...]} mappings.
            if "ids" in arg or "name" in arg:
                ids = coerce_ids(arg.get("ids", default_ids))
                name = str(arg.get("name", bc_type) or bc_type)
                return [(name, ids)]

            normalized: List[Tuple[str, List[int]]] = []
            for name, ids_val in arg.items():
                ids_list = coerce_ids(ids_val)
                if ids_list:
                    normalized.append((str(name), ids_list))
            return normalized

        ids = coerce_ids(arg)
        return [(bc_type, ids)]

    bc_group: Dict[str, Dict[str, Any]] = {}
    grouped_bc_ids: Set[int] = set()

    def add_bc_group(name: str, ids: List[int], bc_type: str) -> None:
        normalized_ids = sorted({int(i) for i in ids if i is not None})
        faces = bc_faces_for(name, normalized_ids, bc_type) if normalized_ids else []
        if not faces:
            return
        bc_group[name] = {
            "type": bc_type,
            "pty_ids": normalized_ids,
            "faces": faces,
        }
        grouped_bc_ids.update(normalized_ids)

    for bc_type, group in [
        ("inlet", normalize_bc_arg(inlet_ids, "inlet", [5])),
        ("outlet", normalize_bc_arg(outlet_ids, "outlet", [6])),
        ("symm_slip", normalize_bc_arg(symm_slip_ids, "symm_slip", [4])),
        ("wall", normalize_bc_arg(wall_ids, "wall", [2])),
    ]:
        for name, ids in group:
            add_bc_group(name, ids, bc_type)

    # add any user-defined boundary groups keyed by name -> list of pty ids
    if custom_bc_ids:
        for name, ids in custom_bc_ids.items():
            add_bc_group(name, ids or [], name)

    # Outer faces: sb2 == -1 (means '0' in file when sb_zero_based_in_file=False; or literal 0 if already zero-based)
    # To cover both cases, treat original token value 0 as "no neighbor".
    # Since we already applied sb_offset, "no neighbor" now appears as -1.
    all_surfaces: List[Dict[str, int]] = []
    ungrouped_surfaces: List[Dict[str, int]] = []
    for p in patches_df.itertuples(index=False):
        if p.sb2 == -1:
            face = face_dict(p.sb1, p.L1i, p.L1j, p.L1k, p.H1i, p.H1j, p.H1k, p.pty)
            all_surfaces.append(face)
            if p.pty not in grouped_bc_ids:
                ungrouped_surfaces.append(face.copy())

    # Periodic faces: pty == 3 (explicit pairs)
    periodic_faces: List[Dict[str, Dict[str, int]]] = []
    for p in patches_df.itertuples(index=False):
        if p.pty == 3 and p.sb2 != -1:
            periodic_faces.append(
                pair_dict(
                    p.sb1, (p.L1i, p.L1j, p.L1k, p.H1i, p.H1j, p.H1k),
                    p.sb2, (p.L2i, p.L2j, p.L2k, p.H2i, p.H2j, p.H2k),
                )
            )
            
    # GIF faces grouped by sf1 for pty in 12..21 or == 1000
    gif_faces: List[Dict[str, int]] = []
    for p in patches_df.itertuples(index=False):
#? This is how we identify gifs inside of a grid pro connectivity file. 
#? If the PTY is between 12 and 21 these are gifs
        if (12 <= p.pty <= 21) or (p.pty == 1000):
            face_temp = face_dict(p.sb1, p.L1i, p.L1j, p.L1k, p.H1i, p.H1j, p.H1k)
            face_temp["id"] = p.pty  # type: ignore
            gif_faces.append(face_temp)

    # Volume zones (unique superblock ptys; odd→fluid, even→solid) with contiguous ids
    volume_zones: List[Dict[str, object]] = []
    prevZoneType = ""
    cid = 0
    for id,v in enumerate(superblock_ptys):
        zone = "fluid" if (v % 2 != 0) else "solid"
        if zone is not prevZoneType:
            cid += 1
            prevZoneType = zone
        volume_zones.append({"block_index": id, "zone_type": zone, "contiguous_index": cid})
        

    return {
        "face_matches": connections,
        "all_surfaces": all_surfaces,
        "ungrouped_surfaces": ungrouped_surfaces,
        "bc_group": bc_group,
        "gif_faces": gif_faces,
        "periodic_faces": periodic_faces,
        "volume_zones": volume_zones,
        "blocksizes": superblock_sizes,
        "patches": patches_df,
    }


def bc_faces_by_type(
    bc_group: Dict[str, Any],
    bc_type: str,
    *,
    unique: bool = False,
) -> List[Dict[str, int]]:
    """Collect faces for a specific boundary-condition type across all named groups.

    Args:
        bc_group: Mapping returned by :func:`read_gridpro_connectivity` (supports
            both the legacy ``{name: [faces]}`` and the new
            ``{name: {"type": ..., "faces": [...]}}`` shapes).
        bc_type: Base BC type (``"inlet"``, ``"outlet"``, ``"wall"``, ``"symm_slip"``).
        unique: When True, keep only the first face per unique ``id``.

    Returns:
        List of face dictionaries belonging to the requested base type.
    """
    faces: List[Dict[str, int]] = []
    seen: Set[int] = set()
    target = (bc_type or "").lower()

    for entry in bc_group.values():
        entry_type = None
        entry_faces: List[Dict[str, int]] = []

        if isinstance(entry, dict) and "faces" in entry:
            entry_faces = entry.get("faces", []) or []
            entry_type = entry.get("type") or entry.get("bc_type")
        elif isinstance(entry, list):
            entry_faces = entry
        else:
            continue

        for face in entry_faces:
            ftype = entry_type or face.get("bc_type") or face.get("bc_name")
            if not isinstance(ftype, str) or ftype.lower() != target:
                continue
            fid = face.get("id")
            if unique and isinstance(fid, int):
                if fid in seen:
                    continue
                seen.add(fid)
            faces.append(face)
    return faces
