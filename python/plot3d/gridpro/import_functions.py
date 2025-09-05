from typing import Dict, List, Tuple, Optional
import numpy as np
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

def read_gridpro_to_blocks(filename: str,encoding: str = "utf-8", comment_prefixes: Tuple[str, ...] = ("#", "//")) -> List[Block]:
    """_summary_

    Args:
        filename (str): _description_
        encoding (str, optional): _description_. Defaults to "utf-8".
        comment_prefixes (Tuple[str, ...], optional): _description_. Defaults to ("#", "//").

    Raises:
        ValueError: _description_

    Returns:
        List[Block]: _description_
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
) -> Dict[str, object]:
    """
        Parse a GridPro connectivity file and return dictionaries compatible with
        plot3d.connectivity_fast-style shapes.

        Expected P-line layout (as provided):
        P pid sb1 sf1 sb2 sf2 fmap L1i L1j L1k H1i H1j H1k L2i L2j L2k H2i H2j H2k pty lbid

        Returns:
        {
            "face_matches":    List[{"block1":{...}, "block2":{...}}],
            "outer_faces":    List[{... face dict ...}],
            "bc_group":       {"inlet":[...], "outlet":[...], "symm_slip":[...], "wall":[...]},
            "gifs":      List[List[face_dict]],  # grouped by sf1
            "periodic_faces": List[{"block1":{...}, "block2":{...}}],
            "volume_zones":   List[{"pty":int, "zone":str, "contiguous_index":int}],
        }
    """
    # ---------------------------- parsing ----------------------------
    superblock_ptys: List[int] = []
    patches: List[Dict[str, object]] = []

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

        patches.append({
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
            elif tag == "P":
                parse_patch(toks)

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
    for p in patches:
        if p["sb2"] != -1 and (p["pty"] in (1, 3)):
            connections.append(
                pair_dict(
                    p["sb1"], (p["L1i"], p["L1j"], p["L1k"], p["H1i"], p["H1j"], p["H1k"]), # type: ignore
                    p["sb2"], (p["L2i"], p["L2j"], p["L2k"], p["H2i"], p["H2j"], p["H2k"]), # type: ignore
                )
            )

    # Outer faces: sb2 == -1 (means '0' in file when sb_zero_based_in_file=False; or literal 0 if already zero-based)
    # To cover both cases, treat original token value 0 as "no neighbor".
    # Since we already applied sb_offset, "no neighbor" now appears as -1.
    outer_faces: List[Dict[str, int]] = []
    pty_exclude = [6,5,4,2]
    for p in patches:
        if p["sb2"] == -1 and (p["pty"] not in pty_exclude):
            outer_faces.append(
                face_dict(p["sb1"], p["L1i"], p["L1j"], p["L1k"], p["H1i"], p["H1j"], p["H1k"], p["pty"]) # type: ignore
            )

    # Boundary-condition groups (by pty)
    def bc_faces_for(pty_val: int) -> List[Dict[str, int]]:
        return [
            face_dict(p["sb1"], p["L1i"], p["L1j"], p["L1k"], p["H1i"], p["H1j"], p["H1k"], p["pty"]) # type: ignore
            for p in patches if p["pty"] == pty_val
        ]

    bc_group = {
        "inlet":     bc_faces_for(5),
        "outlet":    bc_faces_for(6),
        "symm_slip": bc_faces_for(4),
        "wall":      bc_faces_for(2),
    }
    
    # Periodic faces: pty == 3 (explicit pairs)
    periodic_faces: List[Dict[str, Dict[str, int]]] = []
    for p in patches:
        if p["pty"] == 3 and p["sb2"] != -1:
            periodic_faces.append(
                pair_dict(
                    p["sb1"], (p["L1i"], p["L1j"], p["L1k"], p["H1i"], p["H1j"], p["H1k"]), # type: ignore
                    p["sb2"], (p["L2i"], p["L2j"], p["L2k"], p["H2i"], p["H2j"], p["H2k"]), # type: ignore
                )
            )
            
    # GIF faces grouped by sf1 for pty in 12..21 or == 1000
    gif_faces: List[Dict[str, int]] = []
    for p in patches:
        if (12 <= p["pty"] <= 21) or (p["pty"] == 1000): # type: ignore
            face_temp = face_dict(p["sb1"], p["L1i"], p["L1j"], p["L1k"], p["H1i"], p["H1j"], p["H1k"]) # type: ignore
            face_temp["id"] = p["pty"] # type: ignore
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
        "outer_faces": outer_faces,
        "bc_group": bc_group,
        "gif_faces": gif_faces,
        "periodic_faces": periodic_faces,
        "volume_zones": volume_zones,
    }
