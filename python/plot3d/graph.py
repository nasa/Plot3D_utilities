# graph.py  (inside plot3d package)
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional
import inspect
import sys

if sys.platform != "win32":
    try:
        import pymetis  # type: ignore
        HAS_PYMETIS = True
    except Exception:  # pragma: no cover - optional dependency
        pymetis = None  # type: ignore
        HAS_PYMETIS = False
else:
    pymetis = None  # type: ignore
    HAS_PYMETIS = False
from .block import Block  # <-- use your existing Block


# ---------------------------------------------------------------------------
# Build weighted graph from connectivity_fast face_matches
# ---------------------------------------------------------------------------
def build_weighted_graph_from_face_matches(
    face_matches: List[dict],
    n_blocks: int,
    aggregate: str = "sum",
    ignore_self_matches: bool = True,
) -> Tuple[Dict[int, List[int]], Dict[int, Dict[int, int]]]:
    """
    Convert connectivity_fast face_matches into adjacency + weights.
    
    Note:
        Aggregate controls how we combine weights when the same two blocks are connected by multiple faces.
    
    Why this matters:
        connectivity_fast can return more than one face between the same pair of blocks 
        (e.g. a block is split and has two non-contiguous interfaces with its neighbor).
        Each face gives a weight = dI * dJ * dK (number of shared nodes).
        METIS expects one edge per block-pair, with a single weight.
        So if there are multiple faces, we must decide how to merge them. That’s what aggregate does.
        
    Args:
        face_matches (List[dict]): Output from connectivity_fast
        n_blocks (int): number of blocks in a mesh 
        aggregate (str, optional): 'sum'|'max'|'min'. Controls how weights are combined. Defaults to "sum".
        ignore_self_matches (bool, optional): ignores self matching (i==j). Defaults to True.

    Raises:
        ValueError: if aggregate is not one of 'sum','max','min'

    Returns:
        Tuple[Dict[int, List[int]], Dict[int, Dict[int, int]]]: 
            * adj_list (Dict[int, List[int]]): Neighbors for each block.
            * edge_w   (Dict[int, Dict[int, int]]): Edge weights (u->v).
    """
    if aggregate not in {"sum", "max", "min"}:
        raise ValueError("aggregate must be one of {'sum','max','min'}")

    pair_weight: Dict[Tuple[int, int], int] = {}

    for m in face_matches:
        i = int(m["block1"]["block_index"])
        j = int(m["block2"]["block_index"])
        if ignore_self_matches and i == j:
            continue

        IMIN = int(m["block1"]["IMIN"]); JMIN = int(m["block1"]["JMIN"]); KMIN = int(m["block1"]["KMIN"])
        IMAX = int(m["block1"]["IMAX"]); JMAX = int(m["block1"]["JMAX"]); KMAX = int(m["block1"]["KMAX"])
        dI = max(abs(IMAX - IMIN), 1)
        dJ = max(abs(JMAX - JMIN), 1)
        dK = max(abs(KMAX - KMIN), 1)
        w = dI * dJ * dK  # Edge weight = communication cost (number of face nodes)

        a, b = (i, j) if i < j else (j, i)
        if (a, b) not in pair_weight:
            pair_weight[(a, b)] = w
        else:
            if aggregate == "sum":
                pair_weight[(a, b)] += w
            elif aggregate == "max":
                pair_weight[(a, b)] = max(pair_weight[(a, b)], w)
            else:  # "min"
                pair_weight[(a, b)] = min(pair_weight[(a, b)], w)

    # adjacency list and edge weights
    adj_list: Dict[int, List[int]] = {u: [] for u in range(n_blocks)}
    edge_w: Dict[int, Dict[int, int]] = {u: {} for u in range(n_blocks)}

    for (a, b), w in pair_weight.items():
        adj_list[a].append(b)
        adj_list[b].append(a)
        edge_w[a][b] = w
        edge_w[b][a] = w

    for u in range(n_blocks):
        adj_list[u] = sorted(set(adj_list[u]))

    return adj_list, edge_w


def csr_from_adj_and_weights(
    adj_list: Dict[int, List[int]],
    edge_w: Dict[int, Dict[int, int]],
) -> Tuple[List[int], List[int], List[int]]:
    """
    Build CSR arrays (xadj, adjncy, eweights) from adjacency + weights.

    Returns:
        xadj (List[int]): prefix sum of neighbors
        adjncy (List[int]): flattened neighbor list
        eweights (List[int]): edge weights aligned with adjncy
    """
    xadj: List[int] = [0]
    adjncy: List[int] = []
    eweights: List[int] = []
    count = 0
    for u in sorted(adj_list.keys()):
        for v in adj_list[u]:
            adjncy.append(v)
            eweights.append(edge_w[u].get(v, 1))
            count += 1
        xadj.append(count)
    return xadj, adjncy, eweights


# ---------------------------------------------------------------------------
# pymetis compatibility wrapper
# ---------------------------------------------------------------------------
def _metis_part_graph_compat(
    nparts: int,
    xadj: List[int],
    adjncy: List[int],
    vwgt: Optional[List[int]] = None,
    eweights: Optional[List[int]] = None,
):
    """
    Call pymetis.part_graph in a way that's compatible with multiple pymetis versions.
    Tries keyword names first, then falls back to positional args.

    Args:
        nparts (int): number of partitions
        xadj (List[int]): CSR row pointer
        adjncy (List[int]): CSR neighbor list
        vwgt (Optional[List[int]]): vertex weights
        eweights (Optional[List[int]]): edge weights

    Returns:
        (edgecut, parts)
    """
    if not HAS_PYMETIS:
        raise RuntimeError(
            "pymetis is not available. On Windows it is skipped during install; "
            "use Linux/macOS to enable METIS-based partitioning."
        )

    sig_params = set(inspect.signature(pymetis.part_graph).parameters.keys())  # type: ignore[attr-defined]

    # Prefer keyword args when supported
    if {"xadj", "adjncy"}.issubset(sig_params):
        kwargs = {"xadj": xadj, "adjncy": adjncy}
        # Vertex weights
        if "vwgt" in sig_params and vwgt is not None:
            kwargs["vwgt"] = vwgt
        elif "vweights" in sig_params and vwgt is not None:
            kwargs["vweights"] = vwgt  # alternate name
        # Edge weights
        if "eweights" in sig_params and eweights is not None:
            kwargs["eweights"] = eweights
        elif "adjwgt" in sig_params and eweights is not None:
            kwargs["adjwgt"] = eweights  # alternate name

        return pymetis.part_graph(nparts, **kwargs) # type: ignore

    # Fallback: positional signature (older builds)
    return pymetis.part_graph(nparts, xadj, adjncy, vwgt, None, eweights) # type: ignore


# ---------------------------------------------------------------------------
# Partition entrypoint (takes face_matches directly)
# ---------------------------------------------------------------------------
def partition_from_face_matches(
    face_matches: List[dict],
    block_sizes: List[int],
    nparts: int,
    favor_blocksize: bool = True,
    aggregate: str = "sum",
    ignore_self_matches: bool = True,
) -> Tuple[List[int], Dict[int, List[int]], Dict[int, Dict[int, int]]]:
    """
    Partition a graph derived from `face_matches` using pymetis.
    
    Note:
        Aggregate controls how we combine weights when the same two blocks are connected by multiple faces.
    
    Why this matters:
        connectivity_fast can return more than one face between the same pair of blocks 
        (e.g. a block is split and has two non-contiguous interfaces with its neighbor).
        Each face gives a weight = dI * dJ * dK (number of shared nodes).
        METIS expects one edge per block-pair, with a single weight.
        So if there are multiple faces, we must decide how to merge them. That’s what aggregate does.

    Returns
    -------
    parts : List[int]
        Partition id (0-based) for each block.
    adj_list : Dict[int, List[int]]
        Adjacency list used for the partitioning.
    edge_w : Dict[int, Dict[int, int]]
        Edge weights.
    """
    if not HAS_PYMETIS:
        raise RuntimeError(
            "METIS partitioning is disabled because pymetis is unavailable. "
            "Install pymetis (Linux/macOS) or run on a platform where it is supported."
        )

    n_blocks = len(block_sizes)
    adj_list, edge_w = build_weighted_graph_from_face_matches(
        face_matches, n_blocks,
        aggregate=aggregate,
        ignore_self_matches=ignore_self_matches,
    )
    xadj, adjncy, eweights = csr_from_adj_and_weights(adj_list, edge_w)

    vwgt: Optional[List[int]] = None
    if favor_blocksize:
        vwgt = block_sizes

    _edgecut, parts = _metis_part_graph_compat(
        nparts=nparts,
        xadj=xadj,
        adjncy=adjncy,
        vwgt=vwgt,
        eweights=eweights,
    )
    return parts, adj_list, edge_w


# ---------------------------------------------------------------------------
# DDCMP writer
# ---------------------------------------------------------------------------
def write_ddcmp(
    parts: Sequence[int],
    blocksizes: List[int],
    adj_list: Dict[int, List[int]],
    edge_weights: Optional[Dict[int, Dict[int, int]]] = None,
    filename: str = "ddcmp.dat",
) -> None:
    """
    Writes ddcmp.dat and ddcmp_info.txt.

    Notes:
        * parts are 0-based in memory, but written 1-based in the file (to match your C#).
        * edge_weights affects the per-partition 'edge_work' if provided.
    """
    n_proc = (max(parts) + 1) if parts else 0
    n_isp = n_proc
    n_blocks = len(parts)

    out = Path(filename)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write(f"{n_proc}\n{n_isp}\n{n_blocks}\n")
        for b_idx in range(n_blocks):
            f.write(f"{b_idx + 1} {parts[b_idx] + 1}\n")
        for isp in range(n_isp):
            f.write(f"{isp + 1} {isp}\n")

    communication_work = [0] * n_proc
    partition_edge_weights = [0] * n_proc
    volume_nodes = [0] * n_proc

    for b, bsize in enumerate(blocksizes):
        pid = parts[b]
        volume_nodes[pid] += bsize

    ew = edge_weights or {}
    for b in range(len(blocksizes)):
        pid = parts[b]
        for nbr in adj_list.get(b, []):
            nbr_pid = parts[nbr]
            if nbr_pid != pid:
                communication_work[pid] += 1
                partition_edge_weights[pid] += int(ew.get(b, {}).get(nbr, 1))

    info_path = out.parent / "ddcmp_info.txt"
    with info_path.open("w", encoding="utf-8") as f:
        for i in range(n_proc):
            block_count = sum(1 for p in parts if p == i)
            f.write(f"Parition {i:d} has {block_count} blocks\n")
        f.write(f"Number of partitions/processors {n_proc}\n")
        for i in range(n_proc):
            f.write(
                f"Parition or processor {i:d} has communication work {communication_work[i]:d} "
                f"edge_work {partition_edge_weights[i]:d} volume_nodes {volume_nodes[i]:d}\n"
            )
        f.write(
            f"Total communication work {sum(communication_work):d} "
            f"Total edge_work {sum(partition_edge_weights):d}\n"
        )
