"""Graph construction and METIS-based partitioning utilities for Plot3D meshes.

This module converts face-match connectivity data (output of ``connectivity_fast``)
into a weighted block-adjacency graph, partitions that graph with METIS via
``pymetis``, and writes the resulting partition assignment to a DDCMP input file.

Note:
    ``pymetis`` is an optional dependency.  On Windows it is skipped at install
    time; METIS-based partitioning is therefore only available on Linux/macOS.
"""
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
    """Convert ``connectivity_fast`` face matches into a weighted block-adjacency graph.

    ``connectivity_fast`` can return more than one face between the same pair of
    blocks (e.g. when a block has two non-contiguous interfaces with its
    neighbor).  Each face contributes a weight equal to ``dI * dJ * dK`` (the
    number of shared nodes on that face).  Because METIS expects exactly one
    edge per block-pair, multiple faces between the same pair are merged
    according to *aggregate*.

    Args:
        face_matches: Output from ``connectivity_fast``.  Each element is a
            dict containing ``block1`` and ``block2`` sub-dicts with keys
            ``block_index``, ``IMIN``, ``JMIN``, ``KMIN``, ``IMAX``, ``JMAX``,
            and ``KMAX``.
        n_blocks: Total number of blocks in the mesh.
        aggregate: Strategy for merging multiple face weights between the same
            block pair.  Must be one of ``’sum’``, ``’max’``, or ``’min’``.
            Defaults to ``’sum’``.
        ignore_self_matches: If ``True``, face matches where both sides belong
            to the same block are discarded.  Defaults to ``True``.

    Returns:
        A two-element tuple ``(adj_list, edge_w)`` where:

        - **adj_list** (``Dict[int, List[int]]``): Sorted neighbor list for
          each block index.
        - **edge_w** (``Dict[int, Dict[int, int]]``): Mapping
          ``edge_w[u][v]`` -> integer edge weight for the undirected edge
          between blocks *u* and *v*.

    Raises:
        ValueError: If *aggregate* is not one of ``’sum’``, ``’max’``,
            or ``’min’``.
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
    """Build CSR arrays from an adjacency list and edge-weight mapping.

    Produces the three arrays required by METIS/pymetis in Compressed Sparse
    Row (CSR) format.

    Args:
        adj_list: Neighbor list for each block index, keyed by block index.
            Neighbors must be sorted (as returned by
            :func:`build_weighted_graph_from_face_matches`).
        edge_w: Edge weights, where ``edge_w[u][v]`` is the integer weight of
            the undirected edge between blocks *u* and *v*.

    Returns:
        A three-element tuple ``(xadj, adjncy, eweights)`` where:

        - **xadj** (``List[int]``): Row-pointer array of length
          ``n_blocks + 1``; ``xadj[u]`` is the index in *adjncy* where
          block *u*'s neighbors begin.
        - **adjncy** (``List[int]``): Flattened neighbor indices.
        - **eweights** (``List[int]``): Edge weights aligned element-wise
          with *adjncy*.
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
    """Call ``pymetis.part_graph`` compatibly across multiple pymetis versions.

    Different builds of pymetis use different parameter names (e.g. ``vwgt``
    vs. ``vweights``, ``eweights`` vs. ``adjwgt``).  This wrapper inspects the
    live signature and uses keyword arguments when available, falling back to
    the older positional calling convention otherwise.

    Args:
        nparts: Number of desired partitions.
        xadj: CSR row-pointer array (see :func:`csr_from_adj_and_weights`).
        adjncy: CSR neighbor-index array.
        vwgt: Per-vertex (block) weights, or ``None`` to weight all vertices
            equally.
        eweights: Per-edge weights aligned with *adjncy*, or ``None`` to
            weight all edges equally.

    Returns:
        A two-element tuple ``(edgecut, parts)`` as returned by
        ``pymetis.part_graph``, where *edgecut* is the total cut weight and
        *parts* is a list of 0-based partition IDs, one per block.

    Raises:
        RuntimeError: If ``pymetis`` is not installed or is unavailable on
            the current platform.
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
    """Partition a block-connectivity graph into *nparts* parts using METIS.

    Builds a weighted adjacency graph from *face_matches* (see
    :func:`build_weighted_graph_from_face_matches`), converts it to CSR format,
    and calls METIS via :func:`_metis_part_graph_compat`.

    When *favor_blocksize* is ``True``, each block’s vertex weight is set to
    its node count so that METIS balances computational work rather than simply
    counting blocks.

    Args:
        face_matches: Output from ``connectivity_fast``.  Each element is a
            dict with ``block1`` and ``block2`` sub-dicts describing the shared
            face extents and block indices.
        block_sizes: Number of nodes in each block (used as vertex weights when
            *favor_blocksize* is ``True``).
        nparts: Desired number of partitions.
        favor_blocksize: If ``True``, use block node counts as METIS vertex
            weights so that partitions balance computational load.
            Defaults to ``True``.
        aggregate: Strategy for merging multiple face weights between the same
            block pair.  One of ``’sum’``, ``’max’``, or ``’min’``.
            Defaults to ``’sum’``.
        ignore_self_matches: If ``True``, self-matching faces (same block on
            both sides) are discarded before building the graph.
            Defaults to ``True``.

    Returns:
        A three-element tuple ``(parts, adj_list, edge_w)`` where:

        - **parts** (``List[int]``): 0-based partition ID for each block.
        - **adj_list** (``Dict[int, List[int]]``): Adjacency list used for
          the partitioning.
        - **edge_w** (``Dict[int, Dict[int, int]]``): Edge weights used for
          the partitioning.

    Raises:
        RuntimeError: If ``pymetis`` is not installed or is unavailable on
            the current platform.
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
    """Write DDCMP partition files from a METIS partition assignment.

    Produces two output files in the directory of *filename*:

    - **ddcmp.dat** -- main partition assignment file consumed by downstream
      C# / Fortran solvers.  Block and processor indices are written
      **1-based** even though *parts* is 0-based in memory.
    - **ddcmp_info.txt** -- human-readable summary reporting block counts,
      communication work, edge work, and volume-node counts per partition.

    Args:
        parts: 0-based partition ID for each block, as returned by
            :func:`partition_from_face_matches`.
        blocksizes: Number of nodes in each block; used to compute per-partition
            volume-node totals.
        adj_list: Block-adjacency list (keyed by block index) used to identify
            cross-partition edges.
        edge_weights: Optional edge-weight mapping ``edge_weights[u][v]``
            used to compute per-partition edge work.  If ``None``, every
            cross-partition edge counts as weight 1.
        filename: Path to the primary output file.  The info file is written
            to the same directory with the name ``ddcmp_info.txt``.
            Defaults to ``'ddcmp.dat'``.

    Returns:
        None
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
