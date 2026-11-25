"""Structured-grid mesh-quality battery.

Mirrors plot3d-rs ``src/mesh_quality.rs`` (commits ``9ebe7d3`` — the
DEGENERATE-vs-INVERTED cell split; ``ed01fd0`` — the element-type census;
``a6c9d9c`` — ``cell_aspect_ratio`` returns ``inf`` for a zero-length edge
instead of a huge finite artefact; and the threshold-tightening commits
``b6036d2``/``94c4def``). The Rust module is itself a port of the
CFD-readiness checks in ``tgs-py-grc/python/tgs_py/quality/``
(``metrics_3d.py``, ``checks_3d.py``, ``thresholds.py``); the additions on
top of that reference are block handedness (left-handed / negative-Jacobian
detection + fix) and per-cell negative-volume detection.

`untangle_folded_lines` / `untangle_mesh` repair folded grid lines (a node
placed beyond its successor, reversing interior faces) without moving any
block-face node; the ``plot3d-untangle <mesh>`` command runs it on a file.

Departure from the Rust source, deliberate: the per-cell metric kernels
here (``cell_signed_volume``, ``cell_volume_divergence``, ``cell_aspect_ratio``,
``cell_skewness``, ...) are full-field numpy broadcasts over the whole
``(NI,NJ,NK)`` cell grid at once, not per-``(i,j,k)`` scalar functions as
in Rust — this module is meant to run on multi-million-cell meshes, where a
per-cell Python loop (like ``Block.cell_volumes()`` in ``block.py``, which
despite writing into numpy arrays is a plain triple-nested Python loop) would
be far too slow. The reported violation taxonomy, severities, and threshold
values are kept identical to Rust so results are directly comparable.

This module is deliberately kept separate from ``blockfunctions.make_right_handed``
(a differently-scoped, multi-block function that also remaps connectivity
indices) and from ``glennht/validation.py`` (GlennHT-specific BC/connectivity
checks that are informational about handedness rather than diagnostic about
it) — see ``make_block_right_handed`` and each of those modules' own
docstrings for the cross-reference.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

from .block import Block

# =============================================================================
# vec3 helpers (operate on (..., 3) arrays, broadcasting over the leading
# cell-field axes)
# =============================================================================


def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.cross(a, b, axis=-1)


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sum(a * b, axis=-1)


def _norm(a: np.ndarray) -> np.ndarray:
    return np.sqrt(_dot(a, a))


# =============================================================================
# Per-cell metric kernels — full-field numpy broadcasts, port of
# tgs-py-grc `metrics_3d.py` via plot3d-rs `mesh_quality.rs`
# =============================================================================


def cell_dims(block: Block) -> Tuple[int, int, int]:
    """Cell counts ``(NI, NJ, NK)`` for a block.

    Any of them is ``0`` when the corresponding node dimension is ``< 2``
    (the block has no cells in that direction).
    """
    return (
        max(block.IMAX - 1, 0),
        max(block.JMAX - 1, 0),
        max(block.KMAX - 1, 0),
    )


# The 8 corners of cell (i,j,k), in the same (di,dj,dk) order as Rust's
# `cell_volume_divergence`/`cell_has_collapsed_edge`/`cell_distinct_node_count`
# build their `n[0..7]` arrays: di fastest, then dj, then dk. This exact
# order matters for `cell_volume_divergence`'s FACES indices below; the
# other kernels only need it for consistency/reuse.
_CORNER_OFFSETS = [
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
    (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1),
]


def _corner_list(block: Block) -> List[np.ndarray]:
    """The 8 corner-coordinate fields of every cell, each shaped
    ``(NI, NJ, NK, 3)``, built via sliced views (no per-cell Python loop).
    Zero-size axes (a block with < 2 nodes on some axis) flow through as
    zero-size array dimensions automatically.
    """
    NI, NJ, NK = cell_dims(block)
    X, Y, Z = block.X, block.Y, block.Z
    out = []
    for di, dj, dk in _CORNER_OFFSETS:
        sub_x = X[di : di + NI, dj : dj + NJ, dk : dk + NK]
        sub_y = Y[di : di + NI, dj : dj + NJ, dk : dk + NK]
        sub_z = Z[di : di + NI, dj : dj + NJ, dk : dk + NK]
        out.append(np.stack([sub_x, sub_y, sub_z], axis=-1))
    return out


def cell_signed_volume(block: Block) -> np.ndarray:
    """Signed cell volume ``e_i . (e_j x e_k)`` (scalar triple product,
    anchored at corner ``(i,j,k)``) — the handedness / Jacobian indicator:
    ``> 0`` right-handed, ``< 0`` left-handed, ``~= 0`` degenerate.
    """
    n = _corner_list(block)
    p0 = n[0]
    ei, ej, ek = n[1] - p0, n[2] - p0, n[4] - p0
    return _dot(ei, _cross(ej, ek))


# Each face as 4 corner indices into `_corner_list`'s output, wound so the
# area vector points OUT. Mirrors Rust's `FACES` constant exactly.
_FACES = (
    (0, 4, 6, 2),  # i-low
    (1, 3, 7, 5),  # i-high
    (0, 1, 5, 4),  # j-low
    (2, 6, 7, 3),  # j-high
    (0, 2, 3, 1),  # k-low
    (4, 5, 7, 6),  # k-high
)


def cell_volume_divergence(block: Block) -> np.ndarray:
    """True hexahedral cell volume by the divergence theorem,
    ``V = (1/3) sum_faces (r_c . S)``, each quad face's area vector taken
    as ``0.5 (d1 x d2)`` from its diagonals. This is the operator the
    solver actually integrates with, so it is the authority on whether a
    cell has usable volume — contrast `cell_signed_volume`, which anchors
    on one corner and reads exactly 0 on a cell with a collapsed edge at
    that corner even though the cell has real, positive volume.
    """
    n = _corner_list(block)
    v = np.zeros(n[0].shape[:3])
    for ia, ib, ic, idd in _FACES:
        a, b, c, d = n[ia], n[ib], n[ic], n[idd]
        centroid = (a + b + c + d) * 0.25
        s = _cross(c - a, d - b)
        v = v + _dot(centroid, 0.5 * s)
    return v / 3.0


def cell_has_collapsed_edge(block: Block) -> npt.NDArray[np.bool_]:
    """Whether each cell has at least one pair of coincident corner nodes
    — i.e. it sits on a collapsed/pinched grid line (a deliberate, standard
    turbomachinery construction, e.g. an O-grid folded onto a blade-tip
    camber line). Exact (bit-identical) equality — NO tolerance: a
    tolerance would wrongly sweep in merely-close-but-distinct nodes.
    """
    n = _corner_list(block)
    shape3 = n[0].shape[:3]
    collapsed = np.zeros(shape3, dtype=bool)
    for a in range(8):
        for c in range(a + 1, 8):
            collapsed |= np.all(n[a] == n[c], axis=-1)
    return collapsed


def cell_aspect_ratio(block: Block) -> np.ndarray:
    """Per-cell aspect ratio ``max(edge_len) / min(edge_len)`` over the
    three edge vectors. Always ``>= 1``; a perfect cube returns ``1.0``.

    A cell with a zero-length edge has NO aspect ratio, and this returns
    ``inf`` for it rather than a finite number (port of the ``a6c9d9c``
    fix). That cell sits on a collapsed grid line — see
    `cell_has_collapsed_edge` — and the condition belongs to the
    degenerate-cell check, which already reports it by name. Flooring the
    denominator instead (the pre-fix behavior) does not avoid the problem,
    it disguises it: a legitimate short edge over a collapsed one produces
    an astronomically large but precise-looking finite number. ``inf``
    still exceeds any threshold, so nothing is silenced.
    """
    n = _corner_list(block)
    p0 = n[0]
    ei, ej, ek = n[1] - p0, n[2] - p0, n[4] - p0
    li, lj, lk = _norm(ei), _norm(ej), _norm(ek)
    mx = np.maximum(np.maximum(li, lj), lk)
    mn = np.minimum(np.minimum(li, lj), lk)
    with np.errstate(divide="ignore", invalid="ignore"):
        ar = np.where(mn == 0.0, np.inf, mx / mn)
    return ar


def _edge_angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Acute angle in degrees between two edge-vector fields, in [0, 90]."""
    denom = np.maximum(_norm(a) * _norm(b), 1e-30)
    cos_t = np.clip(np.abs(_dot(a, b) / denom), 0.0, 1.0)
    return np.degrees(np.arccos(cos_t))


def cell_skewness(block: Block) -> np.ndarray:
    """Per-cell equiangle skewness in degrees: ``90 - min(orthogonality)``
    over the three edge-vector pairs. ``0`` = perfectly orthogonal, ``90``
    = degenerate.
    """
    n = _corner_list(block)
    p0 = n[0]
    ei, ej, ek = n[1] - p0, n[2] - p0, n[4] - p0
    ortho_min = np.minimum(
        np.minimum(_edge_angle_deg(ei, ej), _edge_angle_deg(ej, ek)),
        _edge_angle_deg(ei, ek),
    )
    return 90.0 - ortho_min


def cell_distinct_node_count(block: Block) -> npt.NDArray[np.int64]:
    """Number of DISTINCT corner nodes (of the 8) per cell, via the same
    exact-equality rule as `cell_has_collapsed_edge` (O-grid pinch nodes
    are bit-identical, so a tolerance would sweep in merely-close nodes).
    """
    n = _corner_list(block)
    shape3 = n[0].shape[:3]
    count = np.zeros(shape3, dtype=np.int64)
    for a in range(8):
        dup = np.zeros(shape3, dtype=bool)
        for c in range(a):
            dup |= np.all(n[a] == n[c], axis=-1)
        count += np.where(dup, 0, 1)
    return count


def cell_centroid(block: Block) -> np.ndarray:
    """Centroid of every cell — the mean of its 8 corner nodes. Shape
    ``(NI, NJ, NK, 3)``.
    """
    n = _corner_list(block)
    total = n[0]
    for arr in n[1:]:
        total = total + arr
    return total / 8.0


# =============================================================================
# Block handedness — left-handed / negative-Jacobian detection + fix
# =============================================================================


class Handedness(Enum):
    """Handedness of a block, from the sign of its median signed cell volume."""

    #: Median signed volume > 0 — the normal, FV-ready orientation.
    RIGHT_HANDED = "RIGHT_HANDED"
    #: Median signed volume < 0 — the block's (i,j,k) indexing is
    #: mirror-flipped; `make_block_right_handed` fixes it by flipping one axis.
    LEFT_HANDED = "LEFT_HANDED"
    #: Median signed volume ~= 0, or the block has no cells — flipping
    #: cannot help; the per-cell negative/degenerate-volume checks flag it.
    DEGENERATE = "DEGENERATE"


def block_handedness(block: Block) -> Handedness:
    """Classify a block's handedness from the median signed cell volume
    (median is robust to a handful of locally-bad cells). A left-handed
    block has essentially every cell at negative signed volume.
    """
    sv = cell_signed_volume(block).ravel()
    if sv.size == 0:
        return Handedness.DEGENERATE
    median = float(np.median(sv))
    # Scale the "~= 0" tolerance by the typical cell size so it works for
    # meshes in any units.
    typical = max(abs(float(sv.max())), abs(float(sv.min())))
    eps = typical * 1e-12 + np.finfo(np.float64).tiny
    if median > eps:
        return Handedness.RIGHT_HANDED
    elif median < -eps:
        return Handedness.LEFT_HANDED
    else:
        return Handedness.DEGENERATE


def make_block_right_handed(block: Block) -> Tuple[Block, Optional[int]]:
    """If `block` is left-handed, return a right-handed copy (one
    structured axis reversed) and the flipped axis index; otherwise return
    a copy and ``None``. Non-mutating: the input block is never modified.

    Flipping one axis is a reflection — it negates every cell's signed
    volume, so a left-handed block becomes right-handed. The geometry
    (cells, volumes, faces) is unchanged, only the (i,j,k) traversal
    direction, so this is physics-neutral relabeling.

    A ``DEGENERATE`` block is returned unchanged (``None``) — flipping
    cannot help; `run_all`'s negative/degenerate-volume checks report it.

    Named `make_block_right_handed`, not `make_right_handed`, to avoid
    colliding with `blockfunctions.make_right_handed` — a differently-scoped
    function that fixes a *list* of blocks in place, on a fixed axis, and
    also remaps connectivity/periodic/outer-face index records. The two are
    deliberately kept separate; this one is a single-block, non-mutating,
    axis-auto-selecting primitive with no connectivity awareness.
    """
    if block_handedness(block) != Handedness.LEFT_HANDED:
        return Block(block.X.copy(), block.Y.copy(), block.Z.copy()), None

    # Flip the first structured axis that actually has cells (i preferred).
    if block.IMAX > 1:
        axis = 0
    elif block.JMAX > 1:
        axis = 1
    else:
        axis = 2

    X = np.flip(block.X, axis=axis).copy()
    Y = np.flip(block.Y, axis=axis).copy()
    Z = np.flip(block.Z, axis=axis).copy()
    return Block(X, Y, Z), axis


# =============================================================================
# Folded grid lines — geometry-only detection and repair
# =============================================================================
#
# A folded grid line is one whose edge between two consecutive nodes points
# opposite to the same edge on the neighbouring parallel lines: a node has
# been placed beyond its successor. The cells either side keep positive
# volume, but the face between them sees its owner and neighbour centres on
# the wrong sides (a.d <= 0), which a finite-volume solver must refuse. The
# repair rebuilds the folded edge from its neighbours, moving only
# block-interior nodes so block surfaces and inter-block point matches stay
# bit-exact. No boundary-condition knowledge is used. Mirrors plot3d-rs
# ``mesh_quality::untangle_folded_lines``.

_AXIS_NAMES = ("i", "j", "k")


@dataclass(frozen=True)
class ReversedFace:
    """An interior face whose area vector ``a`` does not point from its owner
    cell toward its neighbour cell: ``cos(a, d) <= 0``, ``d`` the
    owner-to-neighbour centre vector. ``a`` is ``0.5 (d1 x d2)`` from the
    face diagonals; centres are 8-node means.
    """

    block: int
    #: Face family: 0 = i-face, 1 = j-face, 2 = k-face.
    axis: int
    #: Owner cell ``(i, j, k)``; the neighbour is the next cell along ``axis``.
    owner: Tuple[int, int, int]
    #: ``cos(a, d)``, relative to the block's handedness.
    cos: float


@dataclass(frozen=True)
class NonPositiveCell:
    """A cell whose signed divergence-theorem volume, relative to the block's
    handedness, is not ``> 0``."""

    block: int
    cell: Tuple[int, int, int]
    volume: float


@dataclass(frozen=True)
class FoldedEdge:
    """Edge ``start -> start + 1`` along ``axis`` whose vector ``e`` opposes
    the mean ``m`` of the same edge on the unfolded neighbouring parallel
    lines: ``e.m <= 0``."""

    block: int
    axis: int
    start: Tuple[int, int, int]
    #: ``e.m``.
    dot: float


class RepairConstruction(Enum):
    """How a relocated node ``n`` on a folded edge ``n``-``a`` was placed
    (``a`` the fixed anchor)."""

    #: ``X(n) = X(a) - u``, ``u`` the mean ``n -> a`` edge over the healthy
    #: neighbouring lines.
    NEIGHBOUR_MEAN = "NeighbourMean"
    #: Used when NEIGHBOUR_MEAN would fold the edge on the other side of
    #: ``n``: ``X(n) = X(f) + r (X(a) - X(f))``, ``f`` the next node along
    #: the line, ``r = |w| / (|w| + |u|)``, ``w`` the healthy mean
    #: ``f -> n`` edge.
    CHORD_REDISTRIBUTION = "ChordRedistribution"


@dataclass(frozen=True)
class MovedNode:
    """One node relocated by `untangle_folded_lines`."""

    block: int
    node: Tuple[int, int, int]
    old: Tuple[float, float, float]
    new: Tuple[float, float, float]
    displacement: float
    construction: RepairConstruction


#: Entries printed per list by `UntangleReport.format_report` and
#: `UntangleError` before the rest are only counted.
UNTANGLE_REPORT_LIMIT = 1000


def _fmt_ijk(n) -> str:
    return f"({n[0]},{n[1]},{n[2]})"


def _fmt_xyz(p) -> str:
    return f"({p[0]:.15e}, {p[1]:.15e}, {p[2]:.15e})"


def _step(n, axis):
    m = list(n)
    m[axis] += 1
    return tuple(m)


def _edge_of_cell(axis, start, c) -> bool:
    return start[axis] == c[axis] and all(
        start[q] in (c[q], c[q] + 1) for q in range(3) if q != axis)


class _Bounded:
    """Appends at most `UNTANGLE_REPORT_LIMIT` lines of one kind, counting the rest."""

    def __init__(self):
        self.shown = 0
        self.hidden = 0

    def push(self, out: List[str], line) -> None:
        if self.shown < UNTANGLE_REPORT_LIMIT:
            out.append(line())
            self.shown += 1
        else:
            self.hidden += 1

    def tail(self, out: List[str], what: str) -> None:
        if self.hidden:
            out.append(f"  ... {self.hidden} more {what} not shown "
                       f"(first {UNTANGLE_REPORT_LIMIT} listed)")


@dataclass
class UntangleReport:
    """What `untangle_folded_lines` / `untangle_mesh` found and changed, per
    block. `is_empty` for a mesh that needed no repair. Inside an
    `UntangleError` it holds what was found before the refusal, and the
    moves listed were planned, not applied."""

    blocks_scanned: int = 0
    reversed_faces: List[ReversedFace] = field(default_factory=list)
    nonpositive_cells: List[NonPositiveCell] = field(default_factory=list)
    folded_edges: List[FoldedEdge] = field(default_factory=list)
    moved: List[MovedNode] = field(default_factory=list)
    max_displacement: float = 0.0
    #: Blocks whose repaired geometry was re-checked.
    rechecked_blocks: List[int] = field(default_factory=list)
    #: Offenders left after the repair (non-empty only in a refusal).
    remaining_reversed_faces: List[ReversedFace] = field(default_factory=list)
    remaining_nonpositive_cells: List[NonPositiveCell] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.reversed_faces or self.nonpositive_cells or self.folded_edges
                    or self.moved or self.remaining_reversed_faces
                    or self.remaining_nonpositive_cells)

    def extend(self, other: "UntangleReport") -> None:
        self.blocks_scanned += other.blocks_scanned
        self.reversed_faces.extend(other.reversed_faces)
        self.nonpositive_cells.extend(other.nonpositive_cells)
        self.folded_edges.extend(other.folded_edges)
        self.moved.extend(other.moved)
        self.max_displacement = max(self.max_displacement, other.max_displacement)
        self.rechecked_blocks.extend(other.rechecked_blocks)
        self.remaining_reversed_faces.extend(other.remaining_reversed_faces)
        self.remaining_nonpositive_cells.extend(other.remaining_nonpositive_cells)

    def summary(self) -> str:
        """The summary line."""
        n_blocks = len({f.block for f in self.reversed_faces}
                       | {c.block for c in self.nonpositive_cells})
        return (
            f"untangle folded grid lines: scanned {self.blocks_scanned} block(s); "
            f"{len(self.reversed_faces)} reversed interior face(s) and "
            f"{len(self.nonpositive_cells)} non-positive cell(s) in {n_blocks} block(s), "
            f"localised to {len(self.folded_edges)} folded edge(s); "
            f"{len(self.moved)} node(s) moved, max displacement {self.max_displacement:e}"
        )

    def format_report(self) -> str:
        """Plain-text report: the `summary` line, then per block every
        reversed face (owner -> neighbour cell, direction, cos(a,d), and the
        folded edge and node it was localised to), every non-positive cell,
        every moved node (old -> new, displacement, construction) and the
        re-check result. Each list prints at most `UNTANGLE_REPORT_LIMIT`
        entries and counts the rest. Same layout as plot3d-rs."""
        out = [self.summary()]
        faces, cells, moves, left = _Bounded(), _Bounded(), _Bounded(), _Bounded()
        blocks = sorted({f.block for f in self.reversed_faces}
                        | {c.block for c in self.nonpositive_cells}
                        | {e.block for e in self.folded_edges}
                        | {m.block for m in self.moved}
                        | {f.block for f in self.remaining_reversed_faces}
                        | {c.block for c in self.remaining_nonpositive_cells})
        for b in blocks:
            bf = [f for f in self.reversed_faces if f.block == b]
            bc = [c for c in self.nonpositive_cells if c.block == b]
            be = [e for e in self.folded_edges if e.block == b]
            bm = [m for m in self.moved if m.block == b]
            out.append("")
            out.append(f"block {b}: {len(bf)} reversed face(s), {len(bc)} non-positive "
                       f"cell(s), {len(be)} folded edge(s), {len(bm)} node(s) moved")

            def localised(cs, be=be, bm=bm):
                hits = []
                for e in be:
                    if not any(_edge_of_cell(e.axis, e.start, c) for c in cs):
                        continue
                    end = _step(e.start, e.axis)
                    node = next((m.node for m in bm if m.node in (e.start, end)), None)
                    moved = f", node {_fmt_ijk(node)} moved" if node is not None else ""
                    hits.append(f"folded {_AXIS_NAMES[e.axis]}-edge "
                                f"{_fmt_ijk(e.start)}-{_fmt_ijk(end)}{moved}")
                return "; ".join(hits) if hits else "not localised to a folded edge"

            for f in bf:
                nb = _step(f.owner, f.axis)
                faces.push(out, lambda f=f, nb=nb: (
                    f"  reversed {_AXIS_NAMES[f.axis]}-face: cell {_fmt_ijk(f.owner)} -> "
                    f"{_fmt_ijk(nb)}, cos(a,d) {f.cos:.4e}; {localised([f.owner, nb])}"))
            for c in bc:
                cells.push(out, lambda c=c: (
                    f"  non-positive cell {_fmt_ijk(c.cell)}: volume {c.volume:e}; "
                    f"{localised([c.cell])}"))
            for m in bm:
                moves.push(out, lambda m=m: (
                    f"  moved node {_fmt_ijk(m.node)}: {_fmt_xyz(m.old)} -> {_fmt_xyz(m.new)}, "
                    f"displacement {m.displacement:e} [{m.construction.value}]"))
            rf = [f for f in self.remaining_reversed_faces if f.block == b]
            rc = [c for c in self.remaining_nonpositive_cells if c.block == b]
            if b in self.rechecked_blocks:
                out.append(f"  re-check: {len(rf)} reversed face(s), {len(rc)} non-positive "
                           "cell(s)")
                for f in rf:
                    left.push(out, lambda f=f: (
                        f"    still reversed {_AXIS_NAMES[f.axis]}-face: cell "
                        f"{_fmt_ijk(f.owner)} -> {_fmt_ijk(_step(f.owner, f.axis))}, "
                        f"cos(a,d) {f.cos:.4e}"))
                for c in rc:
                    left.push(out, lambda c=c: (
                        f"    still non-positive cell {_fmt_ijk(c.cell)}: volume {c.volume:e}"))
            elif bf or bc:
                out.append("  re-check: not run (refused before repair)")
        faces.tail(out, "reversed face(s)")
        cells.tail(out, "non-positive cell(s)")
        moves.tail(out, "moved node(s)")
        left.tail(out, "remaining offender(s)")
        return "\n".join(out) + "\n"

    def __str__(self) -> str:
        return self.format_report()


@dataclass(frozen=True)
class Refusal:
    """One reason a block was refused, located by block and node (or cell)
    where possible."""

    block: int
    node: Optional[Tuple[int, int, int]]
    reason: str


class UntangleError(ValueError):
    """`untangle_folded_lines` / `untangle_mesh` refused; nothing was
    changed. ``refusals`` carries every refusal found, not only the first;
    ``report`` what was detected (moves in it planned, not applied)."""

    def __init__(self, refusals: List[Refusal], report: UntangleReport):
        self.refusals = list(refusals)
        self.report = report
        super().__init__(self.format())

    def format(self) -> str:
        out = [
            f"untangle folded grid lines REFUSED at {len(self.refusals)} location(s); "
            "nothing was changed (moves below were planned, not applied)",
            self.report.format_report().rstrip("\n"),
            "",
            "refusals:",
        ]
        shown = _Bounded()
        for r in self.refusals:
            shown.push(out, lambda r=r: (
                f"  block {r.block} {_fmt_ijk(r.node)}: {r.reason}" if r.node is not None
                else f"  block {r.block}: {r.reason}"))
        shown.tail(out, "refusal(s)")
        return "\n".join(out) + "\n"

    def __str__(self) -> str:
        return self.format()


def _orientation_sign(block: Block) -> float:
    """``-1`` for a left-handed block, else ``+1``: face area vectors and
    divergence volumes flip sign with handedness, so both checks judge the
    block as if `make_block_right_handed` had been applied."""
    return -1.0 if block_handedness(block) == Handedness.LEFT_HANDED else 1.0


def _nodes(block: Block) -> np.ndarray:
    return np.stack([block.X, block.Y, block.Z], axis=-1).astype(np.float64)


def _face_area_fields(P: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Area vectors ``0.5 (d1 x d2)`` of every i-, j- and k-face, with the
    corner order of plot3d-rs ``metrics::compute_face_metrics`` (pointing
    toward +i/+j/+k on a right-handed block). Shapes ``(I, NJ, NK, 3)``,
    ``(NI, J, NK, 3)``, ``(NI, NJ, K, 3)``."""
    # i-face: p0=(i,j,k) p1=(i,j+1,k) p2=(i,j+1,k+1) p3=(i,j,k+1)
    si = 0.5 * _cross(P[:, 1:, 1:] - P[:, :-1, :-1], P[:, :-1, 1:] - P[:, 1:, :-1])
    # j-face: p0=(i,j,k) p1=(i,j,k+1) p2=(i+1,j,k+1) p3=(i+1,j,k)
    sj = 0.5 * _cross(P[1:, :, 1:] - P[:-1, :, :-1], P[1:, :, :-1] - P[:-1, :, 1:])
    # k-face: p0=(i,j,k) p1=(i+1,j,k) p2=(i+1,j+1,k) p3=(i,j+1,k)
    sk = 0.5 * _cross(P[1:, 1:, :] - P[:-1, :-1, :], P[:-1, 1:, :] - P[1:, :-1, :])
    return si, sj, sk


def find_inverted_geometry(
    block: Block, block_index: int = 0
) -> Tuple[List[ReversedFace], List[NonPositiveCell]]:
    """Every reversed interior face and non-positive-volume cell of one block.

    A face is reversed when ``cos(a, d) <= 0`` (see `ReversedFace`); a face
    whose cosine is undefined (zero area, coincident centres, non-finite) is
    not counted. A cell is non-positive when its signed `cell_volume_divergence`
    is not ``> 0``. Both are relative to the block's handedness. A block with
    fewer than two nodes on any axis has no cells and reports nothing.
    """
    if min(block.IMAX, block.JMAX, block.KMAX) < 2:
        return [], []
    s = _orientation_sign(block)
    P = _nodes(block)
    vol = s * cell_volume_divergence(block)
    cells = [
        NonPositiveCell(block_index, (int(i), int(j), int(k)), float(vol[i, j, k]))
        for i, j, k in zip(*np.nonzero(~(vol > 0.0)))
    ]
    # Cell order in the report follows plot3d-rs: k outermost, i fastest.
    cells.sort(key=lambda c: (c.cell[2], c.cell[1], c.cell[0]))

    centres = cell_centroid(block)
    si, sj, sk = _face_area_fields(P)
    faces: List[ReversedFace] = []
    for axis, area in enumerate((si, sj, sk)):
        # Interior faces of this family: node planes 1..N-2 along `axis`.
        sl = [slice(None)] * 3
        sl[axis] = slice(1, -1)
        a = area[tuple(sl)]
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis] = slice(0, -1)
        hi[axis] = slice(1, None)
        d = centres[tuple(hi)] - centres[tuple(lo)]
        with np.errstate(divide="ignore", invalid="ignore"):
            nrm = np.sqrt(_dot(a, a) * _dot(d, d))
            cos = s * _dot(a, d) / nrm
        bad = (nrm > 0) & np.isfinite(cos) & (cos <= 0.0)
        for i, j, k in zip(*np.nonzero(bad)):
            faces.append(ReversedFace(block_index, axis, (int(i), int(j), int(k)),
                                      float(cos[i, j, k])))
    faces.sort(key=lambda f: (f.axis, f.owner[2], f.owner[1], f.owner[0]))
    return faces, cells


def _edge(P: np.ndarray, n: Tuple[int, int, int], axis: int) -> np.ndarray:
    m = list(n)
    m[axis] += 1
    return P[tuple(m)] - P[tuple(n)]


def _parallel_starts(dims, n, axis) -> List[Tuple[int, int, int]]:
    out = []
    for q in range(3):
        if q == axis:
            continue
        if n[q] > 0:
            m = list(n)
            m[q] -= 1
            out.append(tuple(m))
        if n[q] + 1 < dims[q]:
            m = list(n)
            m[q] += 1
            out.append(tuple(m))
    return out


def _healthy_parallel_mean(P, dims, folded, n, axis) -> Optional[np.ndarray]:
    """Mean of the same edge over the neighbouring parallel lines not in
    ``folded`` with non-zero length; ``None`` when there is none."""
    acc = np.zeros(3)
    count = 0
    for s in _parallel_starts(dims, n, axis):
        v = _edge(P, s, axis)
        if (axis, s) not in folded and float(v @ v) > 0.0:
            acc = acc + v
            count += 1
    return acc * (1.0 / count) if count else None


def _edge_dot_healthy(P, dims, folded, n, axis) -> Optional[float]:
    e = _edge(P, n, axis)
    if float(e @ e) == 0.0:
        return None
    m = _healthy_parallel_mean(P, dims, folded, n, axis)
    return None if m is None else float(e @ m)


def _classify_folded_edges(P: np.ndarray) -> set:
    """The folded edges ``(axis, (i, j, k))`` of a block. A first pass flags
    every edge with ``e.m <= 0`` against the mean of all its neighbouring
    parallel lines; flagged edges are then re-judged against the mean of
    their unflagged neighbours only, until nothing changes (a large fold
    drags the all-neighbour mean of the healthy lines next to it). An edge
    with no unflagged neighbour stays flagged. The set only shrinks."""
    dims = P.shape[:3]
    folded = set()
    for axis in range(3):
        E = np.diff(P, axis=axis)
        S = np.zeros_like(E)
        C = np.zeros(E.shape[:3])
        for q in range(3):
            if q == axis:
                continue
            dst = [slice(None)] * 3
            src = [slice(None)] * 3
            dst[q], src[q] = slice(1, None), slice(0, -1)
            S[tuple(dst)] += E[tuple(src)]
            C[tuple(dst)] += 1
            dst[q], src[q] = slice(0, -1), slice(1, None)
            S[tuple(dst)] += E[tuple(src)]
            C[tuple(dst)] += 1
        with np.errstate(divide="ignore", invalid="ignore"):
            m = S * (1.0 / C)[..., None]
            dv = _dot(E, m)
        flag = (_dot(E, E) > 0) & (C > 0) & (_dot(m, m) > 0) & np.isfinite(dv) & (dv <= 0)
        for i, j, k in zip(*np.nonzero(flag)):
            folded.add((axis, (int(i), int(j), int(k))))
    while True:
        keep = set()
        for axis, n in folded:
            v = _edge_dot_healthy(P, dims, folded, n, axis)
            if v is None or not (v > 0.0):
                keep.add((axis, n))
        if len(keep) == len(folded):
            return keep
        folded = keep


def _on_block_face(dims, n) -> bool:
    return any(n[q] == 0 or n[q] + 1 == dims[q] for q in range(3))


def _other_incident_edge_folded(dims, folded, n, exclude) -> bool:
    for axis in range(3):
        starts = []
        if n[axis] + 1 < dims[axis]:
            starts.append(tuple(n))
        if n[axis] > 0:
            m = list(n)
            m[axis] -= 1
            starts.append(tuple(m))
        for s in starts:
            if (axis, s) != exclude and (axis, s) in folded:
                return True
    return False


def _place_folded_node(P, dims, folded, moved, anchor, axis, m):
    """Position of ``moved`` on the folded edge to the fixed ``anchor``
    (``m`` the healthy mean of that edge, oriented start -> start+1).
    Returns ``(position, construction, nodes it is built from)`` or raises
    ``ValueError`` with the refusal reason. Mirrors plot3d-rs
    ``place_folded_node``."""
    up = anchor[axis] > moved[axis]
    u = m if up else -m  # healthy n -> a edge
    xa = P[anchor]
    p1 = xa - u
    if up:
        far = None
        if moved[axis] > 0:
            f = list(moved)
            f[axis] -= 1
            far = tuple(f)
    else:
        far = None
        if moved[axis] + 1 < dims[axis]:
            f = list(moved)
            f[axis] += 1
            far = tuple(f)
    if far is None:
        return p1, RepairConstruction.NEIGHBOUR_MEAN, [anchor]
    w_start = far if up else moved
    mw = _healthy_parallel_mean(P, dims, folded, w_start, axis)
    if mw is None:
        return p1, RepairConstruction.NEIGHBOUR_MEAN, [anchor]
    w = mw if up else -mw  # healthy f -> n edge
    xf = P[far]
    if float((p1 - xf) @ w) > 0.0:
        return p1, RepairConstruction.NEIGHBOUR_MEAN, [anchor]
    c = xa - xf
    if not (float(c @ w) > 0.0 and float(c @ u) > 0.0):
        raise ValueError(
            f"anchor ({anchor[0]},{anchor[1]},{anchor[2]}) does not lie beyond node "
            f"({far[0]},{far[1]},{far[2]}) along the line; rebuilding the two edges "
            "between them cannot unfold it"
        )
    lw, lu = float(np.sqrt(w @ w)), float(np.sqrt(u @ u))
    r = lw / (lw + lu)
    return xf + r * c, RepairConstruction.CHORD_REDISTRIBUTION, [anchor, far]


def untangle_folded_lines(block: Block, block_index: int = 0) -> Tuple[Block, UntangleReport]:
    """Repair folded grid lines in one block, using geometry only.

    1. **Detect** reversed interior faces and non-positive cells
       (`find_inverted_geometry`). None found: an exact copy of the block is
       returned with an empty report.
    2. **Localise** to the edges of the affected cells whose vector opposes
       the mean of the same edge on the unfolded neighbouring parallel lines.
    3. **Repair** each folded edge by relocating one endpoint and keeping the
       other (the anchor) fixed. A node on any block face (an index at 0 or
       max) is never moved, so surfaces and inter-block point matches stay
       exact. With one endpoint on a block face the face node is the anchor;
       with both interior, the endpoint whose other incident edges are not
       folded. The relocated node takes the healthy neighbouring lines' mean
       edge (`RepairConstruction.NEIGHBOUR_MEAN`), or, where that would fold
       the next edge along the line, sits on the chord to the next node in
       the neighbours' edge-length ratio
       (`RepairConstruction.CHORD_REDISTRIBUTION`).
    4. **Re-check**: no reversed face and no non-positive cell may remain.

    Raises `UntangleError` (no partial fix, naming every block and node):
    geometry not localisable to a folded edge; a folded edge with both
    endpoints on block faces; an ambiguous anchor; no healthy neighbouring
    line; an anchor that does not lie beyond the next node along the line;
    two different relocations of one node, or one built from a moved node;
    any inverted geometry left after the repair.

    Mirrors plot3d-rs ``mesh_quality::untangle_folded_lines``.
    """
    faces, cells = find_inverted_geometry(block, block_index)
    if not faces and not cells:
        return (Block(block.X.copy(), block.Y.copy(), block.Z.copy()),
                UntangleReport(blocks_scanned=1))

    P = _nodes(block)
    dims = P.shape[:3]
    candidates = set()

    def add_cell_edges(c):
        for axis in range(3):
            p, q = (axis + 1) % 3, (axis + 2) % 3
            for dp in (0, 1):
                for dq in (0, 1):
                    s = list(c)
                    s[p] += dp
                    s[q] += dq
                    candidates.add((axis, tuple(s)))

    for f in faces:
        add_cell_edges(f.owner)
        nb = list(f.owner)
        nb[f.axis] += 1
        add_cell_edges(tuple(nb))
    for c in cells:
        add_cell_edges(c.cell)

    folded = _classify_folded_edges(P)
    # Same order as plot3d-rs's BTreeSet<(axis, [i, j, k])>.
    folded_edges = []
    for axis, start in sorted(candidates & folded):
        dv = _edge_dot_healthy(P, dims, folded, start, axis)
        if dv is None:
            dv = 0.0
        folded_edges.append(FoldedEdge(block_index, axis, start, dv))
    report = UntangleReport(blocks_scanned=1, reversed_faces=faces,
                            nonpositive_cells=cells, folded_edges=folded_edges)

    refusals: List[Refusal] = []
    if not folded_edges:
        for f in faces:
            refusals.append(Refusal(
                block_index, f.owner,
                f"reversed {_AXIS_NAMES[f.axis]}-face (owner cell, cos {f.cos:.4e}) is not "
                "localisable to a folded grid-line edge"))
        for c in cells:
            refusals.append(Refusal(
                block_index, c.cell,
                f"non-positive cell (volume {c.volume:e}) is not localisable to a folded "
                "grid-line edge"))
        raise UntangleError(refusals, report)

    moves = {}  # node -> (position, construction, built-from nodes)
    for fe in folded_edges:
        axis = fe.axis
        n0 = fe.start
        n1 = list(n0)
        n1[axis] += 1
        n1 = tuple(n1)
        m = _healthy_parallel_mean(P, dims, folded, n0, axis)
        if m is None:
            refusals.append(Refusal(
                block_index, n0,
                f"folded {_AXIS_NAMES[axis]}-edge has no healthy neighbouring line to "
                "rebuild it from"))
            continue
        f0, f1 = _on_block_face(dims, n0), _on_block_face(dims, n1)
        if f0 and f1:
            refusals.append(Refusal(
                block_index, n0,
                f"folded {_AXIS_NAMES[axis]}-edge has both endpoints on block faces; "
                "block-face nodes are never moved"))
            continue
        if f0:
            moved, anchor = n1, n0
        elif f1:
            moved, anchor = n0, n1
        else:
            key = (axis, n0)
            o0 = _other_incident_edge_folded(dims, folded, n0, key)
            o1 = _other_incident_edge_folded(dims, folded, n1, key)
            if o0 and not o1:
                moved, anchor = n0, n1
            elif o1 and not o0:
                moved, anchor = n1, n0
            else:
                refusals.append(Refusal(
                    block_index, n0,
                    f"folded {_AXIS_NAMES[axis]}-edge with both endpoints interior has no "
                    "unambiguous anchor (other folded edges at both or neither endpoint)"))
                continue
        try:
            plan = _place_folded_node(P, dims, folded, moved, anchor, axis, m)
        except ValueError as exc:
            refusals.append(Refusal(block_index, moved, str(exc)))
            continue
        if moved in moves:
            if not np.array_equal(moves[moved][0], plan[0]):
                refusals.append(Refusal(
                    block_index, moved,
                    "two folded edges ask for different relocations of this node"))
        else:
            moves[moved] = plan
    for node, (_, _, built_from) in sorted(moves.items()):
        for f in built_from:
            if f in moves:
                refusals.append(Refusal(
                    block_index, node,
                    f"relocation is built from node ({f[0]},{f[1]},{f[2]}), which is "
                    "itself relocated"))
    for node, (new, construction, _) in sorted(moves.items()):
        old = P[node]
        disp = float(np.sqrt(np.sum((new - old) ** 2)))
        report.max_displacement = max(report.max_displacement, disp)
        report.moved.append(MovedNode(
            block_index, node, tuple(float(v) for v in old), tuple(float(v) for v in new),
            disp, construction))
    if refusals:
        raise UntangleError(refusals, report)

    X, Y, Z = block.X.copy(), block.Y.copy(), block.Z.copy()
    for m in report.moved:
        X[m.node], Y[m.node], Z[m.node] = m.new
    repaired = Block(X, Y, Z)

    faces_after, cells_after = find_inverted_geometry(repaired, block_index)
    report.rechecked_blocks.append(block_index)
    if not faces_after and not cells_after:
        return repaired, report
    for f in faces_after:
        refusals.append(Refusal(
            block_index, f.owner,
            f"{_AXIS_NAMES[f.axis]}-face cell {_fmt_ijk(f.owner)} -> "
            f"{_fmt_ijk(_step(f.owner, f.axis))} still reversed after repair, "
            f"cos(a,d) {f.cos:.4e}"))
    for c in cells_after:
        refusals.append(Refusal(
            block_index, c.cell,
            f"cell still non-positive after repair, volume {c.volume:e}"))
    report.remaining_reversed_faces = faces_after
    report.remaining_nonpositive_cells = cells_after
    raise UntangleError(refusals, report)


def untangle_mesh(blocks: Sequence[Block]) -> Tuple[List[Block], UntangleReport]:
    """`untangle_folded_lines` over every block (block index = position).
    All-or-nothing: if any block is refused, `UntangleError` carries every
    block's refusals and the report of every block (moves planned, not
    applied), and no repaired mesh is returned."""
    out: List[Block] = []
    report = UntangleReport()
    refusals: List[Refusal] = []
    for bi, b in enumerate(blocks):
        try:
            rb, rep = untangle_folded_lines(b, bi)
        except UntangleError as exc:
            refusals.extend(exc.refusals)
            report.extend(exc.report)
            continue
        out.append(rb)
        report.extend(rep)
    if refusals:
        raise UntangleError(refusals, report)
    return out, report


def untangle_plot3d_file(filename: str, out_filename: Optional[str] = None,
                         report_filename: Optional[str] = None) -> UntangleReport:
    """Read a Plot3D file, `untangle_mesh` it, and write
    ``<stem>_untangled<ext>`` plus ``<stem>_untangled_report.txt`` next to
    it (or to the given paths), in the input's detected format. Refuses to
    overwrite an existing file. Raises `UntangleError` if the mesh is
    refused; only the report (with the refusal detail) is written then.
    """
    import os.path as osp
    from .read import read_plot3D, _detect_plot3d_format
    from .write import write_plot3D

    stem, ext = osp.splitext(filename)
    out_filename = out_filename or f"{stem}_untangled{ext}"
    report_filename = report_filename or f"{stem}_untangled_report.txt"
    for p in (out_filename, report_filename):
        if osp.exists(p):
            raise FileExistsError(f"refusing to overwrite existing {p}")

    fmt = _detect_plot3d_format(filename)
    blocks = read_plot3D(filename, binary=fmt["binary"], big_endian=fmt["big_endian"],
                         read_double=fmt["read_double"], fortran=fmt["fortran"])
    try:
        repaired, report = untangle_mesh(blocks)
    except UntangleError as exc:
        with open(report_filename, "w") as f:
            f.write(f"input:  {filename}\noutput: none (refused)\n\n{exc}")
        raise
    write_plot3D(out_filename, repaired, binary=fmt["binary"], big_endian=fmt["big_endian"],
                 double_precision=fmt["read_double"], fortran=fmt["fortran"])
    with open(report_filename, "w") as f:
        f.write(f"input:  {filename}\noutput: {out_filename}\n\n")
        f.write(report.format_report())
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point of the ``plot3d-untangle <mesh>`` command: untangles folded
    grid lines in a Plot3D file (see `untangle_plot3d_file`) and prints the
    report. Returns ``1`` on a refusal or an existing output file."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="plot3d-untangle",
        description="Untangle folded grid lines in a Plot3D mesh; writes "
                    "<stem>_untangled<ext> and <stem>_untangled_report.txt.")
    parser.add_argument("mesh", help="input Plot3D file (format auto-detected)")
    parser.add_argument("-o", "--output", default=None, help="output mesh path")
    args = parser.parse_args(argv)
    try:
        report = untangle_plot3d_file(args.mesh, args.output)
    except (UntangleError, FileExistsError) as exc:
        print(exc)
        return 1
    print(report.format_report(), end="")
    return 0


# =============================================================================
# Thresholds — port of tgs-py-grc `thresholds.py` (core fields only)
# =============================================================================


@dataclass(frozen=True)
class Thresholds:
    """CFD-readiness thresholds for the core quality checks. Three named
    presets — `Thresholds.STRICT`, `Thresholds.STANDARD`, `Thresholds.RELAXED`
    — assigned onto the class below (a frozen dataclass cannot construct
    class-level instances of itself inside its own body).
    """

    #: Per-cell skewness (deg): the 99th-percentile cell limit.
    skew_p99_deg: float
    #: Per-cell skewness (deg): the absolute-worst cell limit.
    skew_max_deg: float
    #: Minimum interior-edge orthogonality angle (deg); 90 deg = perfect.
    min_orthogonality_deg: float
    #: Maximum aspect ratio for interior cells.
    max_ar_interior: float
    #: Maximum aspect ratio for wall first-cells (BL legitimately high).
    max_ar_wall: float
    #: Minimum cell volume relative to the block median.
    min_cell_volume_ratio: float
    #: Cell layers dropped from each axis-endpoint before computing the
    #: skewness percentiles (excludes wall first-cells from the stats).
    boundary_drop: int = 2

    @staticmethod
    def from_preset_name(name: str) -> "Thresholds":
        """Resolve a case-insensitive preset name (``"strict"`` /
        ``"standard"`` / ``"relaxed"``). Unknown names fall back to
        ``STANDARD``.
        """
        key = (name or "").upper()
        if key == "STRICT":
            return Thresholds.STRICT
        if key == "RELAXED":
            return Thresholds.RELAXED
        return Thresholds.STANDARD


Thresholds.STRICT = Thresholds(
    skew_p99_deg=50.0,
    skew_max_deg=60.0,
    min_orthogonality_deg=30.0,
    max_ar_interior=100_000.0,
    max_ar_wall=200_000.0,
    min_cell_volume_ratio=1e-5,
    boundary_drop=2,
)
Thresholds.STANDARD = Thresholds(
    skew_p99_deg=80.0,
    skew_max_deg=75.0,
    min_orthogonality_deg=15.0,
    max_ar_interior=5_000.0,
    max_ar_wall=100_000.0,
    min_cell_volume_ratio=1e-6,
    boundary_drop=2,
)
Thresholds.RELAXED = Thresholds(
    skew_p99_deg=87.0,
    skew_max_deg=90.0,
    min_orthogonality_deg=1.0,
    max_ar_interior=30_000.0,
    max_ar_wall=math.inf,
    min_cell_volume_ratio=1e-9,
    boundary_drop=2,
)


# =============================================================================
# Violation model — port of tgs-py-grc `Violation`
# =============================================================================


class Severity(Enum):
    """Severity of a quality violation."""

    #: Advisory — the solver runs, accuracy may suffer.
    WARN = "WARN"
    #: Fatal — breaks the finite-volume discretization (negative /
    #: degenerate cell volume).
    ERROR = "ERROR"


@dataclass(frozen=True)
class CellLocation:
    """Where a violating cell lives — for pointing the user straight at it."""

    block: int
    i: int
    j: int
    k: int
    centroid: Tuple[float, float, float]


@dataclass
class Violation:
    """A single mesh-quality threshold breach."""

    #: Short check name — "negative_volume", "degenerate_cell", "skewness",
    #: "aspect_ratio", "orthogonality", "min_cell_volume".
    check: str
    severity: Severity
    #: The observed value that breached the threshold.
    actual: float
    #: The threshold it breached.
    threshold: float
    #: The worst cell's location, when the check is per-cell.
    location: Optional[CellLocation]
    #: Human-readable one-liner.
    message: str


# =============================================================================
# Rust-order tie-breaking helpers (module-private)
# =============================================================================
#
# Rust's per-cell scan is i-fastest: `(k * ncj + j) * nci + i`. NumPy's
# natural C-order flattening of a (NI,NJ,NK) array is the OPPOSITE — k
# fastest. On an exact tie (which symmetric test fixtures, e.g. a unit
# cube, WILL produce), a naive `np.argmax`/`np.argmin` picks a different
# cell than Rust's scan would, even though the reported value/threshold
# match. Transposing to (NK,NJ,NI) before raveling makes NumPy's C-order
# flatten iterate i fastest, matching Rust's scan order exactly, so ties
# resolve to the same first-encountered cell.


def _rust_order_argmax(
    field: np.ndarray, mask: Optional[npt.NDArray[np.bool_]] = None
) -> Tuple[float, Tuple[int, int, int]]:
    """Max of a (NI,NJ,NK) field and its (i,j,k), breaking ties the way
    Rust's i-fastest scan (`if v > max { ... }`) does. `mask` excludes
    cells from consideration (as if they didn't exist) without needing a
    separate compaction step.
    """
    effective = field if mask is None else np.where(mask, field, -np.inf)
    ft = np.transpose(effective, (2, 1, 0))
    flat_idx = int(np.argmax(ft))
    k, j, i = np.unravel_index(flat_idx, ft.shape)
    return float(field[i, j, k]), (int(i), int(j), int(k))


def _rust_order_argmin(
    field: np.ndarray, mask: Optional[npt.NDArray[np.bool_]] = None
) -> Tuple[float, Tuple[int, int, int]]:
    """Min of a (NI,NJ,NK) field and its (i,j,k), breaking ties the way
    Rust's i-fastest scan (`if v < min { ... }`) does.
    """
    effective = field if mask is None else np.where(mask, field, np.inf)
    ft = np.transpose(effective, (2, 1, 0))
    flat_idx = int(np.argmin(ft))
    k, j, i = np.unravel_index(flat_idx, ft.shape)
    return float(field[i, j, k]), (int(i), int(j), int(k))


# =============================================================================
# Report
# =============================================================================


@dataclass
class MeshQualityReport:
    """The full mesh-quality report for a multi-block mesh."""

    #: Threshold preset name used ("STANDARD" etc.).
    preset: str
    #: Per-block handedness (index b = block b).
    handedness: List[Handedness]
    #: Every threshold breach found, across all blocks.
    violations: List[Violation]

    def n_error(self) -> int:
        """Number of Error-severity violations (negative/degenerate volume)."""
        return sum(1 for v in self.violations if v.severity == Severity.ERROR)

    def n_warn(self) -> int:
        """Number of Warn-severity violations (skewness/AR/orthogonality)."""
        return sum(1 for v in self.violations if v.severity == Severity.WARN)

    def passes(self) -> bool:
        """True when there are no Error violations — the mesh is
        FV-discretizable (warnings allowed).
        """
        return self.n_error() == 0

    def format_report(self) -> str:
        """Human-readable multi-line report — handedness summary, the
        violation counts, then each violation with its location.
        """
        lines = [
            f"Mesh quality ({self.preset} preset): {len(self.handedness)} block(s), "
            f"{self.n_error()} error(s), {self.n_warn()} warning(s)"
        ]
        left = [b for b, h in enumerate(self.handedness) if h == Handedness.LEFT_HANDED]
        degen = [b for b, h in enumerate(self.handedness) if h == Handedness.DEGENERATE]
        if left:
            lines.append(f"  left-handed blocks (need flipping): {left}")
        if degen:
            lines.append(f"  degenerate blocks: {degen}")
        if not left and not degen:
            lines.append("  all blocks right-handed")
        for v in self.violations:
            sev = "ERROR" if v.severity == Severity.ERROR else "warn "
            if v.location is not None:
                loc = v.location
                lines.append(
                    f"  [{sev}] {v.check} — block {loc.block} cell "
                    f"({loc.i},{loc.j},{loc.k}) at "
                    f"({loc.centroid[0]:.4f},{loc.centroid[1]:.4f},{loc.centroid[2]:.4f}): "
                    f"{v.message}"
                )
            else:
                lines.append(f"  [{sev}] {v.check}: {v.message}")
        return "\n".join(lines) + "\n"


# =============================================================================
# Checks — port of tgs-py-grc `checks_3d.py` (core checks) + negative volume
# =============================================================================

#: Max per-cell negative-volume violations listed individually before the
#: report collapses the rest into a "+N more" line.
_MAX_LISTED_NEGATIVE = 64


def run_all(
    blocks: Sequence[Block],
    thresholds: Optional[Thresholds] = None,
    preset_name: str = "STANDARD",
) -> MeshQualityReport:
    """Run the full quality battery on a multi-block mesh.

    Per block (blocks with fewer than 2 nodes on any axis are skipped —
    they have no cells): handedness classification, then the per-cell
    checks — negative signed volume vs. degenerate collapsed-line volume
    (the `9ebe7d3` split), min-cell-volume ratio, equiangle skewness (p99 +
    max, boundary-cropped), minimum orthogonality, and aspect ratio
    (interior vs. wall/overall max — with non-finite values from collapsed
    cells excluded from BOTH before either max is computed, the `a6c9d9c`
    fix). Negative volume -> Error; everything else -> Warn.

    Unlike Rust's `run_all` (which always requires both a `Thresholds` and
    a preset name), `thresholds=None` here resolves via
    `Thresholds.from_preset_name(preset_name)` — a deliberate,
    zero-required-Thresholds-arg convenience.

    Does not mutate the blocks and never raises for a bad mesh — it
    returns the report and the caller decides what is fatal. Apply
    `make_block_right_handed` to fix left-handed blocks first so the
    handedness column reads clean.
    """
    if thresholds is None:
        thresholds = Thresholds.from_preset_name(preset_name)
    t = thresholds

    handedness: List[Handedness] = []
    violations: List[Violation] = []

    for bi, b in enumerate(blocks):
        handedness.append(block_handedness(b))

        NI, NJ, NK = cell_dims(b)
        if NI == 0 or NJ == 0 or NK == 0:
            continue  # no cells — nothing to score

        signed = cell_signed_volume(b)
        true_vol = cell_volume_divergence(b)
        collapsed = cell_has_collapsed_edge(b)
        centroid = cell_centroid(b)

        # --- negative signed volume vs. degenerate collapsed line ---
        #
        # Two DIFFERENT defects hide behind "signed volume <= 0":
        #   INVERTED   — genuinely turned inside out; divergence-theorem
        #                volume is negative. Fatal.
        #   DEGENERATE — sits on a collapsed/pinched grid line; the
        #                corner-anchored triple product is exactly 0 but
        #                the cell has real positive volume. Not fatal.
        neg_or_zero = signed <= 0.0
        degenerate_mask = neg_or_zero & (true_vol > 0.0) & collapsed
        inverted_mask = neg_or_zero & ~degenerate_mask

        neg_count = int(np.count_nonzero(inverted_mask))
        degen_count = int(np.count_nonzero(degenerate_mask))

        if neg_count > 0:
            # argwhere on the (NK,NJ,NI)-transposed mask yields (k,j,i)
            # rows in the same i-fastest order Rust's scan uses.
            coords = np.argwhere(np.transpose(inverted_mask, (2, 1, 0)))
            n_listed = min(neg_count, _MAX_LISTED_NEGATIVE)
            for row in coords[:n_listed]:
                k, j, i = int(row[0]), int(row[1]), int(row[2])
                sv, vd = float(signed[i, j, k]), float(true_vol[i, j, k])
                violations.append(
                    Violation(
                        check="negative_volume",
                        severity=Severity.ERROR,
                        actual=sv,
                        threshold=0.0,
                        location=CellLocation(
                            bi, i, j, k, tuple(centroid[i, j, k].tolist())
                        ),
                        message=(
                            f"signed cell volume {sv:.3e} <= 0 and divergence-theorem "
                            f"volume {vd:.3e} (inverted cell)"
                        ),
                    )
                )
            if neg_count > _MAX_LISTED_NEGATIVE:
                violations.append(
                    Violation(
                        check="negative_volume",
                        severity=Severity.ERROR,
                        actual=float(neg_count),
                        threshold=0.0,
                        location=None,
                        message=(
                            f"block {bi}: {neg_count} cells with non-positive signed "
                            f"volume ({neg_count - _MAX_LISTED_NEGATIVE} more not "
                            f"listed individually)"
                        ),
                    )
                )

        if degen_count > 0:
            coords = np.argwhere(np.transpose(degenerate_mask, (2, 1, 0)))
            k0, j0, i0 = int(coords[0][0]), int(coords[0][1]), int(coords[0][2])
            vd0 = float(true_vol[i0, j0, k0])
            violations.append(
                Violation(
                    check="degenerate_cell",
                    severity=Severity.WARN,
                    actual=float(degen_count),
                    threshold=0.0,
                    location=CellLocation(
                        bi, i0, j0, k0, tuple(centroid[i0, j0, k0].tolist())
                    ),
                    message=(
                        f"block {bi}: {degen_count} cell(s) on a COLLAPSED GRID LINE "
                        f"(coincident corner nodes). These are not inverted — the "
                        f"divergence-theorem volume is positive (e.g. {vd0:.3e} at the "
                        f"first such cell) — so the discretization is well posed. "
                        f"But such cells are orders of magnitude smaller than their "
                        f"neighbours and will throttle an explicit local time step; "
                        f"merge them into a neighbour rather than advancing them as "
                        f"independent control volumes."
                    ),
                )
            )

        # --- degenerate: min |volume| relative to the block median ---
        all_degenerate_are_collapsed = degen_count > 0 and neg_count == 0
        absvol = np.abs(true_vol)
        median = float(np.median(absvol))
        if median <= 0.0:
            violations.append(
                Violation(
                    check="min_cell_volume",
                    severity=Severity.ERROR,
                    actual=0.0,
                    threshold=1.0,
                    location=None,
                    message=(
                        f"block {bi}: median cell volume is non-positive — "
                        f"block is degenerate"
                    ),
                )
            )
        else:
            vmin, (i, j, k) = _rust_order_argmin(absvol)
            ratio = vmin / median
            if ratio < t.min_cell_volume_ratio:
                sev = (
                    Severity.WARN
                    if (all_degenerate_are_collapsed and bool(collapsed[i, j, k]))
                    else Severity.ERROR
                )
                violations.append(
                    Violation(
                        check="min_cell_volume",
                        severity=sev,
                        actual=ratio,
                        threshold=t.min_cell_volume_ratio,
                        location=CellLocation(
                            bi, i, j, k, tuple(centroid[i, j, k].tolist())
                        ),
                        message=(
                            f"min cell volume / median = {ratio:.2e} < "
                            f"{t.min_cell_volume_ratio:.0e} (near-degenerate)"
                        ),
                    )
                )

        # --- skewness (p99 + max, boundary-cropped), orthogonality (global) ---
        skew = cell_skewness(b)
        drop = t.boundary_drop
        can_crop = drop > 0 and NI > 2 * drop and NJ > 2 * drop and NK > 2 * drop
        if can_crop:
            i0c, i1c, j0c, j1c, k0c, k1c = drop, NI - drop, drop, NJ - drop, drop, NK - drop
        else:
            i0c, i1c, j0c, j1c, k0c, k1c = 0, NI, 0, NJ, 0, NK
        cropped = skew[i0c:i1c, j0c:j1c, k0c:k1c]
        if cropped.size == 0:
            # Cropping emptied everything — fall back to the full field.
            cropped = skew
            i0c = j0c = k0c = 0
        skew_p99 = float(np.percentile(cropped, 99))
        skew_max, (li, lj, lk) = _rust_order_argmax(cropped)
        smi, smj, smk = li + i0c, lj + j0c, lk + k0c

        if skew_p99 > t.skew_p99_deg:
            violations.append(
                Violation(
                    check="skewness",
                    severity=Severity.WARN,
                    actual=skew_p99,
                    threshold=t.skew_p99_deg,
                    location=None,
                    message=f"skewness p99 = {skew_p99:.1f}° > {t.skew_p99_deg:.1f}°",
                )
            )
        if skew_max > t.skew_max_deg:
            violations.append(
                Violation(
                    check="skewness",
                    severity=Severity.WARN,
                    actual=skew_max,
                    threshold=t.skew_max_deg,
                    location=CellLocation(
                        bi, smi, smj, smk, tuple(centroid[smi, smj, smk].tolist())
                    ),
                    message=f"skewness {skew_max:.1f}° > {t.skew_max_deg:.1f}°",
                )
            )

        global_skew_max = float(skew.max())
        ortho_min = 90.0 - global_skew_max
        if ortho_min < t.min_orthogonality_deg:
            violations.append(
                Violation(
                    check="orthogonality",
                    severity=Severity.WARN,
                    actual=ortho_min,
                    threshold=t.min_orthogonality_deg,
                    location=None,
                    message=(
                        f"min orthogonality angle {ortho_min:.1f}° < "
                        f"{t.min_orthogonality_deg:.1f}°"
                    ),
                )
            )

        # --- aspect ratio: interior vs. overall(="wall") max ---
        #
        # The finite-value mask is built ONCE, up front, and used for BOTH
        # the wall/overall max and the interior max below — mirroring
        # exactly where Rust's `a6c9d9c` fix places its
        # `if !v.is_finite() { continue; }`, so a collapsed cell's `inf`
        # aspect ratio cannot silently set either max.
        ar = cell_aspect_ratio(b)
        finite_mask = np.isfinite(ar)

        interior_ok = NI > 2 and NJ > 2 and NK > 2
        if interior_ok:
            interior_selector = np.zeros((NI, NJ, NK), dtype=bool)
            interior_selector[1:-1, 1:-1, 1:-1] = True
        else:
            interior_selector = np.ones((NI, NJ, NK), dtype=bool)

        wall_mask = finite_mask
        interior_mask = finite_mask & interior_selector

        if np.any(wall_mask):
            wall_max, (wi, wj, wk) = _rust_order_argmax(ar, mask=wall_mask)
        else:
            wall_max, (wi, wj, wk) = float("-inf"), (0, 0, 0)
        if np.any(interior_mask):
            interior_max, (ii, ij, ik) = _rust_order_argmax(ar, mask=interior_mask)
        else:
            interior_max, (ii, ij, ik) = float("-inf"), (0, 0, 0)

        if interior_max > t.max_ar_interior:
            violations.append(
                Violation(
                    check="aspect_ratio",
                    severity=Severity.WARN,
                    actual=interior_max,
                    threshold=t.max_ar_interior,
                    location=CellLocation(
                        bi, ii, ij, ik, tuple(centroid[ii, ij, ik].tolist())
                    ),
                    message=(
                        f"interior aspect ratio {interior_max:.0f} > "
                        f"{t.max_ar_interior:.0f}"
                    ),
                )
            )
        if wall_max > t.max_ar_wall:
            violations.append(
                Violation(
                    check="aspect_ratio",
                    severity=Severity.WARN,
                    actual=wall_max,
                    threshold=t.max_ar_wall,
                    location=CellLocation(
                        bi, wi, wj, wk, tuple(centroid[wi, wj, wk].tolist())
                    ),
                    message=f"wall aspect ratio {wall_max:.0f} > {t.max_ar_wall:.0f}",
                )
            )

    return MeshQualityReport(preset=preset_name, handedness=handedness, violations=violations)


# =============================================================================
# Element-type inventory — ADS / Code Leo `ELEMTYPE` + `ADVOLAREA` analogue
# =============================================================================
#
# A structured hex cell that sits on an O-grid pinch line collapses to a
# lower element (e.g. a wedge = PRISM). This block reproduces a Code
# Leo-style ELEMTYPE/ADVOLAREA census so a load echoes the reference
# solver's own printout. Diagnostic only: it changes no solver state.


class ElementType(Enum):
    """The standard element a structured hex cell collapses to, by
    distinct corner-node count. Mirrors Code Leo's ELEMTYPE categories.
    """

    #: 8 distinct nodes — a normal hexahedron.
    HEX = "HEX"
    #: 6 — one edge collapsed (a wedge; two coincident node-pairs).
    PRISM = "PRISM"
    #: 5 — one face collapsed toward a point.
    PYRAMID = "PYRAMID"
    #: 4 — a tetrahedron.
    TET = "TET"
    #: 7, or < 4 — a non-standard partial collapse.
    OTHER = "OTHER"

    @staticmethod
    def from_distinct_nodes(n: int) -> "ElementType":
        """Classify by the number of distinct corner nodes."""
        return {
            8: ElementType.HEX,
            6: ElementType.PRISM,
            5: ElementType.PYRAMID,
            4: ElementType.TET,
        }.get(n, ElementType.OTHER)


@dataclass
class BlockElementSummary:
    """Per-block element-type tally + minimum divergence-theorem cell volume."""

    n_hex: int
    n_prism: int
    n_pyramid: int
    n_tet: int
    n_other: int
    #: Minimum (raw, signed) divergence-theorem cell volume in the block.
    min_volume: float
    #: The (i,j,k) of the minimum-volume cell.
    min_volume_cell: Tuple[int, int, int]


@dataclass
class ElementInventory:
    """Whole-mesh element census — the analogue of Code Leo's ELEMTYPE
    (type counts) + ADVOLAREA (per-element volume) load passes.
    """

    per_block: List[BlockElementSummary]
    n_hex: int
    n_prism: int
    n_pyramid: int
    n_tet: int
    n_other: int
    total: int
    global_min_volume: float
    #: (block, i, j, k) of the global-minimum-volume cell.
    global_min_cell: Tuple[int, int, int, int]

    def format_ads_style(self) -> str:
        """Code Leo-style census string (ELEMTYPE counts + ADVOLAREA
        min-vol), so a load echoes the reference solver's own printout.
        """
        lines = [f"Element census (ELEMTYPE analogue) — {self.total} elements:"]
        lines.append(f"  {self.n_hex:>10} ELEMENTS OF HEX     TYPE")
        lines.append(f"  {self.n_prism:>10} ELEMENTS OF PRISM   TYPE")
        lines.append(f"  {self.n_pyramid:>10} ELEMENTS OF PYRAMID TYPE")
        lines.append(f"  {self.n_tet:>10} ELEMENTS OF TET     TYPE")
        if self.n_other > 0:
            lines.append(
                f"  {self.n_other:>10} ELEMENTS OF OTHER   TYPE "
                f"(non-standard partial collapse)"
            )
        gb, gi, gj, gk = self.global_min_cell
        lines.append(
            f"Volume (ADVOLAREA analogue): global MIN VOLUME "
            f"{self.global_min_volume:.6e} at block {gb} cell ({gi},{gj},{gk})"
        )
        for bi, blk in enumerate(self.per_block):
            i, j, k = blk.min_volume_cell
            other = f", {blk.n_other} other" if blk.n_other > 0 else ""
            lines.append(
                f"  block {bi:>2}: MIN VOLUME {blk.min_volume:.4e} at ({i},{j},{k})  "
                f"[{blk.n_hex} hex, {blk.n_prism} prism, {blk.n_pyramid} pyr, "
                f"{blk.n_tet} tet{other}]"
            )
        return "\n".join(lines) + "\n"


def element_type_inventory(blocks: Sequence[Block]) -> ElementInventory:
    """Classify every cell of every block by collapsed-node topology and
    tally per-block volumes. A cheap second load-time pass, independent of
    `run_all`. Vectorized via `cell_distinct_node_count` + boolean masks,
    not per-cell dispatch.
    """
    per_block: List[BlockElementSummary] = []
    t_hex = t_prism = t_pyr = t_tet = t_other = 0
    g_min = math.inf
    g_cell = (0, 0, 0, 0)

    for bi, b in enumerate(blocks):
        NI, NJ, NK = cell_dims(b)
        if NI == 0 or NJ == 0 or NK == 0:
            per_block.append(
                BlockElementSummary(0, 0, 0, 0, 0, math.inf, (0, 0, 0))
            )
            continue

        distinct = cell_distinct_node_count(b)
        n_hex = int(np.count_nonzero(distinct == 8))
        n_prism = int(np.count_nonzero(distinct == 6))
        n_pyr = int(np.count_nonzero(distinct == 5))
        n_tet = int(np.count_nonzero(distinct == 4))
        n_other = NI * NJ * NK - (n_hex + n_prism + n_pyr + n_tet)

        true_vol = cell_volume_divergence(b)
        min_volume, (i, j, k) = _rust_order_argmin(true_vol)

        s = BlockElementSummary(n_hex, n_prism, n_pyr, n_tet, n_other, min_volume, (i, j, k))
        if s.min_volume < g_min:
            g_min = s.min_volume
            g_cell = (bi, i, j, k)

        t_hex += n_hex
        t_prism += n_prism
        t_pyr += n_pyr
        t_tet += n_tet
        t_other += n_other
        per_block.append(s)

    total = t_hex + t_prism + t_pyr + t_tet + t_other
    return ElementInventory(
        per_block=per_block,
        n_hex=t_hex,
        n_prism=t_prism,
        n_pyramid=t_pyr,
        n_tet=t_tet,
        n_other=t_other,
        total=total,
        global_min_volume=(g_min if math.isfinite(g_min) else 0.0),
        global_min_cell=g_cell,
    )
