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
from dataclasses import dataclass
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
