"""Node-for-node face-correspondence certification.

Mirrors plot3d-rs's ``correspondence.rs`` (commit ``0e1b1a1``). A conformal
structured interface is not "four corners agree" -- every node, interior
nodes included, must correspond under exactly one of the 8 structured
permutation mappings, within tolerance:

- zero permutations passing is a definite rejection (:class:`ExceedsTolerance`,
  reporting the best-scoring permutation's worst-node discrepancy for
  diagnostics);
- more than one permutation passing is *also* a rejection
  (:class:`Ambiguous`) -- a well-formed match has exactly one valid
  structured mapping between two patches, so passing under two different
  permutations means the patch geometry is too symmetric/degenerate to
  certify unambiguously, not a match to guess at.

This module depends only on ``.block`` and ``.permutation`` -- not on
``.connectivity``, ``.periodicity``, or ``.verify`` -- so that those modules
can import it without an import cycle. Unlike Rust's per-node early-exit
scan, every candidate permutation here is checked with one vectorized numpy
distance array: face patches are at most a few thousand points, and
simplicity matters more than micro-optimization for this module.
"""

from typing import Callable, List, NamedTuple, Optional, Tuple

import numpy as np

from .block import Block
from .permutation import Patch, apply_permutation, extract_patch_grid


class NodeDiscrepancy(NamedTuple):
    """One node pair's separation.

    ``node_a``/``node_b`` are each block's global ``(i, j, k)`` index (not
    the patch-local ``(u, v)``) -- this mirrors ``correspondence.rs``, whose
    ``NodeDiscrepancy.node_a``/``node_b`` are produced via ``Patch::ijk``
    (global structured indices), so a caller can locate the discrepancy
    directly in the block without re-deriving it from patch bounds.
    """
    distance: float
    node_a: Tuple[int, int, int]
    node_b: Tuple[int, int, int]


class CertifiedMapping(NamedTuple):
    """A structured mapping under which every node pair is within tolerance."""
    permutation_index: int
    plane: str  # 'in-plane' or 'cross-plane' -- see connectivity._orient_vec_to_permutation
    worst: NodeDiscrepancy
    nodes_checked: int


class MappingFailure(Exception):
    """Base class for a failed certification attempt."""


class IncompatibleDimensions(MappingFailure):
    """No permutation makes the two patches' dimensions agree."""

    def __init__(self, dims_a: Tuple[int, int], dims_b: Tuple[int, int]):
        self.dims_a = dims_a
        self.dims_b = dims_b
        super().__init__(
            f"patch dimensions {dims_a} and {dims_b} admit no structured mapping")


class ExceedsTolerance(MappingFailure):
    """Every dimension-admissible permutation has a node pair beyond tolerance.

    ``worst`` is the largest separation of the permutation that came
    closest (``best_permutation``).
    """

    def __init__(self, best_permutation: int, worst: NodeDiscrepancy):
        self.best_permutation = best_permutation
        self.worst = worst
        super().__init__(
            "no structured mapping keeps every node within tolerance; "
            f"permutation {best_permutation} comes closest, failing at "
            f"A{worst.node_a} <-> B{worst.node_b} with separation {worst.distance:e}")


class Ambiguous(MappingFailure):
    """More than one permutation passes: the geometry cannot distinguish them."""

    def __init__(self, permutations: List[int]):
        self.permutations = permutations
        super().__init__(
            f"{len(permutations)} structured mappings {permutations} all fit "
            "within tolerance; the patch geometry cannot distinguish them")


def _uv_axes(const_axis: int) -> Tuple[int, int]:
    """The two varying axes in ascending order, given the constant axis."""
    if const_axis == 0:
        return (1, 2)
    if const_axis == 1:
        return (0, 2)
    return (0, 1)


def _patch_ijk(patch: Patch, u: int, v: int) -> Tuple[int, int, int]:
    """Global ``(i, j, k)`` of patch-local node ``(u, v)``.

    Mirrors Rust's ``Patch::ijk``.
    """
    ua, va = _uv_axes(patch.const_axis())
    idx = list(patch.lo)
    idx[ua] = patch.lo[ua] + u
    idx[va] = patch.lo[va] + v
    return (idx[0], idx[1], idx[2])


def _map_uv(perm: int, u: int, v: int, nu_b: int, nv_b: int) -> Tuple[int, int]:
    """Patch B's local ``(u, v)`` (pre-permutation) that patch A's ``(u, v)``
    corresponds to under ``perm``.

    Mirrors Rust's ``correspondence::map_uv`` exactly. Used only to translate
    an index position on the vectorized distance array back into patch B's
    own raw grid for :class:`NodeDiscrepancy` reporting -- the distance
    values themselves come from :func:`apply_permutation` on the full grid.
    """
    if perm & 4:
        gu, gv = v, u
    else:
        gu, gv = u, v
    if perm & 1:
        gu = nu_b - 1 - gu
    if perm & 2:
        gv = nv_b - 1 - gv
    return gu, gv


def _plane_for(patch_a: Patch, patch_b: Patch) -> str:
    """'in-plane' when both patches share a constant axis, else 'cross-plane'.

    Same convention as ``connectivity._orient_vec_to_permutation``.
    """
    return 'in-plane' if patch_a.const_axis() == patch_b.const_axis() else 'cross-plane'


def _scan_permutation(
    grid_a: np.ndarray, patch_a: Patch,
    grid_b: np.ndarray, patch_b: Patch,
    perm: int,
) -> Optional[Tuple[float, NodeDiscrepancy]]:
    """Max distance and worst node pair for one permutation.

    Returns ``None`` if ``perm``'s permuted dims don't match ``grid_a``'s.
    """
    nu_b, nv_b = grid_b.shape[0], grid_b.shape[1]
    permuted_b = apply_permutation(grid_b, perm)
    if permuted_b.shape[:2] != grid_a.shape[:2]:
        return None

    dist = np.linalg.norm(grid_a - permuted_b, axis=-1)
    u, v = (int(x) for x in np.unravel_index(int(np.argmax(dist)), dist.shape))
    max_dist = float(dist[u, v])
    ub, vb = _map_uv(perm, u, v, nu_b, nv_b)
    worst = NodeDiscrepancy(
        distance=max_dist,
        node_a=_patch_ijk(patch_a, u, v),
        node_b=_patch_ijk(patch_b, ub, vb),
    )
    return max_dist, worst


def _extract_grids(
    block_a: Block, patch_a: Patch, block_b: Block, patch_b: Patch,
    transform: Optional[Callable[[np.ndarray], np.ndarray]],
) -> Tuple[np.ndarray, np.ndarray]:
    grid_a = extract_patch_grid(block_a, patch_a)
    grid_b = extract_patch_grid(block_b, patch_b)
    if transform is not None:
        grid_b = transform(grid_b)
    return grid_a, grid_b


def certify_correspondence(
    block_a: Block, patch_a: Patch, block_b: Block, patch_b: Patch,
    tol: float, transform: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> CertifiedMapping:
    """Certify that ``patch_b`` (optionally ``transform``-ed) corresponds
    node-for-node to ``patch_a`` under exactly one structured permutation.

    Mirrors plot3d-rs ``correspondence::certify_correspondence`` (commit
    ``0e1b1a1``). Every permutation whose dims are compatible with
    ``patch_a`` is scanned (not just the first passing one); zero passing
    raises :class:`ExceedsTolerance`, more than one raises
    :class:`Ambiguous`, and no dimension-compatible permutation at all
    raises :class:`IncompatibleDimensions`.

    ``transform``, if given, is applied to ``patch_b``'s extracted grid
    before comparison (e.g. a periodic rotation bringing block B's points
    into block A's frame) -- it must return an array of the same shape it
    was given.
    """
    grid_a, grid_b = _extract_grids(block_a, patch_a, block_b, patch_b, transform)

    candidates: List[Tuple[int, float, NodeDiscrepancy]] = []
    for perm in range(8):
        result = _scan_permutation(grid_a, patch_a, grid_b, patch_b, perm)
        if result is not None:
            max_dist, worst = result
            candidates.append((perm, max_dist, worst))

    if not candidates:
        raise IncompatibleDimensions(
            (int(grid_a.shape[0]), int(grid_a.shape[1])),
            (int(grid_b.shape[0]), int(grid_b.shape[1])))

    plane = _plane_for(patch_a, patch_b)
    passing = [c for c in candidates if c[1] <= tol]

    if len(passing) == 1:
        perm, _, worst = passing[0]
        return CertifiedMapping(
            permutation_index=perm,
            plane=plane,
            worst=worst,
            nodes_checked=int(grid_a.shape[0] * grid_a.shape[1]))

    if len(passing) == 0:
        best_permutation, _, worst = min(candidates, key=lambda c: c[1])
        raise ExceedsTolerance(best_permutation=best_permutation, worst=worst)

    raise Ambiguous(permutations=sorted(c[0] for c in passing))


def certify_permutation(
    block_a: Block, patch_a: Patch, block_b: Block, patch_b: Patch,
    perm_idx: int, tol: float,
    transform: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> CertifiedMapping:
    """Certify one declared permutation only -- no search over the other 7.

    Mirrors plot3d-rs ``correspondence::certify_permutation`` (commit
    ``0e1b1a1``). Raises :class:`IncompatibleDimensions` if ``perm_idx``'s
    permuted dims don't match ``patch_a``'s, or :class:`ExceedsTolerance`
    (with ``best_permutation=perm_idx``) if the max distance under that one
    permutation exceeds ``tol`` -- there is no fallback search, so a wrong
    declared orientation is a definite failure, never silently replaced by a
    different one found via search. Cannot raise :class:`Ambiguous` (only
    one permutation is ever tried).
    """
    grid_a, grid_b = _extract_grids(block_a, patch_a, block_b, patch_b, transform)

    result = _scan_permutation(grid_a, patch_a, grid_b, patch_b, perm_idx)
    if result is None:
        raise IncompatibleDimensions(
            (int(grid_a.shape[0]), int(grid_a.shape[1])),
            (int(grid_b.shape[0]), int(grid_b.shape[1])))

    max_dist, worst = result
    if max_dist > tol:
        raise ExceedsTolerance(best_permutation=perm_idx, worst=worst)

    return CertifiedMapping(
        permutation_index=perm_idx,
        plane=_plane_for(patch_a, patch_b),
        worst=worst,
        nodes_checked=int(grid_a.shape[0] * grid_a.shape[1]))
