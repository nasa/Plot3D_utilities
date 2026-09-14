"""Shared low-level face-orientation helpers.

This module holds the bottom-level building blocks used by
``connectivity.py``, ``verify.py``, and the ``correspondence.py``
certification module: the 8 canonical permutation matrices, canonical-grid
extraction, permutation application, and the ``Patch`` abstraction for a
single-constant-axis sub-region of a block. It sits below all three of those
modules (depending only on ``.block`` and numpy) specifically so that they
can all import from here without creating an import cycle.

Mirrors the corresponding low-level pieces of plot3d-rs's
``face_record.rs``/``correspondence.rs``.
"""

from typing import NamedTuple, Tuple

import numpy as np

from .block import Block


PERMUTATION_MATRICES = np.array([
    [[ 1,  0], [ 0,  1]],   # 0: identity
    [[-1,  0], [ 0,  1]],   # 1: u reversed
    [[ 1,  0], [ 0, -1]],   # 2: v reversed
    [[-1,  0], [ 0, -1]],   # 3: both reversed
    [[ 0,  1], [ 1,  0]],   # 4: swapped
    [[ 0, -1], [ 1,  0]],   # 5: swap + u reversed
    [[ 0,  1], [-1,  0]],   # 6: swap + v reversed
    [[ 0, -1], [-1,  0]],   # 7: swap + both reversed
], dtype=np.int8)
"""8 canonical 2x2 signed permutation matrices.

Bit encoding: ``index = u_reversed | (v_reversed << 1) | (swapped << 2)``
"""


def _constant_axis(lb: list, ub: list) -> int:
    """Return the index (0/1/2) of the single axis on which ``lb == ub``.

    Returns -1 if no such axis exists (or more than one does, in which
    case the first constant axis found is returned).
    """
    for d in range(3):
        if lb[d] == ub[d]:
            return d
    return -1


def extract_canonical_grid(block: Block, lb: list, ub: list) -> Tuple[np.ndarray, int, int]:
    """Extract face as a canonical 2D grid (nu, nv, 3) in ascending index order.

    Finds the constant axis, then extracts points with the first varying axis
    as the outer loop (u) and the second as the inner loop (v), both ascending.

    Args:
        block: Block to extract from.
        lb: Lower diagonal corner [i, j, k].
        ub: Upper diagonal corner [i, j, k].

    Returns:
        (grid, nu, nv) where grid has shape (nu, nv, 3).

    Raises:
        ValueError: If no constant axis is found.
    """
    lo = [min(lb[d], ub[d]) for d in range(3)]
    hi = [max(lb[d], ub[d]) for d in range(3)]

    const_dim = _constant_axis(lo, hi)
    if const_dim < 0:
        raise ValueError(f"No constant axis found: lo={lo}, hi={hi}")

    vary = [d for d in range(3) if d != const_dim]
    d0, d1 = vary
    nu = hi[d0] - lo[d0] + 1
    nv = hi[d1] - lo[d1] + 1

    grid = np.empty((nu, nv, 3))
    idx = [0, 0, 0]
    idx[const_dim] = lo[const_dim]
    for u in range(nu):
        idx[d0] = lo[d0] + u
        for v in range(nv):
            idx[d1] = lo[d1] + v
            grid[u, v] = [block.X[idx[0], idx[1], idx[2]],
                          block.Y[idx[0], idx[1], idx[2]],
                          block.Z[idx[0], idx[1], idx[2]]]
    return grid, nu, nv


def apply_permutation(grid: np.ndarray, perm_idx: int) -> np.ndarray:
    """Apply a pre-computed permutation matrix to a 2D face grid.

    Uses bit operations on ``perm_idx`` (0-7) to flip and/or transpose the grid.
    The permutation matrix is looked up from ``PERMUTATION_MATRICES``, not recalculated.

    Bit encoding: ``perm_idx = u_reversed | (v_reversed << 1) | (swapped << 2)``

    Args:
        grid: Face grid with shape (nu, nv, 3).
        perm_idx: Permutation index 0-7.

    Returns:
        Permuted grid with shape (out_nu, out_nv, 3).
    """
    g = grid
    if perm_idx & 1:
        g = g[::-1, :, :]    # flip u
    if perm_idx & 2:
        g = g[:, ::-1, :]    # flip v
    if perm_idx & 4:
        g = g.transpose(1, 0, 2)  # swap u, v
    return np.ascontiguousarray(g)


class _PatchFields(NamedTuple):
    block_index: int
    lo: Tuple[int, int, int]
    hi: Tuple[int, int, int]


class Patch(_PatchFields):
    """An ascending-normalized sub-region of a single block's face.

    ``lo``/``hi`` are elementwise ``lo <= hi`` (never a directed lb/ub
    diagonal) and must agree on exactly one axis -- the constant axis
    identifying which of the block's six logical faces (or a sub-region of
    one) this patch lies on.

    Note: subclasses the plain (functional-syntax) ``_PatchFields`` rather
    than inheriting ``NamedTuple`` directly, because ``typing.NamedTuple``
    forbids overriding ``__new__`` in the class body -- this is the
    standard workaround for adding constructor validation to a NamedTuple.
    """
    __slots__ = ()

    def __new__(cls, block_index: int, lo: Tuple[int, int, int], hi: Tuple[int, int, int]):
        lo = tuple(lo)
        hi = tuple(hi)
        if len(lo) != 3 or len(hi) != 3:
            raise ValueError(f"Patch lo/hi must be length-3, got lo={lo}, hi={hi}")
        if any(lo[d] > hi[d] for d in range(3)):
            raise ValueError(f"Patch lo/hi must be ascending-normalized: lo={lo}, hi={hi}")
        n_constant = sum(1 for d in range(3) if lo[d] == hi[d])
        if n_constant != 1:
            raise ValueError(
                f"Patch must have exactly one constant axis, found {n_constant}: "
                f"lo={lo}, hi={hi}")
        return super().__new__(cls, block_index, lo, hi)

    def const_axis(self) -> int:
        """The single axis (0/1/2) on which ``lo == hi``.

        Mirrors plot3d-rs's ``Patch::const_axis`` (``correspondence.rs``).
        Always succeeds -- the constructor already guarantees exactly one
        such axis exists.
        """
        for d in range(3):
            if self.lo[d] == self.hi[d]:
                return d
        raise AssertionError("Patch invariant violated: no constant axis")


def patch_from_bounds(block_index: int, a: Tuple[int, int, int],
                       b: Tuple[int, int, int]) -> Patch:
    """Build a :class:`Patch` from two (possibly non-ascending) corner triples.

    ``a``/``b`` may be given in diagonal (directed lb/ub) convention; this
    takes the elementwise min/max to normalize to ascending order.

    Args:
        block_index: Index of the block the patch lies on.
        a: First corner [i, j, k].
        b: Second corner [i, j, k].

    Returns:
        Patch with ``lo <= hi`` elementwise.
    """
    lo = tuple(min(a[d], b[d]) for d in range(3))
    hi = tuple(max(a[d], b[d]) for d in range(3))
    return Patch(block_index, lo, hi)


def extract_patch_grid(block: Block, patch: Patch) -> np.ndarray:
    """Extract a patch's node coordinates as an (nu, nv, 3) ascending-order array.

    A :class:`Patch` is exactly the ascending-normalized version of what
    :func:`extract_canonical_grid` extracts from a directed lb/ub face
    record, so this is a thin wrapper over it.

    Args:
        block: Block to extract from.
        patch: Patch describing the sub-region to extract.

    Returns:
        Array of shape (nu, nv, 3).
    """
    grid, _, _ = extract_canonical_grid(block, list(patch.lo), list(patch.hi))
    return grid
