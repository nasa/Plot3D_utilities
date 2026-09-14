"""Axisymmetric body-of-revolution -> 2D meridional finite-volume geometry.

A body of revolution stores the same ``(x, r)`` mesh once per theta plane:
every theta plane of a true body of revolution is identical, so a 3D
:class:`~plot3d.block.Block` built by revolving a wall curve about the x-axis
can be collapsed to a single 2D ``(x, r)`` node grid with no loss of
information, then turned into 2D axisymmetric finite-volume geometry (cell
volumes, face weights, face normals) suitable for a 2D solver that reproduces
the 3D body of revolution exactly.

Deliberately distinct from two other, differently-scoped "flatten" concepts
in this package, to avoid confusion:

* :mod:`plot3d.flatmesh` (``flatten_mesh``/``FlatMesh``) flattens a
  *multi-block 3D structured* mesh into an *unstructured* finite-volume graph
  (full 3D hex cells, owner/neighbor adjacency across welded block
  boundaries). It has no axisymmetric/2D-collapse logic and does not require
  a body of revolution.
* :mod:`plot3d.glennht.plot3d_flatten_deck` exports a GlennHT solver "flatten
  deck" file format; it is a boundary-condition/connectivity export, not a
  geometry module.

This module instead collapses a *single*, *axisymmetric* block's redundant
theta dimension and builds the resulting 2D mesh's finite-volume metrics.

Everything here is axisymmetric: a cell's true volume is
``2*pi * integral(r dA)`` and a face's true area is ``2*pi * integral(r dL)``;
the common ``2*pi`` divides out of a finite-volume balance, so
:func:`build_metrics` stores the per-radian quantities
(``vol = area * r_centroid``, face weight ``s = length * r_midpoint``). The
axis is not a singularity: the ``j = 0`` face weight is exactly zero because
``r = 0`` there.
"""
from typing import NamedTuple, Tuple

import numpy as np

from .block import Block


def block_radius(block: Block) -> np.ndarray:
    """Radius about the x-axis at every node of a block.

    Args:
        block (Block): Body-of-revolution block about the x-axis.

    Returns:
        np.ndarray: ``r`` of shape ``(IMAX, JMAX, KMAX)``.
    """
    return np.sqrt(block.Y ** 2 + block.Z ** 2)


def axisymmetry_error(block: Block) -> float:
    """How far the block is from being a true body of revolution.

    Compares the radius at every theta plane against the first plane. For a
    mesh that is a genuine body of revolution this is round-off sized, which
    is the license to keep only one plane.

    Args:
        block (Block): Block to check.

    Returns:
        float: ``max |r(i,j,k) - r(i,j,0)|``.
    """
    r = block_radius(block)
    return float(np.max(np.abs(r - r[:, :, 0:1])))


def flatten_to_meridional(block: Block, k_index: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Extract one constant-theta slice as a 2D ``(x, r)`` node grid.

    Args:
        block (Block): Body-of-revolution block about the x-axis.
        k_index (int, optional): Which theta plane to keep. Defaults to 0.

    Returns:
        Tuple[np.ndarray, np.ndarray]: ``(x2d, r2d)``, each of shape
        ``(IMAX, JMAX)``. Index ``j = 0`` sits on the axis, ``j = JMAX-1`` on
        the wall.
    """
    x2d = np.ascontiguousarray(block.X[:, :, k_index], dtype=float)
    r2d = np.ascontiguousarray(block_radius(block)[:, :, k_index], dtype=float)
    return x2d, r2d


def node_count_reduction(block: Block) -> Tuple[int, int, float]:
    """Nodes before and after flattening.

    Args:
        block (Block): The 3D block.

    Returns:
        Tuple[int, int, float]: ``(nodes_3d, nodes_2d, factor)``.
    """
    n3 = int(block.X.size)
    n2 = int(block.IMAX * block.JMAX)
    return n3, n2, n3 / n2


class MeridionalMetrics(NamedTuple):
    """Finite-volume geometry of a flattened meridional mesh.

    A :class:`typing.NamedTuple` of arrays is automatically a valid JAX
    pytree, so downstream code can pass this straight through ``jax.jit``
    with no registration, without this module depending on JAX itself.

    Index convention for a node grid of shape ``(NI+1, NJ+1)``: cells
    ``(NI, NJ)``; ``j = 0`` touches the axis, ``j = NJ-1`` touches the wall.
    i-faces ``(NI+1, NJ)`` are constant-i faces with normal toward +i.
    j-faces ``(NI, NJ+1)`` are constant-j faces with normal toward +j.
    """
    xc: np.ndarray    # (NI, NJ)    cell centroid x
    rc: np.ndarray    # (NI, NJ)    cell centroid r
    area: np.ndarray  # (NI, NJ)    planar (x, r) cell area
    vol: np.ndarray   # (NI, NJ)    cell volume per radian = area * rc
    si: np.ndarray    # (NI+1, NJ)  i-face weight = length * r_mid
    nix: np.ndarray   # (NI+1, NJ)  i-face unit normal, x component
    nir: np.ndarray   # (NI+1, NJ)  i-face unit normal, r component
    sj: np.ndarray    # (NI, NJ+1)  j-face weight = length * r_mid
    njx: np.ndarray   # (NI, NJ+1)  j-face unit normal, x component
    njr: np.ndarray   # (NI, NJ+1)  j-face unit normal, r component


def build_metrics(x2d: np.ndarray, r2d: np.ndarray) -> MeridionalMetrics:
    """Build cell volumes, face weights and face normals from a node grid.

    Args:
        x2d (np.ndarray): Node x coordinates, shape ``(NI+1, NJ+1)``.
        r2d (np.ndarray): Node r coordinates, shape ``(NI+1, NJ+1)``.

    Returns:
        MeridionalMetrics: Geometry arrays as plain numpy.
    """
    x = np.asarray(x2d, dtype=float)
    r = np.asarray(r2d, dtype=float)

    # --- cells: corners in counter-clockwise order in the (x, r) plane -------
    x0, r0 = x[:-1, :-1], r[:-1, :-1]
    x1, r1 = x[1:, :-1], r[1:, :-1]
    x2, r2 = x[1:, 1:], r[1:, 1:]
    x3, r3 = x[:-1, 1:], r[:-1, 1:]

    # Shoelace formula written as two triangle cross products.
    area = 0.5 * np.abs((x2 - x0) * (r3 - r1) - (x3 - x1) * (r2 - r0))
    xc = 0.25 * (x0 + x1 + x2 + x3)
    rc = 0.25 * (r0 + r1 + r2 + r3)
    vol = area * rc

    # --- i-faces: node (i, j) -> node (i, j+1), tangent points toward +j -----
    dxi = x[:, 1:] - x[:, :-1]
    dri = r[:, 1:] - r[:, :-1]
    li = np.sqrt(dxi ** 2 + dri ** 2)
    # Rotating the tangent by -90 degrees gives the +i-pointing normal.
    nix = dri / li
    nir = -dxi / li
    si = li * 0.5 * (r[:, 1:] + r[:, :-1])

    # --- j-faces: node (i, j) -> node (i+1, j), tangent points toward +i -----
    dxj = x[1:, :] - x[:-1, :]
    drj = r[1:, :] - r[:-1, :]
    lj = np.sqrt(dxj ** 2 + drj ** 2)
    # Rotating the tangent by +90 degrees gives the +j-pointing normal.
    njx = -drj / lj
    njr = dxj / lj
    sj = lj * 0.5 * (r[1:, :] + r[:-1, :])

    return MeridionalMetrics(xc=xc, rc=rc, area=area, vol=vol,
                              si=si, nix=nix, nir=nir,
                              sj=sj, njx=njx, njr=njr)


def enclosed_volume(metrics: MeridionalMetrics) -> float:
    """Total duct volume implied by the metrics, ``2*pi * sum(vol)``.

    Args:
        metrics (MeridionalMetrics): Mesh metrics.

    Returns:
        float: Volume of the full body of revolution.
    """
    return float(2.0 * np.pi * np.sum(np.asarray(metrics.vol)))


def analytic_volume(x: np.ndarray, r_wall: np.ndarray) -> float:
    """Volume of the body of revolution, ``integral(pi R**2) dx``.

    Trapezoidal rule written out by hand (``np.trapezoid`` is numpy >= 2.0
    only and ``np.trapz`` has since been removed, so neither name is safe).

    Args:
        x (np.ndarray): Axial stations.
        r_wall (np.ndarray): Wall radius at those stations.

    Returns:
        float: The reference volume :func:`enclosed_volume` should match.
    """
    x = np.asarray(x, dtype=float)
    a = np.pi * np.asarray(r_wall, dtype=float) ** 2
    return float(np.sum(0.5 * (a[1:] + a[:-1]) * np.diff(x)))
