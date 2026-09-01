"""Stage 4a - Metrics: turn the 2D node grid into finite-volume geometry.

The solver never looks at node coordinates again.  All it needs is, per cell, a
volume and a source area, and per face, an area-like weight and a unit normal.

Everything is **axisymmetric**: the 2D ``(x, r)`` mesh represents a full body of
revolution, so a cell's true volume is ``2*pi * integral(r dA)`` and a face's
true area is ``2*pi * integral(r dL)``.  The common ``2*pi`` divides out of the
finite-volume balance, so we store the per-radian quantities::

    vol   = area * r_centroid        (cell volume  / 2*pi)
    s     = length * r_midpoint      (face area    / 2*pi)

Two things fall out of this for free:

* the axis is not a singularity - ``sj[:, 0]`` is exactly zero because ``r = 0``
  there, so no flux can cross the axis no matter what the ghost state is;
* the duct's real area variation (``A ~ r**2``) is built into the metrics, so
  the 2D solve reproduces the 3D body of revolution exactly.

Index convention for a node grid of shape ``(NI+1, NJ+1)``:

* cells ``(NI, NJ)``; ``j = 0`` touches the axis, ``j = NJ-1`` touches the wall
* i-faces ``(NI+1, NJ)`` - constant-i faces, normal points toward +i
* j-faces ``(NI, NJ+1)`` - constant-j faces, normal points toward +j
"""
from typing import NamedTuple

import numpy as np


class Metrics(NamedTuple):
    """Finite-volume geometry of the flattened meridional mesh.

    A :class:`typing.NamedTuple` of arrays is automatically a valid JAX pytree,
    so this can be passed straight through ``jax.jit`` with no registration.
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


def build_metrics(x2d: np.ndarray, r2d: np.ndarray) -> Metrics:
    """Build cell volumes, face weights and face normals from a node grid.

    Args:
        x2d (np.ndarray): Node x coordinates, shape ``(NI+1, NJ+1)``.
        r2d (np.ndarray): Node r coordinates, shape ``(NI+1, NJ+1)``.

    Returns:
        Metrics: Geometry arrays as plain numpy (convert with :func:`to_jax`).
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

    return Metrics(xc=xc, rc=rc, area=area, vol=vol,
                   si=si, nix=nix, nir=nir,
                   sj=sj, njx=njx, njr=njr)


def to_jax(metrics: Metrics) -> Metrics:
    """Copy a :class:`Metrics` onto the JAX default device as float64 arrays.

    Args:
        metrics (Metrics): Metrics holding numpy arrays.

    Returns:
        Metrics: The same fields as ``jax.numpy`` arrays.
    """
    import jax.numpy as jnp
    return Metrics(*[jnp.asarray(f, dtype=jnp.float64) for f in metrics])


def to_numpy(metrics: Metrics) -> Metrics:
    """Copy a :class:`Metrics` back to plain numpy for the loop-based solver.

    Args:
        metrics (Metrics): Metrics holding numpy or jax arrays.

    Returns:
        Metrics: The same fields as numpy ``float64`` arrays.
    """
    return Metrics(*[np.asarray(f, dtype=float) for f in metrics])


def enclosed_volume(metrics: Metrics) -> float:
    """Total duct volume implied by the metrics, ``2*pi * sum(vol)``.

    Args:
        metrics (Metrics): Mesh metrics.

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


def ghost_cell_centres(xc: np.ndarray, rc: np.ndarray) -> dict:
    """Where the ghost cells of :mod:`euler_bc` would sit, for plotting only.

    The solver never stores ghost *positions* - :mod:`euler_bc` only builds
    ghost *states* (flow variables) on the fly.  This helper reconstructs the
    positions they conceptually occupy: one extra ring of cell centres outside
    each boundary, placed by mirroring the boundary-adjacent centre-to-centre
    spacing outward by one cell.  No corner ghosts, matching the solver, which
    pads the i and j directions separately.

    Args:
        xc (np.ndarray): Cell-centre x, shape ``(NI, NJ)``.
        rc (np.ndarray): Cell-centre r, shape ``(NI, NJ)``.

    Returns:
        dict: ``{"inlet", "outlet", "axis", "wall"} -> (x, r)`` arrays of
        ghost-cell centres.  The axis ring has ``r < 0``.
    """
    xc = np.asarray(xc)
    rc = np.asarray(rc)
    return {
        "inlet": (2 * xc[0] - xc[1], 2 * rc[0] - rc[1]),
        "outlet": (2 * xc[-1] - xc[-2], 2 * rc[-1] - rc[-2]),
        "axis": (2 * xc[:, 0] - xc[:, 1], 2 * rc[:, 0] - rc[:, 1]),
        "wall": (2 * xc[:, -1] - xc[:, -2], 2 * rc[:, -1] - rc[:, -2]),
    }


if __name__ == "__main__":
    from duct_geometry import duct_radius
    from duct_flatten import flatten_to_meridional
    from duct_mesh import revolve_duct

    x, r_wall = duct_radius(201)
    block = revolve_duct(x, r_wall)
    x2d, r2d = flatten_to_meridional(block)
    m = build_metrics(x2d, r2d)

    # Cross check against the analytic volume of the body of revolution.
    analytic = analytic_volume(x, r_wall)
    print(f"cell shape        : {m.vol.shape}")
    print(f"sum(vol) * 2*pi   : {enclosed_volume(m):.8f}")
    print(f"integral(pi R^2)dx: {analytic:.8f}")
    print(f"relative error    : {abs(enclosed_volume(m) / analytic - 1):.3e}")
    print(f"max |sj[:, 0]|    : {np.abs(m.sj[:, 0]).max():.3e}  (axis faces)")
