"""Stage 4a - Metrics: JAX/BC-specific glue on top of plot3d's meridional metrics.

Cell volumes, face weights and face normals now live in
:mod:`plot3d.meridional_flatten` (``build_metrics`` / ``MeridionalMetrics``).
This module keeps only what must not become a library dependency: converting
those metrics to/from JAX arrays, and reconstructing ghost-cell positions for
this example's own Euler boundary-condition scheme.
"""
import numpy as np

from plot3d import MeridionalMetrics


def to_jax(metrics: MeridionalMetrics) -> MeridionalMetrics:
    """Copy a :class:`MeridionalMetrics` onto the JAX default device as float64 arrays.

    Args:
        metrics (MeridionalMetrics): Metrics holding numpy arrays.

    Returns:
        MeridionalMetrics: The same fields as ``jax.numpy`` arrays.
    """
    import jax.numpy as jnp
    return MeridionalMetrics(*[jnp.asarray(f, dtype=jnp.float64) for f in metrics])


def to_numpy(metrics: MeridionalMetrics) -> MeridionalMetrics:
    """Copy a :class:`MeridionalMetrics` back to plain numpy for the loop-based solver.

    Args:
        metrics (MeridionalMetrics): Metrics holding numpy or jax arrays.

    Returns:
        MeridionalMetrics: The same fields as numpy ``float64`` arrays.
    """
    return MeridionalMetrics(*[np.asarray(f, dtype=float) for f in metrics])


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
    from duct_mesh import revolve_duct
    from plot3d import (analytic_volume, build_metrics, enclosed_volume,
                         flatten_to_meridional)

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
