"""Stage 2 - 3D mesh: revolve the wall radius into a Plot3D block.

The duct is a body of revolution about the x-axis, so a 3D mesh is built by
sweeping the meridional ``(x, r)`` node distribution through a range of angles.

We sweep a **wedge** (90 degrees by default) instead of the full 360 degrees.
The two end planes at theta = 0 and theta = wedge are *symmetry planes*, not
physical walls.  A wedge keeps the wireframe plot readable and the mesh file
small, and because the flow is axisymmetric it carries exactly the same
information as a full revolve.

Radial spacing uses a geometric **inflation layer** clustered at the wall, the
way a viscous mesh would be built, so the flatten stage has something visibly
non-uniform to show off.
"""
from typing import Tuple

import numpy as np
from plot3d import Block


def radial_fractions(n_radial: int = 25, growth: float = 1.15) -> np.ndarray:
    """Normalized radial node positions with a wall-clustered inflation layer.

    Cell heights grow geometrically as you march *away* from the wall, so the
    first cell off the wall is the thinnest one.

    Args:
        n_radial (int, optional): Number of radial nodes. Defaults to 25.
        growth (float, optional): Cell-to-cell growth ratio. Defaults to 1.15.
            Use ``1.0`` for uniform spacing.

    Returns:
        np.ndarray: ``s`` of shape ``(n_radial,)`` increasing from ``s[0] = 0``
        on the axis to ``s[-1] = 1`` at the wall.
    """
    h = growth ** np.arange(n_radial - 1, dtype=float)  # index 0 = wall cell
    h = h / h.sum()
    d = np.concatenate([[0.0], np.cumsum(h)])           # distance from wall
    s = 1.0 - d[::-1]                                   # flip: 0 = axis, 1 = wall
    s[0] = 0.0
    s[-1] = 1.0
    return s


def revolve_duct(x: np.ndarray,
                 r_wall: np.ndarray,
                 n_radial: int = 25,
                 n_theta: int = 13,
                 wedge_deg: float = 90.0,
                 growth: float = 1.15) -> Block:
    """Revolve the wall radius curve into a 3D :class:`plot3d.Block`.

    The convention matches :meth:`plot3d.Block.cylindrical`, which recovers
    ``r = sqrt(Y**2 + Z**2)`` and ``theta = arctan2(Y, Z)`` about the x-axis::

        X[i,j,k] = x[i]
        Y[i,j,k] = r[i,j] * sin(theta[k])
        Z[i,j,k] = r[i,j] * cos(theta[k])

    Args:
        x (np.ndarray): Axial stations, shape ``(ni,)``.
        r_wall (np.ndarray): Wall radius at those stations, shape ``(ni,)``.
        n_radial (int, optional): Nodes from axis to wall. Defaults to 25.
        n_theta (int, optional): Nodes across the wedge. Defaults to 13.
        wedge_deg (float, optional): Wedge angle in degrees. Defaults to 90.0.
        growth (float, optional): Inflation-layer growth ratio. Defaults to 1.15.

    Returns:
        Block: Block of shape ``(ni, n_radial, n_theta)``.
    """
    s = radial_fractions(n_radial, growth)
    theta = np.radians(np.linspace(0.0, wedge_deg, n_theta))

    r2d = r_wall[:, None] * s[None, :]                  # (ni, n_radial)

    # Broadcast to (ni, n_radial, n_theta); no Python loops.
    X = np.broadcast_to(x[:, None, None], r2d.shape + (n_theta,)).copy()
    Y = r2d[:, :, None] * np.sin(theta)[None, None, :]
    Z = r2d[:, :, None] * np.cos(theta)[None, None, :]
    return Block(X, np.ascontiguousarray(Y), np.ascontiguousarray(Z))


def inflation_report(n_radial: int = 25,
                     growth: float = 1.15) -> Tuple[float, float, float]:
    """Describe the inflation layer as fractions of the local wall radius.

    The fractions are the same at every axial station (the same ``s`` array is
    scaled by the local radius), so no wall radius is needed here.

    Args:
        n_radial (int, optional): Nodes from axis to wall. Defaults to 25.
        growth (float, optional): Inflation-layer growth ratio. Defaults to 1.15.

    Returns:
        Tuple[float, float, float]: ``(wall_cell_fraction, axis_cell_fraction,
        clustering_ratio)`` where the fractions are cell heights divided by the
        local radius and the ratio is coarsest/finest.
    """
    s = radial_fractions(n_radial, growth)
    ds = np.diff(s)
    wall_frac = float(ds[-1])
    axis_frac = float(ds[0])
    return wall_frac, axis_frac, axis_frac / wall_frac


if __name__ == "__main__":
    from plot3d import plot_blocks

    from duct_geometry import duct_radius

    x, r_wall = duct_radius(201)
    block = revolve_duct(x, r_wall)
    wall_frac, axis_frac, ratio = inflation_report()

    print(f"block shape     : {block.IMAX} x {block.JMAX} x {block.KMAX}")
    print(f"nodes           : {block.X.size}")
    print(f"wall cell       : {100 * wall_frac:.3f}% of local R")
    print(f"axis cell       : {100 * axis_frac:.3f}% of local R")
    print(f"clustering ratio: {ratio:.1f} : 1")

    plot_blocks([block])
