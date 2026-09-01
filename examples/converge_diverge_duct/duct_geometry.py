"""Stage 1 - Geometry: the wall radius of a converging-diverging duct.

The duct wall is a single Bezier curve in the meridional ``(x, r)`` plane.  The
control polygon starts wide, pinches down to a throat, then opens up again, so
the resulting duct is a classic converging-diverging (de Laval) nozzle.

Nothing here knows about meshes or flow solvers - it only answers the question
"what is the wall radius at this axial station?".

Requires ``pyturbo-aero`` (``uv pip install pyturbo-aero`` or
``pip install -e ".[examples]"``).
"""
from typing import Tuple

import numpy as np
from pyturbo.helper import bezier

# Bezier control polygon of the wall in the meridional (x, r) plane.
# The repeated points near x = 0 flatten the curve out into a smooth throat.
CONTROL_X = [-0.90, -0.70, -0.60, -0.05, 0.00, 0.05, 0.75, 0.80, 1.00]
CONTROL_R = [0.25, 0.25, 0.25, 0.10, 0.10, 0.10, 0.50, 0.50, 0.50]


def duct_radius(n_axial: int = 201) -> Tuple[np.ndarray, np.ndarray]:
    """Wall radius sampled at evenly spaced axial stations.

    The Bezier curve is parametric in ``t``, which does **not** give evenly
    spaced ``x``.  We therefore sample it densely and interpolate onto a uniform
    axial grid, which is what the mesh generator wants.

    Args:
        n_axial (int, optional): Number of axial stations. Defaults to 201.

    Returns:
        Tuple[np.ndarray, np.ndarray]: ``(x, R)``, each of shape ``(n_axial,)``.
        ``x`` is uniformly spaced, ``R`` is the wall radius at those stations.
    """
    curve = bezier(CONTROL_X, CONTROL_R)
    t = np.linspace(0.0, 1.0, max(20 * n_axial, 2000))
    dense_x, dense_r = curve.get_point(t)
    dense_x = np.asarray(dense_x, dtype=float).ravel()
    dense_r = np.asarray(dense_r, dtype=float).ravel()

    # x(t) is monotonically increasing because the control x values are, so a
    # plain 1D interpolation is all that is needed to re-parameterize by x.
    x = np.linspace(dense_x[0], dense_x[-1], n_axial)
    r = np.interp(x, dense_x, dense_r)
    return x, r


def throat(x: np.ndarray, r: np.ndarray) -> Tuple[int, float, float]:
    """Locate the narrowest station of the duct.

    Args:
        x (np.ndarray): Axial stations.
        r (np.ndarray): Wall radius at those stations.

    Returns:
        Tuple[int, float, float]: ``(index, x_throat, r_throat)``.
    """
    i = int(np.argmin(r))
    return i, float(x[i]), float(r[i])


def area_ratio(r: np.ndarray) -> np.ndarray:
    """Local flow area divided by the throat area, ``A(x)/A_throat``.

    Circular cross sections, so the area ratio is simply ``(r/r_throat)**2``.

    Args:
        r (np.ndarray): Wall radius at each axial station.

    Returns:
        np.ndarray: ``A/A_throat`` at each station (minimum value 1.0).
    """
    return (r / r.min()) ** 2


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    x, r = duct_radius()
    i_t, x_t, r_t = throat(x, r)
    ar = area_ratio(r)

    print(f"axial stations : {x.size}  from x={x[0]:.3f} to x={x[-1]:.3f}")
    print(f"throat         : x={x_t:.4f}  r={r_t:.4f}  (index {i_t})")
    print(f"inlet  A/At    : {ar[0]:.3f}")
    print(f"exit   A/At    : {ar[-1]:.3f}")

    plt.figure(figsize=(8, 4))
    plt.plot(x, r, "C0", label="wall")
    plt.plot(x, -r, "C0")
    plt.plot([x_t], [r_t], "ro", label=f"throat (x={x_t:.2f})")
    plt.axhline(0.0, color="k", lw=0.5, ls="--")
    plt.xlabel("x")
    plt.ylabel("r")
    plt.title("Converging-diverging duct wall")
    plt.legend()
    plt.axis("equal")
    plt.tight_layout()
    plt.show()
