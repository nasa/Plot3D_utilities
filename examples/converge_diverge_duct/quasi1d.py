"""Analytic reference: quasi-one-dimensional nozzle flow with a normal shock.

Classic compressible-flow textbook material, used here purely as an independent
check on the CFD result.  The model assumes the flow is uniform across each
cross section and depends only on the local area - a good approximation for a
slender duct, and a useful sanity bound even for a steep one.

Given a choked throat and a back pressure that is too high for a fully
supersonic exit but too low for a fully subsonic duct, a normal shock parks
itself in the diverging section at exactly the station that makes the
post-shock subsonic expansion arrive at the back pressure.  This module finds
that station.
"""
from typing import NamedTuple

import numpy as np
from scipy.optimize import brentq

GAMMA = 1.4


class Quasi1DSolution(NamedTuple):
    """Analytic solution sampled at the requested axial stations."""
    x: np.ndarray          #: axial stations
    mach: np.ndarray       #: Mach number
    p: np.ndarray          #: static pressure, Pa
    x_shock: float         #: shock location
    mach_shock: float      #: Mach number just upstream of the shock
    mach_exit: float       #: exit Mach number
    p_exit: float          #: exit static pressure, Pa
    x_throat: float        #: throat location
    area_ratio_exit: float #: A_exit / A_throat


def area_ratio_from_mach(m: float) -> float:
    """Isentropic area ratio ``A/A*`` at a given Mach number.

    Args:
        m (float): Mach number (> 0).

    Returns:
        float: ``A/A*``.
    """
    e = (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0))
    return (1.0 / m) * ((2.0 / (GAMMA + 1.0))
                        * (1.0 + 0.5 * (GAMMA - 1.0) * m ** 2)) ** e


def mach_from_area_ratio(ar: float, supersonic: bool = False) -> float:
    """Invert the area-Mach relation.

    Two Mach numbers give the same area ratio, one either side of sonic; pick
    with ``supersonic``.

    Args:
        ar (float): ``A/A*``, must be >= 1.
        supersonic (bool, optional): Take the supersonic branch. Defaults to False.

    Returns:
        float: Mach number.
    """
    ar = max(float(ar), 1.0)
    if ar == 1.0:
        return 1.0
    bracket = (1.0, 50.0) if supersonic else (1.0e-6, 1.0)
    return float(brentq(lambda m: area_ratio_from_mach(m) - ar, *bracket))


def p_over_p0(m):
    """Isentropic static/stagnation pressure ratio.

    Args:
        m: Mach number.

    Returns:
        ``p/p0``.
    """
    return (1.0 + 0.5 * (GAMMA - 1.0) * np.asarray(m) ** 2) ** (-GAMMA / (GAMMA - 1.0))


def normal_shock(m1: float):
    """Normal-shock jump relations.

    Args:
        m1 (float): Upstream Mach number (>= 1).

    Returns:
        tuple: ``(m2, p2_over_p1, p02_over_p01)``.
    """
    g = GAMMA
    m2 = np.sqrt((1.0 + 0.5 * (g - 1.0) * m1 ** 2) / (g * m1 ** 2 - 0.5 * (g - 1.0)))
    p2_p1 = 1.0 + 2.0 * g / (g + 1.0) * (m1 ** 2 - 1.0)
    t2_t1 = (1.0 + 0.5 * (g - 1.0) * m1 ** 2) / (1.0 + 0.5 * (g - 1.0) * m2 ** 2)
    rho2_rho1 = p2_p1 / t2_t1
    p02_p01 = (rho2_rho1 ** (g / (g - 1.0))) * (p2_p1 ** (-1.0 / (g - 1.0)))
    return float(m2), float(p2_p1), float(p02_p01)


def solve(x: np.ndarray, r_wall: np.ndarray,
          p0: float = 2.0e5, t0: float = 500.0,
          p_back: float = 1.0e5) -> Quasi1DSolution:
    """Quasi-1D solution for a choked duct with a normal shock.

    Args:
        x (np.ndarray): Axial stations.
        r_wall (np.ndarray): Wall radius at those stations.
        p0 (float, optional): Inlet stagnation pressure, Pa. Defaults to 2.0e5.
        t0 (float, optional): Inlet stagnation temperature, K. Defaults to 500.0.
        p_back (float, optional): Back pressure, Pa. Defaults to 1.0e5.

    Returns:
        Quasi1DSolution: Mach and pressure distributions plus the shock station.

    Raises:
        ValueError: If the back pressure does not place a shock inside the duct.
    """
    x = np.asarray(x, dtype=float)
    r_wall = np.asarray(r_wall, dtype=float)
    i_throat = int(np.argmin(r_wall))
    ar = (r_wall / r_wall[i_throat]) ** 2          # A/A_throat = A/A*
    ar_exit = float(ar[-1])

    # Diverging section only; A must increase monotonically for the shock search.
    x_div = x[i_throat:]
    ar_div = np.maximum.accumulate(ar[i_throat:])

    def exit_pressure(x_shock: float) -> float:
        """Exit static pressure produced by a shock at ``x_shock``."""
        ar_s = float(np.interp(x_shock, x_div, ar_div))
        m1 = mach_from_area_ratio(ar_s, supersonic=True)
        _, _, p0_ratio = normal_shock(m1)
        # The post-shock flow has a larger sonic area because p0 dropped.
        ar_exit_2 = ar_exit * p0_ratio
        m_e = mach_from_area_ratio(ar_exit_2, supersonic=False)
        return p0 * p0_ratio * float(p_over_p0(m_e))

    lo, hi = float(x_div[0]), float(x_div[-1])
    if (exit_pressure(lo) - p_back) * (exit_pressure(hi) - p_back) > 0.0:
        raise ValueError("back pressure does not place a normal shock inside "
                         "the diverging section")
    x_shock = float(brentq(lambda xs: exit_pressure(xs) - p_back, lo, hi))

    ar_s = float(np.interp(x_shock, x_div, ar_div))
    m1 = mach_from_area_ratio(ar_s, supersonic=True)
    _, _, p0_ratio = normal_shock(m1)

    # Sample the whole duct.
    mach = np.empty_like(x)
    p = np.empty_like(x)
    for i, (xi, ari) in enumerate(zip(x, ar)):
        if xi <= x[i_throat]:                      # converging: subsonic
            m = mach_from_area_ratio(ari, supersonic=False)
            p[i] = p0 * p_over_p0(m)
        elif xi < x_shock:                         # diverging, pre-shock
            m = mach_from_area_ratio(ari, supersonic=True)
            p[i] = p0 * p_over_p0(m)
        else:                                      # post-shock: new p0, subsonic
            m = mach_from_area_ratio(ari * p0_ratio, supersonic=False)
            p[i] = p0 * p0_ratio * p_over_p0(m)
        mach[i] = m

    return Quasi1DSolution(x=x, mach=mach, p=p,
                           x_shock=x_shock, mach_shock=m1,
                           mach_exit=float(mach[-1]), p_exit=float(p[-1]),
                           x_throat=float(x[i_throat]),
                           area_ratio_exit=ar_exit)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from duct_geometry import duct_radius

    x, r_wall = duct_radius(401)
    sol = solve(x, r_wall)

    print(f"throat        : x = {sol.x_throat:.4f}")
    print(f"A_exit/A_throat: {sol.area_ratio_exit:.3f}")
    print(f"shock         : x = {sol.x_shock:.4f}, M1 = {sol.mach_shock:.4f}")
    print(f"exit          : M = {sol.mach_exit:.4f}, p = {sol.p_exit:.1f} Pa")

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    ax[0].plot(sol.x, sol.mach)
    ax[0].axvline(sol.x_shock, color="r", ls="--")
    ax[0].set_ylabel("Mach")
    ax[1].plot(sol.x, sol.p / 1000.0)
    ax[1].axvline(sol.x_shock, color="r", ls="--")
    ax[1].set_ylabel("p [kPa]")
    ax[1].set_xlabel("x")
    fig.suptitle("Quasi-1D reference solution")
    plt.tight_layout()
    plt.show()
