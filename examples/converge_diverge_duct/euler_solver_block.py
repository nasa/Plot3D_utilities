"""The same solver, written the slow way: nested Python loops over plain floats.

This module is a **deliberate duplicate** of :mod:`euler_solver_flatten`.  The
governing equations, the Rusanov flux, the ghost-cell boundary conditions, the
local time step and the 4-stage Runge-Kutta update are all identical - they are
simply re-derived here for a single cell at a time, with ``math.sqrt`` and
builtin ``max``/``abs``, inside ``for i: for j:`` loops.

Nothing here imports JAX.  That is the point: run both solvers on the same mesh
from the same initial state and they agree to round-off, but one of them is two
or three orders of magnitude slower.  Flattening the mesh is what makes the flow
field a dense rectangular array in the first place; writing the solver as array
operations on that array is what lets a compiler vectorize it.

It would have been shorter to call :mod:`euler_physics` with scalar arguments,
but that would have measured JAX's per-operation dispatch cost rather than the
cost of not vectorizing, so the formulas are spelled out again here.
"""
import math
from typing import List, Optional, Tuple

import numpy as np

from euler_metrics import Metrics

GAMMA = 1.4
CP = 1005.0
R_GAS = CP * (GAMMA - 1.0) / GAMMA
RHO_FLOOR = 1.0e-4
P_FLOOR = 10.0

# Runge-Kutta stage coefficients, same as the vectorized solver.
RK_ALPHA = (0.25, 1.0 / 3.0, 0.5, 1.0)


# --------------------------------------------------------------------------
# Scalar physics
# --------------------------------------------------------------------------

def _primitives(U):
    """Primitive variables of one cell.

    Args:
        U: Sequence of 4 floats, ``[rho, rho*u, rho*v, rho*E]``.

    Returns:
        tuple: ``(rho, u, v, p)`` as floats.
    """
    rho = max(float(U[0]), RHO_FLOOR)
    u = float(U[1]) / rho
    v = float(U[2]) / rho
    p = max((GAMMA - 1.0) * (float(U[3]) - 0.5 * rho * (u * u + v * v)), P_FLOOR)
    return rho, u, v, p


def _conservative(rho, u, v, p):
    """Conservative state of one cell.

    Args:
        rho: Density.
        u: Axial velocity.
        v: Radial velocity.
        p: Static pressure.

    Returns:
        tuple: 4 floats.
    """
    return (rho, rho * u, rho * v, p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v))


def _flux(U, nx, nr):
    """Physical Euler flux of one cell projected on a normal.

    Args:
        U: Sequence of 4 floats.
        nx: Normal, axial component.
        nr: Normal, radial component.

    Returns:
        tuple: 4 floats.
    """
    rho, u, v, p = _primitives(U)
    qn = u * nx + v * nr
    return (rho * qn,
            rho * u * qn + p * nx,
            rho * v * qn + p * nr,
            (float(U[3]) + p) * qn)


def _rusanov(UL, UR, nx, nr):
    """Rusanov numerical flux across one face.

    Args:
        UL: State on the -side, 4 floats.
        UR: State on the +side, 4 floats.
        nx: Normal, axial component.
        nr: Normal, radial component.

    Returns:
        tuple: 4 floats.
    """
    rhoL, uL, vL, pL = _primitives(UL)
    rhoR, uR, vR, pR = _primitives(UR)
    aL = math.sqrt(GAMMA * pL / rhoL)
    aR = math.sqrt(GAMMA * pR / rhoR)
    smax = max(abs(uL * nx + vL * nr) + aL, abs(uR * nx + vR * nr) + aR)
    fL = _flux(UL, nx, nr)
    fR = _flux(UR, nx, nr)
    return tuple(0.5 * (fL[c] + fR[c]) - 0.5 * smax * (float(UR[c]) - float(UL[c]))
                 for c in range(4))


# --------------------------------------------------------------------------
# Scalar boundary conditions
# --------------------------------------------------------------------------

def _inlet(U, p0, t0):
    """Subsonic inflow ghost cell from a stagnation state.

    Args:
        U: Interior state, 4 floats.
        p0: Stagnation pressure, Pa.
        t0: Stagnation temperature, K.

    Returns:
        tuple: Ghost state, 4 floats.
    """
    rho, u, v, p = _primitives(U)
    a = math.sqrt(GAMMA * p / rho)
    a0_sq = GAMMA * R_GAS * t0
    r_minus = u - 2.0 * a / (GAMMA - 1.0)

    A = (GAMMA + 1.0) / (GAMMA - 1.0)
    B = 2.0 * r_minus
    C = 0.5 * (GAMMA - 1.0) * r_minus * r_minus - a0_sq
    a_b = (-B + math.sqrt(max(B * B - 4.0 * A * C, 0.0))) / (2.0 * A)

    u_b = r_minus + 2.0 * a_b / (GAMMA - 1.0)
    t_b = a_b * a_b / (GAMMA * R_GAS)
    m_b = u_b / a_b
    p_b = p0 * (1.0 + 0.5 * (GAMMA - 1.0) * m_b * m_b) ** (-GAMMA / (GAMMA - 1.0))
    return _conservative(p_b / (R_GAS * t_b), u_b, 0.0, p_b)


def _outlet(U, p_back):
    """Subsonic outflow ghost cell at a fixed back pressure.

    Args:
        U: Interior state, 4 floats.
        p_back: Back pressure, Pa.

    Returns:
        tuple: Ghost state, 4 floats.
    """
    rho, u, v, p = _primitives(U)
    a = math.sqrt(GAMMA * p / rho)
    p_b = p if abs(u) / a >= 1.0 else p_back
    return _conservative(rho, u, v, p_b)


def _reflect(U, nx, nr):
    """Slip-wall / axis ghost cell.

    Args:
        U: Interior state, 4 floats.
        nx: Face normal, axial component.
        nr: Face normal, radial component.

    Returns:
        tuple: Ghost state, 4 floats.
    """
    rho, u, v, p = _primitives(U)
    qn = u * nx + v * nr
    return _conservative(rho, u - 2.0 * qn * nx, v - 2.0 * qn * nr, p)


# --------------------------------------------------------------------------
# Residual, time step, Runge-Kutta - all as explicit loops
# --------------------------------------------------------------------------

def residual(U: np.ndarray, m: Metrics, cfg) -> np.ndarray:
    """Net flux imbalance of every cell, one face at a time.

    Args:
        U (np.ndarray): Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics as numpy arrays.
        cfg: Object with ``p0``, ``t0``, ``p_back`` attributes.

    Returns:
        np.ndarray: Residual, shape ``(NI, NJ, 4)``.
    """
    ni, nj = m.vol.shape
    R = np.zeros((ni, nj, 4))

    # i-direction faces (inlet at i = 0, outlet at i = ni)
    for i in range(ni + 1):
        for j in range(nj):
            if i == 0:
                UR = U[0, j]
                UL = _inlet(UR, cfg.p0, cfg.t0)
            elif i == ni:
                UL = U[ni - 1, j]
                UR = _outlet(UL, cfg.p_back)
            else:
                UL = U[i - 1, j]
                UR = U[i, j]
            f = _rusanov(UL, UR, m.nix[i, j], m.nir[i, j])
            s = m.si[i, j]
            for c in range(4):
                if i > 0:
                    R[i - 1, j, c] += f[c] * s
                if i < ni:
                    R[i, j, c] -= f[c] * s

    # j-direction faces (axis at j = 0, wall at j = nj)
    for i in range(ni):
        for j in range(nj + 1):
            nx, nr = m.njx[i, j], m.njr[i, j]
            if j == 0:
                UR = U[i, 0]
                UL = _reflect(UR, nx, nr)
            elif j == nj:
                UL = U[i, nj - 1]
                UR = _reflect(UL, nx, nr)
            else:
                UL = U[i, j - 1]
                UR = U[i, j]
            f = _rusanov(UL, UR, nx, nr)
            s = m.sj[i, j]
            for c in range(4):
                if j > 0:
                    R[i, j - 1, c] += f[c] * s
                if j < nj:
                    R[i, j, c] -= f[c] * s

    # Axisymmetric pressure source in the radial momentum equation.
    for i in range(ni):
        for j in range(nj):
            R[i, j, 2] -= _primitives(U[i, j])[3] * m.area[i, j]
    return R


def local_time_step(U: np.ndarray, m: Metrics, cfg) -> np.ndarray:
    """Per-cell stable time step.

    Args:
        U (np.ndarray): Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics as numpy arrays.
        cfg: Object with a ``cfl`` attribute.

    Returns:
        np.ndarray: Time step, shape ``(NI, NJ)``.
    """
    ni, nj = m.vol.shape
    dt = np.zeros((ni, nj))
    for i in range(ni):
        for j in range(nj):
            rho, u, v, p = _primitives(U[i, j])
            a = math.sqrt(GAMMA * p / rho)
            total = 0.0
            for nx, nr, s in ((m.nix[i, j], m.nir[i, j], m.si[i, j]),
                              (m.nix[i + 1, j], m.nir[i + 1, j], m.si[i + 1, j]),
                              (m.njx[i, j], m.njr[i, j], m.sj[i, j]),
                              (m.njx[i, j + 1], m.njr[i, j + 1], m.sj[i, j + 1])):
                total += (abs(u * nx + v * nr) + a) * s
            dt[i, j] = cfg.cfl * m.vol[i, j] / total
    return dt


def step(U: np.ndarray, m: Metrics, cfg) -> Tuple[np.ndarray, float]:
    """One 4-stage Runge-Kutta step with a frozen local time step.

    Args:
        U (np.ndarray): Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics as numpy arrays.
        cfg: Boundary-condition / CFL configuration.

    Returns:
        tuple: ``(U_new, rms_density_residual)``.
    """
    ni, nj = m.vol.shape
    dt = local_time_step(U, m, cfg)

    U0 = U
    Uk = U
    res0 = None
    for alpha in RK_ALPHA:
        R = residual(Uk, m, cfg)
        if res0 is None:
            res0 = R
        Uk = np.empty_like(U0)
        for i in range(ni):
            for j in range(nj):
                scale = alpha * dt[i, j] / m.vol[i, j]
                for c in range(4):
                    Uk[i, j, c] = U0[i, j, c] - scale * R[i, j, c]

    rms = 0.0
    for i in range(ni):
        for j in range(nj):
            d = res0[i, j, 0] / m.vol[i, j]
            rms += d * d
    return Uk, math.sqrt(rms / (ni * nj))


def solve(metrics: Metrics,
          cfg,
          n_iter: int = 50,
          tol: float = 0.0,
          report_every: int = 10,
          U0: Optional[np.ndarray] = None,
          verbose: bool = False) -> Tuple[np.ndarray, List[Tuple[int, float]]]:
    """March in time - same signature as :func:`euler_solver_flatten.solve`.

    Args:
        metrics (Metrics): Mesh metrics as numpy arrays (see
            :func:`euler_metrics.to_numpy`).
        cfg: Boundary-condition / CFL configuration.
        n_iter (int, optional): Number of iterations. Defaults to 50.
        tol (float, optional): Relative residual stop criterion, ``0`` to run
            all ``n_iter``. Defaults to 0.0.
        report_every (int, optional): History interval. Defaults to 10.
        U0 (np.ndarray, optional): Starting state. Required - this solver is
            only ever used for the like-for-like comparison, so it takes its
            initial state from the vectorized solver.
        verbose (bool, optional): Print the residual as it goes. Defaults to False.

    Returns:
        tuple: ``(U, history)`` where ``history`` is a list of
        ``(iteration, relative_residual)``.

    Raises:
        ValueError: If ``U0`` is not supplied.
    """
    if U0 is None:
        raise ValueError("euler_solver_block.solve needs an explicit U0 "
                         "(use euler_solver_flatten.initial_condition)")
    U = np.array(U0, dtype=float)
    history: List[Tuple[int, float]] = []
    res0: Optional[float] = None

    for it in range(1, n_iter + 1):
        U, res = step(U, metrics, cfg)
        if it == 1 or it % report_every == 0 or it == n_iter:
            if res0 is None:
                res0 = max(res, 1e-300)
            rel = res / res0
            history.append((it, rel))
            if verbose:
                print(f"    iter {it:6d}   residual {rel:.3e}")
            if tol > 0.0 and rel < tol:
                break
    return U, history


if __name__ == "__main__":
    import time

    from duct_flatten import flatten_to_meridional
    from duct_geometry import duct_radius
    from duct_mesh import revolve_duct
    from euler_metrics import build_metrics
    from euler_solver_flatten import Config, initial_condition
    from euler_metrics import to_jax

    x, r_wall = duct_radius(41)
    block = revolve_duct(x, r_wall, n_radial=13, n_theta=7)
    metrics = build_metrics(*flatten_to_meridional(block))
    cfg = Config()
    U0 = np.asarray(initial_condition(to_jax(metrics), cfg), dtype=float)

    t = time.perf_counter()
    U, hist = solve(metrics, cfg, n_iter=10, report_every=5, U0=U0, verbose=True)
    elapsed = time.perf_counter() - t
    print(f"10 iterations on {metrics.vol.size} cells: {elapsed:.2f} s "
          f"({1000 * elapsed / 10:.1f} ms/iteration)")
