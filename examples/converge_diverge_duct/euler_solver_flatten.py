"""Stage 4d - The vectorized (JAX) solver on the flattened mesh.

This is the solver you would actually use.  The whole flow field is a single
``(NI, NJ, 4)`` array and every operation - fluxes, boundary conditions, the
Runge-Kutta update - is a whole-array expression.  There is no per-cell Python
loop anywhere, so ``jax.jit`` can fuse the lot into one compiled kernel that
runs on CPU, GPU or TPU unchanged.

See :mod:`euler_solver_block` for the same physics written as nested Python
loops, and ``main.py`` for the side-by-side timing.

Scheme
------
* Finite volume, cell centred, on the axisymmetric metrics of
  :mod:`euler_metrics`.
* First-order Rusanov flux (:mod:`euler_physics`).
* Ghost-cell boundary conditions (:mod:`euler_bc`).
* 4-stage Runge-Kutta with ``alpha = (1/4, 1/3, 1/2, 1)`` and a **local** time
  step frozen across the stages - this is a steady-state solve, so each cell is
  allowed to march at its own maximum stable rate.

Notes on JAX
------------
* ``jax_enable_x64`` is switched on below.  In float32 the density residual
  stalls around 1e-4 and never reaches the tolerance; this solve genuinely
  needs float64.
* There is exactly one ``jax.jit``, applied to a single time step.  The
  iteration loop stays an ordinary Python ``for`` loop so it is readable.
  ``jax.lax.fori_loop`` / ``scan`` would remove the per-iteration dispatch
  overhead, at the cost of hiding the loop - not worth it for a tutorial.
* The residual is pulled back to the host only every ``report_every``
  iterations.  Every ``float(device_array)`` forces a device synchronization,
  and doing that each iteration is a classic way to throw away most of JAX's
  speed.
"""
from typing import List, NamedTuple, Optional, Tuple

import jax

jax.config.update("jax_enable_x64", True)  # must happen before any array work

import jax.numpy as jnp  # noqa: E402

from euler_bc import inlet_ghost, outlet_ghost, reflect_ghost  # noqa: E402
from euler_metrics import Metrics  # noqa: E402
from euler_physics import (GAMMA, R_GAS, conservative, mach_from_p_ratio,  # noqa: E402
                           primitives, rusanov, sound_speed, t_over_t0)


class Config(NamedTuple):
    """Boundary conditions and CFL number for a run.

    Like :class:`euler_metrics.Metrics`, a ``NamedTuple`` is already a JAX
    pytree, so this passes through ``jax.jit`` without any registration.
    """
    p0: float = 2.0e5      #: inlet stagnation pressure, Pa
    t0: float = 500.0      #: inlet stagnation temperature, K
    p_back: float = 1.0e5  #: outlet static pressure, Pa
    cfl: float = 1.5       #: CFL number (drop to ~0.5 if a case goes unstable)


def initial_condition(metrics: Metrics, cfg: Config):
    """Isentropic guess from a linear static-pressure ramp along the duct.

    Starting from a plausible field rather than a uniform one roughly halves the
    iteration count and keeps the startup transient mild.

    Args:
        metrics (Metrics): Mesh metrics (JAX arrays).
        cfg (Config): Boundary conditions.

    Returns:
        Conservative state of shape ``(NI, NJ, 4)``.
    """
    x = metrics.xc
    frac = (x - x.min()) / (x.max() - x.min())
    p = cfg.p0 + (cfg.p_back - cfg.p0) * frac
    m = mach_from_p_ratio(p, cfg.p0)
    t = cfg.t0 * t_over_t0(m)
    rho = p / (R_GAS * t)
    u = m * jnp.sqrt(GAMMA * R_GAS * t)
    return conservative(rho, u, jnp.zeros_like(u), p)


def _with_ghosts_i(U, cfg: Config):
    """Add inlet/outlet ghost columns, giving shape ``(NI+2, NJ, 4)``."""
    g_in = inlet_ghost(U[0], cfg.p0, cfg.t0)[None, :, :]
    g_out = outlet_ghost(U[-1], cfg.p_back)[None, :, :]
    return jnp.concatenate([g_in, U, g_out], axis=0)


def _with_ghosts_j(U, m: Metrics):
    """Add axis/wall ghost rows, giving shape ``(NI, NJ+2, 4)``."""
    g_axis = reflect_ghost(U[:, 0], m.njx[:, 0], m.njr[:, 0])[:, None, :]
    g_wall = reflect_ghost(U[:, -1], m.njx[:, -1], m.njr[:, -1])[:, None, :]
    return jnp.concatenate([g_axis, U, g_wall], axis=1)


def residual(U, m: Metrics, cfg: Config):
    """Net flux imbalance of every cell, ``sum(F.n S) - source``.

    At steady state this is zero.  ``dU/dt = -residual / vol``.

    Args:
        U: Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics.
        cfg (Config): Boundary conditions.

    Returns:
        Array of shape ``(NI, NJ, 4)``.
    """
    Ui = _with_ghosts_i(U, cfg)
    Uj = _with_ghosts_j(U, m)

    fi = rusanov(Ui[:-1], Ui[1:], m.nix, m.nir) * m.si[..., None]
    fj = rusanov(Uj[:, :-1], Uj[:, 1:], m.njx, m.njr) * m.sj[..., None]

    # Axisymmetric source: the r-momentum equation picks up p * dA because the
    # cell's own pressure acts on the (shrinking) circumferential faces.
    _, _, _, p = primitives(U)
    zero = jnp.zeros_like(p)
    source = jnp.stack([zero, zero, p * m.area, zero], axis=-1)

    return (fi[1:] - fi[:-1]) + (fj[:, 1:] - fj[:, :-1]) - source


def local_time_step(U, m: Metrics, cfg: Config):
    """Per-cell stable time step, ``CFL * vol / sum((|qn| + a) * S)``.

    Args:
        U: Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics.
        cfg (Config): Boundary conditions (uses ``cfl``).

    Returns:
        Array of shape ``(NI, NJ)``.
    """
    rho, u, v, p = primitives(U)
    a = sound_speed(rho, p)

    def spectral(nx, nr, s):
        return (jnp.abs(u * nx + v * nr) + a) * s

    total = (spectral(m.nix[:-1], m.nir[:-1], m.si[:-1])
             + spectral(m.nix[1:], m.nir[1:], m.si[1:])
             + spectral(m.njx[:, :-1], m.njr[:, :-1], m.sj[:, :-1])
             + spectral(m.njx[:, 1:], m.njr[:, 1:], m.sj[:, 1:]))
    return cfg.cfl * m.vol / total


def _step(U, m: Metrics, cfg: Config):
    """One 4-stage Runge-Kutta step with a frozen local time step.

    Args:
        U: Conservative state, shape ``(NI, NJ, 4)``.
        m (Metrics): Mesh metrics.
        cfg (Config): Boundary conditions.

    Returns:
        tuple: ``(U_new, rms_density_residual)``.
    """
    dt = local_time_step(U, m, cfg)
    scale = (dt / m.vol)[..., None]

    U0 = U
    res0 = None
    for alpha in (0.25, 1.0 / 3.0, 0.5, 1.0):
        res = residual(U, m, cfg)
        if res0 is None:
            res0 = res
        U = U0 - alpha * scale * res

    rms = jnp.sqrt(jnp.mean((res0[..., 0] / m.vol) ** 2))
    return U, rms


#: The one and only jit in this tutorial.  Re-traced whenever the mesh shape
#: changes, which ``main.py``'s coarse comparison mesh deliberately triggers.
step = jax.jit(_step)


def solve(metrics: Metrics,
          cfg: Config,
          n_iter: int = 20000,
          tol: float = 1e-6,
          report_every: int = 100,
          U0=None,
          verbose: bool = False) -> Tuple[jnp.ndarray, List[Tuple[int, float]]]:
    """March to steady state.

    Args:
        metrics (Metrics): Mesh metrics as JAX arrays (see
            :func:`euler_metrics.to_jax`).
        cfg (Config): Boundary conditions and CFL.
        n_iter (int, optional): Maximum iterations. Defaults to 20000.
        tol (float, optional): Stop when the density residual has dropped by
            this factor from its initial value. ``0.0`` disables the check.
            Defaults to 1e-6.
        report_every (int, optional): Host sync / history interval. Defaults to 100.
        U0 (optional): Starting state. Defaults to
            :func:`initial_condition`.
        verbose (bool, optional): Print the residual as it goes. Defaults to False.

    Returns:
        tuple: ``(U, history)`` where ``history`` is a list of
        ``(iteration, relative_residual)``.

    Note:
        The residual typically plateaus somewhere between 1e-6 and 1e-8 once the
        shock has stopped moving to within a cell width.  That is normal
        first-order shock-capturing behaviour, not a bug.
    """
    U = initial_condition(metrics, cfg) if U0 is None else jnp.asarray(U0)
    history: List[Tuple[int, float]] = []
    res0: Optional[float] = None

    for it in range(1, n_iter + 1):
        U, res = step(U, metrics, cfg)
        if it == 1 or it % report_every == 0 or it == n_iter:
            r = float(res)
            if res0 is None:
                res0 = max(r, 1e-300)
            rel = r / res0
            history.append((it, rel))
            if verbose:
                print(f"    iter {it:6d}   residual {rel:.3e}")
            if tol > 0.0 and rel < tol:
                break
    return U, history


# --------------------------------------------------------------------------
# Unit check: uniform flow down a straight pipe must not move at all
# --------------------------------------------------------------------------

def straight_pipe_residual(ni: int = 41, nj: int = 13,
                           mach_number: float = 0.3) -> Tuple[float, float]:
    """Residual of an exactly uniform flow in a constant-radius pipe.

    This is the single most valuable check in the whole solver.  A uniform state
    in a straight pipe is an exact solution, so the residual must be machine
    zero.  Getting that requires the axisymmetric source term, every face
    normal orientation, and the wall reflection to all be right at once - one
    sign error anywhere and this check fails loudly.

    Args:
        ni (int, optional): Axial cells + 1. Defaults to 41.
        nj (int, optional): Radial cells + 1. Defaults to 13.
        mach_number (float, optional): Uniform Mach number. Defaults to 0.3.

    Returns:
        tuple: ``(max_abs_residual, reference_flux_scale)``.  Their ratio is the
        number that should be near machine epsilon.
    """
    import numpy as np

    from euler_metrics import build_metrics, to_jax

    x = np.linspace(0.0, 1.0, ni)
    r = np.linspace(0.0, 0.1, nj)
    metrics = to_jax(build_metrics(*np.meshgrid(x, r, indexing="ij")))

    # A uniform state, and the stagnation/back conditions consistent with it, so
    # the boundary conditions reproduce the interior exactly.
    t, p = 300.0, 1.0e5
    rho = p / (R_GAS * t)
    u = mach_number * float(jnp.sqrt(GAMMA * R_GAS * t))
    t0 = t * (1.0 + 0.5 * (GAMMA - 1.0) * mach_number ** 2)
    p0 = p * (t0 / t) ** (GAMMA / (GAMMA - 1.0))
    cfg = Config(p0=p0, t0=t0, p_back=p)

    shape = metrics.vol.shape
    U = conservative(jnp.full(shape, rho), jnp.full(shape, u),
                     jnp.zeros(shape), jnp.full(shape, p))

    res = residual(U, metrics, cfg)
    scale = float(jnp.max(jnp.abs(rusanov(U, U, 1.0, 0.0))) * jnp.max(metrics.si))
    return float(jnp.max(jnp.abs(res))), scale


if __name__ == "__main__":
    res, scale = straight_pipe_residual()
    print("straight pipe, uniform flow")
    print(f"  max |residual| : {res:.3e}")
    print(f"  flux scale     : {scale:.3e}")
    print(f"  ratio          : {res / scale:.3e}  (want ~1e-16)")
