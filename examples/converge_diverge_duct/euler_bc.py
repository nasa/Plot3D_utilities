"""Stage 4c - Boundary conditions as ghost cells.

Every boundary is handled the same way: given the first interior cell, produce a
*ghost* cell just outside the domain such that the Rusanov flux between them is
the boundary flux we want.  The solver then treats boundary faces exactly like
interior faces, which keeps the residual loop free of special cases.

Three conditions are needed for this duct:

* :func:`inlet_ghost` - subsonic inflow at fixed stagnation state, using the
  outgoing Riemann invariant so the boundary does not over-specify the flow;
* :func:`outlet_ghost` - subsonic outflow at fixed back pressure;
* :func:`reflect_ghost` - a slip wall, reused unchanged on the axis.
"""
import jax.numpy as jnp

from euler_physics import GAMMA, R_GAS, conservative, primitives, sound_speed


def inlet_ghost(U_in, p0, t0):
    """Subsonic inflow ghost state from a fixed stagnation condition.

    Only one piece of information leaves the domain through a subsonic inlet:
    the left-running Riemann invariant ``R- = u - 2a/(gamma-1)``.  Combining it
    with the stagnation enthalpy implied by ``T0`` gives a quadratic for the
    boundary speed of sound::

        A*a_b**2 + B*a_b + C = 0
        A = (gamma+1)/(gamma-1),  B = 2*R-,  C = (gamma-1)/2 * R-**2 - a0**2

    Args:
        U_in: Conservative state of the first interior cell column, shape
            ``(NJ, 4)``.
        p0: Inlet stagnation pressure, Pa.
        t0: Inlet stagnation temperature, K.

    Returns:
        Ghost state of shape ``(NJ, 4)``, flow aligned with the axis.
    """
    rho, u, v, p = primitives(U_in)
    a = sound_speed(rho, p)

    a0_sq = GAMMA * R_GAS * t0
    r_minus = u - 2.0 * a / (GAMMA - 1.0)

    A = (GAMMA + 1.0) / (GAMMA - 1.0)
    B = 2.0 * r_minus
    C = 0.5 * (GAMMA - 1.0) * r_minus ** 2 - a0_sq
    disc = jnp.maximum(B ** 2 - 4.0 * A * C, 0.0)
    a_b = (-B + jnp.sqrt(disc)) / (2.0 * A)          # positive root

    u_b = r_minus + 2.0 * a_b / (GAMMA - 1.0)
    t_b = a_b ** 2 / (GAMMA * R_GAS)
    m_b = u_b / a_b
    p_b = p0 * (1.0 + 0.5 * (GAMMA - 1.0) * m_b ** 2) ** (-GAMMA / (GAMMA - 1.0))
    rho_b = p_b / (R_GAS * t_b)
    return conservative(rho_b, u_b, jnp.zeros_like(u_b), p_b)


def outlet_ghost(U_out, p_back):
    """Subsonic outflow ghost state at a fixed back pressure.

    Density and both velocity components are extrapolated from the interior;
    only pressure is imposed.  If the exit happens to go supersonic (which it
    does during startup, before the shock settles) nothing may be imposed at
    all, so pressure is extrapolated instead.

    Args:
        U_out: Conservative state of the last interior cell column, shape
            ``(NJ, 4)``.
        p_back: Static back pressure, Pa.

    Returns:
        Ghost state of shape ``(NJ, 4)``.
    """
    rho, u, v, p = primitives(U_out)
    supersonic = jnp.abs(u) / sound_speed(rho, p) >= 1.0
    p_b = jnp.where(supersonic, p, p_back)
    return conservative(rho, u, v, p_b)


def reflect_ghost(U_wall, nx, nr):
    """Slip-wall (and axis) ghost state: mirror the velocity about the face.

    The normal velocity is reversed and the tangential velocity is kept, so the
    flux through the face carries pressure only.  The formula is invariant to
    the sign of the normal, so the same function serves the outer wall
    (``j = NJ-1``) and the axis (``j = 0``).

    On the axis the face weight ``sj[:, 0]`` is exactly zero anyway, so this is
    belt-and-braces there rather than load bearing.

    Args:
        U_wall: Conservative state of the adjacent interior cells.
        nx: Face unit normal, axial component.
        nr: Face unit normal, radial component.

    Returns:
        Ghost state with the same shape as ``U_wall``.
    """
    rho, u, v, p = primitives(U_wall)
    qn = u * nx + v * nr
    return conservative(rho, u - 2.0 * qn * nx, v - 2.0 * qn * nr, p)


if __name__ == "__main__":
    import jax

    from euler_physics import mach

    # Same precision as the solver (euler_solver_flatten.py sets this too);
    # without it this standalone demo would silently run in float32.
    jax.config.update("jax_enable_x64", True)

    # A wall reflection must exactly cancel the normal velocity.
    U = conservative(jnp.array([1.2]), jnp.array([100.0]), jnp.array([30.0]),
                     jnp.array([1.0e5]))
    G = reflect_ghost(U, 0.0, 1.0)
    rho, u, v, p = primitives(0.5 * (U + G))
    print(f"mirrored normal velocity at face: {float(v[0]):.3e}  (want 0)")

    # A stagnant interior cell at the inlet stagnation state must stay stagnant.
    t0, p0 = 500.0, 2.0e5
    rho0 = p0 / (R_GAS * t0)
    U0 = conservative(jnp.array([rho0]), jnp.array([0.0]), jnp.array([0.0]),
                      jnp.array([p0]))
    print(f"inlet ghost Mach from rest     : {float(mach(inlet_ghost(U0, p0, t0))[0]):.3e}")

    # A slower interior cell must be accelerated by the stagnation condition.
    U1 = conservative(jnp.array([1.3]), jnp.array([80.0]), jnp.array([0.0]),
                      jnp.array([1.8e5]))
    print(f"inlet ghost Mach (moving cell) : {float(mach(inlet_ghost(U1, p0, t0))[0]):.4f}")
    print(f"outlet ghost p                 : "
          f"{float(primitives(outlet_ghost(U1, 1.0e5))[3][0]):.1f} Pa")
