"""Stage 4b - Physics: gas model, state conversions and the numerical flux.

Everything here is written with ``jax.numpy`` and is shape agnostic: the same
functions work on a single state, a row of states, or the whole ``(NI, NJ, 4)``
field.  That shape-agnosticism is what makes the vectorized solver possible.

The conservative state vector, with the last axis holding the four components::

    U = [rho, rho*u, rho*v, rho*E]

for the 2D axisymmetric Euler equations, closed with a calorically perfect gas::

    p = (gamma - 1) * (rho*E - 0.5*rho*(u**2 + v**2))
    a = sqrt(gamma * p / rho)
"""
import jax.numpy as jnp

#: Ratio of specific heats for air.
GAMMA = 1.4
#: Specific heat at constant pressure, J/(kg K).
CP = 1005.0
#: Specific gas constant, J/(kg K).  R = Cp * (gamma - 1) / gamma.
R_GAS = CP * (GAMMA - 1.0) / GAMMA

#: Floors applied inside :func:`primitives`.  These are a NaN guard for bad
#: startup transients only - they must not be active at convergence.
RHO_FLOOR = 1.0e-4
P_FLOOR = 10.0


def conservative(rho, u, v, p):
    """Pack primitive variables into the conservative state vector.

    Args:
        rho: Density.
        u: Axial velocity.
        v: Radial velocity.
        p: Static pressure.

    Returns:
        Array with a trailing axis of size 4: ``[rho, rho*u, rho*v, rho*E]``.
    """
    rho_e = p / (GAMMA - 1.0) + 0.5 * rho * (u ** 2 + v ** 2)
    return jnp.stack([rho, rho * u, rho * v, rho_e], axis=-1)


def primitives(U):
    """Unpack the conservative state vector into primitive variables.

    Args:
        U: Conservative state with a trailing axis of size 4.

    Returns:
        tuple: ``(rho, u, v, p)``, each with ``U``'s leading shape.
    """
    rho = jnp.maximum(U[..., 0], RHO_FLOOR)
    u = U[..., 1] / rho
    v = U[..., 2] / rho
    p = jnp.maximum((GAMMA - 1.0) * (U[..., 3] - 0.5 * rho * (u ** 2 + v ** 2)),
                    P_FLOOR)
    return rho, u, v, p


def sound_speed(rho, p):
    """Speed of sound.

    Args:
        rho: Density.
        p: Static pressure.

    Returns:
        Speed of sound ``sqrt(gamma p / rho)``.
    """
    return jnp.sqrt(GAMMA * p / rho)


def mach(U):
    """Local Mach number of a conservative state.

    Args:
        U: Conservative state with a trailing axis of size 4.

    Returns:
        Mach number ``sqrt(u**2 + v**2) / a``.
    """
    rho, u, v, p = primitives(U)
    return jnp.sqrt(u ** 2 + v ** 2) / sound_speed(rho, p)


def flux(U, nx, nr):
    """Physical Euler flux projected onto a face normal.

    Args:
        U: Conservative state with a trailing axis of size 4.
        nx: Unit normal, axial component.
        nr: Unit normal, radial component.

    Returns:
        Array with a trailing axis of size 4, the flux ``F*nx + G*nr``.
    """
    rho, u, v, p = primitives(U)
    qn = u * nx + v * nr
    return jnp.stack([rho * qn,
                      rho * u * qn + p * nx,
                      rho * v * qn + p * nr,
                      (U[..., 3] + p) * qn], axis=-1)


def rusanov(UL, UR, nx, nr):
    """Rusanov (local Lax-Friedrichs) numerical flux.

    First order and slightly dissipative, but it has no tunable constants, needs
    only one ghost layer, and stays robust through the strong shock this duct
    produces - the right trade for a tutorial.

    Args:
        UL: State on the -side of the face.
        UR: State on the +side of the face.
        nx: Unit normal, axial component (pointing from L to R).
        nr: Unit normal, radial component.

    Returns:
        Numerical flux with a trailing axis of size 4.
    """
    rhoL, uL, vL, pL = primitives(UL)
    rhoR, uR, vR, pR = primitives(UR)
    qnL = uL * nx + vL * nr
    qnR = uR * nx + vR * nr
    smax = jnp.maximum(jnp.abs(qnL) + sound_speed(rhoL, pL),
                       jnp.abs(qnR) + sound_speed(rhoR, pR))
    return 0.5 * (flux(UL, nx, nr) + flux(UR, nx, nr)) \
        - 0.5 * smax[..., None] * (UR - UL)


# --------------------------------------------------------------------------
# Isentropic relations (used by the inlet BC and the initial condition)
# --------------------------------------------------------------------------

def p_over_p0(m):
    """Static-to-stagnation pressure ratio at Mach number ``m``.

    Args:
        m: Mach number.

    Returns:
        ``p/p0``.
    """
    return (1.0 + 0.5 * (GAMMA - 1.0) * m ** 2) ** (-GAMMA / (GAMMA - 1.0))


def t_over_t0(m):
    """Static-to-stagnation temperature ratio at Mach number ``m``.

    Args:
        m: Mach number.

    Returns:
        ``T/T0``.
    """
    return 1.0 / (1.0 + 0.5 * (GAMMA - 1.0) * m ** 2)


def mach_from_p_ratio(p, p0):
    """Isentropic Mach number implied by a static/stagnation pressure pair.

    Args:
        p: Static pressure.
        p0: Stagnation pressure.

    Returns:
        Mach number (0 where ``p >= p0``).
    """
    ratio = jnp.maximum((p0 / p) ** ((GAMMA - 1.0) / GAMMA) - 1.0, 0.0)
    return jnp.sqrt(2.0 / (GAMMA - 1.0) * ratio)


if __name__ == "__main__":
    import jax

    # Same precision as the solver (euler_solver_flatten.py sets this too);
    # without it this standalone demo would silently run in float32.
    jax.config.update("jax_enable_x64", True)

    print(f"gamma = {GAMMA}, Cp = {CP}, R_gas = {R_GAS:.4f} J/(kg K)")

    # Round trip: primitives -> conservative -> primitives.
    U = conservative(1.2, 150.0, 5.0, 1.0e5)
    print("round trip     :", [float(q) for q in primitives(U)])
    print("Mach           :", float(mach(U)))

    # A uniform state must produce the exact physical flux (no dissipation).
    n = (1.0, 0.0)
    f_exact = flux(U, *n)
    f_rus = rusanov(U, U, *n)
    print("uniform flux err:", float(jnp.max(jnp.abs(f_exact - f_rus))))
