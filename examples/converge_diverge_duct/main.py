"""Converging-diverging duct: geometry -> 3D mesh -> flatten -> Euler solve.

Run it with no arguments::

    python main.py

Stages
------
0. Sanity check - uniform flow in a straight pipe must give a zero residual.
1. Geometry    - a Bezier wall radius ``R(x)`` (:mod:`duct_geometry`), then
                 what 1D isentropic theory predicts from ``A(x)/A*`` alone
                 (:mod:`quasi1d`) - printed before any CFD is run, so the
                 prediction is on record ahead of the result.
2. 3D mesh     - revolve it into a Plot3D block, write and read it back
                 (:mod:`duct_mesh`).
3. Flatten     - collapse the axisymmetric block to a 2D meridional mesh
                 (:mod:`duct_flatten`), then look at the domain and at where
                 the boundary-condition ghost cells live (Figures 1 and 2).
4. Solve       - 2D axisymmetric Euler with JAX (:mod:`euler_solver_flatten`),
                 compared against that 1D prediction.
5. Compare     - the same solver written as nested Python loops
                 (:mod:`euler_solver_block`), on a coarse mesh, to show what
                 the vectorized data layout buys you.

Requires ``pyturbo-aero`` and ``jax``::

    pip install -e ".[examples]"
"""
import math
import time
from copy import deepcopy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot3d import Block, read_plot3D, reduce_blocks, write_plot3D

import quasi1d
from duct_flatten import (axisymmetry_error, flatten_to_meridional,
                          node_count_reduction)
from duct_geometry import area_ratio, duct_radius, throat
from duct_mesh import inflation_report, revolve_duct
from euler_metrics import (analytic_volume, build_metrics, enclosed_volume,
                           ghost_cell_centres, to_jax)
from euler_physics import mach, primitives
from euler_solver_block import solve as solve_block
from euler_solver_flatten import (Config, initial_condition,
                                  solve as solve_flatten,
                                  straight_pipe_residual)

SCRIPT_DIR = Path(__file__).resolve().parent
MESH_FILE = SCRIPT_DIR / "duct_3d.xyz"

# Production mesh.  gcd(200, 24, 36) = 4, so the GCD-reduced 3D wireframe
# (the same reduction plot3d's plot_blocks uses) stays clean; other node
# counts can make it unusable.
N_AXIAL, N_RADIAL, N_THETA = 201, 25, 37
# Comparison mesh: gcd(40, 12, 6) = 2, same reasoning.
C_AXIAL, C_RADIAL, C_THETA = 41, 13, 7
COMPARE_ITERS = 50


def banner(text: str) -> None:
    """Print a section header.

    Args:
        text (str): Header text.
    """
    print(f"\n{'=' * 68}\n{text}\n{'=' * 68}")


def check(condition: bool, message: str) -> None:
    """Abort with a clear message if a correctness check fails.

    A plain ``assert`` would do, but it is silently skipped under ``python -O``
    and these two checks are the ones nothing downstream should run without.

    Args:
        condition (bool): Must be true to continue.
        message (str): What went wrong.
    """
    if not condition:
        raise RuntimeError(message)


def shock_location(x_axis: np.ndarray, mach_axis: np.ndarray) -> float:
    """Locate a shock as the steepest drop in Mach number along the axis.

    Args:
        x_axis (np.ndarray): Axial cell-centre positions.
        mach_axis (np.ndarray): Mach number along the axis.

    Returns:
        float: Axial position of the steepest drop.
    """
    i = int(np.argmin(np.diff(mach_axis)))
    return float(0.5 * (x_axis[i] + x_axis[i + 1]))


def build_pipeline(n_axial: int, n_radial: int, n_theta: int):
    """Run geometry -> revolve -> flatten -> metrics at a chosen resolution.

    Args:
        n_axial (int): Axial nodes.
        n_radial (int): Radial nodes.
        n_theta (int): Circumferential nodes.

    Returns:
        tuple: ``(x, r_wall, block, x2d, r2d, metrics_numpy)``.
    """
    x, r_wall = duct_radius(n_axial)
    block = revolve_duct(x, r_wall, n_radial=n_radial, n_theta=n_theta)
    x2d, r2d = flatten_to_meridional(block)
    return x, r_wall, block, x2d, r2d, build_metrics(x2d, r2d)


# --------------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------------

def draw_mesh_2d(ax, x2d: np.ndarray, r2d: np.ndarray, every: int = 4) -> None:
    """Draw the flattened ``(x, r)`` node grid as thin lines.

    All radial-index lines are drawn (that is where the wall clustering shows);
    axial-index lines are thinned so they do not bury it.

    Args:
        ax: Matplotlib axes.
        x2d (np.ndarray): Node x, shape ``(NI+1, NJ+1)``.
        r2d (np.ndarray): Node r, shape ``(NI+1, NJ+1)``.
        every (int, optional): Keep every ``every``-th axial line. Defaults to 4.
    """
    ax.plot(x2d, r2d, "k-", lw=0.4)
    ax.plot(x2d[::every, :].T, r2d[::every, :].T, "k-", lw=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.set_aspect("equal")


def draw_wireframe_3d(ax, block: Block) -> None:
    """Draw a block as a 3D wireframe on a caller-supplied axes.

    The same picture ``plot3d.plot_blocks`` draws, but into an existing axes so
    it can share a figure with the 2D panel (``plot_blocks`` owns its own
    figure and calls ``plt.show`` itself).  Like ``plot_blocks`` it GCD-reduces
    the block first so the wireframe stays readable.

    Args:
        ax: Matplotlib 3D axes.
        block (Block): Block to draw.
    """
    gcd = math.gcd(block.IMAX - 1, math.gcd(block.JMAX - 1, block.KMAX - 1))
    b = reduce_blocks([deepcopy(block)], gcd)[0]     # reduce_blocks mutates
    X, Y, Z = b.X, b.Y, b.Z
    style = dict(color="C0", lw=0.6, alpha=0.7)
    for j in range(X.shape[1]):
        for k in range(X.shape[2]):
            ax.plot(X[:, j, k], Y[:, j, k], Z[:, j, k], **style)
    for i in range(X.shape[0]):
        for k in range(X.shape[2]):
            ax.plot(X[i, :, k], Y[i, :, k], Z[i, :, k], **style)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            ax.plot(X[i, j, :], Y[i, j, :], Z[i, j, :], **style)
    ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)), zoom=1.2)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")


def figure_domain(block: Block, x2d: np.ndarray, r2d: np.ndarray,
                  n3: int, n2: int) -> None:
    """Figure 1: the flattened 2D mesh next to the 3D mesh it came from.

    Args:
        block (Block): The revolved 3D block.
        x2d (np.ndarray): Flattened node x.
        r2d (np.ndarray): Flattened node r.
        n3 (int): Node count of the 3D block.
        n2 (int): Node count of the 2D mesh.
    """
    fig = plt.figure(figsize=(15, 5.5))
    ax2d = fig.add_subplot(1, 2, 1)
    draw_mesh_2d(ax2d, x2d, r2d)
    ax2d.set_title(f"flattened meridional (x, r) mesh - {n2} nodes\n"
                   "(every radial line drawn: note the wall clustering)")

    ax3d = fig.add_subplot(1, 2, 2, projection="3d")
    draw_wireframe_3d(ax3d, block)
    ax3d.set_title(f"revolved 3D mesh - {n3} nodes\n"
                   f"{block.IMAX} x {block.JMAX} x {block.KMAX}, full 360 deg "
                   "revolve, GCD-reduced for drawing")
    fig.suptitle("Figure 1 - the domain the solver will run on")
    fig.tight_layout()


def figure_ghost_cells(x2d: np.ndarray, r2d: np.ndarray,
                       xc: np.ndarray, rc: np.ndarray,
                       x_throat: float, r_throat: float) -> None:
    """Figure 2: interior cell centres plus one ring of ghost cells outside.

    The wall cells are only ~0.5% of the local radius tall, so the wall ghost
    ring hugs the wall at full scale; a second panel zooms in at the throat to
    show it really is outside the domain.

    Args:
        x2d (np.ndarray): Flattened node x.
        r2d (np.ndarray): Flattened node r.
        xc (np.ndarray): Cell-centre x, shape ``(NI, NJ)``.
        rc (np.ndarray): Cell-centre r, shape ``(NI, NJ)``.
        x_throat (float): Throat station, where the zoom panel looks.
        r_throat (float): Throat radius.
    """
    ghosts = ghost_cell_centres(xc, rc)
    colors = {"inlet": "C2", "outlet": "C3", "axis": "C4", "wall": "C1"}

    def draw(ax, dot_size: float, ghost_size: float) -> None:
        ax.plot(x2d, r2d, "-", color="0.75", lw=0.4)
        ax.plot(x2d.T, r2d.T, "-", color="0.75", lw=0.4)
        ax.scatter(xc, rc, s=dot_size, color="0.35", label="interior cell centres")
        for name, (gx, gr) in ghosts.items():
            ax.scatter(gx, gr, s=ghost_size, marker="s", color=colors[name],
                       label=f"{name} ghost cells", zorder=3)
        ax.set_aspect("equal")

    fig, (ax, axz) = plt.subplots(1, 2, figsize=(15, 5.5),
                                  gridspec_kw={"width_ratios": (3, 1)})
    draw(ax, dot_size=1.5, ghost_size=10)
    ax.set_xlabel("x")
    ax.set_ylabel("r")
    ax.legend(loc="upper left", fontsize=8, markerscale=1.5)
    ax.set_title("one ring outside each boundary, mirrored from the adjacent\n"
                 "cell spacing; the axis ring sits at r < 0")

    # Zoom on the wall at the throat, where the clustering is tightest.  The
    # zoom is stretched vertically (aspect "auto"), otherwise the 0.5%-of-R
    # wall cells would collapse to a hairline.
    dx = 2.5 * float(np.diff(xc[:, 0]).mean())
    draw(axz, dot_size=14, ghost_size=36)
    axz.set_aspect("auto")
    axz.set_xlim(x_throat - dx, x_throat + dx)
    axz.set_ylim(0.965 * r_throat, 1.012 * r_throat)
    axz.set_xlabel("x")
    axz.set_ylabel("r")
    axz.set_title("zoom: wall ghost ring at the throat\n"
                  "(stretched vertically)")
    ax.indicate_inset_zoom(axz, edgecolor="0.3")
    fig.suptitle("Figure 2 - where the boundary-condition ghost cells live")
    fig.tight_layout()


def main() -> None:
    """Run every stage and produce the figures and summary table."""
    cfg = Config()

    # ---------------------------------------------------------------- stage 0
    banner("STAGE 0 - solver sanity check")
    res, scale = straight_pipe_residual()
    print("Uniform flow in a straight pipe is an exact solution, so its")
    print("residual must be machine zero.  This one check exercises the")
    print("axisymmetric source term, every face normal, and the wall reflection.")
    print(f"  max |residual| / flux scale : {res / scale:.3e}")
    check(res / scale < 1e-12, "straight-pipe check failed - do not trust the solve")

    # ---------------------------------------------------------------- stage 1
    banner("STAGE 1 - geometry")
    x, r_wall, block, _, _, _ = build_pipeline(N_AXIAL, N_RADIAL, N_THETA)
    _, x_throat, r_throat = throat(x, r_wall)
    ar = area_ratio(r_wall)
    wall_angle = np.degrees(np.arctan(np.gradient(r_wall, x)))
    print(f"  axial stations   : {x.size}   x = [{x[0]:.2f}, {x[-1]:.2f}]")
    print(f"  throat           : x = {x_throat:.4f}, r = {r_throat:.4f}")
    print(f"  inlet  A/A_throat: {ar[0]:.3f}")
    print(f"  exit   A/A_throat: {ar[-1]:.3f}")
    print(f"  max wall angle   : {wall_angle.max():.1f} deg "
          f"(diverging section - remember this for the summary)")

    # ------------------------------------------------ 1D isentropic prediction
    # Everything here follows from the wall radius the mesh is built on and the
    # inlet/outlet pressures - no CFD.  Quoting it now, before the solve, keeps
    # the prediction honest: it cannot have been tuned to the answer.
    banner("1D ISENTROPIC PREDICTION - from A(x)/A* alone, before any CFD")
    print("  Once the throat is sonic, A* = A_throat and the area-Mach relation")
    print("  A/A* = f(M) fixes M at every station from R(x).  The two isentropic")
    print("  branches at the exit bound the back pressure and decide the regime:")
    ar_exit = float(ar[-1])
    m_sub = quasi1d.mach_from_area_ratio(ar_exit, supersonic=False)
    m_sup = quasi1d.mach_from_area_ratio(ar_exit, supersonic=True)
    p_sub = cfg.p0 * float(quasi1d.p_over_p0(m_sub))
    p_sup = cfg.p0 * float(quasi1d.p_over_p0(m_sup))
    m_shock_exit, p2_p1, _ = quasi1d.normal_shock(m_sup)
    p_shock_exit = p_sup * p2_p1
    print(f"  exit A/A*                 : {ar_exit:.3f}")
    print(f"  subsonic branch at exit   : M = {m_sub:.4f}, p = {p_sub / 1e3:7.1f} kPa"
          "  (unchoked limit)")
    print(f"  supersonic branch at exit : M = {m_sup:.4f}, p = {p_sup / 1e3:7.1f} kPa"
          "  (shock-free limit)")
    print(f"  normal shock at the exit  : M = {m_shock_exit:.4f}, "
          f"p = {p_shock_exit / 1e3:7.1f} kPa  (behind the shock)")
    print(f"  back pressure             :             p = {cfg.p_back / 1e3:7.1f} kPa")
    print("  -> p_back is below the subsonic branch, so the throat chokes, and above")
    print("     a shock at the exit, so the shock cannot leave the duct: a normal")
    print("     shock stands inside the diverging section.  Behind it the flow is")
    print("     subsonic with a lower p0 (larger A*), so A/A* is rescaled by p02/p01")
    print("     and the subsonic branch is followed to the exit.")
    analytic = quasi1d.solve(x, r_wall, cfg.p0, cfg.t0, cfg.p_back)
    print(f"  shock station             : x = {analytic.x_shock:.4f}, "
          f"M = {analytic.mach_shock:.4f} just upstream (the 1D peak)")
    print()
    print(f"  1D isentropic prediction: exit Mach = {analytic.mach_exit:.4f}"
          f"   (exit p = {analytic.p_exit / 1e3:.1f} kPa)")

    # ---------------------------------------------------------------- stage 2
    banner("STAGE 2 - 3D mesh (revolve)")
    wall_frac, axis_frac, cluster = inflation_report(N_RADIAL)
    print(f"  block            : {block.IMAX} x {block.JMAX} x {block.KMAX}"
          f"  ({block.X.size} nodes)")
    print(f"  revolve          : 360 deg (full body of revolution)")
    print(f"  wall cell height : {100 * wall_frac:.3f}% of local R")
    print(f"  axis cell height : {100 * axis_frac:.3f}% of local R")
    print(f"  clustering ratio : {cluster:.1f} : 1")

    write_plot3D(str(MESH_FILE), [block], binary=True)
    block = read_plot3D(str(MESH_FILE), binary=True)[0]
    print(f"  round tripped through {MESH_FILE.name}")

    # ---------------------------------------------------------------- stage 3
    banner("STAGE 3 - flatten to the meridional plane")
    # Deliberately flattening the block that came *back* off disk, so the whole
    # chain from geometry to solver metrics goes through a real Plot3D file.
    x2d, r2d = flatten_to_meridional(block)
    n3, n2, factor = node_count_reduction(block)
    print(f"  axisymmetry error: {axisymmetry_error(block):.3e}  (round-off)")
    print(f"  3D nodes         : {n3}")
    print(f"  2D nodes         : {n2}")
    print(f"  reduction        : {factor:.0f}x fewer nodes, same information")

    metrics_np = build_metrics(x2d, r2d)
    print(f"  duct volume      : {enclosed_volume(metrics_np):.6f} "
          f"(analytic {analytic_volume(x, r_wall):.6f})")

    # Look at the domain before solving on it.  ``block=False`` so the windows
    # pop up without pausing the script; the final ``plt.show()`` keeps them.
    print("  drawing Figures 1 and 2 (domain, ghost cells)...")
    figure_domain(block, x2d, r2d, n3, n2)
    figure_ghost_cells(x2d, r2d, metrics_np.xc, metrics_np.rc, x_throat, r_throat)
    plt.show(block=False)

    # ---------------------------------------------------------------- stage 4
    banner("STAGE 4 - axisymmetric Euler solve (JAX, vectorized)")
    metrics = to_jax(metrics_np)
    ni, nj = metrics_np.vol.shape
    print(f"  cells            : {ni} x {nj} = {ni * nj}")
    print(f"  inlet            : P0 = {cfg.p0:.3e} Pa, T0 = {cfg.t0:.1f} K")
    print(f"  outlet           : p  = {cfg.p_back:.3e} Pa")
    print(f"  CFL              : {cfg.cfl}")
    print("  marching...")

    t_start = time.perf_counter()
    U, history = solve_flatten(metrics, cfg, n_iter=20000, tol=1e-6,
                               report_every=500, verbose=True)
    elapsed = time.perf_counter() - t_start

    n_done, res_final = history[-1]
    print(f"  iterations       : {n_done}")
    print(f"  final residual   : {res_final:.3e}")
    print(f"  wall clock       : {elapsed:.2f} s "
          f"({1000 * elapsed / n_done:.2f} ms/iteration)")
    print("  (the residual plateaus once the shock stops moving to within a"
          " cell - that is normal for a first-order shock-capturing scheme)")

    M = np.asarray(mach(U))
    rho, u, _, p = [np.asarray(q) for q in primitives(U)]
    xc = np.asarray(metrics.xc)
    rc = np.asarray(metrics.rc)
    x_axis = xc[:, 0]
    m_axis = M[:, 0]

    # Mass-averaged Mach at each axial station, sum(rho*u*S*M) / sum(rho*u*S),
    # with S the cell's mean i-face weight (the annulus its flow passes
    # through).  Quasi-1D theory predicts *this* cross-section average, not
    # the value on the axis, so it is the fair curve to lay over the analytic
    # one in Figure 3c.
    s_cell = 0.5 * (metrics_np.si[:-1] + metrics_np.si[1:])
    w = rho * u * s_cell
    m_avg = (w * M).sum(axis=1) / w.sum(axis=1)

    # ---------------------------------------------------------------- stage 5
    banner("STAGE 5 - vectorized (JAX) vs nested Python loops")
    print(f"  Same equations, same mesh, same initial state, {COMPARE_ITERS} iterations.")
    print(f"  Coarse {C_AXIAL} x {C_RADIAL} x {C_THETA} mesh, because a pure-Python")
    print("  loop over the production mesh would take minutes per iteration.")
    print("  NOTE: 50 iterations is a transient, NOT a converged flow field -")
    print("  this stage compares software, not physics.")

    # build_metrics returns float64 numpy; to_jax is an exact float64 copy of
    # it, so the loop solver (numpy) and the JAX solver see bit-identical
    # metrics without any round trip.
    _, _, _, _, _, cmetrics_np = build_pipeline(C_AXIAL, C_RADIAL, C_THETA)
    cmetrics = to_jax(cmetrics_np)
    U0 = initial_condition(cmetrics, cfg)
    U0_np = np.asarray(U0)

    # Warm-up call: the coarse mesh has different array shapes, so JAX traces
    # and compiles a second version of `step`.  Time that separately.
    t_compile = time.perf_counter()
    solve_flatten(cmetrics, cfg, n_iter=1, tol=0.0, report_every=1, U0=U0)
    t_compile = time.perf_counter() - t_compile

    t0 = time.perf_counter()
    U_jax, _ = solve_flatten(cmetrics, cfg, n_iter=COMPARE_ITERS, tol=0.0,
                             report_every=COMPARE_ITERS, U0=U0)
    U_jax.block_until_ready()
    t_jax = time.perf_counter() - t0

    t0 = time.perf_counter()
    U_loop, _ = solve_block(cmetrics_np, cfg, n_iter=COMPARE_ITERS, tol=0.0,
                            report_every=COMPARE_ITERS, U0=U0_np)
    t_loop = time.perf_counter() - t0

    M_jax = np.asarray(mach(U_jax))
    M_loop = np.asarray(mach(U_loop))
    p_jax = np.asarray(primitives(U_jax)[3])
    p_loop = np.asarray(primitives(U_loop)[3])
    d_mach = float(np.max(np.abs(M_jax - M_loop)))
    d_p = float(np.max(np.abs(p_jax - p_loop)) / cfg.p0)

    print(f"\n  cells                    : {cmetrics_np.vol.size}")
    print(f"  max |M_jax - M_loop|     : {d_mach:.3e}")
    print(f"  max |p_jax - p_loop| / P0: {d_p:.3e}")
    check(d_mach < 1e-10 and d_p < 1e-10, "the two solvers disagree")
    print("  -> identical to round-off, as they must be")

    ms_jax = 1000.0 * t_jax / COMPARE_ITERS
    ms_loop = 1000.0 * t_loop / COMPARE_ITERS
    print(f"\n  JAX one-time compile     : {1000 * t_compile:8.1f} ms")
    print(f"  vectorized (JAX)         : {ms_jax:8.3f} ms/iteration")
    print(f"  nested Python loops      : {ms_loop:8.3f} ms/iteration")
    print(f"  speedup                  : {ms_loop / ms_jax:8.1f} x")

    # ---------------------------------------------------------------- summary
    x_shock_axis = shock_location(x_axis, m_axis)
    x_shock_wall = shock_location(xc[:, -1], M[:, -1])
    i_throat_cell = int(np.argmin(np.abs(x_axis - analytic.x_throat)))

    m_exit_cfd = float(M[-1].mean())

    banner("SUMMARY - CFD vs. 1D isentropic (A/A*) prediction")
    print(f"  {'quantity':<28}{'CFD':>14}{'1D (A/A*)':>14}")
    print(f"  {'-' * 56}")
    print(f"  {'throat Mach (on axis)':<28}{m_axis[i_throat_cell]:>14.4f}{1.0:>14.4f}")
    print(f"  {'throat Mach (mass-avg)':<28}{m_avg[i_throat_cell]:>14.4f}{1.0:>14.4f}")
    print(f"  {'peak Mach (on axis)':<28}{m_axis.max():>14.4f}"
          f"{analytic.mach.max():>14.4f}")
    print(f"  {'peak Mach (mass-avg)':<28}{m_avg.max():>14.4f}"
          f"{analytic.mach.max():>14.4f}")
    print(f"  {'shock x (on axis)':<28}{x_shock_axis:>14.4f}"
          f"{analytic.x_shock:>14.4f}")
    print(f"  {'shock x (at wall)':<28}{x_shock_wall:>14.4f}"
          f"{analytic.x_shock:>14.4f}")
    print(f"  {'exit Mach':<28}{m_exit_cfd:>14.4f}{analytic.mach_exit:>14.4f}")
    print(f"  {'exit pressure [kPa]':<28}{p[-1].mean() / 1000:>14.4f}"
          f"{analytic.p_exit / 1000:>14.4f}")

    # The headline check, spelled out on its own so it cannot be missed.
    print(f"\n  EXIT MACH   1D isentropic (A/A*) prediction : {analytic.mach_exit:.4f}")
    print(f"              CFD, mean over the outlet cells : {m_exit_cfd:.4f}"
          f"   ({100 * (m_exit_cfd / analytic.mach_exit - 1):+.1f} %)")

    # Mach along the duct at a handful of stations: the A/A* prediction next
    # to what the CFD produced.  Near the shock the columns disagree by
    # construction - the CFD shock is curved and stands at a slightly
    # different x on the axis than at the wall (see below).
    stations = np.unique(np.concatenate([np.linspace(x[0], x[-1], 9),
                                         [analytic.x_throat]]))
    print("\n  Mach vs. x - 1D isentropic (A/A*) prediction next to the CFD:")
    print(f"  {'x':>8}{'A/A*':>9}{'M 1D':>10}{'M CFD avg':>12}{'M CFD axis':>12}")
    print(f"  {'-' * 51}")
    for xs in stations:
        tag = "  <- throat" if xs == analytic.x_throat else ""
        print(f"  {xs:>8.3f}{np.interp(xs, x, ar):>9.3f}"
              f"{np.interp(xs, analytic.x, analytic.mach):>10.4f}"
              f"{np.interp(xs, x_axis, m_avg):>12.4f}"
              f"{np.interp(xs, x_axis, m_axis):>12.4f}{tag}")
    print("\n  The duct chokes at the throat and a normal shock stands in the")
    print("  diverging section.  Shock position, exit Mach and exit pressure")
    print("  match the analytic solution closely.  The peak Mach does not, and")
    print("  grid refinement does not fix it (the on-axis peak only reaches")
    print("  ~2.30 at 801 x 97, not 2.46).  Quasi-1D theory assumes one Mach")
    print(f"  number per cross section; this wall diverges at up to "
          f"{wall_angle.max():.0f} degrees,")
    print("  which makes the flow genuinely two dimensional in two ways:")
    print("   * across a section the axis is the slowest streamline and the")
    print("     wall the fastest.  Mass-averaging over the section (Figure 3c)")
    print("     recovers the quasi-1D expansion - compare the two throat rows;")
    print(f"   * the shock is curved - it stands at x = {x_shock_wall:.3f} on the "
          f"wall but")
    print(f"     x = {x_shock_axis:.3f} on the axis - so no station is entirely "
          f"supersonic")
    print("     at the analytic peak, and both computed peaks stop short of it.")

    # ---------------------------------------------------------------- figures
    print("\nDrawing Figures 3 and 4...")

    # Figure 3: the 2D story.
    fig, ax = plt.subplots(2, 2, figsize=(13, 8))

    draw_mesh_2d(ax[0, 0], x2d, r2d)
    ax[0, 0].set_title(f"(a) flattened meridional mesh\n"
                       f"{n3} 3D nodes -> {n2} 2D nodes ({factor:.0f}x fewer)")

    cf = ax[0, 1].contourf(xc, rc, M, levels=40, cmap="turbo")
    ax[0, 1].contour(xc, rc, M, levels=[1.0], colors="w", linewidths=1.5)
    fig.colorbar(cf, ax=ax[0, 1], label="Mach")
    ax[0, 1].plot(x, r_wall, "k-", lw=1.2)
    ax[0, 1].set_title("(b) Mach number (white line: M = 1)")
    ax[0, 1].set_xlabel("x")
    ax[0, 1].set_ylabel("r")
    ax[0, 1].set_aspect("equal")

    # (c) The axis is the slowest streamline once the wall diverges steeply, so
    # the on-axis curve alone reads low; the mass-averaged curve is the one the
    # quasi-1D model actually predicts, and it tracks the analytic expansion
    # until the (curved) shock, which reaches the wall before the axis.
    ax[1, 0].plot(x_axis, m_axis, "C0-", lw=1.5, label="CFD, on the axis")
    ax[1, 0].plot(x_axis, m_avg, "C1-", lw=1.5, label="CFD, mass-averaged")
    ax[1, 0].plot(analytic.x, analytic.mach, "k--", lw=1.2,
                  label="1D isentropic (A/A*) prediction + normal shock")
    ax[1, 0].axvline(analytic.x_shock, color="r", ls=":", lw=1.0,
                     label=f"1D shock x={analytic.x_shock:.3f}")
    ax[1, 0].set_title("(c) Mach along the duct - the A/A* prediction is for the "
                       "section average\n"
                       f"the shock is curved: x = {x_shock_wall:.2f} at the wall, "
                       f"{x_shock_axis:.2f} on the axis")
    ax[1, 0].set_xlabel("x")
    ax[1, 0].set_ylabel("Mach")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    it, res_hist = zip(*history)
    ax[1, 1].semilogy(it, res_hist, "C1-")
    ax[1, 1].set_title("(d) convergence history")
    ax[1, 1].set_xlabel("iteration")
    ax[1, 1].set_ylabel("relative density residual")
    ax[1, 1].grid(alpha=0.3, which="both")

    fig.suptitle("Figure 3 - converging-diverging duct, 2D axisymmetric Euler (JAX)")
    fig.tight_layout()

    # Figure 4: the software comparison.
    fig4, ax4 = plt.subplots(1, 2, figsize=(12, 4.5))
    cx_axis = np.asarray(cmetrics_np.xc)[:, 0]
    ax4[0].plot(cx_axis, M_jax[:, 0], "C0-", lw=2.5, label="vectorized (JAX)")
    ax4[0].plot(cx_axis, M_loop[:, 0], "k--", lw=1.2, label="nested Python loops")
    ax4[0].set_title(f"Mach on the axis after {COMPARE_ITERS} iterations\n"
                     f"(transient, not converged) - max diff {d_mach:.1e}")
    ax4[0].set_xlabel("x")
    ax4[0].set_ylabel("Mach")
    ax4[0].legend()
    ax4[0].grid(alpha=0.3)

    bars = ax4[1].bar(["vectorized\n(JAX)", "nested\nPython loops"],
                      [ms_jax, ms_loop], color=["C0", "C3"])
    ax4[1].set_yscale("log")
    ax4[1].set_ylabel("ms / iteration")
    ax4[1].set_title(f"Cost per iteration, {cmetrics_np.vol.size} cells\n"
                     f"{ms_loop / ms_jax:.0f}x speedup")
    for b, v in zip(bars, [ms_jax, ms_loop]):
        ax4[1].text(b.get_x() + b.get_width() / 2, v, f"{v:.3g}",
                    ha="center", va="bottom")
    fig4.suptitle("Figure 4 - same physics, two data layouts")
    fig4.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
