# Converging-diverging duct: geometry → 3D mesh → flatten → Euler solve

A complete, self-contained tutorial pipeline that takes a parametric duct wall,
turns it into a real Plot3D block, flattens it back down to a 2D meridional
mesh, and solves the axisymmetric Euler equations on it with JAX.

The flow chokes at the throat and forms a normal shock in the diverging
section — a genuinely interesting result to check against the classic quasi-1D
analytic solution.

## Run it

```bash
pip install -e ".[examples]"      # from the repository root: pyturbo-aero + jax
cd examples/converge_diverge_duct
python main.py
```

No arguments, no input files. Takes about 10 seconds on a laptop CPU.

## What each file does

| File | Responsibility |
|---|---|
| `duct_geometry.py` | Bézier wall radius `R(x)` in the meridional plane. |
| `duct_mesh.py` | Revolve `(x, R)` about the x-axis into a `plot3d.Block`, with a wall-clustered inflation layer. |
| `duct_flatten.py` | Collapse the axisymmetric 3D block back to a 2D `(x, r)` node grid. |
| `euler_metrics.py` | Node grid → axisymmetric cell volumes, face areas and face normals. |
| `euler_physics.py` | Gas model, state conversions, Rusanov flux, isentropic relations. |
| `euler_bc.py` | Ghost-cell boundary conditions: subsonic inlet, back-pressure outlet, slip wall / axis. |
| `euler_solver_flatten.py` | The vectorized solver, JIT-compiled with JAX. |
| `euler_solver_block.py` | The same physics in nested Python `for` loops — the slow reference, used only for the comparison. |
| `quasi1d.py` | Analytic quasi-1D solution with a normal shock. |
| `main.py` | Runs everything and draws the figures. |

Every module is importable with no side effects, and every one of them can also
be run on its own (`python duct_mesh.py`, `python euler_metrics.py`, …) to see
just that stage.

## The two payoffs of flattening

**Less data.** The 201 × 25 × 37 revolved block has 185,925 nodes. Because the
duct is a body of revolution, every one of the 37 constant-θ planes carries
identical `(x, r)` information, so 5,025 nodes say exactly the same thing —
37× less data, verified by an `axisymmetry_error` of ~1e-16.

**Faster code.** Once the mesh is a dense 2D array, the whole flow field is a
single `(NI, NJ, 4)` array and the solver becomes array expressions that
`jax.jit` can fuse and vectorize. `main.py` runs the identical physics both ways
on the same coarse mesh from the same initial state:

* they agree to ~1e-16 (they are the same equations, after all);
* the vectorized version is roughly **200× faster per iteration**.

## Physics and numerics

2D axisymmetric Euler, written in the `r`-weighted form so the 2D solve is
exactly equivalent to the 3D body of revolution:

```
∂(rU)/∂t + ∂(rF)/∂x + ∂(rG)/∂r = [0, 0, p, 0]
```

* First-order Rusanov (local Lax-Friedrichs) flux — no tunable constants, one
  ghost layer, robust through a strong shock.
* Ghost-cell BCs: Riemann-invariant subsonic inlet (`P0 = 2.0e5 Pa`,
  `T0 = 500 K`), fixed back pressure outlet (`p = 1.0e5 Pa`), reflective slip
  wall. The axis needs no special treatment — its face area is `∝ r`, which is
  exactly zero there.
* 4-stage Runge-Kutta, `α = (1/4, 1/3, 1/2, 1)`, local time step at `CFL = 1.5`.
* `float64` throughout (`jax_enable_x64`). In `float32` the residual stalls
  around 1e-4.

## Expected output

Stage 0 prints a straight-pipe unit check. Uniform flow in a constant-radius
pipe is an exact solution, so the residual must be machine zero
(~1e-16 relative). If that number is not tiny, nothing downstream is
trustworthy — it is the single check that exercises the axisymmetric source
term, every face-normal orientation and the wall reflection at once.

Right after the geometry stage — before any CFD runs — `main.py` prints what 1D
isentropic theory predicts from the duct's own `A(x)/A*`: the subsonic and
supersonic branches at the exit, why the back pressure forces a normal shock
inside the diverging section, where that shock stands, and the headline
prediction

```
  1D isentropic prediction: exit Mach = 0.1283   (exit p = 100.0 kPa)
```

The main solve then converges in roughly 6,000 iterations (~3.5 s) to a
relative density residual below 1e-6, and prints a table like:

```
  quantity                               CFD     1D (A/A*)
  --------------------------------------------------------
  throat Mach (on axis)               0.9640        1.0000
  throat Mach (mass-avg)              0.9929        1.0000
  peak Mach (on axis)                 2.0867        2.4592
  peak Mach (mass-avg)                2.0837        2.4592
  shock x (on axis)                   0.3255        0.3139
  shock x (at wall)                   0.2685        0.3139
  exit Mach                           0.1250        0.1283
  exit pressure [kPa]               100.0664      100.0000
```

followed by the exit Mach spelled out on its own (prediction 0.1283, CFD
0.1250) and a short station-by-station table of Mach vs. `x` — the `A/A*`
prediction next to the mass-averaged and on-axis CFD values.

Shock location, exit Mach and exit pressure all land close to the analytic
values. The peak Mach does not, and **grid refinement does not fix it** — the
on-axis peak reaches only ~2.30 even at 801 × 97, never the analytic 2.46. The
reason is not numerical diffusion. Quasi-1D theory assumes one Mach number per
cross section, and this duct's wall diverges at up to 24°, which makes the flow
genuinely two-dimensional in two ways:

1. **the Mach number varies across a section** — the axis is the slowest
   streamline, the wall the fastest (≈1.97 on the axis vs ≈2.20 at the wall at
   `x ≈ 0.23`). Averaging over the section weighted by mass flow,
   `Σ ρu·S·M / Σ ρu·S`, recovers the quasi-1D expansion curve: compare the two
   throat rows above, and the orange vs blue curves in figure 3(c);
2. **the shock is curved** — it stands at `x ≈ 0.27` on the wall but
   `x ≈ 0.33` on the axis, so no station is entirely supersonic at the analytic
   peak, and both computed peaks stop short of it.

So the on-axis Mach is an inherently poor proxy for the quasi-1D value here.
The mass-averaged curve is the like-for-like comparison, and the remaining
gap is the 2D shock, which no quasi-1D model can represent.

Four figures are produced:

1. **domain overview**, before solving — the flattened 2D `(x, r)` mesh (wall
   clustering visible) beside the revolved 3D mesh it came from;
2. **ghost cells** — the 2D domain with one ring of ghost-cell positions
   outside each boundary (inlet, outlet, axis, wall), plus a zoom at the
   throat wall where the clustering is tightest. The solver never stores these
   positions (`euler_bc.py` only builds ghost *states*); they are reconstructed
   for the picture by `euler_metrics.ghost_cell_centres`;
3. a 2 × 2 panel: the flattened mesh, the Mach field with the sonic line, Mach
   along the duct (on-axis and mass-averaged) against the quasi-1D solution,
   and the convergence history;
4. the JAX-vs-loop comparison: overlaid Mach profiles (they coincide) and a
   log-scale cost-per-iteration bar chart.

The 50-iteration comparison in figure 4 is a **transient, not a converged
answer** — that stage compares software, not flow accuracy.

## Notes

* The mesh is a full 360° body of revolution — the first and last θ planes
  sit at the same physical location, closing the mesh rather than leaving it
  as a wedge slice.
* Node counts are chosen so `gcd(200, 24, 36) = 4` (and `gcd(40, 12, 6) = 2`
  for the coarse mesh); the 3D wireframe is GCD-reduced before drawing, the
  same way `plot3d.plot_blocks` does it, and other counts can make it
  unusable.
* `main.py` writes `duct_3d.xyz` and reads it straight back, so the solver
  metrics really are built from a round-tripped Plot3D file. `*.xyz` is already
  gitignored.
* The coarse comparison mesh has different array shapes from the production
  mesh, so JAX traces and compiles `step` a second time. That compile cost is
  timed and reported separately — a cheap demonstration of JAX's
  shape-dependent recompilation.
