# low-Mach validation: Nek5000 `lowMach_test` (method of manufactured solutions)

This case validates Neko's low-Mach solver against Nek5000's
`short_tests/lowMach_test` — a linearized manufactured solution from
Tomboulides, Lee & Orszag (JCP 146, 1998). It is the case that exposed (and now
verifies the fix for) three structural bugs in the low-Mach path.

## The manufactured solution

Domain: Nek5000's own mesh (`lowMach_test.re2`, x ∈ [-1,1], y ∈ [0,1], converted
to `lowMach_nek.nmsh` with `rea2nbin`; Neko extrudes the 2-D mesh to a thin slab).
With δ = 0.2, P0 = R = cp = μ = λ = 1:

```
u = T = 0.5 (3 + tanh(x/δ)),   v = w = 0
ρ   = P0/(R T) = 1/T                         (ideal-gas EOS)
q   = sech²(x/δ)/δ · (0.5 + tanh(x/δ)/δ)     (manufactured heat source)
Q_T = ∇·u = du/dx = 0.5/δ (1 - tanh²(x/δ))   (thermal divergence = QTL)
```

This `(u, T, p)` is the exact steady state of the intended low-Mach equations
```
ρ (u·∇)u = -∇p + ∇·[μ(∇u + ∇uᵀ) - (2/3)μ(∇·u)I]
ρ cp (u·∇)T = ∇·(λ∇T) + q
∇·u = Q_T
```
so the driver sets the IC = the analytic state, pins the x-faces with it, and
`user%compute` reports the MAX error of `u`, `T`, `Q_T` vs the closed forms each
step (mirroring Nek5000 `USERCHK`).

## Result

Neko **reproduces and holds** the manufactured solution to discretization
accuracy. Over a t = 2.0 run (2000 steps, dt = 1e-3 = Nek5000's timestep),
`Normal end`, 0 non-convergences, errors flat the whole way:

| quantity | max error |
|----------|-----------|
| `u` (velocity)            | 2.0e-5 |
| `T` (temperature)         | 3.0e-7 |
| `Q_T` (thermal divergence)| 7.2e-3 |

(`Q_T` ~7e-3 is the spectral-element discretization of the sharp δ=0.2 tanh
second derivative; it converges with resolution.) See `viz_profiles.png`
(every GLL point on the analytic curve) and `viz_contour.png`.

## Run it

```bash
export LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib
../../install/bin/makeneko lowmach_validation.f90
./neko lowmach_validation.case
# visual comparison vs analytic:
/home/hochi/projects/matrixenv/bin/python render_validation.py   # -> viz_*.png
```

Config notes: low-Mach requires `full_stress_formulation: true` +
`velocity_solver.type: coupled_cg`. At μ=1 the velocity Helmholtz is
stiffness-dominated, so use `dt ≤ 1e-3`. The Nek5000 mesh is well conditioned;
`genmeshbox` boxes stall the variable-density pressure solve.

## The three bugs this case found and fixed

All three were needed for the manufactured solution to be held; each was isolated
by a targeted diagnostic (frozen-T field, dt/order/operator-invariance tests, and
a 7-agent adversarial term-by-term audit of Neko vs Nek5000 `plan4`/`makeq`).

1. **`Q_T` conduction operator** (`fluid_lowmach_pnpn.f90`). `∇·(λ∇T)` was built
   with `cdtp` — the weak/transpose operator `Dᵀ` — which does **not** compose
   with `Binv` into a strong divergence. It returned ±900 oscillatory garbage
   even for `T = 2 + x²` (whose Laplacian is the constant 2). Fixed by using
   `dudxyz` (the strong collocation derivative), the same primitive Neko's `div()`
   operator uses. Verified: `T=2+x² → Q_T = 2` to 1e-11.

2. **missing source in `Q_T`** (`fluid_lowmach_pnpn.f90`). `Q_T` was conduction
   only; the volumetric heat release belongs in the thermal divergence,
   `Q_T = [∇·(λ∇T) + q]/(ρ cp T)` (Nek5000 `userqtl_scig`). Added an optional
   registry field `<temperature>_qdot` that is added to `Q_T`; off by default so
   existing cases are bit-for-bit unchanged. `Q_T_field` is now registered
   (`'Q_T'`) so user hooks / field writers can read it.

3. **energy-source density mis-weighting** (`scalar_pnpn.f90`) — the structural
   cause of the steady drift. In the low-Mach `ρcp` path the user source `q` was
   bundled into the RHS with the advection and the whole RHS multiplied by
   `ρcp(x) = 1/T`, so the solver enforced `… + ρcp·q` instead of `… + q`, leaving
   a `(ρcp - 1)·q` steady residual everywhere `T≠1`. The source (and weak-BC
   flux) is now snapshotted raw **before** the `ρcp` advection weighting and
   re-added un-weighted — mirroring Nek5000's split of `makeuq` (raw source) from
   `convab` (vtrans-weighted advection). **This affects every low-Mach case with a
   scalar source** (heat release, species reactions).

The adversarial audit also confirmed the **momentum** side is structurally exact
at the manufactured solution (field-ρ advection, the `(4/3)μ u''` viscous trace
split is not double-counted in `ax_helm_full`, the `-(2/3)∇(μQ_T)` dilatation,
`∇μ = 0` since μ=1, BDF cancellation) — it is not a source of error.

## Paper-style convergence test: linear initial condition

The hold test above starts AT the analytic steady state. Tomboulides & Orszag
instead converge TO it, which additionally exercises the scheme during a strong
density transient. The user file supports this via the optional case key
`"linear_ic": true` (read in `startup`): the IC becomes the straight line
between the same x = ±1 boundary values, u = T = 1.5 + 0.5·tanh(1/δ)·x, and the
run must find the tanh front on its own.

```bash
cd run_linear_ic && LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib \
  ../neko linear_ic.case   # end_time 4.0 (≈6 transient time constants)
```

Result (dt = 1e-3, t = 4.0, `Normal end`): errors decay exponentially from
O(1) and land on the SAME discretization floor as the hold test —

| quantity | linear IC, t=4 | hold test |
|----------|---------------|-----------|
| `u`   | 2.03e-5 | 2.0e-5 |
| `T`   | 4.4e-7  | 3.0e-7 |
| `Q_T` | 7.17e-3 | 7.2e-3 |

### Bug 4 (found by this test): pressure projection is unsafe under a ρ transient

With `pressure_solver.projection_space_size: 4` (the hold-test setting) the
transient run **blows up at t ≈ 0.45**: clean decay for ~470 steps, then a
high-wavenumber velocity disturbance grows ~35%/step to NaN, while every KSP
solve still converges to tolerance. Discriminating experiments
(`run_linic_dthalf/`, `run_linic_noproj/`):

- dt 1e-3, projection 4 → NaN at t ≈ 0.47
- dt 5e-4, projection 4 → NaN at t ≈ 0.42 (same PHYSICAL time — not a CFL/dt issue;
  log shows `WARNING: New vector not linearly indepependent!` at onset)
- dt 1e-3, projection 0 → stable, converges to the floor above

Cause: the projection basis stores `A·x_i` products for the pressure operator
`A = ∇·(1/ρ ∇)`, but in low-Mach ρ(x,t) evolves every step, so the stored
products go stale; at onset the projection step AMPLIFIED the pre-projection
residual 25× (6.4e-3 → 0.16) instead of reducing it, and the injected pressure
error feeds back through u → T → ρ. The hold test never sees this because ρ is
time-constant from step one.

**Fix (verified, recommended): keep projection but set
`pressure_solver.projection_reorthogonalize_basis: true`** — the existing flag
the 2026-07-01 audit predicted would be needed (its mechanism 4). Verified in
`run_linic_reorth/`: dt 1e-3, projection 4 + reorthogonalize runs the full
transient to `Normal end` with final errors BIT-FOR-BIT identical to the
projection-off run (u 2.027e-5, T 4.421e-7, Q_T 7.173e-3), zero
linear-dependence warnings, AND keeps the acceleration: ~2 pressure
iterations/step late-run vs ~40 with projection off. Use this in any low-Mach
case with a density transient (e.g. heated-case startups);
`projection_space_size: 0` is the fallback. Nek5000's projection has the same
structural staleness.
