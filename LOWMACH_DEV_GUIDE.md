# Low-Mach (variable-density) solver for Neko — developer guide

Branch: `feature/lowmach`.  Scheme name: `"lowmach"` → `fluid_lowmach_pnpn_t`
(`src/fluid/fluid_lowmach_pnpn.f90`), registered in
`src/fluid/fluid_base_fctry.f90`.  CPU backend only.

This is the detailed design document. For the short version see
`LOWMACH_IMPLEMENTATION.md`. Audience: a developer who knows Neko's
incompressible `fluid_pnpn` solver and wants to understand/extend this.

--------------------------------------------------------------------------------
## 0. What problem this solves

The incompressible solver assumes constant density. We need the **low-Mach
number** regime: a flow that is acoustically incompressible (Ma → 0, no acoustic
waves) but has **large density variation from heating** (e.g. an air channel,
T = 300 → 1600 K, ρ varies ~5×). Boussinesq is invalid at that density ratio;
the fully compressible Euler solver is the wrong regime (high-speed, inviscid,
acoustically stiff). The low-Mach model filters acoustics, keeps a **non-zero
velocity divergence** set by thermal expansion, and couples to the temperature
(scalar) equation through an ideal-gas equation of state.

--------------------------------------------------------------------------------
## 1. Governing equations

Pressure is split into a uniform thermodynamic part `P0(t)` and a hydrodynamic
part `p(x,t)`. Open domain ⇒ `P0 = const`, `dP0/dt = 0`.

**EOS (ideal gas):**            `ρ = P0 / (R T)`

**Continuity / divergence constraint** (from `Dρ/Dt = −ρ∇·u` and the EOS):

    ∇·u = Q_T ,    Q_T = ∇·(λ∇T) / (ρ cp T)

(`Q_T` is the *thermal divergence*; derived using the energy equation
`ρ cp DT/Dt = ∇·(λ∇T)` to eliminate `DT/Dt`.)

**Momentum** (full deviatoric stress, because `∇·u ≠ 0`):

    ρ Du/Dt = −∇p + ∇·τ + f ,
    τ = μ(∇u + ∇uᵀ) − (2/3) μ (∇·u) I

**Energy** — solved by the Neko **scalar** named `temperature`:

    ρ cp (∂T/∂t + u·∇T) = ∇·(λ∇T)

Consistency requirement (critical): `ρ`, `λ`, `cp` used on the momentum/`Q_T`
side must be the **same fields** the scalar equation uses. Otherwise the imposed
`∇·u = Q_T` does not match the temperature the scalar transports and the scheme
blows up (this was a real bug — see §6).

--------------------------------------------------------------------------------
## 2. Architecture / class layout

```
fluid_scheme_incompressible_t           (existing)
        └── fluid_pnpn_t                 (existing Pn-Pn projection solver)
                └── fluid_lowmach_pnpn_t (NEW; overrides init + step)
```

`fluid_lowmach_pnpn_t` (`fluid_lowmach_pnpn.f90`) reuses *everything* from
`fluid_pnpn` (mesh, dofmap, gs, BCs, projection, the `Ax_vel`/`Ax_prs`
operators, the krylov solvers, the AB/BDF time scheme) and adds:

| member | meaning |
|---|---|
| `low_mach_enabled` | master switch; if false, `step` runs vanilla pnpn |
| `P0, R_specific` | EOS constants (`case.fluid.low_mach.P0 / .R_specific`) |
| `T_eos_min, T_eos_max` | clamp on T used in the EOS (robustness, §3) |
| `k_cond, cp` | fallback constant conductivity / cp (used only if the scalar fields are absent) |
| `T_field_name, T_ptr` | name + lazily-resolved pointer to the temperature field |
| `lambda_ptr, cp_ptr` | lazily-resolved pointers to the scalar's `λ(T)`, `cp` fields (§5) |
| `Q_T_field` | the thermal divergence field |
| `lm_prs_res, lm_vel_res` | low-Mach residual operators (abstract `lowmach_*_res_t`) |

Residual operators:
- `src/fluid/lowmach_residual.f90` — abstract `lowmach_prs_res_t` /
  `lowmach_vel_res_t` (interfaces include `Q_T` and a `rho` field).
- `src/fluid/lowmach_res_fctry.f90` — factory.
- `src/fluid/bcknd/cpu/lowmach_res_cpu.f90` — CPU implementations
  `lowmach_prs_res_cpu_compute` / `lowmach_vel_res_cpu_compute`.

(No device backend yet; `init` errors out on `NEKO_BCKND_DEVICE`.)

--------------------------------------------------------------------------------
## 3. EOS density update — `lowmach_update_density`

Called at the top of `step`. Recomputes `ρ` from the current `T`:

    rho(i) = P0 / (R_specific * clamp(T(i), T_eos_min, T_eos_max))

`this%rho` is the registry field `fluid_rho`, *shared* with the scalar scheme
(passed in as the scalar's density), so updating it propagates to the energy
equation. The **clamp** exists because a temperature over/undershoot (Gibbs
wiggle near a sharp gradient) would otherwise send `ρ = P0/(R T) → ∞` or
**negative**, which destroys both the physics and the pressure conditioning.
Set `case.fluid.low_mach.T_eos_min/T_eos_max` to the physical T range (generous).

--------------------------------------------------------------------------------
## 4. Momentum residual — `lowmach_vel_res_cpu_compute`

Goal: assemble `∇·τ` for variable `ρ, μ` and non-zero `∇·u`, **without changing
`ax_helm`'s interface** (it has 6 backends and is performance-critical).

Decomposition:
- `∇·[μ(∇u+∇uᵀ)]` → the **existing coupled full-stress operator** `ax_helm_full`
  via `Ax%compute_vector(...)`, with `μ` carried in `coef%h1`, `h2 = ρ·bd/dt`.
  Verified (kernel `ax_helm_stress_lx`) that it computes `μ(∇u+∇uᵀ)` with **no**
  `−2/3` trace term — so it is purely the deformation part.
- dilatation `−(2/3)∇(μ Q_T)` → assembled **explicitly** with `opgrad` on the
  field `μ·Q_T` (Q_T is known this step) and added to the residual exactly like
  the `−∇p` and `+f` terms.

Therefore low-Mach **requires `full_stress_formulation = true`** so that
`Ax_vel` is the coupled operator (enforced in `init`), which in turn requires a
coupled velocity solver (`velocity_solver.type = "coupled_cg"`).

This is the same "implicit operator + explicit correction" pattern Neko's
*incompressible* stress formulation (`pnpn_res_stress_cpu`) already uses; nothing
new is threaded into `ax_helm`.

--------------------------------------------------------------------------------
## 5. Thermal divergence `Q_T` — `lowmach_update_Q_T`

Computed once per step from `T`. **This is the part that took the most care; it
matches Nek5000 `userqtl_scig` (`Nek5000/core/plan4.f`).**

Discretization (SEM weak→strong reconstruction), identical structure to Nek5000:

    g = opgrad(T)                  ! weak grad
    gs_op(g, ADD); g *= Binv       ! assemble -> strong grad
    g *= λ(x)                      ! multiply by conductivity FIELD
    Q = cdtp_x(gx)+cdtp_y(gy)+cdtp_z(gz)   ! weak divergence
    gs_op(Q, ADD); Q *= Binv       ! assemble -> strong
    Q_T = Q / (ρ cp T)             ! divide by rho*cp*T

The decisive detail: `λ` and `cp` are the scalar's **spatially varying fields**,
resolved lazily from the registry as `"<temperature>_lambda_tot"` (fallback
`"_lambda"`) and `"_cp"` — i.e. Nek5000's `vdiff`/`vtrans` of the temperature
field. They are filled by the user `material_properties` hook (e.g. Sutherland
`λ(T)`). Using a **constant** `k`/`cp` instead (the original code) makes `Q_T`
inconsistent with the energy balance and was the cause of the heated-channel
blow-up. If the scalar fields are absent, it falls back to `k_cond`/`cp`.

The weak `cdtp` drops the boundary surface integral — **so does Nek5000's
`opdiv`**; that is not a bug, the volumetric `∇·(λ∇T)` carries the wall flux when
the gradient is resolved.

--------------------------------------------------------------------------------
## 6. Pressure residual — `lowmach_prs_res_cpu_compute`

Builds the RHS of the pressure-Poisson with three additions over incompressible:

1. **Variable-density operator**: `coef%h1 = 1/ρ` (via `invers2`), so the
   pressure operator is `∇·((1/ρ)∇p)` (vs the scalar `1/ρ_const` the
   incompressible solver uses, `pnpn_res_cpu.f90:66`).
2. **Variable-property viscous force** in rotational form (mirrors the validated
   incompressible stress residual): `wa = μ(∇×∇×u) − Sᵀ(2∇μ)`, plus the
   low-Mach corrections obtained by substituting `∇·u = Q_T` into the rotational
   identity:

        Δwa = −(4/3) μ ∇Q_T + (2/3) Q_T ∇μ

   (`∇μ` is reused from the `Sᵀ∇μ` term; `∇Q_T` via `dudxyz`.)
3. **Continuity source**: `p_res += (bd/dt) Q_T B` — the `div u = Q_T`
   constraint (Nek5000 `plan4` `admcol3` convention).

--------------------------------------------------------------------------------
## 7. Variable-density pressure SOLVE — `lowmach_pressure_solve` (the hard part)

### The problem
The pressure operator now has a *spatially varying* coefficient `h1 = 1/ρ(x)`.
Neko's multigrid preconditioners (`hsmg`, `phmg`) are built/tuned for the
**constant-coefficient** Laplacian (the incompressible solver literally uses a
scalar density). With variable `1/ρ`, GMRES **stagnates** (≈0 % residual
reduction) and the run blows up — *only* once density varies.

Diagnosis was by elimination (all in `examples/heated_channel/`): incompressible
`pnpn` on the same mesh/BCs converges; low-Mach stalls even at ΔT = 10 K and
with `Q_T` off but variable ρ on; `phmg` also fails. ⇒ it is the variable-ρ
operator, not BCs/mesh/`Q_T`.

### The fix — defect correction
Solve the true variable-ρ system `A_var p = b` by an outer iteration whose inner
solve uses a **constant reference operator** `A_ref` (hsmg's home turf):

    p ← p + A_ref⁻¹ ( b − A_var p )

- `A_var` (the true residual `b − A_var p`, BCs and all) is re-evaluated each
  outer pass by calling `lm_prs_res%compute` — the exact, MMS-verified residual.
- `A_ref` is the constant-coefficient operator `h1 = c_ref`, inverted by the
  existing `ksp_prs` + `hsmg` (~40–100 iters).
- **`c_ref = ½(min(1/ρ) + max(1/ρ))`** — the midpoint of the `1/ρ` range. This
  makes the fixed-point contraction factor `(b−a)/(b+a) < 1` for *any* density
  ratio (≈0.33–0.4 for a 2× ratio ⇒ ~14 outer passes).

Implementation notes:
- `c_Xh%h1` is swapped to `c_ref` only around the inner solve; `ax_helm` is never
  modified.
- Inner solve tolerance must be **tight** (`pressure_solver.absolute_tolerance
  ≈ 1e-6`); a loose inner tol gives an inaccurate correction and slow outer
  convergence.
- Convergence: `outer_rtol = 1e-4`, `outer_atol = 1e-8`, `max_outer = 60`. Sets
  `ksp_result%converged` so `fluid_step_info` labels it correctly (this flag was
  the source of the spurious "did not converge for Pressure" log line).
- Why the EOS clamp matters here too: a Gibbs node with `T → 0` makes
  `1/ρ → ∞`, so the `1/ρ` range → `[~0, …]`, `c_ref` is wrecked and the
  contraction → 1. The clamp keeps the range bounded.

Possible future optimisation: replace the fixed-point outer loop with FGMRES
(use `A_ref` as a flexible preconditioner) to cut the outer count for stiff
density ratios; cache the `p`-independent part of the residual.

--------------------------------------------------------------------------------
## 8. Variable-density RHS makers

`lm_rhs_ext_field_rho`, `lm_rhs_bdf_field_rho`, `lm_rhs_oifs_field_rho` are
variable-ρ counterparts of `rhs_maker_*_cpu`: they assemble the AB-extrapolated
forcing history and the BDF lagged-velocity term weighting **per node by ρ(x)**
instead of a scalar ρ. Consequence for source terms: a user `source_term`
provides an **acceleration**; the rhs-maker multiplies it by ρ to get the body
force. (So gravity is passed as `(ρ_ref/ρ − 1) g` → ×ρ → `−(ρ−ρ_ref) g`, the
reduced buoyancy.)

--------------------------------------------------------------------------------
## 9. Time step — `fluid_lowmach_pnpn_step`

Per step (low-Mach branch):
1. `lowmach_update_density` — ρ from EOS (clamped).
2. `lowmach_update_Q_T` — Q_T from T, λ(T), cp.
3. extrapolate velocity (`sumab`), compute source term, advection.
4. variable-ρ RHS assembly (`lm_rhs_*_field_rho`).
5. `update_material_properties` (user hook sets μ(T); for the scalar, λ(T), cp).
6. **pressure**: `lowmach_pressure_solve` (defect correction, §7).
7. **velocity**: `lm_vel_res%compute` (full stress + dilatation, §4), then the
   coupled velocity solve (`solve_coupled`).
8. optional `vol_flow%adjust` if `flow_rate_force` is set.

Property/coupling timing: the scalar is solved after the fluid in the simulation
loop, so at step *n* the fluid sees `T`, `λ`, `ρ` from step *n−1* (first-order
explicit coupling — standard, fine).

--------------------------------------------------------------------------------
## 10. Material properties / coupling (user side)

Variable `μ(T)`, `λ(T)`, `cp` are supplied through the Neko `material_properties`
user hook:
- `scheme_name == 'fluid'`  → `properties = [rho, mu]`; set `mu = μ(T)`
  (Sutherland). `ax_helm_full` then uses it (variable-μ stress). Do **not**
  overwrite `rho` (the EOS owns it).
- `scheme_name == 'temperature'` → `properties = [cp, lambda]`; set
  `lambda = λ(T)`, `cp`. `Q_T` reads the same `λ`/`cp` fields ⇒ consistent.

EOS `P0`, `R_specific` come from `case.fluid.low_mach`. Gravity/forcing come from
a user `source_term`.

--------------------------------------------------------------------------------
## 11. Required configuration

```jsonc
"fluid": {
  "scheme": "lowmach",
  "full_stress_formulation": true,          // REQUIRED (coupled stress operator)
  "velocity_solver": { "type": "coupled_cg", ... },   // REQUIRED with full stress
  "pressure_solver": { "type": "gmres", "preconditioner": {"type":"hsmg"},
                       "absolute_tolerance": 1e-6 },   // tight inner tol
  "low_mach": {
     "enabled": true,
     "P0": 101325.0, "R_specific": 287.0,   // EOS (air); default 1.0 (non-dim)
     "cp": 1005.0, "k_conductivity": 0.06,  // fallback consts if no scalar fields
     "T_eos_min": 250.0, "T_eos_max": 700.0,
     "temperature_field": "temperature"
  }
},
"scalar": { "name": "temperature", "enabled": true, ... }   // the energy equation
```
Also: a `scalar` with `λ`, `cp` set consistently (material_properties hook or
`startup` add); do not run the scalar in its non-dimensional `cp=ρ=1` mode if you
want dimensional physics.

Build (this worktree — NOT the April `/home/hochi/neko-install`):
`make -C src && make install`, then `install/bin/makeneko user.f90`.

--------------------------------------------------------------------------------
## 12. Verification

- **Operator-level MMS** — `examples/mms_lowmach/` (`./run.sh`): manufactured
  fields on a periodic box; the discrete `∇·τ` (velocity residual) and
  `∇·(∇·τ)` (pressure residual) converge **spectrally** to the analytic value
  with **variable μ**. Confirms signs/coefficients of the stress + dilatation
  terms. (Calls the shipped residual routines.)
- **Coupled time loop** — `examples/rayleigh_benard/rayleigh_lowmach*` runs to
  Normal end.
- **Heated channel (dimensional, Sutherland)** — `examples/heated_channel/`:
  runs to Normal end with gravity + Q_T + μ(T) + variable ρ + inlet T-ramp once
  the consistent-`Q_T` fix is in (`run_fullfix.log`).
- **Turbulent channel (Re_τ=180, DNS, order 7)** — `examples/turb_channel/`:
  set up and runs stably (`flow_rate_force`, hot/cold walls, Reichardt+
  perturbation IC). Reaching developed turbulence is cluster-scale.

Still TODO: quantitative validation vs a reference (Nusselt, mean/RMS profiles).

--------------------------------------------------------------------------------
## 13. Known limitations / future work

- **CPU only** — no CUDA/HIP/OpenCL backend (`init` errors on device).
- **`dP0/dt` not implemented** — open domains only; closed/confined heated
  cavities (rising background pressure) unsupported.
- **Pressure-solve cost** — defect correction is ~7–14 inner solves/step;
  FGMRES acceleration + residual caching would speed it up.
- **Q_T at no-slip heated walls** — `∇·u = Q_T` vs `u = 0` is the usual low-Mach
  wall subtlety; stable here, worth scrutiny for accuracy.
- **Strong heating transients** — a temporal heating ramp + small `dt` may be
  needed when imposing very large ΔT instantly.

--------------------------------------------------------------------------------
## 14. File-by-file change map

| File | Role |
|---|---|
| `src/fluid/fluid_lowmach_pnpn.f90` | scheme: init (reads JSON, enforces full-stress, allocates Q_T, builds residuals), EOS+clamp (`lowmach_update_density`), consistent Q_T (`lowmach_update_Q_T`), variable-ρ rhs makers, `step`, defect-correction pressure (`lowmach_pressure_solve`) |
| `src/fluid/bcknd/cpu/lowmach_res_cpu.f90` | momentum residual (full stress + dilatation), pressure residual (variable-ρ operator + rotational stress + Δwa + Q_T source) |
| `src/fluid/lowmach_residual.f90`, `lowmach_res_fctry.f90` | abstract residual types + factory (interfaces carry `Q_T`, `rho`) |
| `src/fluid/fluid_base_fctry.f90` | registers `"lowmach"` scheme |
| `examples/mms_lowmach/` | MMS verification (driver + run.sh + README) |
| `examples/heated_channel/` | dimensional Sutherland heated channel (mesh gen, case, user.f90: μ/k(T), gravity, wall-T ramp) |
| `examples/turb_channel/lowmach_turb.{f90,case}` | Re_τ=180 turbulent channel (DNS, order 7) |

--------------------------------------------------------------------------------
## 15. Debugging lessons (so the next person doesn't repeat them)

1. **Consistency beats magnitude.** The heated-channel blow-up was not a residual
   sign error (MMS was clean) nor gravity nor the inlet corner — it was `Q_T`
   using a constant `k` inconsistent with the scalar's `λ(T)`. Cross-checking
   Nek5000 `userqtl_scig` showed the structure was already right; only the
   property fields were wrong.
2. **Isolate by toggling one physics term at a time** (Q_T off, gravity off,
   incompressible baseline on the same mesh). That chain pinned each cause
   unambiguously.
3. **The pressure preconditioner is constant-coefficient.** Any variable-density
   elliptic operator in Neko needs the defect-correction (or a coefficient-aware
   preconditioner); do not expect `hsmg` to handle `1/ρ(x)` directly.
4. **Always set `ksp_monitor_t%converged`** in a custom solve, or
   `fluid_step_info` will mislabel it.
5. **`build`/`install` target the worktree, not the stale system install.**
