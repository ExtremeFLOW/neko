# Low-Mach (variable-density) solver for Neko — implementation notes

Branch: `feature/lowmach`. Scheme: `fluid.scheme = "lowmach"`
(`fluid_lowmach_pnpn_t`, extends `fluid_pnpn_t`).

Target application: an **open channel with a strongly heated wall** (air,
T ≈ 300 → 1600 K), low Mach number. Density varies ~5×, and μ, k, cp all vary
with T. This is the low-Mach variable-density regime (acoustically filtered,
non-zero velocity divergence from thermal expansion) — *not* Boussinesq, which
is invalid at this density ratio.

This document records the physics, the bugs found, the Nek5000 cross-check, the
fixes, and how everything was validated. It is meant for a developer review.

---

## 1. Governing equations (what the code solves)

Open domain ⇒ thermodynamic pressure `P0 = const` (no `dP0/dt`). Ideal gas.

- **EOS:** `ρ = P0 / (R T)`
- **Continuity / divergence constraint:** `∇·u = Q_T`, with
  `Q_T = ∇·(λ∇T) / (ρ cp T)` (the thermal divergence; from `Dρ/Dt` + the energy
  equation `ρ cp DT/Dt = ∇·(λ∇T)`).
- **Momentum:** `ρ Du/Dt = −∇p + ∇·τ`, with the full deviatoric stress
  `τ = μ(∇u + ∇uᵀ) − (2/3) μ (∇·u) I`. (+ body force `f`, e.g. gravity.)
- **Energy:** solved by the Neko scalar (`temperature`):
  `ρ cp (∂T/∂t + u·∇T) = ∇·(λ∇T)`.

`λ = λ(T)`, `μ = μ(T)` (Sutherland), `cp` ≈ const for air. `ρ`, `λ`, `cp` MUST
be the same fields in the momentum/Q_T side and the energy side — consistency is
essential (see §5).

---

## 2. Files changed

| File | What |
|---|---|
| `src/fluid/bcknd/cpu/lowmach_res_cpu.f90` | momentum residual = full variable-property viscous stress + dilatation; pressure residual = stress form + low-Mach `∇Q_T`/`Sᵀ∇μ` corrections + `Q_T` continuity source |
| `src/fluid/fluid_lowmach_pnpn.f90` | EOS + **EOS temperature clamp**; **consistent variable-λ(T) Q_T** (Nek-matching); **defect-correction variable-density pressure solve**; `full_stress_formulation` requirement |
| `examples/mms_lowmach/` | method-of-manufactured-solutions verification of the viscous residuals |
| `examples/heated_channel/` | the dimensional heated-channel case (mesh gen, `.case`, user `.f90` with Sutherland μ/k, gravity, streamwise wall-T ramp) |

CPU backend only. `ax_helm`'s interface was deliberately **not** changed (see §4).

---

## 3. Viscous stress in the residuals (momentum + pressure)

The deviatoric stress `∇·τ` splits into a part the existing operators give and an
explicit part:

- `∇·[μ(∇u+∇uᵀ)]` → the **existing** coupled full-stress operator `ax_helm_full`
  (`Ax%compute_vector`, μ carried in `coef%h1`). No new operator arguments.
- dilatation `−(2/3)∇(μ Q_T)` → assembled explicitly with `opgrad` (Q_T known).

`lowmach_vel_res_cpu_compute`: `h1=μ, h2=ρ·bd/dt`, `Ax%compute_vector`, then add
`−(2/3)·opgrad(μ Q_T)`. Requires `full_stress_formulation = true` (so `Ax_vel`
is the coupled operator) ⇒ enforced in `init`.

`lowmach_prs_res_cpu_compute`: variable-density pressure operator (`h1 = 1/ρ`),
rotational viscous form `μ(∇×∇×u) − Sᵀ(2∇μ)` (mirrors Neko's *incompressible*
stress residual) **plus** the low-Mach corrections derived from substituting
`∇·u = Q_T` into the rotational identity:
`Δ = −(4/3) μ ∇Q_T + (2/3) Q_T ∇μ`, plus the continuity source
`p_res += (bd/dt) Q_T B` (Nek5000 `plan4` `admcol3` convention).

**Validation — `examples/mms_lowmach/`** (`./run.sh`): manufactured fields on a
periodic box, spectral convergence of the discrete `∇·τ` (and `∇·∇·τ` for the
pressure side) to the analytic value, with **variable μ**. L∞ error drops ~2
decades per +2 polynomial order → signs and coefficients are correct.

Note (checked in `ax_helm_full_cpu.f90`): the full-stress kernel computes
`s11 = 2μ ∂u/∂x …`, i.e. `μ(∇u+∇uᵀ)` with **no** `−2/3` trace term — so the
explicit dilatation is not double-counted.

---

## 4. The variable-density pressure solve (the hard numerical bug)

**Symptom:** the pressure GMRES stagnated (≈0 % residual reduction) and the run
blew up, but *only once density varied* (step 1 with uniform T converged fine).

**Diagnosis (by elimination, `examples/heated_channel/`):**
- incompressible `pnpn` on the *same mesh + BCs* converges (`run_inc.log`).
- low-Mach stalls even at ΔT=10 K, and with `Q_T` off but variable ρ on
  (`run_nograv.log`) → it is the variable-ρ **operator**, not Q_T/BCs/mesh.
- `phmg` also fails.

**Root cause:** the low-Mach pressure Poisson has coefficient `h1 = 1/ρ(x)`
(spatially varying). The incompressible solver uses a *scalar* density
(`pnpn_res_cpu.f90:66`, `rho_val = rho%x(1)`), so its operator is
constant-coefficient — which is what Neko's `hsmg`/`phmg` multigrid is built for.
With variable `1/ρ` the preconditioner is ineffective.

**Fix — defect correction** (`lowmach_pressure_solve` in `fluid_lowmach_pnpn.f90`):
solve the *true* variable-ρ system by iterating
`p ← p + A_ref⁻¹ (b − A_var p)`, where
- `A_var` (true, `h1=1/ρ`) is re-evaluated each pass by `lm_prs_res` (this is the
  exact residual, BCs and all), and
- `A_ref` (constant `h1 = c_ref`) is inverted by the existing `ksp`/`hsmg` — its
  home turf, ~50 iters.
- `c_ref = ½(min(1/ρ) + max(1/ρ))` (midpoint of the `1/ρ` range) makes the fixed
  point contract with factor `(b−a)/(b+a) < 1` for any density ratio.

Converges ~0.4/iter (≈14 outer passes) for a 2× ratio. `ax_helm` is never
modified; only `c_Xh%h1` is swapped around the existing solve.

---

## 5. The Q_T thermal divergence (the physics bug) — Nek5000 cross-check

**Symptom:** with the pressure solver fixed, the *isolated* variable-ρ case ran
to Normal end (`run_final.log`), but the **full case (gravity + Q_T)** still blew
up at step ~2 (pressure RHS → 1e3 → NaN). Isolation showed it was **Q_T**, not
gravity (`run_qt_only.log` blows up identically with gravity off), and not the
inlet corner (mesh refinement + a smooth wall-T ramp did not help, `run_ramp.log`).

**Nek5000 reference** — `userqtl_scig` in `Nek5000/core/plan4.f`:
```
opgrad(t)          ! weak grad T
opdssum            ! assemble (gs add)
opcolv(.,binvm1)   ! ×1/B  -> strong grad
opcolv(.,vdiff(...,2))   ! × λ  (the VARIABLE temperature conductivity field)
opdiv              ! weak divergence
dssum; col2(.,binvm1)    ! assemble; ×1/B -> strong
... / (vtrans(...,2) * t) ! ÷ (ρ cp · T)   [vtrans of T = ρ cp]
```
i.e. `Q_T = ∇·(λ∇T)/(ρ cp T)`, using the **field** `λ` (`vdiff`) and `ρ cp`
(`vtrans`) of the temperature field, *not* constants. The boundary surface term
is dropped (weak `opdiv`) — same as Neko's `cdtp`; that was never the issue.

**Root cause in Neko:** `lowmach_update_Q_T` had the *identical structure* to
`userqtl_scig` but multiplied by a **constant `k_cond`** and divided by a
**constant `cp`**. With `k_const ≠ λ(T)` the imposed `∇·u = Q_T` is inconsistent
with the temperature the scalar actually transports → wrong, large source → blow-up.

**Fix:** `lowmach_update_Q_T` now resolves the scalar's fields from the registry
(`"<temperature>_lambda_tot"`/`"_lambda"`, `"_cp"`) and uses `λ(T)` and `cp` per
node (falling back to the constants if absent). One-line-per-node change; the
opgrad→cdtp structure is unchanged (it already matched Nek5000). After this the
pressure RHS is ~2e-5 (was ~1e3) and the full case is stable.

The scalar's `λ(T)`, `μ(T)` are supplied by the user `material_properties` hook
(Sutherland) in `examples/heated_channel/heated_channel.f90`; `Q_T` reads the
same `λ` field ⇒ consistent.

---

## 6. EOS temperature clamp

`ρ = P0/(R T)` blows up (or goes negative) if `T` over/undershoots near a sharp
gradient (Gibbs). That single bad node makes the `1/ρ` range `~[0,1.7]`, which
also wrecks the defect-correction contraction. `lowmach_update_density` now
clamps `T` to `[T_eos_min, T_eos_max]` (JSON `case.fluid.low_mach.T_eos_min/max`)
before the EOS. Defensive guard; with the ramp + refined mesh it rarely binds.

---

## 7. The heated-channel case (`examples/heated_channel/`)

- Mesh: `genmeshbox 0 0.2 0 0.02 0 0.04 24 16 4 .false. .false. .true.`
  (streamwise x inflow→outflow, wall-normal y, periodic z).
- Fluid: `scheme=lowmach`, `full_stress_formulation=true`,
  `velocity_solver=coupled_cg` (required with full stress),
  `pressure_solver=gmres`/`hsmg`, `low_mach{P0=101325,R_specific=287,
  k_conductivity,cp,T_eos_min=250,T_eos_max=700,temperature_field=temperature}`.
- BCs: inlet `velocity_value [1,0,0]`; outlet `outflow`; walls `no_slip`;
  scalar: inlet `dirichlet 300`, outlet `neumann 0`, **hot wall `user`** (the
  ramp), cold wall `dirichlet 300`.
- User `.f90`: `material_properties` sets Sutherland μ(T) (fluid) and λ(T)=k(T)
  (scalar); `source_term` adds reduced buoyancy `−(ρ−ρ_ref)g` (the rhs-maker
  multiplies the source by ρ, so the hook passes the acceleration form);
  `dirichlet_conditions` ramps the hot-wall T smoothly 300→600 over 3 cm
  (`smoothstep`) to keep the inlet corner continuous.

**Result:** `run_fullfix.log` → **Normal end, 0 NaN**, `field0.f00000` written.
Pressure converges ~1e-8 every step.

---

## 8. Config requirements (low-Mach)

- `case.fluid.scheme = "lowmach"`, `case.fluid.low_mach.enabled = true`.
- `case.fluid.full_stress_formulation = true` (else `init` errors).
- `case.fluid.velocity_solver.type = "coupled_cg"` (full stress ⇒ coupled solve).
- a `scalar` named to match `low_mach.temperature_field`, with `λ`, `cp` set
  consistently (e.g. via a `material_properties` hook); for variable properties,
  do NOT run the scalar in its non-dimensional `cp=ρ=1` mode.

Build (this worktree, not the April `/home/hochi/neko-install`):
`make -C src && make install`, then `install/bin/makeneko user.f90`.

---

## 9. Validation status

- **Viscous residuals:** MMS, spectral convergence, variable μ (`examples/mms_lowmach`). ✅
- **Variable-density pressure solver:** isolated case runs to Normal end. ✅
- **Full heated channel:** runs to Normal end, stable, output written. ✅
- **Physical/quantitative validation against a reference (Nusselt, profiles):**
  NOT yet done — this only establishes the solver runs and is internally
  consistent.

## 10. Known limitations / future work

- **CPU only** — no device (CUDA/HIP) backend; `init` errors on device.
- **`dP0/dt` not implemented** — open domains only (fine for the channel).
- **Defect-correction cost** — ~7–14 inner solves/step; could be Krylov-
  accelerated (FGMRES with `A_ref` as preconditioner) for fewer outer passes.
  ~23 cosmetic "defect correction not converged" notices appear when it hits the
  outer cap at `res ≈ 1e-8`; loosen `outer_atol` in `lowmach_pressure_solve` to
  silence.
- **`Q_T` near no-slip heated walls** — the constraint `∇·u=Q_T` vs `u=0` is the
  usual low-Mach wall subtlety; behaves here but worth scrutiny for accuracy.
- **Quantitative validation** and a temporal heating ramp for very strong ΔT
  (e.g. the full 1600 K target) are the next steps.
