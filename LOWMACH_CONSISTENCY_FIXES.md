# Low-Mach ↔ Nek5000 consistency review & fixes

Date: 2026-05-30. Branch: `feature/lowmach`.

This document records a physical-consistency review of the Neko low-Mach
implementation against the established Nek5000 low-Mach (Tomboulides PnPn)
solver, and the fixes applied as a result. It is meant to be read alongside
`LOWMACH_IMPLEMENTATION.md` (the original dev writeup) and the project status
memory.

## TL;DR

The **fluid** side (momentum / pressure / variable-property viscous stress) was
found **correct and consistent with Nek5000** — including the parts that look
unusual (the rotational-form stress and the `−4/3 / +2/3` low-Mach coefficients,
the explicit `−(2/3)∇(μQ_T)` dilatation, the `(bd/dt)Q_T·B` pressure source, the
variable-density `1/ρ` Poisson). No changes were needed there.

The real inconsistencies were in the **coupling to the scalar (energy) equation**
and in two property/EOS details. The 600 K test case ran, but these would break
the 300→1600 K target. Four fixes were applied (one critical, three high). A
fifth item (a one-step ordering lag) is **documented but deliberately not
changed** because the safe remedy has cross-case implications — see below.

All changes are gated/no-op for non-low-Mach cases, so the rest of the Neko
test suite is unaffected.

---

## What was verified correct (no change)

Confirmed by derivation and by cross-reading Nek5000 `plan4.f` / `navier1.f`:

- **Rotational-form viscous stress** in the pressure residual
  (`lowmach_res_cpu.f90`): `wa = μ(∇×∇×u) − Sᵀ(2∇μ)` plus the correction
  `−(4/3)μ∇Q_T + (2/3)Q_T∇μ` equals the exact `−div(τ)` for
  `τ = μ(∇u+∇uᵀ) − (2/3)μ(∇·u)I` with `∇·u = Q_T`. The `4/3` arises as
  `−2 + 2/3` (the rotational identity's missing `−2μ∇(∇·u)` plus the dilatation's
  `+(2/3)μ∇Q_T`). The coefficients are right.
- **Velocity residual** full-stress operator + explicit `−(2/3)∇(μQ_T)`: the
  coupled `ax_helm_full` computes only `div[μ(∇u+∇uᵀ)]` (no `−2/3` trace), so the
  explicit dilatation is **not** double-counted.
- **`(bd/dt)·Q_T·B`** pressure source matches Nek `admcol3(respr,QTL,bm1,dtbd)`.
- **Constant P0 / no dP0/dt** is the physically correct branch for the *open*
  heated channel (matches Nek default `ifdp0dt=.false.`); for the closed
  Rayleigh–Bénard case `ortho()` removes the incompatible pressure mode exactly
  as Nek does. Only relevant if a *closed* heated domain is run later.

---

## Fixes applied

### Fix 1 (CRITICAL) — energy equation now uses the full ρ·cp(x) field

**Problem.** The fluid imposes `div u = Q_T = div(λ∇T)/(ρ(x)·cp(x)·T(x))` with the
fully spatially-varying `ρcp`. But the scalar (energy) solver transported `T`
with a **constant point-1** `ρcp`:
- `scalar_residual_cpu.f90` took `rhocp` as a scalar `real` (`cfill(h2, rhocp·bd/dt)`),
- `scalar_pnpn.f90` passed `rho%x(1,1,1,1)*cp%x(1,1,1,1)` to the residual and
  `rho%x(1,1,1,1)` to the EXT/BDF makers (no `cp`, no spatial variation).

So the temperature the scalar produced did **not** satisfy the energy balance
`Q_T` assumes. Over 300→1600 K, `ρcp` varies > 5×, making `div u = Q_T`
inconsistent at leading order — the same failure class that produced the
heated-channel blow-up. Nek5000 avoids this by carrying `vtrans(:,:,:,2)=ρcp` as
a **full field** in both the T-Helmholtz (`h2 = vtrans·bd/dt`) and `userqtl_scig`.

**Change.** A gated variable-`ρcp` path in `scalar_pnpn.f90`, enabled by the new
case key `case.scalar.low_mach.enabled` (default `.false.` → stock path is
bit-for-bit unchanged):
- New member `scalar_pnpn_t%low_mach_rhocp`, read in `scalar_pnpn_init`.
- In `scalar_pnpn_step`: build a scratch field `ρcp(x) = rho(x)·cp(x)` and use it
  to density-weight the EXT/BDF/OIFS RHS via new helpers
  `lm_scalar_rhs_ext_field`, `lm_scalar_rhs_bdf_field`, `lm_scalar_rhs_oifs_field`
  (counterparts of the fluid's `lm_rhs_*_field_rho`).
- The residual is computed inline with `h1 = lambda_tot`,
  `h2 = ρcp(x)·bd/dt` (a **field**, via `cmult2`), replicating
  `scalar_residual_cpu` but with the variable coefficient. `ρcp` is rebuilt with
  the freshly-updated `cp` after `update_material_properties`.

`rho` here is the shared EOS density (`scalar_scheme%rho => rho`, the same field
the fluid updates), so the scalar now sees the variable density it always had
access to but previously collapsed to point 1.

> Note: this also fixes a *latent* stock inconsistency — even for constant
> properties, the stock scalar uses `ρcp` on the LHS mass term but only `ρ` on
> the BDF/EXT RHS, which disagree when `cp ≠ 1`. The low-Mach path uses `ρcp`
> consistently everywhere.

Files: `src/scalar/scalar_pnpn.f90`, `examples/heated_channel/heated_channel.case`.

### Fix 2 (HIGH) — EOS temperature clamp made self-consistent

**Problem.** `lowmach_update_density` forms `rho = P0/(R·T_clamp)` with
`T_clamp = min(max(T,T_eos_min),T_eos_max)`, but `lowmach_update_Q_T` divided by
the **unclamped** `T`. When `T` leaves the clamp window, `ρ·T ≠ P0/R`, breaking
the ideal-gas relation the `Q_T = div(λ∇T)/(ρcpT)` form relies on — precisely in
the hottest cells. The case also set `T_eos_max = 700 K`, below the 1600 K
target, so the clamp would activate over a large hot region.

**Change.**
- `lowmach_update_Q_T` now divides by the **same clamped** temperature, so
  `ρ·T_clamp = P0/R` holds identically (`src/fluid/fluid_lowmach_pnpn.f90`).
- `heated_channel.case`: `T_eos_max 700 → 2000` so the clamp is a true overshoot
  guard, never active in the physical 300–1600 K range.

### Fix 3 (HIGH) — `lambda_tot` and `mu_tot` refreshed in laminar runs

**Problem.** `scalar_scheme_update_material_properties` and
`fluid_scheme_update_material_properties` only wrote `lambda_tot` / `mu_tot`
inside their turbulence branches. With no SGS model (the laminar heated channel),
those `*_tot` fields stayed frozen at their init value, so a user-written
Sutherland `k(T)` / `μ(T)` (into `lambda` / `mu`) **never reached** the diffusion
operator, the viscous stress, or `Q_T` — silently defeating the project's
temperature-dependent properties.

**Change.** An unconditional `field_copy(lambda_tot, lambda)` /
`field_copy(mu_tot, mu)` right after the user material-properties hook, before the
turbulence branches add `ν_t` on top. No-op for constant properties; correct for
variable. Files: `src/scalar/scalar_scheme.f90`,
`src/fluid/fluid_scheme_incompressible.f90`.

---

## Deferred (documented, not changed): one-step coupling lag

**Issue.** Neko steps the fluid before the scalar (`simulation.f90:157` then
`:166`), and the low-Mach fluid updates `rho` and `Q_T` at the top of its step
from the temperature in the registry — which is still `Tⁿ` at that instant. So
`rho` and `Q_T` lag the temperature by one step. Nek5000 does the opposite order
(`heat → setprop → qthermal → fluid`, `drive1.f:284–292`), evaluating them from
`Tⁿ⁺¹`.

**Why it is not changed here.** Both orders are legitimate first-order (O(dt))
operator splits (confirmed by adversarial review), and the consistency fixes
above remove the actual blow-up mechanism. The clean remedy — advancing the
scalar before the fluid — is a **global** change to `simulation.f90` that alters
the time-splitting for *every* coupled Neko case (buoyancy/Boussinesq, etc.) and
their reference results. That is out of scope for a low-Mach-only fix and should
be a deliberate, separately-validated decision.

**If you want it:** the minimal, lowest-risk option is to recompute `rho` and
`Q_T` from `Tⁿ⁺¹` *after* the scalar step (e.g. re-call
`lowmach_update_density` / `lowmach_update_Q_T` at the end of the coupled step so
they are ready for the next fluid solve), or reorder scalar-before-fluid behind a
flag. Ask and this can be added.

---

## How to build, enable, and run

Build (this worktree — not the April `/home/hochi/neko-install`):

```
LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib make -C src
make install
install/bin/makeneko examples/heated_channel/heated_channel.f90
```

Enable the variable-`ρcp` energy coupling in the case file:

```json
"fluid":  { "scheme": "lowmach", "full_stress_formulation": true,
            "low_mach": { "enabled": true, "T_eos_max": 2000.0, ... } },
"scalar": { "low_mach": { "enabled": true }, ... }
```

`case.scalar.low_mach.enabled = true` turns on the field-`ρcp` path; omitting it
(or `false`) leaves the stock constant-property scalar untouched.

When enabled, the log prints
`Scalar low-Mach: variable rho*cp(x) weighting enabled ...` at init — use that to
confirm the flag was parsed.

## Files changed

| File | Fix |
|------|-----|
| `src/scalar/scalar_pnpn.f90` | 1 — gated field-`ρcp` energy equation + helpers |
| `src/scalar/scalar_scheme.f90` | 3 — laminar `lambda_tot` refresh |
| `src/fluid/fluid_scheme_incompressible.f90` | 3 — laminar `mu_tot` refresh |
| `src/fluid/fluid_lowmach_pnpn.f90` | 2 — `Q_T` uses the clamped `T` |
| `examples/heated_channel/heated_channel.case` | 2 — `T_eos_max→2000`; scalar `low_mach.enabled` |

No new source files, factory, or build-system changes; nothing outside the
low-Mach path changes behavior.
