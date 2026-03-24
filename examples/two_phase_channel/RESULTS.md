# Two-phase turbulent channel — run results

## Case overview

The goal is to validate CDI (Olsson–Kreiss interface sharpening) and CSF (Brackbill
surface tension) in a turbulent channel flow. The turbulence is the stress-test
environment; the interface capturing method is what is being validated.

Weber number: **We = ρ U_b² R / σ = R / σ** (U_b=1, ρ=1).

| Run | Case file | We | σ | R | y_c | IC | Status | Purpose |
|-----|-----------|:--:|---|---|-----|----|--------|---------|
| `channel_test_v4` | `_v4.case` | 730 | 4.1×10⁻⁴ | 0.3 | 0 | Turbulent Reichardt | **Completed** t=0–5 | High-We reference |
| `channel_test_laminar` | `_laminar.case` | 1 | 0.3 | 0.3 | 0 | Laminar Poiseuille | **Planned** | Ground-truth CDI/CSF baseline |
| `channel_test_we10` | `_we10.case` | 10 | 0.03 | 0.3 | 0 | `fluid00004.chkp` + drop | **Running** t=20→25 | First turbulent CDI/CSF test |
| `channel_test_we1` | `_we1.case` | 1 | 0.3 | 0.3 | 0 | `fluid00004.chkp` + drop | **Planned** | Primary validation (after we10) |
| `channel_single_phase` | `_single_phase.case` | — | — | — | — | Turbulent Reichardt | **Completed** t=0–25 | Fluid spin-up; checkpoint at t=20 |
| `channel_test_restart` | `_restart.case` | 1.33 | 0.3 | 0.4 | 0 | `fluid00004.chkp` + drop | **Blown up** (v1, v2) | We=1 blow-up reference data |
| `channel_test_restart_off` | `_restart_off.case` | 1.33 | 0.3 | 0.4 | 0.3 | `fluid00004.chkp` + drop | **Ready** | Larger off-centre drop; log-law region shear |

Run directories: `/lscratch/sieburgh/simulations/<run_name>/`

---

## Success criteria

For a correctly functioning CDI/CSF implementation:

| Diagnostic | Correct | Failure |
|------------|---------|---------|
| κ_rms | Stable near 2/R = 6.67 (sphere); or slowly varying with known deformation | Unbounded growth |
| φ_max | Stable near 0.995 (IC value) | Monotonic decline toward 0.5 → dissolving |
| φ_min | \|φ_min\| ≪ 0.01 | Large negatives → numerical noise |
| E_kin | Consistent with imposed flow rate | Divergence → CSF sign/magnitude error |

The laminar We=1 case is the ground truth: κ_rms = 6.67 throughout, no deformation.
Any deviation there is a method error, not a physical effect.

---

## channel_single_phase — turbulent spin-up (COMPLETED t=0–25)

**Purpose:** Run the single-phase channel to statistical stationarity, then checkpoint.
Checkpoint `fluid00004.chkp` (t=20) is the fluid IC for all restart cases.

**Setup:** 81×18×27 mesh, N=7, Re_b=2800, 16 MPI ranks, Reichardt IC + perturbations,
end_time=25, checkpoints every 5 TU (`fluid00000.chkp`=t=0, …, `fluid00004.chkp`=t=20,
`fluid00005.chkp`=t=25).

**Turbulence indicator:** u_max drops from ~1.52 (Poiseuille-like IC) to ~1.35–1.45
by t≈10 TU and fluctuates there through t=25 — confirmed turbulent. The instantaneous
maximum includes turbulent fluctuations u'≈0.15–0.20 on top of the mean centreline
U_cl≈1.15; note that u_max ≠ U_cl.

**Mean velocity profile:** Time-averaged over t=20–25 (6 snapshots) matches the
Reichardt model throughout the channel, confirming Re_τ=180 turbulence is established.
See `postprocess_single_phase.py` → `ekin_single_phase.png`, `meanprofile_single_phase.png`.

| t | E_kin | u_max | Notes |
|---|-------|-------|-------|
| 0 | 0.507 | 1.541 | Reichardt IC + perturbations |
| 5 | 0.535 | 1.465 | Transitioning |
| 10 | 0.533 | 1.373 | Approaching turbulent state |
| 15 | 0.531 | 1.363 | Statistically stationary |
| 20 | 0.532 | 1.344 | **checkpoint fluid00004.chkp** |
| 25 | 0.531 | 1.400 | Final state |

---

## channel_test_v4 — high-We reference (We=730, COMPLETED t=0–5)

**Purpose:** High-We reference case. Surface tension negligible; drop deforms freely.
Documents CDI performance under strong interface straining.

**Setup:** ε=0.05, γ=0.015, σ=4.1×10⁻⁴, R=0.3, Re_b=2800, turbulent IC, T=5 TU.

| t | φ_max | φ_min | κ_rms | u_max | E_kin |
|---|-------|-------|-------|-------|-------|
| 0.000 | 0.996 | 0.000  |  6.37 | 1.541 | 0.507 |
| 0.513 | 0.997 | −0.001 | 34.88 | 1.569 | 0.537 |
| 0.770 | 0.998 | −0.002 | 55.12 | 1.569 | 0.537 |
| 1.026 | 1.001 | −0.004 | 55.46 | 1.567 | 0.537 |
| 1.539 | 1.003 | −0.005 | 50.93 | 1.568 | 0.537 |
| 2.309 | 0.999 | −0.006 | 46.77 | 1.526 | 0.538 |
| 3.079 | 0.998 | −0.006 | 43.53 | 1.489 | 0.538 |
| 3.849 | 0.986 | −0.005 | 43.57 | 1.462 | 0.538 |
| 4.490 | 0.971 | −0.007 | 44.73 | 1.459 | 0.538 |
| 4.875 | 0.966 | −0.008 | 46.06 | 1.465 | 0.538 |

**Key observations:**

- κ_rms jumped from 6.37 → ~55 within the first TU due to IC perturbations
  (τ_IC/τ_CDI = 1.5 — CDI barely faster than startup burst). After the burst,
  κ_rms slowly decays toward ~43 but remains well above equilibrium at t=5.
- φ_max drifted down from 0.996 to 0.966 by t=5 — the drop is deforming at
  We=730 where surface tension plays no role. Not a method failure.
- φ_min reached −0.008: acceptable undershoot.
- u_max still declining at t=5 (1.465, not yet at turbulent stationary ~1.17).
  The flow has not reached statistical stationarity in this run.

This case serves as the CDI stress test at very high We, not as the primary
CDI/CSF validation target.

![channel_test_v4 diagnostics](diagnostics_ekin_channel_test_v4.png)

---

## channel_test_laminar — CDI/CSF ground truth (We=1, PLANNED)

**Purpose:** Cleanest possible CDI/CSF test. No turbulence, no startup perturbations,
strong surface tension. The drop barely deforms (mean shear at y=0 is zero by symmetry).
κ_rms = 6.67 throughout is the analytical expectation; any deviation is a method error.

**Setup:** Laminar Poiseuille IC (no perturbations). ε=0.07, γ=0.05,
σ=0.3 (We=1), R=0.3, Re_b=2800, end_time=10.

| t | φ_max | φ_min | κ_rms | u_max |
|---|-------|-------|-------|-------|
| — | — | — | — | — |

---

## channel_test_we1 — primary validation (We=1, PLANNED)

**Purpose:** CDI/CSF under turbulent straining with strong surface tension. The CSF
force is large and measurable: Δp = 2σ/R = 2.0 ≫ ρU_b²/2 = 0.5. Any curvature
error amplifies into a visible spurious velocity. Comparison with the laminar case
isolates what turbulence adds.

**Setup:** ε=0.07, γ=0.05, σ=0.3 (We=1), R=0.3, Re_b=2800, restart from
`fluid00004.chkp`, end_time=25 (runs t=20→25).

| t | φ_max | φ_min | κ_rms | u_max |
|---|-------|-------|-------|-------|
| — | — | — | — | — |

---

## channel_test_we10 — moderate deformation (We=10, RUNNING)

**Purpose:** First turbulent CDI/CSF test with stable explicit timestep.
Δt/Δt_cap ≈ 0.45 — 2.2× inside capillary stability boundary.
Tests CDI under sustained strain at moderate deformation.

**Setup:** ε=0.07, γ=0.05, σ=0.03 (We=10), R=0.3, Re_b=2800, restart from
`fluid00004.chkp`, end_time=25 (runs t=20→25). 16 MPI ranks.

| t | φ_max | φ_min | κ_rms | u_max |
|---|-------|-------|-------|-------|
| — | — | — | — | — |

---

## channel_test_restart — CDI/CSF without startup transient (We=1, READY v2)

**Purpose:** Drop injected analytically into statistically stationary turbulence from
`fluid00004.chkp` (t=20 of single-phase run). The IC perturbations have long decayed —
eliminates the τ_IC/τ_CDI startup burst entirely. Cleanest test of CDI/CSF
performance under sustained turbulent straining, without startup contamination.

**Parameters (v2 — corrected ε, γ, R, CFL):**

| Parameter | Value | Notes |
|-----------|-------|-------|
| restart_file | `fluid00004.chkp` | Single-phase checkpoint at t=20 |
| end_time | 25.0 | **Absolute time.** Neko restores t=20 from checkpoint and runs to t=25 (5 TU) |
| σ (We=1) | 0.3 | Strong surface tension |
| ε | 0.07 | ε/Δ_GLL_x,z=3.17 ✓ (was 0.05 → under-resolved in x,z) |
| γ | 0.05 | γ/u_max=0.033 ✓; τ_CDI=0.098 TU (was 0.015 → too weak) |
| R | 0.4 | Larger drop (D=0.8h, φ(0)=0.997); was R=0.3 |
| target_cfl | 0.2 | Capillary stability: Δt_cap≈0.0015 TU at We=1 (was 0.4 → too large) |
| N | 7 | Polynomial order |
| output_value | 0.1 TU | Field files at t=20.1, 20.2, …, 25.0 |
| output_checkpoints | true | Saves checkpoint at t=25 (for potential extension) |

**Important — end_time is absolute:** Neko reads the checkpoint, restores `t = 20.0`,
and checks `is_done = (t >= end_time)` each step. With `end_time: 25.0` the simulation
runs from t=20 to t=25. Using `end_time: 5.0` would stop immediately (20 ≥ 5).

**How restart works in Neko:** `restart_file` causes Neko to restore u, v, w, p,
and BDF history from the checkpoint. The scalar has no data in the single-phase
checkpoint, so Neko falls through to the user scalar IC (scheme='phase'), which
places the drop analytically via the tanh profile. The fluid IC (scheme='fluid')
is NOT called — the velocity field comes directly from the checkpoint.

**Workflow:**
```bash
cd /lscratch/sieburgh/simulations/channel_test_restart
cp ~/code/neko-multiphase-channel/examples/two_phase_channel/turb_channel_two_phase_restart.case .
mpirun -np 16 ./neko turb_channel_two_phase_restart.case
```

**MPI rank count:** Must match the spin-up run — **16 ranks** for `fluid00004.chkp`.

**Verification at startup (check neko.log):**
- `WARNING: Checkpoint has no scalar field; scalar IC will be applied by user code`
- `Step = 1  t = 0.2000471E+02` — confirms restart from t=20
- First ekin.csv row at t≈20.002: φ_max≈0.997, κ_rms≈2/R=5.0, u_max≈1.34 ✓

**Blow-up trace (v1, ε=0.05, γ=0.015, R=0.3, target_cfl=0.4, We=1.0):**

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 20.002 | 0.996 | 0.000 | 6.37 | 1.341 | Drop injected; κ_rms≈2/R=6.67 ✓ |
| 21.0 | — | — | 8.63 | 1.656 | Gradual deformation, 1.0 TU stable |
| 21.3 | — | — | 19.21 | 2.007 | Rapid growth begins |
| 21.5 | — | — | 66.21 | 7.603 | Runaway |
| 21.53 | — | — | 130 | 24 | Blown up |
| 21.55 | 1.53 | −0.52 | 169 | 42 | Fully diverged |

**Blow-up trace (v2, ε=0.07, γ=0.05, R=0.4, target_cfl=0.2, We=1.33):**

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 20.002 | 0.9954 | 0.000 | 4.758 | 1.341 | Drop injected; κ_rms≈2/R=5.0 ✓ |
| 20.067 | 0.9963 | 0.000 | 4.858 | 1.341 | Growing; φ_max stable |
| 20.132 | 0.9965 | 0.000 | 4.998 | 1.341 | κ_rms still near 5.0 |
| 20.197 | 0.9968 | 0.000 | 5.767 | 1.342 | +21% growth in 0.2 TU |
| 20.263 | 0.9969 | 0.000 | 7.685 | 1.347 | +33% more; interface still intact |
| 20.317 | 0.9971 | 0.000 | 9.577 | 1.746 | u_max elevated; deformation driven |
| 20.356 | 0.9974 | 0.000 | 19.395 | 2.327 | κ_rms doubled; runaway begins |
| 20.378 | 0.9977 | 0.000 | 42.144 | 4.531 | CSF force large; velocity spikes |
| 20.387 | 0.9981 | 0.000 | 101.24 | 9.530 | Fully explosive; κ_rms×5 in 0.03 TU |
| 20.392 | 0.9986 | 0.000 | 170.77 | 15.71 | Interface still intact (φ_max<1!) |
| 20.395 | 1.0497 | 0.000 | 180.17 | 26.66 | φ_max > 1: CDI overloaded |
| 20.398 | 1.2195 | −0.244 | 177.74 | 37.89 | Diverged |
| 20.400 | 1.5418 | −0.436 | 173.41 | 51.40 | Fully diverged |

14 v1 field snapshots (t=20.1–21.5) available. Animation: `python3 animate_two_phase_channel.py --run channel_test_restart`.

---

## channel_test_restart_off — off-centre drop, larger radius (READY)

**Purpose:** Same restart approach as `channel_test_restart` but with a larger drop
placed off-centre in the log-law region. Tests CDI/CSF under asymmetric shear: the drop
experiences non-zero mean shear dU/dy ≠ 0 (unlike y_c=0), stronger velocity gradients,
and turbulent structures closer to the wall. Also exercises a slightly higher We.

**Parameters:**

| Parameter | Value | Notes |
|-----------|-------|-------|
| restart_file | `fluid00004.chkp` | Same turbulent IC as `channel_test_restart` |
| end_time | 25.0 | Runs t=20→25 (5 TU) |
| R | 0.4 | Larger drop (D=0.8h); minimum wall clearance = (1−0.3)−0.4 = 0.3h = 54 wall units |
| y_c | 0.3 | Wall-normal offset; y_c⁺ = 54 wall units (log-law region) |
| σ | 0.3 | We = ρU_b²R/σ = 0.4/0.3 ≈ 1.33 |
| ε | 0.07 | ε/Δ_GLL_x,z=3.17 ✓; R/2ε=2.86 → φ(0)=0.997 ✓ |
| γ | 0.05 | γ/u_max=0.033; τ_CDI=0.098 TU |
| target_cfl | 0.2 | Capillary stability at σ=0.3 |
| N | 7 | Polynomial order |

**Resolution check for R=0.4:** R/2ε = 0.4/(2×0.07) = 2.86 → φ(0) = 0.997 ✓.
Minimum wall clearance 0.3h = 54 wall units — the interface (±2ε = ±0.1h = ±18 wall
units) does not reach the viscous sublayer.

**What to look for vs `channel_test_restart`:**
- κ_rms: may deviate more from 2/R=5.0 (sphere with R=0.4) due to stronger local shear
- φ_max: should remain near 0.999 (tighter IC for R=0.4)
- Drop migration: asymmetric turbophoretic drift possible (drop in shear-dominated region)
- u_max: same turbulent background — should stay ~1.35–1.45

**Workflow:**
```bash
cd /lscratch/sieburgh/simulations/channel_test_restart_off
mpirun -np 16 ./neko turb_channel_two_phase_restart_off.case
```

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| — | — | — | — | — | Not yet run |

---

## CDI/CSF instability analysis — We=1 restart cases (v1 and v2)

Both restart attempts with We=1 (σ=0.3) in turbulent channel flow have blown up
despite different parameter sets. This section documents the two cases in detail,
explains what the diagnostics reveal about the failure mechanism, and describes
what needs to change for stable operation.

### Parameter comparison

| Parameter | v1 | v2 | Analytical target |
|-----------|----|----|-------------------|
| ε | 0.05 | 0.07 | ε/Δ_GLL_x,z ≥ 3: needs ≥0.066 |
| γ | 0.015 | 0.05 | γ/u_max ∈ [0.01, 0.1] |
| R | 0.3 | 0.4 | — |
| We | 1.0 (σ=0.3/R=0.3) | 1.33 (σ=0.3/R=0.4) | — |
| target_cfl | 0.4 | 0.2 | Δt < Δt_cap ≈ 0.0015 TU |
| Δt (nominal) | ≈ 0.0043 TU | ≈ 0.0013 TU | — |
| τ_CDI = ε²/γ_c | 0.167 TU | 0.065 TU | — |
| ε/Δ_GLL_x,z | 2.26 | 3.17 | ≥ 3 ✓ |
| γ/u_max | 0.010 | 0.033 | ≥ 0.01 ✓ |
| Blow-up time | t ≈ 21.55 (1.55 TU) | t ≈ 20.40 (0.40 TU) | — |

v2 had corrected ε and γ, yet blew up **four times faster** than v1. This is
diagnostic: the parameter corrections helped individual CDI metrics (ε/Δ, γ/u_max)
but did not address the root cause.

### What the diagnostics reveal

**φ_max behaviour — what CDI is doing**

The most important observation from v2 is that **φ_max remains below 1.0 throughout
the κ_rms growth phase**. At t=20.392, κ_rms=170 and u_max=15.7, yet φ_max=0.9986.
φ_max does not exceed 1.0 until t=20.395, by which point u_max=26.7 and the blow-up
is already irreversible. This tells us:

- CDI is **functioning correctly** throughout the explosive κ growth
- The interface **thickness is maintained** (φ_max near 1 → tanh profile intact)
- CDI failure (interface dissolving) is **not the primary cause** of the blow-up

This is critical: the blow-up is not due to CDI being unable to maintain the interface
profile. CDI is doing its job. The failure lies elsewhere.

**κ_rms growth rate — drop deformation, not dissolution**

In v2, κ_rms grows monotonically from the very first diagnostic interval:

```
t=20.002 → 20.067: κ_rms 4.758 → 4.858   (+2%   in 0.065 TU)
t=20.067 → 20.132: κ_rms 4.858 → 4.998   (+3%   in 0.065 TU)
t=20.132 → 20.197: κ_rms 4.998 → 5.767   (+15%  in 0.065 TU)
t=20.197 → 20.263: κ_rms 5.767 → 7.685   (+33%  in 0.065 TU)
t=20.263 → 20.317: κ_rms 7.685 → 9.577   (+25%  in 0.054 TU)
t=20.317 → 20.356: κ_rms 9.577 → 19.395  (+103% in 0.039 TU — doubling time ~0.039 TU)
t=20.356 → 20.378: κ_rms 19.39 → 42.14   (+117% in 0.022 TU — doubling time ~0.019 TU)
```

The doubling time accelerates from ~0.5 TU early on to ~0.02 TU during the explosive
phase. For comparison, τ_CDI = 0.065 TU. Once the doubling time falls below τ_CDI,
CDI cannot respond fast enough — but note this only happens at κ_rms>10, after the
runaway has already begun.

The growth of κ_rms while φ_max is stable means the **drop is deforming** — developing
regions of high curvature (sharp protrusions, tips, dimples) under turbulent straining
— rather than dissolving. For a sphere of radius R=0.4, κ = 2/R = 5.0 everywhere.
Any deviation from a sphere increases κ_rms above this value. The question is why
surface tension (which should resist deformation at We=1) is not providing sufficient
restoring force.

**u_max spike and the CSF feedback loop**

In v2, u_max remained stable at 1.341–1.347 through t=20.263 (three diagnostic intervals
of near-steady flow). Then between t=20.263 and t=20.317 it jumped to 1.746. This is
still within the turbulent fluctuation range, but it coincides with κ_rms accelerating
from 7.7 to 9.6. After this point, the growth becomes explosive.

The CSF source term adds a force F_ST = σ κ ∇φ to the momentum RHS. At the interface,
|∇φ| ≈ 1/(2ε) = 7.14 (for ε=0.07). The force magnitude at each diagnostic interval:

```
t=20.263: F_ST ≈ 0.3 × 7.7 × 7.14  ≈  16.5  (moderate)
t=20.317: F_ST ≈ 0.3 × 9.6 × 7.14  ≈  20.6
t=20.356: F_ST ≈ 0.3 × 19.4 × 7.14 ≈  41.5  (already large)
t=20.378: F_ST ≈ 0.3 × 42.1 × 7.14 ≈  90.2  (enormous)
t=20.387: F_ST ≈ 0.3 × 101 × 7.14  ≈ 216    (catastrophic)
```

Each step, the CSF force is applied explicitly — it adds a velocity increment
Δu_CSF ≈ F_ST × Δt to the RHS:

```
t=20.263: Δu ≈ 16.5 × 0.0013 ≈ 0.021  (small)
t=20.356: Δu ≈ 41.5 × 0.00075 ≈ 0.031 (still manageable)
t=20.378: Δu ≈ 90  × 0.00043 ≈ 0.039  (beginning to drive deformation)
t=20.387: Δu ≈ 216 × 0.00017 ≈ 0.037  (comparable to turbulent u')
```

Note that as u_max grows, target_cfl=0.2 causes Δt to shrink proportionally (Δt ∝ 1/u_max).
This partially compensates — the velocity increment Δu_CSF ∝ F_ST/u_max. However, F_ST
grows as κ grows, while u_max lags behind. The feedback loop is:

**κ grows → F_ST grows → velocity spike at interface → flow strains interface → κ grows further**

Surface tension should close this loop (larger κ → larger restoring F_ST toward sphere),
but here F_ST is driving the instability rather than damping it. The reason: the
**direction of the CSF force**. CSF acts as σκ∇φ — it drives flow toward regions of
lower pressure (inside the drop at any convex point). For a nearly spherical drop with
large σ, this correctly creates the Laplace pressure jump. But when the interface
develops sharp tips (high local κ), the force at those tips is enormous and directed
inward — creating jets that drive further deformation of neighboring interface regions.
With explicit time integration, these jets are applied over Δt without any implicit
correction, potentially overshooting the equilibrium and amplifying the deformation.

**κ_rms plateau at ~170–180**

Both v1 and v2 show κ_rms plateauing near 170–180 before starting to decrease during
the diverged phase. This is not a recovery — it is a numerical artifact. Once κ_rms
exceeds ~150, the implied local radius of curvature is:

```
r_local = 2/κ_rms ≈ 2/170 ≈ 0.012
```

This is less than ε=0.07 — the local curvature features are sub-interface in scale and
cannot be resolved. The spectral element representation saturates: κ_rms cannot grow
further because the numerical scheme has lost track of the actual interface shape. The
plateau marks the point where the problem becomes fully unresolved, not a stabilisation.

**E_kin during the blow-up**

A notable feature of v2: E_kin remains near 0.532 (the single-phase turbulent value)
throughout the stable phase and rises only slightly (0.532 → 0.570) through most of
the blow-up. It does not spike until very late (after u_max>25), indicating that the
**bulk flow is not disrupted** until the very end. The CSF instability is localised at
the drop interface and only couples to the bulk velocity field once the velocity spikes
at the interface become large enough to disturb the surrounding flow significantly.
This is consistent with a localised capillary instability, not a global flow breakdown.

### Why v2 blew up faster than v1

Despite having better CDI parameters, v2 blew up in 0.4 TU vs 1.5 TU for v1.
Three contributing factors:

**1. Larger drop (R=0.4 vs R=0.3)**

The drop occupies a larger fraction of the channel cross-section. With R=0.4 and
half-channel height h=1.0, the drop extends over 80% of the half-height. More turbulent
eddies interact with the interface simultaneously, providing a broader turbulent straining
field. The mean strain rate experienced by the interface scales with the number of
interacting turbulent structures, which grows roughly as R²/l_t² where l_t is the
integral length scale. With R=0.4 vs R=0.3, this ratio increases by (0.4/0.3)² ≈ 1.78.

**2. Higher Weber number (We=1.33 vs We=1.0)**

With σ=0.3 fixed and R=0.4, We=R/σ=1.33 rather than 1.0. The surface tension force
per unit interface area is σ/R=0.75, compared to 1.0 for R=0.3. The restoring
acceleration is 25% weaker per unit curvature. While both cases are at low We and
surface tension strongly dominates inertia, the 25% reduction in restoring force means
the feedback loop amplifies faster.

**3. Larger initial curvature sensitivity**

For R=0.4, 2/R=5.0. For R=0.3, 2/R=6.67. A given fractional perturbation to the
interface shape (say 10% deviation from spherical) produces a larger absolute increase
in κ_rms for R=0.3 (larger base curvature means more variation). But the restoring
force scales as σκ — and since σ is the same, the R=0.3 drop actually has stronger
absolute restoring force at equivalent relative deformation. This partially explains
why v1 remained stable for longer.

**Summary: why parameter corrections alone were insufficient**

The parameter corrections (ε, γ, target_cfl) addressed:
- Resolution of the interface profile (ε) ✓
- CDI resharpening rate relative to advection (γ) ✓
- Velocity-based timestep constraint (target_cfl) ✓

But they did not address the fundamental stability issue:

- The **capillary wave stability constraint** is Δt_cap = sqrt(ρΔx³/(2πσ)). For
  σ=0.3 and Δx=0.016, Δt_cap ≈ 0.0015 TU. With target_cfl=0.2 and u_max=1.34,
  Δt ≈ 0.0013 TU — barely inside the limit. The margin is only 13%.
- During the first explosive phase (t=20.27–20.32), u_max rises to 1.75 before
  the CFL controller has reacted. At u_max=1.75 and old Δt=0.0013 TU:
  actual CFL ≈ 0.2 × (1.75/1.34) = 0.26, which exceeds target_cfl=0.2. Meanwhile
  Δt_cap is unchanged. The actual ratio Δt/Δt_cap ≈ 0.0013/0.0015 = 0.87 —
  we were only 13% inside the stability boundary throughout.
- More fundamentally: **the capillary stability limit is a local condition**. The
  minimum Δx in a spectral element mesh is not Δ_elem/N but the GLL spacing near
  the element boundary, which is approximately (Δ_elem/N²) in the worst case.
  For N=7, the minimum GLL spacing near element edges is roughly 0.016/7 × 0.2 ≈ 0.0005.
  The implied capillary limit would then be Δt_cap ≈ sqrt(0.0005³/(2π×0.3)) ≈ 9×10⁻⁶ TU.
  This is far too small to be practical — suggesting that at We=1, explicit CSF may be
  fundamentally incompatible with the SEM GLL point distribution at N=7.

### Diagnostic summary — what went right and what failed

| Diagnostic | v1 | v2 | Conclusion |
|------------|----|----|------------|
| φ_max stability | Degraded late (t>21.5) | Stable throughout blow-up | CDI working ✓ |
| κ_rms growth from t=0 | Gradual (1 TU stable) | Immediate | CSF instability, not CDI |
| u_max at blow-up onset | ~1.7 | ~1.75 | Same trigger level |
| E_kin at blow-up onset | ~0.534 | ~0.532 | Bulk flow undisturbed |
| Time to blow-up | 1.55 TU | 0.40 TU | v2 worse due to R, We |
| κ_rms plateau | ~170 | ~180 | Same: resolution limit |

**The instability is in the explicit CSF treatment, not CDI.** CDI is maintaining the
interface; CSF is creating velocity spikes that the explicit integrator cannot damp.

### Path forward

Three options, in order of increasing We=1 difficulty:

**Option 1: We=10, σ=0.03 (recommended first run)**

σ=0.03 reduces the CSF force by 10× relative to We=1. The capillary stability limit
becomes Δt_cap = sqrt(0.016³/(2π×0.03)) ≈ 0.005 TU — more than 3× larger than the
velocity CFL Δt at target_cfl=0.4. CSF at We=10 should be stable with current parameters.
This is the physically interesting case for CDI validation under moderate deformation.

**Option 2: laminar We=1 (ground-truth validation)**

The laminar case has no turbulent straining — the drop barely deforms (mean shear at
y=0 is zero by symmetry). Without the turbulent amplification, any CSF instability
grows much slower and the laminar flow provides a clean analytical baseline (κ_rms=5.0
throughout). Running the laminar case first would confirm whether CDI+CSF work correctly
in the absence of turbulence before diagnosing the turbulent instability.

**Option 3: We=1 turbulent with target_cfl=0.05**

A 4× reduction in target_cfl (from 0.2 to 0.05) would give Δt ≈ 0.00033 TU —
roughly 4× inside the estimated capillary stability limit even accounting for CFL
controller lag. This would make the We=1 turbulent run ≈4× more expensive (4× more
timesteps), but would definitively test whether the blow-up is timestep-driven or
physical. If the blow-up persists at target_cfl=0.05, the instability is physical
(turbulence genuinely destroys the drop at We=1 with R=0.4); if it disappears,
the explicit CSF timestep constraint was the limiting factor.

---

## Next steps

- [x] `channel_single_phase` completed — `fluid00004.chkp` at t=20 available
- [x] `channel_test_restart` v1 — blew up at t=21.55 (1.55 TU); 14 field files available
- [x] `channel_test_restart` v2 — blew up at t=20.40 (0.40 TU); corrected CDI params not sufficient
- [~] **`channel_test_we10`** (We=10, σ=0.03) — **running** t=20→25; expected stable (Δt/Δt_cap=0.45)
- [ ] **`channel_test_laminar`** (We=1, laminar IC) — validate CDI+CSF without turbulence
- [ ] `channel_test_we1` (We=1, restart) — primary validation; blocked pending stable timestep strategy
- [ ] `channel_test_restart_off` (We=1.33, R=0.4, y_c=0.3) — after stable parameters found
- [ ] Compare κ_rms, φ_max across laminar / We=10 / We=1 to quantify CSF stiffness regime
