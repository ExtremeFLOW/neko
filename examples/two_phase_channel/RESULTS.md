# Two-phase turbulent channel — run results

## Case overview

The goal is to validate CDI (Olsson–Kreiss interface sharpening) and CSF (Brackbill
surface tension) in a turbulent channel flow. The turbulence is the stress-test
environment; the interface capturing method is what is being validated.

Weber number: **We = ρ U_b² R / σ = R / σ** (U_b=1, ρ=1).

| Run | Case file | We | σ | R | y_c | IC | Status | Purpose |
|-----|-----------|:--:|---|---|-----|----|--------|---------|
| `channel_test_v4` | `_v4.case` | 730 | 4.1×10⁻⁴ | 0.3 | 0 | Turbulent Reichardt | **Completed** t=0–5 | High-We reference |
| `channel_test_laminar` | `_laminar.case` | 1 | 0.3 | 0.3 | 0 | Laminar | **Planned** | Ground-truth CDI/CSF baseline |
| `channel_test_we1` | `_we1.case` | 1 | 0.3 | 0.3 | 0 | Turbulent Reichardt | **Planned** | Primary validation |
| `channel_test_we10` | `_we10.case` | 10 | 0.03 | 0.3 | 0 | Turbulent Reichardt | **Planned** | Moderate deformation |
| `channel_single_phase` | `_single_phase.case` | — | — | — | — | Turbulent Reichardt | **Completed** t=0–25 | Fluid spin-up; checkpoint at t=20 |
| `channel_test_restart` | `_restart.case` | 1 | 0.3 | 0.4 | 0 | `fluid00004.chkp` + drop | **Ready** (v2) | CDI/CSF without startup transient |
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

**Setup:** Identical to we1.case except `turbulent_ic: false`. ε=0.05, γ=0.015,
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

**Setup:** ε=0.05, γ=0.015, σ=0.3 (We=1), R=0.3, Re_b=2800, turbulent IC, T=5 TU.

| t | φ_max | φ_min | κ_rms | u_max |
|---|-------|-------|-------|-------|
| — | — | — | — | — |

---

## channel_test_we10 — moderate deformation (We=10, PLANNED)

**Purpose:** CDI/CSF with moderate interface deformation. Inertia and surface tension
comparable. Tests CDI under sustained strain without the free-deformation limit of v4.

**Setup:** ε=0.05, γ=0.015, σ=0.03 (We=10), R=0.3, Re_b=2800, turbulent IC, T=5 TU.

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

**Root cause of v1 blow-up (t=21.5):** Three compounding issues:
1. **ε too small** (0.05): ε/Δ_GLL_x,z=2.26 — under-resolved in streamwise/spanwise;
   inaccurate curvature → CSF force errors.
2. **γ too weak** (0.015): γ/u_max=0.010, at the lower bound. During turbulent
   spikes u_max→2–3: γ/u_max→0.005, below criterion. CDI could not keep up.
3. **target_cfl=0.4 too large for We=1**: capillary stability requires
   Δt < sqrt(Δx³/(2πσ)) ≈ 0.0015 TU; velocity CFL gave Δt≈0.004 TU (3× too large).
   Capillary waves grew unbounded → κ_rms 6→169, u_max 1.3→42 in 1.5 TU.

**Blow-up trace (v1, ε=0.05, γ=0.015, R=0.3, target_cfl=0.4):**

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 20.002 | 0.996 | 0.000 | 6.37 | 1.341 | Drop injected ✓ |
| 21.0 | — | — | 8.63 | 1.656 | Gradual deformation |
| 21.3 | — | — | 19.21 | 2.007 | Rapid growth begins |
| 21.5 | — | — | 66.21 | 7.603 | Runaway |
| 21.53 | — | — | 130 | 24 | Blown up |
| 21.55 | 1.53 | −0.52 | 169 | 42 | Fully diverged |

14 field snapshots (t=20.1–21.5) available. Animation: `python3 animate_two_phase_channel.py --run channel_test_restart`.

**Important — end_time is absolute:** Neko reads the checkpoint, restores `t = 20.0`,
and checks `is_done = (t >= end_time)` each step. With `end_time: 25.0` the simulation
runs from t=20 to t=25. Using `end_time: 5.0` would stop immediately (20 ≥ 5).

**How restart works in Neko:** `restart_file` causes Neko to restore u, v, w, p,
and BDF history from the checkpoint. The scalar has no data in the single-phase
checkpoint, so Neko falls through to the user scalar IC (scheme='phase'), which
places the drop analytically via the tanh profile. The fluid IC (scheme='fluid')
is NOT called — the velocity field comes directly from the checkpoint.

**Workflow (v2):**
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

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| — | — | — | — | — | v2 not yet run |

---

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

## Next steps

- [x] `channel_single_phase` completed — `fluid00004.chkp` at t=20 available
- [ ] `channel_test_restart` v2 (We=1, R=0.4, ε=0.07, γ=0.05, CFL=0.2) — **ready to run**
  - v1 (ε=0.05, γ=0.015, CFL=0.4) blew up at t=21.5; 14 field files available for animation
- [ ] `channel_test_restart_off` (We=1.33, R=0.4, y_c=0.3, same corrected params) — **ready**, run after restart
- [ ] `channel_test_laminar` (We=1) — ground-truth CDI/CSF baseline (no turbulence)
- [ ] `channel_test_we1` (We=1, turbulent Reichardt IC) — primary validation
- [ ] `channel_test_we10` (We=10, turbulent) — moderate deformation
- [ ] Compare κ_rms, φ_max across restart / restart_off / laminar to quantify turbulence effect on CDI
