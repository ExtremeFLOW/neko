# Two-phase turbulent channel — run results

## Case overview

The goal is to validate CDI (Olsson–Kreiss interface sharpening) and CSF (Brackbill
surface tension) in a turbulent channel flow. The turbulence is the stress-test
environment; the interface capturing method is what is being validated.

Weber number: **We = ρ U_b² R / σ = R / σ** (U_b=1, ρ=1).

| Run | Case file | We | σ | R | y_c | IC | Status | Purpose |
|-----|-----------|:--:|---|---|-----|----|--------|---------|
| `channel_single_phase` | `_single_phase.case` | — | — | — | — | Turbulent Reichardt | **Completed** t=0–25 | Fluid spin-up; checkpoint at t=20 |
| `channel_test_v4` | `_v4.case` | 730 | 4.1×10⁻⁴ | 0.3 | 0 | Turbulent Reichardt | **Completed** t=0–5 | High-We reference (σ≈0, no CSF instability) |
| `channel_test_restart` | `_restart.case` | 1.33 | 0.3 | 0.4 | 0 | `fluid00004.chkp` + drop | **Blown up** t=20.40 (v2) | We=1.33 blow-up reference; explicit CSF instability |
| `channel_test_we10` | `_we10.case` | 10 | 0.03 | 0.3 | 0 | `fluid00004.chkp` + drop | **Blown up** t=20.44 | Δt/Δt_cap=0.69 — below 1 is not sufficient |
| `channel_test_we10_diag` | `_we10_diag.case` | 10 | 0.03 | 0.3 | 0 | `fluid00004.chkp` + drop | **Blown up** t=20.447 | Dense snapshots (0.02 TU); same signature; animation produced |
| `channel_test_laminar` | `_laminar.case` | 1 | 0.3 | 0.3 | 0 | Reichardt IC, no perturbations | **Blown up** t=0.90 | Seeding not necessary; same CSF instability |
| `channel_test_sigma0` | `_sigma0.case` | — | 0 | 0.3 | 0 | `fluid00004.chkp` + drop | **Terminated** t=21.24 | CDI-only: κ_rms 6→64 (0.4 TU); CDI intact; SEM element-face artifact dominant |
| `channel_test_sigma0_diag` | `_sigma0_diag.case` | — | 0 | 0.3 | 0 | `fluid00004.chkp` + drop | **Completed** t=20→20.5 | 20 snapshots at 0.02 TU; **confirmed** κ grows with element faces crossed |
| `channel_test_sigma0_gamma025` | `_sigma0_gamma025.case` | — | 0 | 0.3 | 0 | `fluid00004.chkp` + drop | **Completed** t=20→20.26 | γ=0.25 parametric: stronger CDI makes κ 9× worse (sharper interface → larger C0-kink) |
| `channel_test_we1` | `_we1.case` | 1 | 0.3 | 0.3 | 0 | `fluid00004.chkp` + drop | **Deferred** | Superseded by Phase 2 (finer mesh) |
| `channel_test_restart_off` | `_restart_off.case` | 1.33 | 0.3 | 0.4 | 0.3 | `fluid00004.chkp` + drop | **Deferred** | Superseded by Phase 2 (finer mesh) |
| `channel_p2_single_phase` | `_single_phase_p2.case` (108×18×36) | — | — | — | — | Turbulent Reichardt | **Running** job 18985538 | New-mesh spin-up; produces `fluid00004.chkp` for all P2 restarts |
| `channel_p2_sigma0` | `_p2_sigma0.case` | — | 0 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P2 spin-up | Phase 2: 3.1 elem/interface, original code; benchmark κ_rms vs P1 |
| `channel_p2_we10` | `_p2_we10.case` | 10 | 0.04 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P2 σ0 result | Phase 2: We=10 first σ>0 test |
| `channel_p2_we1` | `_p2_we1.case` | 1 | 0.4 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P2 σ0 result | Phase 2: We=1 primary production case |
| `channel_p3_single_phase` | `_single_phase_p3.case` (144×18×48) | — | — | — | — | Turbulent Reichardt | **Running** job 18986273 | Phase 3 spin-up (1 node, 128 ranks, 8h); produces P3 `fluid00004.chkp` |
| `channel_p3_sigma0` | `_p3_sigma0.case` | — | 0 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P3 spin-up | Phase 3: 4.1 elem/interface; mesh convergence point vs P1/P2 |
| `channel_p3_we10` | `_p3_we10.case` | 10 | 0.04 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P3 σ0 result | Phase 3: We=10 on finest mesh |
| `channel_p3_we1` | `_p3_we1.case` | 1 | 0.4 | 0.4 | 0 | `fluid00004.chkp` + drop | **Pending** P3 σ0 result | Phase 3: We=1 on finest mesh |

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

The laminar We=1 case was intended as a ground truth (κ_rms = 6.67 throughout), but
blew up at t=0.90 TU due to the same CSF capillary instability. A static (zero-flow)
or semi-implicit CSF case is needed for a clean baseline.

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

## Single-phase spin-up summary — all meshes

All two-phase restart cases require a turbulent checkpoint at t=20 on the same mesh
and with the same MPI rank count as the restart. One spin-up is needed per mesh.

| Mesh | Case file | Run dir | MPI ranks | Checkpoint | Status |
|------|-----------|---------|-----------|------------|--------|
| 81×18×27 (P1) | `_single_phase.case` | `channel_single_phase/` | 16 | `fluid00004.chkp` (t=20) | **Completed** |
| 108×18×36 (P2) | `_single_phase_p2.case` | `channel_p2_single_phase/` | 128 | `fluid00004.chkp` (t=20) | **Running** job 18985538 |
| 144×18×48 (P3) | `_single_phase_p3.case` | `channel_p3_single_phase/` | 128 | `fluid00004.chkp` (t=20) | **Running** job 18986273 |

**Turbulence indicator:** u_max settles to 1.35–1.45 by t≈10 and fluctuates there.
**Verification:** check `ekin.csv` — u_max fluctuating, no trend, no blow-up at t=20.
**MPI rank count for restarts must match spin-up** (mesh partitioning compatibility).

---

## Mesh convergence overview — Phase 1/2/3

The three meshes form a geometric series (factor 1.33× in each xz direction per step),
giving clean convergence data for interface resolution. All phases use identical physics
(same code, ε=0.09 for P2/P3, R=0.4) to isolate the mesh variable.

| Phase | Mesh | Δxz | Δy | Elements | 4ε/Δxz | Δt_cap (We=10) | Δt/Δt_cap |
|-------|------|-----|----|----------|---------|----------------|-----------|
| P1 | 81×18×27 | 0.155 | 0.111 | 39,366 | 1.8 (ε=0.07) | 0.082 | ~0.016 |
| P2 | 108×18×36 | 0.1164 | 0.111 | 69,984 | 3.1 (ε=0.09) | 0.079 | ~0.010 |
| P3 | 144×18×48 | 0.0873 | 0.111 | 124,416 | **4.1** (ε=0.09) | 0.061 | ~0.007 |

**Scientific question:** with the original SEM code (no Fortran fix), does κ_rms error
scale with interface coverage (4ε/Δ), or is the element-face artifact mesh-independent?
If κ_rms decreases from P1 to P3, the finer mesh partially mitigates the problem.
If κ_rms stays ~64 across all meshes, the Fortran fix (extra GS pass) is essential.

**genmeshbox commands:**
```bash
# All meshes: generate on egidius, then transfer to Dardel via dardel-ftn
# P1 (done 2026-03)
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 81 18 27 .true. .false. .true.
# P2 (done 2026-03-26, 16MB on Dardel)
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 108 18 36 .true. .false. .true.
# P3 (done 2026-03-27, 28MB on Dardel)
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 144 18 48 .true. .false. .true.
mv box.nmsh box_phys_144x18x48.nmsh
```

---

## channel_p2_single_phase — Phase 2 turbulent spin-up (RUNNING job 18985538)

**Purpose:** Turbulent spin-up on the 108×18×36 mesh. Produces `fluid00004.chkp` at
t=20 for all Phase 2 two-phase restart cases.

**Setup:** 108×18×36 mesh, N=7, Re_b=2800, 128 MPI ranks (Dardel, 1 node), Reichardt IC
+ perturbations, end_time=25, checkpoints every 5 TU, field output every 5 TU.
Job 18985538 submitted 2026-03-27.

**Expected outputs** in `$SCRATCH_DIR/channel_p2_single_phase/`:
- `ekin.csv` — u_max should settle to 1.35–1.45 by t=10, fluctuate through t=20
- `fluid00004.chkp` at t=20 → copy to each Phase 2 two-phase run directory
- Field snapshots at t=5, 10, 15, 20, 25

| t | E_kin | u_max | Notes |
|---|-------|-------|-------|
| — | — | — | Not yet completed |

---

## channel_p3_single_phase — Phase 3 turbulent spin-up (RUNNING job 18986273)

**Purpose:** Turbulent spin-up on the 144×18×48 mesh. Produces `fluid00004.chkp` at
t=20 for all Phase 3 two-phase restart cases.

**Setup:** 144×18×48 mesh, N=7, Re_b=2800, 128 MPI ranks (Dardel, 1 node, 8h), Reichardt IC
+ perturbations, end_time=25, checkpoints every 5 TU, field output every 5 TU.

**Mesh:** `box_phys_144x18x48.nmsh` generated on egidius and transferred to Dardel
(2026-03-27). Already at `$SRC/box_phys_144x18x48.nmsh`.

**Submitted** as job 18986273, 2026-03-27. Script: `cluster/job_channel_p3_single_phase.sh`.

| t | E_kin | u_max | Notes |
|---|-------|-------|-------|
| — | — | — | Not yet run |

---

## Phase 2 two-phase cases — 108×18×36 mesh (PENDING spin-up)

All three cases restart from `channel_p2_single_phase/fluid00004.chkp` (128 ranks).
Run in order: sigma0 first (diagnostic, 1 TU), then we10, then we1.
Job scripts: `cluster/job_channel_p2_sigma0.sh`, `_we10.sh`, `_we1.sh`.

### channel_p2_sigma0 — σ=0 CDI diagnostic (PENDING)

**Purpose:** Measure κ_rms on the P2 mesh (3.1 elements across interface). Compare with
Phase 1 baseline (κ_rms ≈ 64 on 81×18×27). If element-face artifact scales with Δ,
expect lower κ_rms here. If κ_rms stays ~64, the Fortran fix is essential.

**Setup:** ε=0.09, γ=0.05, σ=0, R=0.4, end_time=21.0 (1 TU). Target κ_rms ≈ 2/R = 5.0 (spherical).

| t | φ_max | κ_rms | Notes |
|---|-------|-------|-------|
| — | — | — | Not yet run |

### channel_p2_we10 — We=10 (PENDING)

**Setup:** ε=0.09, γ=0.05, σ=0.04, R=0.4, We=10, end_time=25.0.
Δt_cap≈0.079, Δt/Δt_cap≈0.051. Phase 1 baseline blew up at t≈20.44.

| t | φ_max | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|
| — | — | — | — | Not yet run |

### channel_p2_we1 — We=1 (PENDING)

**Setup:** ε=0.09, γ=0.05, σ=0.4, R=0.4, We=1, end_time=25.0.
Δt_cap≈0.025, Δt/Δt_cap≈0.16. Target production case.

| t | φ_max | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|
| — | — | — | — | Not yet run |

---

## Phase 3 two-phase cases — 144×18×48 mesh (PENDING spin-up)

All three cases restart from `channel_p3_single_phase/fluid00004.chkp` (128 ranks).
Job scripts: `cluster/job_channel_p3_sigma0.sh`, `_we10.sh`, `_we1.sh`.

### channel_p3_sigma0 — σ=0 CDI diagnostic (PENDING)

**Purpose:** κ_rms convergence point at finest mesh (4.1 elements across interface).
Compare with P1 (κ_rms≈64, 1.8 elem) and P2 result.

**Setup:** ε=0.09, γ=0.05, σ=0, R=0.4, end_time=21.0 (1 TU).

| t | φ_max | κ_rms | Notes |
|---|-------|-------|-------|
| — | — | — | Not yet run |

### channel_p3_we10 — We=10 (PENDING)

**Setup:** ε=0.09, γ=0.05, σ=0.04, R=0.4, We=10, end_time=25.0.
Δt_cap≈0.061 (at Δx=0.0873), Δt/Δt_cap≈0.042.

| t | φ_max | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|
| — | — | — | — | Not yet run |

### channel_p3_we1 — We=1 (PENDING)

**Setup:** ε=0.09, γ=0.05, σ=0.4, R=0.4, We=1, end_time=25.0.
Δt_cap≈0.019 (at Δx=0.0873), Δt/Δt_cap≈0.13.

| t | φ_max | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|
| — | — | — | — | Not yet run |

---

## channel_test_v4 — high-We reference (We=730, COMPLETED t=0–5)

**Purpose:** High-We reference case. Surface tension negligible; drop deforms freely.
Documents CDI performance under strong interface straining.

**Setup:** ε=0.05, γ=0.015, σ=4.1×10⁻⁴, R=0.3, Re_b=2800, turbulent IC, T=5 TU.

**Capillary timestep stability:** $\Delta t_{\mathrm{cap}} = \sqrt{0.016^3/(2\pi \times 4.1\!\times\!10^{-4})} \approx 0.040$ TU. With `target_cfl=0.4` and $u_{\max}\approx1.5$, $\Delta t \approx 0.003$ TU → $\Delta t/\Delta t_{\mathrm{cap}} \approx 0.08$ — well inside boundary. CSF is negligible at We=730 so capillary waves are not a concern.

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

![channel_test_v4 diagnostics](figures/diagnostics_ekin_channel_test_v4.png)

---

## channel_test_laminar — seeding test (We=1, BLOWN UP t=0.90)

**Purpose:** Test whether turbulent seeding is necessary for the CSF instability.
Same parameters as we1 but Reichardt IC with no velocity perturbations.

**Setup:** Reichardt IC (no perturbations, turbulent\_ic=false). ε=0.07, γ=0.05,
σ=0.3 (We=1), R=0.3, Re_b=2800, end_time=10.

**Capillary timestep stability:** $\Delta t_{\mathrm{cap}} \approx 0.00059$ TU at $\Delta x_{\mathrm{eff}}=0.0087$. With u_max≈1.15 (Reichardt, no perturbations), $\Delta t \approx 0.0015$ TU → $\Delta t/\Delta t_{\mathrm{cap}} \approx 2.9$ — worse than turbulent v2 (ratio=2.2), yet blow-up took twice as long (0.90 TU vs 0.44 TU). The longer survival directly measures the seeding effect of turbulent fluctuations.

**Blow-up trace:**

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 0.000 | 0.981 | 0.000 | 6.097 | 1.149 | IC; κ_rms ≈ 2/R ✓ |
| 0.173 | 0.985 | 0.000 | 6.105 | 1.207 | Flat — no turbulent seeding |
| 0.433 | 0.987 | 0.000 | 6.648 | 1.227 | Slow growth above 6.67 begins |
| 0.606 | 0.988 | 0.000 | 8.089 | 1.246 | Growth accelerating |
| 0.751 | 0.989 | 0.000 | 10.51 | 1.507 | Flow relaminarised to Poiseuille |
| 0.840 | 0.990 | 0.000 | 22.02 | 2.180 | Runaway |
| 0.868 | 0.991 | 0.000 | 40.92 | 3.285 | Explosive |
| 0.893 | 0.992 | 0.000 | 158.8 | 14.70 | Plateau; φ_max still < 1 — CDI intact ✓ |
| 0.896 | 1.061 | −0.192 | 179.4 | 31.17 | φ_max > 1: diverged |

**Key observations:**
- φ_max < 1 throughout the κ_rms runaway — CDI intact, same as all turbulent blow-ups
- Slow incubation (~0.5 TU flat, vs ~0.26 TU for turbulent v2): numerical round-off seeds more slowly than turbulent fluctuations
- After seeding, the growth pattern is identical: same plateau (~155–180), same mechanism
- Flow relaminarised from Reichardt (u_max≈1.15) to Poiseuille (u_max≈1.5) around t=0.75 — this accelerated the blow-up slightly but was not the trigger

**Conclusion:** Turbulent seeding is not necessary — numerical round-off alone triggers the instability. Turbulent fluctuations roughly halve the stable duration. The exact stability boundary is not yet established; We=10 (ratio $= 0.69 < 1$) also blew up.

---

## channel_test_we1 — primary validation (We=1, PLANNED)

**Purpose:** CDI/CSF under turbulent straining with strong surface tension. The CSF
force is large and measurable: Δp = 2σ/R = 2.0 ≫ ρU_b²/2 = 0.5. Any curvature
error amplifies into a visible spurious velocity. Comparison with the laminar case
isolates what turbulence adds.

**Setup:** ε=0.07, γ=0.05, σ=0.3 (We=1), R=0.3, Re_b=2800, restart from
`fluid00004.chkp`, end_time=25 (runs t=20→25).

**Capillary timestep stability:** $\Delta t_{\mathrm{cap}} \approx 0.00059$ TU ($\Delta x_{\mathrm{eff}} = 0.0087$); $\Delta t \approx 0.00130$ TU → $\Delta t/\Delta t_{\mathrm{cap}} \approx 2.2$ — well above 1. Two prior restart attempts (v1, v2, different R) both blew up. See channel_test_restart blow-up analysis and ANALYSIS.md §4.4 for details.

| t | φ_max | φ_min | κ_rms | u_max |
|---|-------|-------|-------|-------|
| — | — | — | — | — |

---

## channel_test_we10 — moderate deformation (We=10, BLOWN UP t=20.44)

**Purpose:** First turbulent CDI/CSF test with a We=10 capillary margin. Blew up at
t≈20.44 TU (~0.44 TU after injection) with the same κ_rms runaway as the We=1 cases.

**Setup:** ε=0.07, γ=0.05, σ=0.03 (We=10), R=0.3, Re_b=2800, restart from
`fluid00004.chkp`, end_time=25 (runs t=20→25). 16 MPI ranks.

**Capillary timestep stability** (see ANALYSIS.md §4.4):

| Quantity | Value |
|----------|-------|
| $\Delta t_{\mathrm{cap}} = \sqrt{\Delta x_{\mathrm{eff}}^3 / (2\pi\sigma)}$, $\Delta x_{\mathrm{eff}} = 0.0087$ | **0.00188 TU** |
| $\Delta t_{\mathrm{cap}}$ (element-average, $\Delta x = 0.016$, naive estimate) | 0.00466 TU |
| $\Delta t$ (observed, `target_cfl=0.2`) | 0.00130 TU |
| $\Delta t / \Delta t_{\mathrm{cap}}$ (corrected, $\Delta x_{\mathrm{eff}}$) | **0.69** |
| $\Delta t / \Delta t_{\mathrm{cap}}$ (naive, element-average) | 0.28 — incorrectly predicted stable |

The naive element-average estimate predicted ratio=0.28 (stable), but the correct
$\Delta x_{\mathrm{eff}}$-based ratio is 0.69. The ratio is below 1 but the run still
blew up. The stability boundary is not yet established.

**Blow-up trace:**

| t | φ_max | φ_min | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 20.002 | 0.981 | 0.000 | 6.10 | 1.341 | Drop injected; κ_rms ≈ 2/R = 6.67 ✓ |
| 20.132 | 0.986 | 0.000 | 6.74 | 1.341 | Stable; φ_max intact |
| 20.263 | 0.987 | 0.000 | 16.54 | 1.426 | Growth onset |
| 20.317 | 0.988 | 0.000 | 27.64 | 1.711 | Crosses 2/ε=28.6 threshold |
| 20.358 | 0.989 | 0.000 | 52.0 | 2.203 | Runaway |
| 20.389 | 0.991 | 0.000 | 93.4 | 3.262 | Explosive |
| 20.423 | 0.996 | 0.000 | 153.8 | 4.667 | Plateau / resolution saturation |
| 20.435 | 0.998 | 0.000 | 158.1 | 4.907 | φ_max still < 1 — CDI intact ✓ |

**Key observation:** φ_max < 1 throughout the blow-up — CDI is functioning correctly.
The instability is in the explicit CSF treatment, not CDI. Identical signature to the
We=1 blow-up cases.

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

## channel_test_we10_diag — dense blow-up diagnostic (We=10, BLOWN UP t=20.447)

**Purpose:** Dense snapshots (every 0.02 TU) of the We=10 turbulent blow-up, for
animation and supervisor presentation. Confirms the same blow-up sequence with full
temporal resolution.

**Setup:** Identical to `channel_test_we10` except `output_value: 0.02` and
`end_time: 20.5`. 20 field files, t=20.02→20.50. Run blew up at t≈20.447 and
reached `end_time` while diverged.

**Blow-up trace:**

| t | φ_max | κ_rms | u_max | Notes |
|---|-------|-------|-------|-------|
| 20.002 | 0.981 | 6.10 | 1.341 | Drop injected; κ_rms ≈ 2/R ✓ |
| 20.132 | 0.986 | 6.74 | 1.341 | Stable |
| 20.263 | 0.987 | 16.5 | 1.426 | Growth onset |
| 20.317 | 0.988 | 27.6 | 1.711 | Crossing 2/ε threshold |
| 20.358 | 0.989 | 52.0 | 2.203 | Runaway |
| 20.435 | 0.998 | 158 | 4.907 | Plateau; φ_max < 1 — CDI intact ✓ |
| 20.447 | **1.002** | 158 | 7.095 | **φ_max > 1: CDI overloaded** |
| 20.458 | 1.015 | 162 | 6.144 | Fully diverged |

Identical to original `channel_test_we10` blow-up sequence. Animation:
`blowup_channel_test_we10_diag.gif` (20 frames, 3 fps, 3 panels: φ/κ/|u|).

---

## channel_test_sigma0 — CDI quality test (σ=0, TERMINATED t=21.24)

**Purpose:** Isolate CDI from CSF. With σ=0 the momentum equation has no surface tension
force. The scalar evolves under CDI only. Tests whether CDI maintains interface sharpness
and whether κ is accurate under turbulent straining, independently of CSF timestep issues.

**Motivation:** φ_max < 1 in all blow-up cases confirms CDI maintains interface amplitude,
but does not confirm accurate normals $\hat{\mathbf{n}} = \nabla\varphi/|\nabla\varphi|$.
Inaccurate normals give inaccurate $\kappa = -\nabla\cdot\hat{\mathbf{n}}$, which would
amplify the CSF force even if the timestep constraint is the primary cause.

**Setup:** ε=0.07, γ=0.05, σ=0.0, R=0.3, Re_b=2800, restart from `fluid00004.chkp`,
end_time=25. 16 MPI ranks. Terminated early at t=21.24 (1.24 TU after injection).

**Diagnostic trace:**

| t | φ_max | κ_rms | κ_max | u_max | Notes |
|---|-------|-------|-------|-------|-------|
| 20.002 | 0.981 | 6.1 | 125 | 1.341 | Drop injected; κ_rms ≈ 2/R ✓ |
| 20.067 | 0.985 | 6.3 | 254 | 1.341 | κ_max already elevated |
| 20.132 | 0.986 | 7.8 | 558 | 1.341 | κ_rms above spherical value |
| 20.197 | 0.987 | 19.9 | 756 | 1.342 | Rapid growth |
| 20.263 | 0.987 | 38.7 | 867 | 1.342 | Continues rising |
| 20.393 | 0.988 | 59.0 | 808 | 1.344 | Approaching plateau |
| 20.588 | 0.989 | 63.9 | 858 | 1.352 | **Peak; ~9.6× spherical value** |
| 20.784 | 0.988 | 61.0 | 816 | 1.359 | Slow decline |
| 21.044 | 0.989 | 55.5 | 812 | 1.356 | Continues declining |
| 21.240 | 0.986 | 51.3 | 859 | 1.355 | Terminated |

**Key observations:**

- **u_max and Ekin unaffected throughout** — no blow-up without CSF. CDI alone does not destabilise the flow.
- **φ_max stable at 0.986–0.989** — CDI maintains interface amplitude.
- **κ_rms spikes from 6.1 to ~64 in 0.4 TU**, then slowly declines. Expected spherical value: 6.67.
- **κ_max reaches 800+** — large point-wise curvature values at highly deformed interface regions.
- The peak κ_rms ≈ 64 is ~9.6× the spherical reference. This reflects both genuine drop deformation (no surface tension to restore shape) and potentially inaccurate normal computation in stretched regions.

**Field snapshot analysis (March 2026, `analyze_sigma0_normals.py`):**

Element-local postprocessing of field0.f00041 (t=20.501) and field0.f00042 (t=21.000):

| Metric | t=20.501 | t=21.000 | Reference |
|--------|----------|----------|-----------|
| Drop centroid displacement | 0.605 (streamwise) | 1.210 (streamwise) | 0 = injection point |
| $\|\nabla\varphi\|_{\rm mean}$ (interface) | 3.927 | 4.372 | 3.571 (CDI theory) |
| Sphericity $r_{\rm std}$ on $\varphi=0.5$ | 0.0290 | 0.0523 | 0 (perfect sphere) |
| Volume $\int\varphi\,dV$ | 0.20965 | 0.20962 | constant |
| $\hat{\mathbf{n}}$ alignment angle | 25.4° | 34.2° | 0° (radial) |
| $\kappa_{\rm rms}$ (postprocess, element-local) | 25.8 | 24.2 | 6.67 (sphere) |
| $\kappa_{\rm rms}$ (Neko, GS-averaged) | ~62 | ~57 | 6.67 (sphere) |

**Interpretation:**

1. **CDI is functioning correctly.** Interface gradient is above (not below) the CDI theory
   value, volume is conserved to 4 significant figures, and φ_max is stable. The CDI is
   not under-performing.

2. **Drop advects at ~U_cl.** Centroid displaced 0.605 units streamwise in 0.5 TU → ~U_cl =
   1.21 U_b (physically correct for a centreline drop). No wall-normal or spanwise drift.

3. **Drop is nearly spherical at t=20.5** (r_std = 0.029 = 9.7%R). Mild deformation by
   t=21.0 (r_std = 0.052 = 17%R). Genuine σ=0 deformation is present but too small to
   account for κ_rms = 62 (which is 9.3× the spherical value).

4. **Element-boundary artifact is the dominant contribution to κ_rms.** Evidence:
   - κ_max = 124 at t=20.002 (before any turbulent deformation), rising to 800+ by t=20.067
   - Element-local postprocessing gives κ_rms = 26, GS-averaged Neko gives κ_rms = 62
     (factor ~2.4): the inter-element gradient jumps substantially boost the Neko estimate
   - The rising κ_rms (6 → 64 over 0.4 TU) likely reflects the drop sweeping through more
     element faces as it advects streamwise at U_cl

5. **Normal misalignment of 25–34° relative to drop centroid** is large for a nearly
   spherical drop and points to element-boundary gradient errors in n̂.

**What this means for the CSF blow-ups:**

The CSF force is driven by the GS-averaged κ, which includes both genuine curvature and
element-boundary artifacts (factor ~2.4 above element-local). At We=10 (σ=0.03):
$F_{\mathrm{ST}} \approx 0.03 \times 62 / 0.14 \approx 13$ — much larger than the spherical
estimate of $0.03 \times 6.67 / 0.14 \approx 1.4$. Both the genuine σ=0 drop deformation
and the SEM element-boundary artifact feed into the capillary instability. The stability
analysis in §7.4–7.6 used κ ≈ 6.67 (spherical); the actual driving κ is ~9× larger.

**High-frequency diagnostic run:** `channel_test_sigma0_diag/` — 20 snapshots at 0.02 TU intervals, t=20.12→20.50. See next section.

---

## channel_test_sigma0_diag — diagnostic confirmation (COMPLETED t=20.12→20.50)

**Purpose:** 20 snapshots at 0.02 TU to track κ growth from the first element-face crossing.

**Key time-series (from `analyze_sigma0_normals.py`):**

| $t$ | centroid disp. | $\hat{\mathbf{n}}$ misalign | $\kappa_\text{rms}$ (Python) | tanh L2 | width |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 20.12 | 0.14 | 4.3° | 6.0 | 0.014 | 0.304 |
| 20.22 | 0.27 | 8.0° | 11.9 | 0.018 | 0.301 |
| 20.32 | 0.41 | 16.0° | 20.1 | 0.025 | 0.298 |
| 20.42 | 0.53 | 23.2° | 26.7 | 0.031 | 0.296 |
| 20.50 | 0.61 | 25.4° | 25.8 | 0.036 | 0.293 |

**Confirmed findings:**
- At $t=20.12$: minimal advection (0.14 units), n̂ misalignment only 4.3°, κ_rms=6.0 ≈ spherical. Very clean start.
- κ_rms grows **linearly with centroid displacement**, not with drop deformation (r_std grows from 0.010 to 0.029 — small)
- Interface width starts 8.6% too wide (0.304 vs theory 0.280) — CDI with γ=0.05 is slow (τ_CDI=1.4 TU)
- tanh L2 residual grows monotonically — CDI cannot re-sharpen as fast as turbulence strains the interface

**Element-face hypothesis confirmed.** κ rises as the drop sweeps through more element faces at U_cl ≈ 1.2 U_b.

---

## channel_test_sigma0_gamma025 — CDI parameter test (COMPLETED t=20.00→20.26)

**Purpose:** Test whether increasing γ from 0.05 to 0.25 (τ_CDI = 0.07/0.25 = 0.28 TU) improves κ accuracy.

**Result: κ is immediately 9× worse.**

| $t$ | $\kappa_\text{rms}$ ($\gamma=0.05$) | $\kappa_\text{rms}$ ($\gamma=0.25$) |
|:---:|:---:|:---:|
| 20.002 | 6.1 | 6.1 |
| 20.07 | 7.8 | **71.3** |
| 20.13 | 7.8 | **85.1** |
| 20.20 | 19.9 | **82.0** |
| 20.26 | 19.9 | **79.5** |

The jump occurs at the **first output step** — before turbulence has time to act.

**Reason:** stronger CDI compression drives the interface to a sharper steady-state (higher |∇φ|). The C0-kink amplitude at element faces scales with |∇φ|, so the artifact is amplified. φ_max recovers faster (CDI genuinely works better), but κ accuracy is catastrophically worse.

**Conclusion:** γ tuning trades CDI profile quality for κ accuracy. With the current κ computation, there is no γ that satisfies both. The only fix is the extra GS pass on n̂ in the Fortran code.

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

Both We=10 (`channel_test_we10`) and laminar We=1 (`channel_test_laminar`) have since
been run — both blew up. See their respective sections above. The explicit CSF
instability occurs across all We values tested so far and in the absence of turbulence.

---

## Next steps

**Phase 1 — confirmed (current mesh, egidius):**
- [x] All CSF cases blow up from explicit capillary instability (We=1, 10, laminar, we10_diag)
- [x] CDI functioning correctly throughout all blow-ups (φ_max < 1, volume conserved)
- [x] σ=0 diagnostic: κ inaccuracy dominated by SEM element-boundary artifact (GS-averaged n̂ → C0 kink → D[N,N]=14 amplification → κ_face~23)
- [x] γ tuning (0.05→0.25): CDI profile improves but κ immediately 9× worse — not the path forward
- [x] Animations produced: `blowup_channel_test_laminar.gif`, `blowup_channel_test_sigma0_diag.gif`, `blowup_channel_test_we10_diag.gif`

**Phase 1 — remaining:**
- [ ] **Fortran fix:** extra GS pass on n̂ after normalisation in `turb_channel_two_phase.f90` (~5 lines). Smooths C1 kink before div(n̂).
- [ ] Verify fix: σ=0 run should give κ_rms ≈ 6.67 (spherical), no growth with advection
- [ ] If fix verified: run We=10 with σ>0 on current mesh

**Phase 2 — in progress (108×18×36, Dardel):**
- [x] Dardel setup: clone, build `neko-channel`, generate mesh
- [x] Case files created: `p2_sigma0`, `p2_we10`, `p2_we1`
- [ ] P2 spin-up complete (job 18985538 running) → check `fluid00004.chkp` at t=20
- [ ] Submit `channel_p2_sigma0` — key question: does κ_rms improve from P1 baseline (~64)?
- [ ] Submit `channel_p2_we10` and `channel_p2_we1` pending σ=0 result

**Phase 3 — planned (144×18×48, Dardel):**
- [x] Case files created: `p3_sigma0`, `p3_we10`, `p3_we1`, `p3_single_phase`
- [x] Job script: `cluster/job_channel_p3_single_phase.sh` (2 nodes, 256 ranks)
- [ ] Generate mesh: `genmeshbox ... 144 18 48 ...` on Dardel login node
- [ ] Verify 2-node job policy on Dardel (naiss2025-3-39 account)
- [ ] Submit P3 spin-up → produces `fluid00004.chkp` for P3 restarts
- [ ] Submit P3 two-phase cases

**Key open question (P1/P2/P3 σ=0 comparison):**
Does κ_rms follow: P1≈64 → P2≈? → P3≈? (decreasing → mesh helps) or stays ~64 (Fortran fix essential)?
