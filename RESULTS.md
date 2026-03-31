# Simulation Results

Convergence series: isotropic meshes L1/L2/L3/L4 with constant ε/Δxz≈0.457.
All runs on Dardel (NAISS account `naiss2025-3-39`), results synced to egidius at
`/lscratch/sieburgh/simulations/<run_name>/`.

See CLAUDE.md for the naming convention mapping (L1/L2/L3 ↔ p2/p3/l3 run prefixes).

## Current status (as of 2026-03-31)

| Category | Status |
|----------|--------|
| L1/L2/L3 single-phase spin-ups | ✓ COMPLETED — turbulence validated |
| L1/L2/L3 σ=0 CDI diagnostic | ✓ COMPLETED and postprocessed |
| L4 single-phase spin-up | NOT STARTED (planned) |
| L4 σ=0 | NOT STARTED (planned, after L4 spin-up) |
| L1/L2/L3 We=10, We=1 | **BLOCKED — awaiting CDI normal fix** |

**Immediate next steps (in order):**
1. **γ sensitivity test**: run L2 σ=0 with γ=0.2 (Γ*=0.14) — rules out low γ as root cause.
   Single case file change, ~3h on Dardel, no code modifications needed.
2. **2D Couette CDI subproject**: clean 2D testbed for diagnosing and fixing normal computation.
3. Return to L1/L2/L3 We cases once CDI normals are resolved.

---

## Single-phase spin-ups

All spin-ups run t=0→25, checkpoint `fluid00004.chkp` written at t=20.
Validation: u_max (last 5 rows of ekin.csv) fluctuating ~1.35–1.45 → TURBULENT.

| Mesh | Run name | Job | Nodes/Ranks | Status | u_max (last 5) | Notes |
|------|----------|-----|-------------|--------|----------------|-------|
| 81×18×27 (baseline) | `channel_single_phase` | — | 16 ranks | COMPLETED | ~1.38 | egidius local run |
| L1: 108×18×36 | `channel_p2_single_phase` | 18985538 | 1/128 | COMPLETED | ~1.38 | |
| L2: 144×24×48 | `channel_p3_single_phase` | 19002586 | 2/256 | COMPLETED | 1.37±0.01 | 2 nodes (1 node OOM) |
| L3: 192×32×64 | `channel_l3_single_phase` | 19003203 | 4/512 | COMPLETED | 1.34±0.00 | 4 nodes (2 nodes OOM) |
| L4: 288×48×96 | `channel_l4_single_phase` | — | ~14/1792 | NOT STARTED | — | planned; ~14 nodes est. |

---

## L1 (108×18×36) — σ=0 CDI diagnostic, COMPLETED and postprocessed

ε=0.053, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25). 128 ranks, 1 node.

| Case | Run name | Job | Status | κ_rms plateau | Postprocessed |
|------|----------|-----|--------|---------------|---------------|
| σ=0 ε=0.053 | `channel_p2_sigma0_eps053` | 19003682 | COMPLETED t=20→25 | ~56 (t≈22) | ✓ time-series, snapshots |
| σ=0 ε=0.09 | `channel_p2_sigma0` | 19003683 | COMPLETED t=20→25 | ~50 (t≈22) | ✓ time-series only |

**Note on L1 We cases:** Early We=10 and We=1 runs (`channel_p2_we10`, `channel_p2_we1`,
jobs 19003684–19003685) were submitted 2026-03-28 but are considered stale:
`channel_p2_we1` blew up immediately at t≈20.22 (capillary instability, same as baseline
mesh); `channel_p2_we10` was killed at t≈22.3. These runs are de-prioritised — the
canonical We cases will be run at L2/L3 once the CDI kink artifact is resolved.

---

## L2 (144×24×48) — σ=0 CDI diagnostic, COMPLETED and postprocessed

ε=0.04, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25). 256 ranks, 2 nodes.

| Case | Run name | Job | Status | κ_rms plateau | Postprocessed |
|------|----------|-----|--------|---------------|---------------|
| σ=0 | `channel_p3_sigma0` | 19050210 | COMPLETED t=20→25 (3h35) | ~54–56 (t≈23) | ✓ time-series, snapshots, normals |
| σ=0 _p2 | `channel_p3_sigma0_p2` | 19050211 | CANCELLED | — | — |
| We=10 | `channel_p3_we10` | — | NOT SUBMITTED | — | Awaiting CDI fix |
| We=1 | `channel_p3_we1` | — | NOT SUBMITTED | — | Awaiting CDI fix |

Note: `channel_p3_sigma0_p2` was cancelled after diagnosing the _p2 extra GS as a no-op
(see σ=0 diagnosis section below).

---

## L3 (192×32×64) — σ=0 CDI diagnostic, COMPLETED and postprocessed

ε=0.03, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25). 512 ranks, 4 nodes.

| Case | Run name | Job | Status | κ_rms plateau | Postprocessed |
|------|----------|-----|--------|---------------|---------------|
| σ=0 | `channel_l3_sigma0` | 19050366 | COMPLETED t=20→25 (5h51) | ~67–68 (t≈24) | ✓ time-series, snapshots, normals |
| We=10 | `channel_l3_we10` | — | NOT SUBMITTED | — | Awaiting CDI fix |
| We=1 | `channel_l3_we1` | — | NOT SUBMITTED | — | Awaiting CDI fix |

---

## L4 (288×48×96) — planned, not yet started

ε=0.02, R=0.4, γ=0.05. 1,327,104 elements, ε/Δxz=0.458, ε/Δ_GLL=3.2, R/ε=20.
Factor 3/2 from L3 (the 4/3 progression gives nz=85.3, non-integer, so 3/2 used instead).
Estimated ~14 nodes / 1792 ranks.

Workflow (when ready):
1. Generate mesh on egidius: `genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 288 48 96 .true. .false. .true.`
2. Write spin-up case file and job script (follow L3 pattern)
3. Run spin-up, validate turbulence
4. Submit σ=0 case; defer We cases until CDI kink is resolved

| Case | Run name | Job | Status |
|------|----------|-----|--------|
| Single-phase spin-up | `channel_l4_single_phase` | — | NOT STARTED |
| σ=0 | `channel_l4_sigma0` | — | NOT STARTED |
| We=10 | `channel_l4_we10` | — | NOT SUBMITTED — awaiting CDI fix |
| We=1 | `channel_l4_we1` | — | NOT SUBMITTED — awaiting CDI fix |

---

## σ=0 CDI diagnostic — findings

σ=0 cases use CDI only (no surface tension). κ_rms measures the curvature residual
from the SEM discretisation. Expected reference value: κ_sphere = 2/R = 5.

### Non-dimensional CDI parameter

The non-dimensional compression parameter is Γ* = γ/u_max. With γ=0.05 and
u_max≈1.38 (instantaneous domain maximum from ekin.csv): **Γ* ≈ 0.036** across
all mesh levels (same γ, same flow). This is less than 1, meaning CDI compression
is ~27× slower than the peak local velocity. The interface is maintained because
the CDI acts on the interface thickness scale ε, not the bulk scale — but Γ*
sets the margin. It is the same value for L1, L2, L3 since γ is not rescaled
with ε in the current setup.

### κ_rms convergence results

| Mesh | ε | κ_rms plateau | κ_rms/(2/R) | Plateau reached | Notes |
|------|---|---------------|-------------|-----------------|-------|
| L1 (ε=0.09) | 0.09 | ~50 | ~10× | t≈22 | Large ε, informative only |
| L1 (ε=0.053) | 0.053 | ~56 | ~11.2× | t≈22 | Convergence series point |
| L2 (ε=0.04) | 0.040 | ~54–56 | ~10.8–11.2× | t≈23 | Slight late-time creep |
| L3 (ε=0.03) | 0.030 | ~67–68 | ~13.4–13.6× | t≈24 | Late-time creep similar to L2 |

**Key finding:** κ_rms does NOT decrease with mesh refinement at constant ε/Δxz.
L3 (finer mesh) reaches a HIGHER plateau than L2. Both L2 and L3 show slow upward
creep after the plateau — consistent with turbulence gradually deforming the interface
and increasing kink amplitude in dynamic equilibrium.

### Kink root cause

The artifact is **intra-element**, not inter-element. After GS(∇φ)+normalize, n̂ at
shared face nodes is already identical between elements — so the extra GS on n̂ (the
_p2 "fix") is a no-op: it sums identical values and divides back.

The actual kink is in the Lagrange polynomial connecting the averaged face-node n̂ to
the element-local first-interior n̂. The SEM endpoint derivative amplification
D[N,N]≈14 (N=7) amplifies this jump. Since the kink curvature scales as ~1/Δx ∝ 1/ε
at constant ε/Δx, **mesh refinement at constant ε/Δxz makes κ_rms worse, not better**.

CDI convergence in the literature refers to interface shape accuracy (φ), not the
computed curvature via div(n̂). The SEM derivative amplification means div(n̂)-based
κ will not converge to κ_sphere regardless of mesh refinement without a different
curvature scheme (height-function, parabolic fit, or pre-smoothed n̂).

### Potential fixes (not yet implemented)

- Repeated GS on ∇φ before normalisation — propagates the face average inward one
  node-ring per pass (requires multiple passes per time step)
- Helmholtz-type smoother on n̂ — filters the intra-element kink
- Fundamentally different curvature scheme — height-function or parabolic fit

### Mesh isotropy analysis

The turbulent channel mesh has a systematic 4.51% anisotropy in the wall-normal
direction: Δy/Δxz ≈ 0.9549 across all mesh levels (consequence of the domain
aspect ratio 2/(4π/3 / (NZ/3)) — see table below).

| Level | Δx=Δz | Δy | Δy/Δxz |
|-------|-------|----|--------|
| L1 | 0.1164 | 0.1111 | 0.9549 |
| L2 | 0.0873 | 0.0833 | 0.9549 |
| L3 | 0.0654 | 0.0625 | 0.9549 |

For a 45° interface (worst case), this anisotropy rotates the computed normal by
**~1.3°** from the true direction. This is ~7% of the early-time kink artifact
(δ_mean≈16–18° for the spherical drop). The anisotropy is a real but minor
contributor — not the dominant source of normal error.

**Implication for new subprojects:** the 2D Couette CDI study should use a perfectly
isotropic mesh (Δx = Δy) to eliminate this confounding factor entirely.

### Normal field diagnostics (generated 2026-03-31)

Quantitative diagnostics run via `postprocess_normals.py` on L2 and L3 σ=0 runs.

**Diagnostic 1 — Angular deviation δ = arccos(n̂_computed · n̂_ideal)**

n̂_computed = ∇φ/|∇φ| from element-local gradient. n̂_ideal = inward radial from
drop centroid. **This reference is only geometrically correct for a spherical drop.**
Once the drop deforms, the radial direction from the centroid is no longer the true
interface normal, so δ at later times conflates two effects: (a) genuine SEM kink
error and (b) drop non-sphericity. Only the t=20.5 result (near-spherical drop) is
an unambiguous measure of the SEM normal error.

| Run | t=20.5 TU (spherical — valid) | t=23.0 TU (deforming — mixed) | t=25.0 TU (deformed — confounded) |
|-----|-------------------------------|-------------------------------|-----------------------------------|
| L2 ε=0.04 | mean=18.4°, p90=41° | mean=39.7°, p90=85° | mean=62.1°, p90=115° |
| L3 ε=0.03 | mean=15.5°, p90=33° | mean=38.0°, p90=76° | mean=57.6°, p90=104° |

**What t=20.5 tells us:** For the nearly spherical drop immediately after restart,
the SEM already produces a mean normal error of ~16–18° (p90≈33–41°). This is the
baseline kink artifact before any interface deformation. A perfectly working CDI would
give ~0°. The 16–18° is **not** caused by the drop being non-spherical — it is the
genuine SEM endpoint derivative kink.

**What the late-time data tells us:** The late-time growth of δ cannot be cleanly
attributed to worsening normals alone (the drop shape changes significantly). However,
the φ profile width (Diagnostic 2) does not have this ambiguity.

**Diagnostic 2 — φ interface profile width (10-90% width, vertical cut through drop top)**

Expected ideal width = 4ε·arctanh(0.8) ≈ 4.4ε. This diagnostic is robust to drop
shape — it only checks whether the interface stays sharp, regardless of geometry.

| Run | ε | t=20.5 TU | t=23.0 TU | t=25.0 TU |
|-----|---|-----------|-----------|-----------|
| L2 ε=0.04 | 0.04 | 4.17ε | 4.17ε | **6.04ε** |
| L3 ε=0.03 | 0.03 | 3.77ε | 3.33ε | **10.24ε** |

**Key finding:** Interface broadens significantly by t=25 on both meshes — more
severely on L3. This is consistent with CDI sharpening failing. However, this
diagnostic alone cannot separate low γ from bad n̂ as the cause.

**What the two diagnostics together confirm:**
- The SEM kink artifact produces ~16° mean normal error even for a spherical drop (D1, t=20.5).
- The φ profile broadens significantly by the end of the run (D2).
- The diagnostic is insufficient to definitively rule out γ being a contributing cause.

### Open questions and recommended next steps

**Question 1 — Is γ too low?**

Γ* = γ/u_max ≈ 0.036 is quite low. A γ sensitivity test (single L2 run with γ=0.2
or 0.5, same case file) would directly answer this: if the profile stays sharp, γ is
the bottleneck; if it still broadens, n̂ is the root cause. This is the cheapest
possible test (no code changes, ~3h on Dardel). **This should be done before
investing in more complex normal smoothing strategies.**

**Question 2 — Can we get a clean normal-quality diagnostic for a non-spherical drop?**

The angular deviation diagnostic breaks down once the drop deforms. A 2D Couette flow
case (drop in linear shear, steady deformation, no turbulence) would provide a
controlled geometry where the true interface shape is known analytically at steady
state. This is the motivation for the planned CDI normals 2D subproject.

**Question 3 — Can the SEM kink be fixed?**

Possible strategies: repeated GS on ∇φ, Helmholtz smoother on n̂, height-function
curvature. The 2D Couette subproject is the right testbed for these fixes before
applying them to the turbulent channel.

**Recommended order:**
1. γ sensitivity test on L2 (cheap, answers open question 1).
2. 2D Couette CDI normals subproject (clean diagnostic, fix candidates).
3. Return to turbulent channel We>0 cases once CDI normals are resolved.

### Postprocessing figures (generated 2026-03-31)

| Figure | Content |
|--------|---------|
| `sigma0_timeseries_channel_{p2_sigma0,p2_sigma0_eps053,p3_sigma0,l3_sigma0}.png` | κ_rms(t), φ extrema, E_kin(t) |
| `sigma0_snapshots_channel_{p2_sigma0_eps053,p3_sigma0,l3_sigma0}.png` | φ/\|u\|/n̂_y/κ rows, 3 snapshots |
| `sigma0_normals_channel_{p3_sigma0,l3_sigma0}.png` | Zoom+quiver: n̂ field around drop top |
| `blowup_channel_{p2_sigma0_eps053,p3_sigma0,l3_sigma0}.gif` | Animated φ/κ/\|u\|, 10 frames |
| `three_meshes_sigma0.gif` | L1/L2/L3 φ stacked, element grid overlaid, 5 frames |
| `normals_angular_dev_channel_{p3_sigma0,l3_sigma0}.png` | Angular deviation δ, all 3 snapshots on one panel |
| `normals_phi_profile_channel_{p3_sigma0,l3_sigma0}.png` | φ profile width over time; ideal tanh reference |
