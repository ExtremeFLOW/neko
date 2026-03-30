# Simulation Results

Convergence series: isotropic meshes L1/L2/L3/L4 with constant ε/Δxz≈0.457.
All runs on Dardel (NAISS account `naiss2025-3-39`), results synced to egidius at
`/lscratch/sieburgh/simulations/<run_name>/`.

---

## Single-phase spin-ups

All spin-ups run t=0→25, checkpoint `fluid00004.chkp` written at t=20.
Validation: u_max (last 5 rows of ekin.csv) < 1.45 → TURBULENT.

| Mesh | Run name | Job | Nodes/Ranks | Status | u_max (last 5) | Notes |
|------|----------|-----|-------------|--------|----------------|-------|
| 81×18×27 (baseline) | `channel_single_phase` | — | 16 ranks | COMPLETED | ~1.38 | egidius local run |
| L1: 108×18×36 | `channel_p2_single_phase` | 18985538 | 1/128 | COMPLETED | ~1.38 | |
| L2: 144×24×48 | `channel_p3_single_phase` | 19002586 | 2/256 | COMPLETED | 1.37±0.01 | 2 nodes (1 node OOM) |
| L3: 192×32×64 | `channel_l3_single_phase` | 19003203 | 4/512 | COMPLETED | 1.34±0.00 | 4 nodes (2 nodes OOM) |
| L4: 288×48×96 | `channel_l4_single_phase` | — | ~14/1792 | NOT STARTED | — | planned; ~14 nodes est. |

---

## L1 (108×18×36) — two-phase cases, jobs submitted 2026-03-28

ε=0.053 (convergence series), R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25).
128 ranks, 1 node.

| Case | Run name | Job | Status | κ_rms at t=25 | Notes |
|------|----------|-----|--------|---------------|-------|
| σ=0 ε=0.053 | `channel_p2_sigma0_eps053` | 19003682 | COMPLETED t=20→25 | ~56 | Convergence series point |
| σ=0 ε=0.09 | `channel_p2_sigma0` | 19003683 | COMPLETED t=20→25 | ~50 | Informative (large ε) |
| We=10 | `channel_p2_we10` | 19003684 | synced | — | |
| We=1 | `channel_p2_we1` | 19003685 | synced | — | |

---

## L2 (144×24×48) — two-phase cases

ε=0.04, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25).
256 ranks, 2 nodes.

| Case | Run name | Job | Status | κ_rms at t=25 | Notes |
|------|----------|-----|--------|---------------|-------|
| σ=0 | `channel_p3_sigma0` | 19050210 | COMPLETED t=20→25 (3h35) | ~55.9 | Convergence series point |
| σ=0 _p2 | `channel_p3_sigma0_p2` | 19050211 | CANCELLED | — | Extra GS on n̂ is a no-op |
| We=10 | `channel_p3_we10` | TBD | PENDING | — | |
| We=1 | `channel_p3_we1` | TBD | PENDING | — | |

---

## L3 (192×32×64) — two-phase cases

ε=0.03, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25).
512 ranks, 4 nodes.

| Case | Run name | Job | Status | κ_rms at t=25 | Notes |
|------|----------|-----|--------|---------------|-------|
| σ=0 | `channel_l3_sigma0` | 19050366 | COMPLETED t=20→25 (5h51) | ~67 | Convergence series point |
| We=10 | `channel_l3_we10` | TBD | PENDING | — | |
| We=1 | `channel_l3_we1` | TBD | PENDING | — | |

---

## L4 (288×48×96) — planned, not yet started

ε=0.02, R=0.4, γ=0.05. 1,327,104 elements, ε/Δxz=0.458, ε/Δ_GLL=3.2, R/ε=20.
Factor 3/2 from L3 (the 4/3 progression gives nz=85.3, non-integer, so 3/2 used instead).
Estimated ~14 nodes / 1792 ranks. Intended for high-resolution visualisation and as the
fourth convergence series point. May need to go higher still for publication quality.

Workflow (when ready):
1. Generate mesh on egidius: `genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 288 48 96 .true. .false. .true.`
2. Write spin-up case file (`turb_channel_single_phase_l4.case`) and job script
3. Run spin-up (~14 nodes, est. wall time TBD)
4. Submit σ=0 and We cases as for L1/L2/L3

| Case | Run name | Job | Status | Notes |
|------|----------|-----|--------|-------|
| Single-phase spin-up | `channel_l4_single_phase` | — | NOT STARTED | ~14 nodes / 1792 ranks |
| σ=0 | `channel_l4_sigma0` | — | NOT STARTED | |
| We=10 | `channel_l4_we10` | — | NOT STARTED | |
| We=1 | `channel_l4_we1` | — | NOT STARTED | |

---

## σ=0 CDI diagnostic — interpretation and kink diagnosis

σ=0 cases use CDI only, no surface tension. κ_rms measures the curvature residual
from the SEM discretisation. Expected baseline: 2/R = 5.

**Kink root cause (2026-03-29):** The artifact is INTRA-ELEMENT, not inter-element.
After GS(∇φ)+normalize, n̂ at shared face nodes is already identical between elements.
An extra GS on n̂ (the _p2 "fix") is therefore a no-op: it sums identical values
and divides back. The actual kink is in the Lagrange polynomial connecting the
averaged face-node n̂ to the element-local first-interior n̂; D[N,N]≈14 (N=7)
amplifies this jump. GS cannot reach interior nodes.

A true fix requires: (a) repeated GS on ∇φ before normalisation (propagates
average inward one node-ring per pass), (b) a Helmholtz-type smoother on n̂,
or (c) a fundamentally different curvature scheme.

**Diagnostics to add in postprocessing:** analyse n̂ directly from φ snapshots
(recomputable via SEM derivative operators), visualise face-node vs first-interior
jump in n̂ and the raw div(n̂) field.

| Mesh | ε | κ_rms plateau | κ_rms / (2/R) | Notes |
|------|---|---------------|---------------|-------|
| L1 (ε=0.09) | 0.09 | ~50 | ~10× | Large ε, informative |
| L1 (ε=0.053) | 0.053 | ~56 | ~11.2× | Plateau reached t≈22 |
| L2 (ε=0.04) | 0.04 | ~54–56 | ~10.8–11.2× | Plateau reached t≈23; slight late-time creep |
| L3 (ε=0.03) | 0.03 | ~67–68 | ~13.4–13.6× | Plateau reached t≈24; late-time creep similar to L2 |

**Key finding (2026-03-30):** κ_rms does NOT decrease with mesh refinement at constant ε/Δxz.
L3 (finer) converges to a HIGHER plateau than L2. Both L2 and L3 show a slow upward creep
after the plateau — consistent with turbulence gradually deforming the interface and increasing
kink amplitude in dynamic equilibrium.

Interpretation: the kink curvature scales as ~1/ε (curvature of the intra-element kink in n̂
is proportional to 1/Δx ∝ 1/ε at constant ε/Δx). CDI convergence in the literature refers
to the INTERFACE SHAPE (φ accuracy), not the COMPUTED curvature via div(n̂). The SEM
derivative amplification D[N,N]≈14 means the kink in n̂ dominates κ_rms regardless of mesh
refinement. A different curvature computation scheme (height-function, parabolic fit, or
pre-smoothed n̂) would be needed to achieve κ_rms → 2/R with mesh refinement.
