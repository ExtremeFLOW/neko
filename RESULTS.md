# Simulation Results

Convergence series: isotropic meshes L1/L2/L3 with constant ε/Δxz≈0.457.
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

## L2 (144×24×48) — two-phase cases (pending submission 2026-03-29)

ε=0.04, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25).
256 ranks, 2 nodes. Build: `neko_two_phase` (job 19050141 running 2026-03-29).

| Case | Run name | Job | Status | Notes |
|------|----------|-----|--------|-------|
| σ=0 standard | `channel_p3_sigma0` | 19050210 | RUNNING | Standard `turb_channel_two_phase.f90` |
| σ=0 _p2 | `channel_p3_sigma0_p2` | 19050211 | CANCELLED | Extra GS on n̂ is a no-op — see note below |
| We=10 | `channel_p3_we10` | TBD | PENDING | |
| We=1 | `channel_p3_we1` | TBD | PENDING | |

_p2 build: job 19050141 (`build_neko_two_phase_p2.sh`), submitted 2026-03-29.

---

## L3 (192×32×64) — two-phase cases (pending submission)

ε=0.03, R=0.4, γ=0.05. Restart from `fluid00004.chkp` (t=20→25).
512 ranks, 4 nodes. Requires `neko_two_phase` binary (same as L2).

| Case | Run name | Job | Status | Notes |
|------|----------|-----|--------|-------|
| σ=0 | `channel_l3_sigma0` | TBD | PENDING | |
| We=10 | `channel_l3_we10` | TBD | PENDING | |
| We=1 | `channel_l3_we1` | TBD | PENDING | |

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

| Mesh | ε | κ_rms (t=25) | κ_rms / (2/R) | Notes |
|------|---|-------------|---------------|-------|
| L1 (ε=0.09) | 0.09 | ~50 | ~10× | Large ε, informative |
| L1 (ε=0.053) | 0.053 | ~56 | ~11× | Convergence series point |
| L2 (ε=0.04) | 0.04 | TBD | TBD | Standard vs _p2 comparison |
| L3 (ε=0.03) | 0.03 | TBD | TBD | |
