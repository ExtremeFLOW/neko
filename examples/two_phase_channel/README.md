# Two-Phase Turbulent Channel Flow

Turbulent channel flow at Re_τ=180 with a droplet, using CDI phase-field
interface sharpening and CSF surface tension. Builds on the validated CSF
implementation from the `spurious_currents` example (curvature sign fix,
commit 637fc27).

## Physics

- Domain: 4π × 2 × 4π/3, walls at y=±1, periodic in x and z
- Re_b = 2800 (Re_τ ≈ 180), matched fluids (ρ=μ=const)
- Driving: `flow_rate_force` in x, U_b = 1
- Single drop, R=0.3–0.4, position configurable via `drop_center_y` (default: channel centre y=0)

## Parameter selection

ε and γ are not free constants — they are derived from the mesh and flow speed.

**Interface width ε** must satisfy two constraints:

1. **Resolution:** ε/Δ_GLL ≥ 2–3, where Δ_GLL ≈ Δ_elem/N is the average GLL
   spacing within an element. The interface must span at least 2–3 GLL points
   or the tanh profile is under-resolved and kappa will be inaccurate.

2. **Drop integrity:** R/2ε > 2, so that φ(0) = ½(1+tanh(R/2ε)) > 0.964 at
   the drop centre. Below this the drop "dissolves" in the IC.

For the 81×18×27 mesh at N=7: Δ_GLL_y = (2/18)/7 = 0.016, Δ_GLL_x,z = 0.022.
The binding constraint is the streamwise/spanwise spacing: ε ≥ 3 × 0.022 = 0.066.
Choose ε=0.07 → ε/Δ_GLL_x,z=3.17 ✓, ε/Δ_GLL_y=4.41 (fine).
For R=0.3: R/2ε=2.14 → φ(0)=0.986 ✓. For R=0.4: R/2ε=2.86 → φ(0)=0.997 ✓.

Note: ε=0.05 was used in early runs and gives ε/Δ_GLL_x,z≈2.3, which is
under-resolved in streamwise/spanwise — this is why ε was increased to 0.07.

**CDI mobility γ** controls interface sharpening strength. The criterion is:

    0.01 ≤ γ / u_max ≤ 0.1

where u_max is the peak velocity (~1.5 for Re_τ=180 with U_b=1). The lower bound
ensures CDI resharpens faster than it is advected away; the upper bound keeps the
CDI flux an order of magnitude below the advective flux so it does not distort
the interface. During turbulent events u_max can spike to 2–3, so γ should keep
γ/u_max above 0.01 even at elevated velocities. The scalar diffusivity is
λ = ε·γ; the CDI resharpening timescale is τ_CDI = ε²/γ.

γ=0.05 gives γ/u_max=0.033 at u_max=1.5 and 0.017 at u_max=3 ✓.
τ_CDI = 0.07²/0.05 = 0.098 TU. λ = 0.07×0.05 = 3.5×10⁻³.

**Capillary timestep (We=1 with explicit CSF):** For σ=0.3 and Δ_GLL=0.016,
the capillary stability limit is Δt_cap ≈ sqrt(Δx³/(2πσ)) ≈ 0.0015 TU.
With U_b=1 and target_cfl=0.4 the advective Δt≈0.004 — 3× too large.
Use target_cfl=0.2 for all two-phase cases with σ≥0.03 (We≤10).

**Surface tension σ** is set by the target Weber number:

    We = ρ U_b² R / σ,   U_b = 1 (bulk velocity, reference scale for method validation)

σ = ρ U_b² R / We = 0.3 / We. Primary validation target: We=1 → σ=0.3.

## Key parameters (primary validation: We=1)

| Parameter | Value | How to choose |
|-----------|-------|---------------|
| ε | **0.07** (Legacy) / **0.053** (L1) / **0.040** (L2) / **0.030** (L3) | ε/Δxz=0.457=const; R/2ε > 2 |
| γ | **0.05** | 0.01 ≤ γ/u_max ≤ 0.1; γ/u_max=0.033 at u_max=1.5 |
| target_cfl | **0.2** | Capillary stability: Δt_cap ≈ 0.0015 TU at We=1 |
| σ | 0.3 (We=1) | We = ρ U_b² R/σ; σ = R/We |
| R | 0.3 (P1) / 0.4 (P2/P3) | D=0.6–0.8h; clearance to wall ≥ 0.3h |
| y_c | 0 (centre) | `drop_center_y` in case file; 0=centreline, 0.3=log-law region |
| N | 7 | Polynomial order |
| Mesh | 81×18×27 (Legacy), 108×18×36 (L1), 144×24×48 (L2), 192×32×64 (L3) | 4/3 refinement per level; ε/Δxz=const |

## Files

**Source modules** (`src/`):
| File | Purpose |
|------|---------|
| `src/turb_channel_two_phase.f90` | Two-phase user module (IC, CDI/CSF source terms, diagnostics) |
| `src/turb_channel_two_phase_p2.f90` | Same with extra GS pass on n̂ — deferred, not yet used |
| `src/turb_channel_single_phase.f90` | Fluid-only module for single-phase spin-up |

**Case files** (`cases/<mesh>/`):
| Directory | Mesh | Cases |
|-----------|------|-------|
| `cases/81x18x27/` | 81×18×27 (Δ=0.155) | Baseline Phase 1 runs (single-phase, v4, sigma0, we10, laminar, …) |
| `cases/108x18x36/` | 108×18×36 (Δ=0.116, ε=0.053) | L1 runs (single-phase, sigma0, we10, we1) |
| `cases/144x24x48/` | 144×24×48 (Δ=0.087, ε=0.040) | L2 runs (single-phase, sigma0, we10, we1) |
| `cases/192x32x64/` | 192×32×64 (Δ=0.065, ε=0.030) | L3 runs (single-phase, sigma0, we10, we1) |

**Postprocessing** (`postprocess/`):
| File | Purpose |
|------|---------|
| `postprocess/postprocess_single_phase.py` | ekin plot + mean velocity profile |
| `postprocess/animate_two_phase_channel.py` | Animate φ and \|u\| field snapshots |
| `postprocess/animate_blowup.py` | φ/κ/\|u\| animation for blow-up diagnostics |
| `postprocess/postprocess_sigma0.py` | σ=0 diagnostics: κ_rms/φ time-series + field snapshots. Args: `--run`, `--R`, `--eps`, `--mesh`, `--no-snapshots` |

**Meshes** (gitignored, generated with `genmeshbox`):
- `box_phys_81x18x27.nmsh`, `box_phys_108x18x36.nmsh`, `box_phys_144x24x48.nmsh`, `box_phys_192x32x64.nmsh`

## Running

```bash
cd examples/two_phase_channel
source ../../setup-env-channel.sh --egidius   # or --local / --cluster
genmeshbox 0 12.5664 -1.0 1.0 0 4.1888 81 18 27 .true. .false. .true.
mv box.nmsh box_phys_81x18x27.nmsh
makeneko src/turb_channel_two_phase.f90
mpirun -np 16 ./neko cases/81x18x27/turb_channel_two_phase_v4.case
```

All turbulent two-phase cases require a single-phase spin-up first (produces
`fluid00004.chkp` at t=20). See `CLAUDE.md` for the full workflow. For
production runs on Dardel use the job scripts in `cluster/`.

## Turbulence: initial condition and when the flow is turbulent

### Initial condition

The turbulent IC uses the **Reichardt profile** as the mean velocity plus two deterministic
perturbations that seed the turbulent self-sustaining cycle:

```
u(y) = (Re_tau/Re_b) * [ (1/k) ln(1 + k y+)
         + (C - ln(k)/k)(1 - exp(-y+/11) - (y+/11) exp(-y+/3)) ]
```

with Re_τ=180, k=0.41, C=5.17. This is a smooth approximation to the turbulent mean
profile that satisfies the no-slip boundary condition exactly; the Reichardt formula
matches the viscous sublayer (y⁺ < 5) and the log-law (y⁺ > 30) continuously.

The perturbations have two parts:

1. **Large-scale sinusoidal** (`eps=0.05`, wavenumbers kx=3, kz=4): a divergence-free
   (u,v,w) wave that injects energy at the integral scale (~1/3 domain length). This
   seeds the large-scale vortical structures (streaks and rolls) that drive the
   self-sustaining turbulence cycle.

2. **Small-scale sinusoidal** (`eps=0.005`, kx=17, kz=13): same form at smaller scales,
   filling the energy cascade more quickly.

3. **Random v-perturbation** (`eps=0.001`): a deterministic pseudo-random function of
   (x,y,z) that breaks symmetry and prevents artificial persistence of the initial
   large-scale modes.

The Reichardt profile is NOT Poiseuille flow (which would be parabolic). It is designed
to match the turbulent mean — so the flow begins near its statistically stationary state
from t=0. This is why transition to sustained turbulence is fast (O(10) convective units)
rather than slow (O(100+) as from a laminar IC).

### Laminar vs turbulent: key diagnostic — u_max

The most direct indicator that the flow is turbulent is the centerline (maximum) velocity:

| Flow state | u_max | Why |
|------------|-------|-----|
| Laminar Poiseuille | 1.5 | Parabolic profile: u_max = 1.5 U_b |
| Turbulent at Re_τ=180 | ~1.35–1.45 | Flat log-law mean (~1.15 U_cl/U_b) plus turbulent fluctuations u'≈0.2 |

Note: u_max in Neko's ekin.csv is the instantaneous maximum over the full domain, not
the time-averaged centreline velocity. At Re_τ=180, turbulent fluctuations add ~0.15–0.20
on top of the mean centreline (~1.15), so the instantaneous maximum fluctuates in the
range 1.35–1.45. This is confirmed by the actual single-phase spin-up run (see
`postprocess_single_phase.py`), which also shows the mean velocity profile matches the
Reichardt model throughout the channel.

The flow is statistically turbulent once u_max has dropped from ~1.5 (initial) and is
fluctuating around a flat mean. A secular upward trend toward 1.5 would indicate
relaminarisation (not expected at Re_b=2800).

The instantaneous u_max fluctuates in turbulent flow; what matters is the **time-mean
being ~1.15–1.20 with no secular trend**. If u_max is drifting upward toward 1.5, the
turbulence has relaminarised (does not happen at Re_b=2800).

### Ekin and statistical stationarity

With `flow_rate_force` maintaining U_b=1, the volume-averaged kinetic energy is:

    Ekin = (1/2V) ∫ u_i u_i dV  ≈  0.5 * ( U_b² * shape_factor + u'² )

For turbulent channel flow at Re_τ=180, Ekin ≈ 0.5–0.6 (depends on the velocity profile
shape factor and turbulent fluctuation intensity). The flow is statistically stationary
when Ekin fluctuates around a flat mean with no trend. Typical time to reach stationarity
from a Reichardt IC: **t ≈ 10–15 convective time units** (t × U_b/h, U_b=h=1 here).

### Checkpoint selection for two-phase restart

**All turbulent two-phase cases restart from `fluid00004.chkp`.** This is the standard
workflow: run the single-phase spin-up first, then inject the drop analytically via the
scalar IC. The laminar case is the only exception (Poiseuille IC, no checkpoint).

The single-phase case (`turb_channel_single_phase.case`) runs to t=25 with checkpoints
at t=5, 10, 15, 20, 25. The recommended restart checkpoint is **`fluid00004.chkp` (t=20)**:
- t=5–10: turbulence establishing, statistics not fully stationary
- t=15–20: statistically stationary; use either checkpoint
- t=25: also valid, but t=20 is sufficient

To verify before restarting: inspect `ekin.csv` from the single-phase run. If u_max
is fluctuating around 1.15–1.20 with no trend for at least 5–10 time units, the flow
is ready.

## Results

See [RESULTS.md](RESULTS.md).
