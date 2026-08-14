# Setup notes — heated low-Mach cases in Neko

Lessons from getting `hotzone_trip` / `hotzone_species` to run cleanly
(2026-08-12/13). Read this before setting up another heated case; most of it
cost hours to learn and is invisible until it bites.

---

## 1. Never impose a hot wall as a hard step — in space OR in time

In a low-Mach solver the wall temperature is a **forcing term for the pressure
field**, because the velocity divergence is set by the heating:

```
div(u)  ∝  DT/Dt
```

So a discontinuity in the wall condition is handed to the solver as an impulse.
There are two independent discontinuities and each needs its own ramp.

### 1a. Spatial — `L_RAMP` (a DISTANCE, metres)

Where an adiabatic wall meets a heated wall, the BC jumps from `neumann 0` to
`dirichlet T_hot` in zero distance. That is a genuine singularity of the
continuous problem: heat flux → ∞, thermal BL thickness → 0. **No mesh resolves
it**, so a spectral element sitting on it rings forever.

Measured (run 4637223, hard step, T_hot = 3):

| t | min T | points T<0.9 |
|---|---|---|
| 0.05 | 0.944 | 0 |
| 0.25 | 0.355 | 2,107 |
| 0.50 | **−0.508** | 5,459 |

Monotonic, *growing*, permanently localised at z = 0.251…0.278 within
0.6–1.2 mm of the floor — i.e. glued to the leading edge. It never heals,
because the singularity keeps manufacturing it.

**Sizing:** span 3–4 elements. Here elements near the edge are ~6.8 mm, so
`L_RAMP = 0.025` m ≈ 3.7 elements ≈ 26 GLL points at N=6.

### 1b. Temporal — `T_RAMP` (a DURATION, sim-time units)

Gas starts at T_cold; setting walls to T_hot at t=0 asks the near-wall gas to
jump instantly. `div(u) ∝ DT/Dt` then spikes, velocities blow up, and the CFL
controller collapses dt.

Measured (run 4637228, spatial ramp only, no temporal ramp):

| | with temporal ramp | without |
|---|---|---|
| dt | 1.29e-5 → 1.86e-5 (grows) | 1.29e-5 → **4.86e-7** (26× collapse) |
| pressure iters | 12–17 | 60–90 |
| step cost | 0.45 s | 1.58 s |
| useful output | ~1.8 flow-through/job | ~0.05 |

**≈30× difference in output per job.** The pressure solver "struggling" was a
*symptom*; tuning its tolerance would have treated the symptom and left the dt
penalty in place.

---

## 2. How to size `T_RAMP` — and the mistake I made

The near-wall thermal film thickens as `δ ≈ sqrt(α t)`, `α = λ/(ρ cp)`.
A polynomial cannot represent a feature thinner than its own element, so while
δ < element height you get ringing no matter what.

For this case (α = 2.20e-5, first element 0.824 mm):

| t | δ | elements |
|---|---|---|
| 0.01 | 0.47 mm | 0.6 |
| 0.03 | 0.81 mm | **1.0** |
| 0.12 | 1.63 mm | 2.0 |
| 0.28 | 2.48 mm | 3.0 |

**I set `T_RAMP = 0.03`, matched to δ = ONE element. That was ~4x too short.**
One element is not "resolved" — you want the film spanning **2–3 elements**
before the wall reaches full temperature, i.e. `T_RAMP ≈ 0.12–0.28` here.

Consequence: a transient undershoot to min T = 0.28 formed around t ≈ 0.05,
convected downstream and flushed out by t = 0.41. Bounded and self-healing
(unlike 1a), but avoidable.

**Rule for next time:**

```
T_RAMP  ≈  (2..3 × first_element_height)^2 / alpha        alpha = lambda/(rho*cp)
L_RAMP  ≈  3..4 × local streamwise element size
```

### Telling the two failure modes apart

| symptom | cause |
|---|---|
| undershoot **grows**, pinned at a fixed z, hugging the wall | spatial step (1a) |
| undershoot **convects downstream and decays**, spread across the duct interior, symmetric floor/ceiling | startup, film thinner than mesh (1b/2) |
| **dt collapses** in the first ~25 steps, pressure iters 3–5× normal | missing temporal ramp |

---

## 3. Use smoothstep, not linear

```fortran
s = min(max(arg, 0.0_rp), 1.0_rp)
s = s*s*(3.0_rp - 2.0_rp*s)          ! C1: S'=0 at both ends
```

A linear ramp has dT/dt **jumping** at both ends, and since `div(u) ∝ DT/Dt`
each jump is a step in the pressure source. Smoothstep's derivative vanishes at
both ends. Cost: peak slope is 1.5× the linear rate, which is fine — a bounded
peak beats an instantaneous one. If still too harsh, quintic smootherstep
`s³(s(6s−15)+10)` is C² (peak 1.875×).

**Fortran gotcha:** `pure function Twall(z, t) result(T)` does not compile —
Fortran is case-insensitive, so dummy `t` collides with result `T`. Name the
time argument `tnow`.

---

## 4. Things that must move together (silent-failure list)

- **Both heated zones must ramp.** Ramping `wallHot` but leaving `substrate` a
  hard dirichlet just moves the discontinuity to their junction. Mark them
  together: `"zone_indices": [4, 5]`. One `Twall(z,t)` serves both, because the
  substrate lies below `Z_HEAT − L_RAMP` so its spatial factor is already 1.
- **`T_eos_max` must exceed `T_hot`.** Raising T_hot 3→6 with `T_eos_max` left
  at 5.0 silently clamps the density EOS near every hot wall. No error, just
  wrong density.
- **Changing Re moves four numbers**, not one: `mu_ref`, the temperature
  `lambda_ref` (= mu_ref/Pr), every species `lambda` (= mu_ref/Sc), and the
  derived `g_nd` in the user file.
- **Checkpoint/output cadence must match the dt you actually get**, not the one
  you hoped for. Run 4637228 would have run 16 h and written **no checkpoint**
  because `checkpoint_value` was sized for a dt 10× larger than reality.
  Size it after seeing dt settle, or set it small.

See also `../neko_bug/README.md` for two upstream Neko issues found along the
way (a crash in the scalar BC error path, and `case.scalar` silently shadowing
`case.scalars`).
