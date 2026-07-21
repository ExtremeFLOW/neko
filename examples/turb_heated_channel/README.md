# Turbulent heated channel (low-Mach, Re_τ ≈ 180, 4× density ratio)

A clean, wall-resolved **turbulent** channel with **strong wall heating**, solved with
the low-Mach (variable-density) Neko scheme. It is the `turb_channel` low-Mach case
turned up from a 2:1 to a **4:1 temperature ratio**, on a mesh refined to resolve the
near-wall gradients.

## Physics

| Quantity | Value |
|---|---|
| Domain (δ = half-height = 1) | x∈[0,2π] streamwise (periodic), y∈[0,2] wall-normal (no-slip), z∈[0,π] spanwise (periodic) |
| Bulk Reynolds number | Re_b = 2800 (→ Re_τ ≈ 180), bulk velocity U_b = 1 |
| Molecular viscosity | μ = 1/Re_b = 3.571×10⁻⁴ (constant) |
| Prandtl number | Pr = 0.71 → λ = μ/Pr |
| **Heating** | hot wall y=0 (zone 3) **T = 4**, cold wall y=2 (zone 4) **T = 1** |
| EOS | ideal gas ρ = P0/(R·T) = 1/T → **ρ ∈ [0.25, 1.0], a 4× density ratio** |
| Streamwise drive | constant volumetric flow rate (`flow_rate_force`), U_b = 1 |
| Buoyancy | none (clean forced variable-density test) |

The density varies 4× through the EOS; that is the dominant low-Mach effect. Because
the walls heat/cool the fluid, ρ, and the divergence constraint ∇·u = Q_T, vary
strongly near the walls — which is what the mesh below is built to resolve.

## Mesh (`make_mesh.sh` + `ydist.csv`)

`genmeshbox` makes a uniform box, but accepts a per-direction vertex-distribution
file. We use **16 × 16 × 12 elements at polynomial order 7**, with the wall-normal
direction **clustered toward both walls** (tanh, γ=2; the 17 vertices are in
`ydist.csv`). x and z stay uniform.

Resulting resolution at Re_τ = 180 (worst/cold wall):

| | value | DNS target |
|---|---|---|
| first off-wall point | **y⁺₁ ≈ 0.27** (0.14 at the hot wall) | < 1 |
| streamwise GLL spacing | **Δx⁺ ≈ 10** | 8–15 |
| spanwise GLL spacing | **Δz⁺ ≈ 6.7** | 4–10 |

The hot-wall conduction sublayer (y⁺_T ≈ Pr^(−1/3) ≈ 1.1) is resolved with margin
(y⁺₁ = 0.14 there).

## Build & run

```bash
cd examples/turb_heated_channel
./make_mesh.sh                                   # -> box.nmsh
LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib \
    /home/hochi/neko-lowmach/install/bin/makeneko turb_heated_channel.f90   # -> ./neko
LD_LIBRARY_PATH=/home/hochi/json-fortran-install/lib \
    ./neko turb_heated_channel.case
```

## Output

- `field0.nek5000` — **u, v, w, p, T** (correctly labelled), every 5 time units.
- `props0.nek5000` — **ρ, μ, k, T**. ⚠️ Neko's `field_writer` packs a field *list*
  into fld slots, so the ParaView labels are scrambled: **`pressure` = ρ,
  `velocity-X` = μ, `velocity-Y` = k, `velocity-Z` = T**. The data is correct.

## What to expect

- The Reichardt-profile IC + perturbations trip transition; bulk velocity stays at
  ~1 (forced), centerline/bulk ≈ 1.16 once turbulent (≈1.5 if it relaminarises).
- The two walls have **different friction velocities** (the lighter hot wall carries
  less drag), so "Re_τ = 180" is nominal (cold-wall anchored); the mean profile is
  asymmetric.
- **Cost / known bottleneck:** the variable-density pressure Poisson is solved by a
  defect-correction loop whose *inner* gmres is deliberately **inexact**
  (`absolute_tolerance 1e-3`, `max_iterations 150`) — a tight tolerance hangs step 1
  at this density ratio. The pressure solve dominates wall-time; full turbulence
  statistics (`end_time` ~200, tens of thousands of steps) is expensive. For a quick
  check, lower `end_time` to ~6 (≈1 flow-through, t_ft = 2π) and confirm the flow
  stays bounded and the perturbations grow.

## Notes / extensions

- `case.scalar.low_mach.enabled = true` activates the variable ρ·cp(x) energy
  coupling — **required** for a consistent heated low-Mach run.
- To make μ(T)/k(T) temperature-dependent (e.g. μ = μ_ref·T^0.7, k = cp·μ/Pr), add a
  `material_properties` user hook — the per-step `mu_tot`/`lambda_tot` refresh already
  propagates it to the solver and to Q_T.
- To add buoyancy (hot-bottom/cold-top → mixed convection), add a `source_term` hook
  and a gravity term; note this introduces unstable stratification.
