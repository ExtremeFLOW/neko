# heated_hump_full — full-physics production configuration

Everything this campaign validated, in one case:

| ingredient | setting | pedigree |
|---|---|---|
| mesh | `hump_mid.nmsh`, 138,240 elem (wake-refined) | built 2026-07-22, zones verified |
| order | 6 (47M GLL pts) | ILES + GJP posture |
| Re_in | 2000 | tripped -> sustained turbulence expected |
| trip | random volume-force trip line at x=2 | `trip.f90` |
| heating | T_hot = 4 (4x density ratio), wall x > 11 | amplitude scaling of validated physics |
| viscosity | power_law mu(T), mu_ref = 1/2000, exponent 0.7, Pr = 0.71 | channel muT soak (t=200 flawless) |
| gravity | do_gravity = .true., g_nd = 0.15 reduced buoyancy | channel soaks |
| scalar stabilization | GJP (0.8/4.0) on temperature | Sigma GJP restart (violations -5x..-100x) |
| outflow | Dong (energy-stable) | backflow fix, validated t=70->100 |
| solver | heat-first + reorth projection + h1=1 + coarse cg-100 | the whole week |

## Build (Sigma) — NOTE: TWO source files, trip module FIRST
```bash
module load buildenv-gcc/2022a-eb
export LD_LIBRARY_PATH=$HOME/json-fortran-install/jsonfortran-gnu-9.2.0/lib:$LD_LIBRARY_PATH
~/neko-lowmach/install/bin/makeneko trip.f90 heated_hump_full.f90
```
(makeneko requires each custom module in its own file, compiled before its user.)

## Run
```bash
sbatch heated_hump_full.sbatch      # 10 nodes/320 ranks: ~150k pts/rank, ~2-4 s/step
                                    # 20 nodes recommended if queue allows
```
end_time 50, fields every 1 tu, checkpoints every 10.

## Recommended: 1-2h shakedown FIRST (validated physics on the new mesh)
```bash
sed -e 's/T_hot  = 4.0_rp/T_hot  = 2.0_rp/' -e 's/do_gravity = .true./do_gravity = .false./' \
    heated_hump_full.f90 > shakedown.f90
# + remove the property_model block from a copy of the case, makeneko, run ~5-10 tu
```
If the shakedown holds its vitals (prs iters, dt, trip signature at x=2), launch the full run.

## What to watch
- startup: `Scalar` before `Fluid`, `Property : power_law`, order 6, user ICs
- prs iterations: expect settled ~20-60 on this mesh with cg-100 (tamg block is the fallback)
- T in [1,4] with GJP (clamp floor 0.8 should stay idle)
- t ~ 8-15: tripped transition; t > 15: turbulent heated buoyant wake
