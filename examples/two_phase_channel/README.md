# Two-Phase Turbulent Channel Flow

Turbulent channel flow at Re_tau=180 with a droplet, using CDI phase-field
interface sharpening and CSF surface tension. Builds on the validated CSF
implementation from the `spurious_currents` example (sign fix, commit 637fc27).

## Physics

- Domain: 4π × 2 × 4/3π, walls at y=±1, periodic in x and z
- Re_b = 2800 (Re_tau ≈ 180), matched fluids (ρ=μ=const)
- Driving: `flow_rate_force` in x, U_b = 1
- Surface tension: We_tau ≈ 2 at σ=4.1e-4, R=0.2, u_tau=0.0643
- Single drop at channel centre (local/validation); two drops at y=±0.5 (production)

## Files

| File | Purpose |
|------|---------|
| `turb_channel_two_phase.f90` | User module |
| `turb_channel_two_phase.case` | Local test case (N=5, R=0.3, end_time=5) |
| `box.nmsh` | 16×18×12 wall-clustered mesh (gitignored; copy from `../turb_channel/`) |

## Local quick-start (macOS)

```bash
cp ../turb_channel/box.nmsh .
source ../../setup-env-channel.sh --local
makeneko turb_channel_two_phase.f90
mpirun -np 4 ./neko turb_channel_two_phase.case
```

## Egidius

```bash
# On egidius (after cloning repo and building Neko):
source ../../setup-env-channel.sh --egidius
makeneko turb_channel_two_phase.f90
mpirun -np <N> ./neko turb_channel_two_phase.case
```

## Dardel (production)

Generate a finer mesh first:
```bash
./contrib/genmeshbox/genmeshbox 0 4 -1 1 0 1.5 24 20 16 .false. .false. .false.
mv box.nmsh examples/two_phase_channel/
```

Then update the case file: `polynomial_order: 7`, `epsilon: 0.10`,
`drop_radius: 0.2`, `end_time: 200`, and submit via `/submit-job`.

## Next steps

- [x] Add `--egidius` option to `setup-env-channel.sh` and `build-neko-channel.sh`
- [x] Clone repo and build Neko on egidius (`/lscratch/sieburgh/local/neko-channel/`)
- [x] `makeneko` compiles cleanly on egidius
- [ ] Run short validation (10–100 steps) on egidius
- [ ] Generate production mesh with `genmeshbox` (24×20×16)
- [ ] Update case file for production (N=7, eps=0.10, two drops, end_time=200)
- [ ] Submit production run to Dardel
