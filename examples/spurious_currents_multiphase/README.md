# Spurious Currents Multiphase Test Case

Stationary circular drop (D=0.4) in a 1x1 box with slip boundaries.
Equal densities (rho=300) and viscosities (mu=0.1), surface tension sigma
variable (swept over 0.05-1.0), giving Laplace number La = sigma*rho*D/mu^2.

Non-dimensional time: t* = sigma*t/(mu*D). Each sweep run targets t*=250.

## Quick start

Generate the mesh:

```bash
genmeshbox 0 1 0 1 -0.1 0 10 10 1 .false. .false. .true.
```

Single run:

```bash
makeneko spurious_currents.f90
mpirun -np 4 ./neko spurious_currents.case
```

Sigma sweep (on Dardel):

```bash
./run_sigma_sweep.sh --cluster
```

## Files

| File | Purpose |
|------|---------|
| `spurious_currents.f90` | User module with CSF surface tension and diagnostics |
| `spurious_currents.case` | Case file (sigma, end_time patched by sweep script) |
| `run_sigma_sweep.sh` | Sweep over sigma values (local or Dardel) |
| `dardel_job.sh` | SLURM template for Dardel |
| `analyze_spurious_currents.ipynb` | Main analysis notebook |
| `sign_convention_diagram.png` | Diagram of the curvature sign bug and fix |
| `FINDINGS.md` | Documentation of the sign bug investigation |

## Bug fix history

See [FINDINGS.md](FINDINGS.md) for details on the curvature sign bug
(kappa = div(n) -> kappa = -div(n)) fixed in commit `637fc27`.
