# Shipping the heated_hump low-Mach case to HPC

Everything you need to build the (CPU) low-Mach Neko fork and run the
flow-over-a-hump case on a cluster.

## 0. Key facts before you start
- **CPU-only.** The low-Mach scheme errors out on GPU
  (`fluid_lowmach_pnpn.f90`: "device backend not yet implemented"). Build and run
  on **CPU compute nodes**, many MPI ranks. Do **not** configure `--with-cuda/hip`.
- **Dependencies:** an MPI compiler (`mpif90`/`mpicc`), BLAS/LAPACK (OpenBLAS or MKL),
  and **json-fortran**. `gslib` is **not** required (the recycling inflow uses Neko's
  native `findpts`).
- The fixes this case depends on (Q_T conduction, energy-source weighting, etc.) are
  **committed on the `feature/lowmach` branch** — ship the repo at that commit.

## 1. What to copy to the cluster
Two things:
1. The **Neko fork** `~/neko-lowmach` (the whole repo, `feature/lowmach` branch).
   Either `git push` to a remote and `git clone` on HPC, or tar it:
   ```bash
   tar czf neko-lowmach.tar.gz -C ~ neko-lowmach \
       --exclude='neko-lowmach/install' --exclude='*.o' --exclude='*.mod' \
       --exclude='examples/*/field0.*' --exclude='examples/*/*.chkp'
   ```
2. The **case files** (already inside `examples/heated_hump/`):
   `heated_hump.f90`, `heated_hump.case`, `hump.nmsh` (the 53k-element mesh, ~13 MB),
   `hump.geo` (to regenerate the mesh if needed), and these scripts.

`hump.nmsh` is self-contained — you do **not** need gmsh/gmsh2nek on the cluster
unless you want to change the mesh (then copy `hump.geo` and re-run the steps in it).

## 2. Build (see `build_neko_hpc.sh`)
On a login/build node, after loading your MPI + BLAS + (CMake for json-fortran) modules:
```bash
bash build_neko_hpc.sh            # builds json-fortran + Neko into ~/neko-lowmach/install
```
Then build the case driver:
```bash
export LD_LIBRARY_PATH=$JSONFORTRAN/lib:$LD_LIBRARY_PATH
~/neko-lowmach/install/bin/makeneko heated_hump.f90    # -> ./neko
```

## 3. Resolution & resources
- Full mesh: **53,760 elements**. At `polynomial_order = 7` that's **27.5 M grid points**
  (the `.case` is already set to order 7). Order 5 = 11.6 M if you want a cheaper run.
- This is a real LES — budget **~256-1024 CPU cores** (a handful of nodes). Rule of
  thumb ~20-50k points/core; tune for your machine. Total RAM should be comfortably
  > ~100 GB across the job (the pressure GMRES Krylov basis is the main consumer).
- `mpirun -n <ranks>` — Neko partitions the mesh automatically (ParMETIS).

## 4. Run (see `job.slurm`)
Edit `job.slurm` for your scheduler/account/partition, then:
```bash
sbatch job.slurm
```

## 5. Tuning the physics for a real turbulent run
All knobs are at the top of `heated_hump.f90` (rebuild with `makeneko` after editing):
- `Re_in` — currently **500** (stable on the coarse local mesh). On the full mesh you
  can raise it (e.g. **5000-10000**) for genuine turbulence; the fine grid resolves it.
- `g_nd` (gravity), `trip_A` (trip), `x_heat=11`, `T_hot=2` — as requested.
- GJP stays on in the `.case` (`gradient_jump_penalty`, factor 0.8, exp 4) for LES stability.
- Start from a checkpoint to restart long runs (`output_checkpoints: true` is set;
  add `"restart_file": "fluidXXXXX.chkp"` under `case` to continue).

## 6. Sanity check on the cluster first
Before the big run, do the 5-step smoke test (cheap) to confirm the build + mesh load:
```bash
mpirun -n 8 ./neko smoke.case      # uniform inflow, ~10 steps, should converge cleanly
```
Then launch the full job.
