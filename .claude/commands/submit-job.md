Compile a Neko user module and submit a simulation job to Dardel.

Arguments: example directory name (under `examples/`); optionally nodes and walltime.
Example: `/submit-job spurious_currents_multiphase`
Example: `/submit-job oscillating_droplet 2 04:00:00`

For generic (non-Neko) Dardel job submission, use `/dardel-submit` instead.

## Steps

1. Parse arguments:
   - `CASE` = name of the directory under `examples/` on the cluster
   - `NODES` = integer, default 1 (1 node = 128 cores)
   - `WALLTIME` = HH:MM:SS, default 02:00:00

2. Sync local changes to the cluster via the transfer node:
   ```bash
   ./cluster/sync.sh up examples/<CASE>
   ```

3. SSH to `dardel` (login node) and:
   ```bash
   source ~/.bashrc
   cd $KTHMECH_PROJECT/software/neko-src/examples/<CASE>
   ```

4. Compile the user module if a matching `.f90` exists:
   ```bash
   makeneko <CASE>.f90
   ```

5. Verify the `.case` JSON file has `output_directory` pointing to scratch:
   `/cfs/klemming/scratch/e/eriksie/<run_identifier>`
   Update it if needed.

6. Check or create `job.sh` using the Neko template:
   ```bash
   #!/bin/bash
   #SBATCH -A naiss2025-3-39
   #SBATCH -t <WALLTIME>
   #SBATCH -N <NODES>
   #SBATCH -J <CASE>
   #SBATCH -p main

   cd $KTHMECH_PROJECT/software/neko-src/examples/<CASE>
   srun -u -n $((NODES * 128)) ./neko <CASE>.case
   ```

7. Submit: `sbatch job.sh` — report the job ID.

8. Tell the user to monitor with: `squeue -u eriksie`

## Paths on cluster
- Neko source + examples: `$KTHMECH_PROJECT/software/neko-src/` = `/cfs/klemming/projects/supr/kthmech/eriksie/software/neko-src/`
- Output should always go to: `/cfs/klemming/scratch/e/eriksie/<run_id>/`
