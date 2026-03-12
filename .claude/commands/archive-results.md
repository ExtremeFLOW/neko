Archive Neko simulation results from Dardel scratch to project storage and optionally transfer to local machine.

Arguments: scratch subdirectory name; optionally a local destination path.
Example: `/archive-results spurious_currents_sigma_sweep`
Example: `/archive-results oscillating_droplet examples/oscillating_droplet/results`

For archiving non-Neko results, use `/dardel-archive` instead.

## Steps

1. SSH to `dardel` and check the scratch directory:
   ```bash
   du -sh /cfs/klemming/scratch/e/eriksie/<RUN_ID>
   ls /cfs/klemming/scratch/e/eriksie/<RUN_ID>/
   ```
   Report size and contents to the user.

2. Copy to project results storage:
   ```bash
   mkdir -p $KTHMECH_PROJECT/results/<RUN_ID>
   cp -r /cfs/klemming/scratch/e/eriksie/<RUN_ID>/* $KTHMECH_PROJECT/results/<RUN_ID>/
   ```
   (`$KTHMECH_PROJECT` = `/cfs/klemming/projects/supr/kthmech/eriksie`)

3. Verify file counts match. Report success.

4. Transfer to local machine via `dardel-ftn`:
   - If the user gave a local destination, rsync there:
     ```bash
     rsync -avz --progress dardel-ftn:$KTHMECH_PROJECT/results/<RUN_ID>/ <LOCAL_DEST>/
     ```
   - Otherwise ask whether they want a local copy and where.

5. Ask before deleting scratch. Only delete with explicit confirmation:
   ```bash
   rm -rf /cfs/klemming/scratch/e/eriksie/<RUN_ID>
   ```

## What data matters for post-processing
For spurious currents simulations, the key files are:
- `ekin.csv` — kinetic energy and diagnostics time series (main analysis input)
- `field0.f0*` — Nek5000-format field snapshots (for visualisation, large)
- `*.log` — job log with solver output
- `*.case` — case configuration used for this run

For oscillating droplet: all of the above plus the droplet shape history if written.
