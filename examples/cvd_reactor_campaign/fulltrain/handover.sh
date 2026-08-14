#!/bin/sh
# Handover run-1 -> run-2: wait until run-1 (job 4637305) finishes writing the
# t=0.4 checkpoint, then cancel it and submit fulltrain2.sbatch (14 nodes,
# restart from that checkpoint). Runs on the login node under nohup; all
# progress in handover.log. Safe against partial writes: the checkpoint must
# be >3 GB AND unchanged for 90 s before it is used.
BASE=/proj/liu-compute-2026-4/users/x_adigh/neko_low_mach_solver
OLDJOB=4637305
cd $BASE/full_train_run_2
echo "handover watcher started $(date)" >> handover.log
while true; do
    size=$(stat -c %s $BASE/fulltrain/fluid00002.chkp 2>/dev/null || echo 0)
    if [ "$size" -gt 3000000000 ]; then
        echo "checkpoint appeared ($size bytes) $(date)" >> handover.log
        sleep 90
        size2=$(stat -c %s $BASE/fulltrain/fluid00002.chkp 2>/dev/null || echo 0)
        if [ "$size2" = "$size" ]; then
            echo "checkpoint stable, handing over $(date)" >> handover.log
            cp $BASE/fulltrain/fluid00002.chkp . || { echo "COPY FAILED" >> handover.log; exit 1; }
            scancel $OLDJOB
            sbatch fulltrain2.sbatch >> handover.log 2>&1
            echo "HANDOVER DONE $(date)" >> handover.log
            exit 0
        fi
        echo "still growing ($size -> $size2), waiting" >> handover.log
    else
        state=$(squeue -h -j $OLDJOB -o %T 2>/dev/null)
        if [ -z "$state" ]; then
            echo "ALERT: job $OLDJOB gone without t=0.4 checkpoint $(date)" >> handover.log
            exit 1
        fi
    fi
    sleep 60
done
