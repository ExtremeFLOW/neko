#!/bin/sh
#SBATCH -A naiss2025-22-129
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH -J epsilon_gamma_sweep
#SBATCH -p main

### MODULES & ENV ####
# module unload PDC
# ml PrgEnv-cray/8.5.0
# ml cpe/23.09
# source ./.env
######################

# nproc # = 256



### RUN ###
cd $SLURM_SUBMIT_DIR/zalesak_sweep/gamma_0.01_epsilon_0.02 
srun -u -n 64 ./neko zalesak.case