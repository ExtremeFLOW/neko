#!/bin/bash -l 

#SBATCH -A snic2022-3-25 
#SBATCH -J testrun-rayleigh-probes
#SBATCH -p gpu
#SBATCH -t 00:01:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL

MPI_GPU_SUPPORT_ENABLED=1
srun -u ./neko rayleigh.case > logfile.rayleigh
