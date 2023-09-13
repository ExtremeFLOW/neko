#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A snic2022-3-25


# The name of the script is myjob
#SBATCH -J myjob
#SBATCH --output=myjob.o%j # Name of stdout output file
#SBATCH --error=myjob.e%j  # Name of stderr error file

# The partition
#SBATCH -p gpu

# 1 hour wall-clock time will be given to this job
#SBATCH -t 00:00:15

# Number of nodes
#SBATCH --nodes=1

# Number of MPI processes per node
#SBATCH --ntasks-per-node=1

#SBATCH --mail-type=all         # Send email at begin and end of job

# Load a ROCm module and set the accelerator target
ml PrgEnv-cray
ml craype-accel-amd-gfx90a
ml rocm

cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

CPU_BIND="map_cpu:48,56,16,24,1,8,32,40"

export MPICH_GPU_SUPPORT_ENABLED=1

srun -u --cpu-bind=${CPU_BIND} ./select_gpu ./neko rayleigh.case
rm -rf ./select_gpu
