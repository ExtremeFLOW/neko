# Cylinder at $Re_D=200$, N=5, 1824 elements
Enabled simcomps:
- fluid_statistics, basic, 2D
- lift/drag calculation
- lambda 2

Included files:
- Gmsh script to generate mesh
- neko .nmsh file (converted from gmsh-re2-nmsh)
- case file

## Insitu data processing

For this example, we have also added some files, named `insitu_*` that can be used
for insitu data processing, i.e., treat the data while the simulation runs.

In this particular case we perform proper orthogonal decomposition (POD). For this purpose,
the user would first need to run a simulation with the standard provided files until the steady
state is reached. Thereafter, the `insitu` files can be executed.

The files related to insitu execution are:
- `insitu_turb_pipe.f90` is a new user file that includes the data streaming in the user module.
- `turb_pipe.case` includes a parameter `istream` that controls the data-stream frequency.
- `insitu_task.py` is the python script that performs the insitu POD.

There are some dependecies for this to work:
- Neko must be installed with `Adios2` support such that data can be transferred to python.
- PySEMTools must be installed, as the python script relies heavily on it.

### Executing neko in mpmd
The insitu run can be performed by executing neko and python in multiple program multiple data (mpmd) mode.

On personal computers with mpi, one can execute, for example:
```bash
mpirun -n 1 ./neko insitu_turb_pipe.case : -n 4 python3 insitu_task.py
```
provided that the visible devices are set up correctly.

On supercomputers, for the particular case of slurm, one can run in mpmd mode by executing a comand such as:
```bash
srun --multi-prog mpmd.conf
```
`mpmd.conf`, in its simplest form, is a file that indicates what application each rank executes. It can be
generated, for example, as:

```bash
cat > mpmd.conf << 'EOF'
0 ./select_gpu ./neko turb_pipe.case
1 python3 insitu_task.py
2 python3 insitu_task.py
3 python3 insitu_task.py
4 python3 insitu_task.py
EOF
```

This example has been succesfully executed on supercomputers such as lumi, leveraging the typically idling
CPUs on a GPU compute node, to execute POD while the simulation runs on GPUs.
