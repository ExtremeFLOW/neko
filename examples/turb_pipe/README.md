# Pipeflow at Re_b=5300 and Re_tau = 180.

## Mesh generation

There is a pre-generated nmsh file that can be used. To regenerate the mesh with
a different resolution or domain size, one needs gmsh and gmsh2nek, the latter
available as a tool in Nek5000. To modify the mesh parameters, you should edit
the .geo file, which serves as input to gmsh.

Once you are happy generate the gmsh mesh with
```
gmsh -3 ./turb_pipe.geo -order 2
```

The next stage is to convert the mesh to the .re2 format of Nek5000. Run
gmsh2nek and the following input when prompted:

```
Enter mesh dimension: 3
Input fluid .msh file name: turb_pipe
Do you have solid mesh ? (0 for no, 1 for yes) 0
 Enter number of periodic boundary surface pairs:
1
 input surface 1 and  surface 2  Boundary ID
1 2
 please give re2 file name:
turb_pipe

```
This will generate a turb_pipe.re2 file, which can be converted to Neko's format
with

```
rea2nbin turb_pipe.re2
```

## Running

To perform a DNS of the case set the polynomial order to 7 in the case file..
should be used. We have two case files, one called pipe.case using
"flow_rate_force" to drive the flow,  and one called pipe_source.case that
instead uses a constant source term. The source term values is based on
computing the shear stress as 2 * (2 * Re_tau/Re_B)^2, which assumes that the
bulk velocity is 1.

We use quite forgiving tolerances, and this case should be possible to run
efficiently on consumer GPUs in single precision, specify '--enable-real=sp'
when configuring.

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
