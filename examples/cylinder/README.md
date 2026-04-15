# Cylinder at $Re_D=200$, N=5, 1824 elements
Enabled simcomps:
- fluid_statistics, basic, 2D
- lift/drag calculation
- lambda 2

Included files:
- Gmsh script to generate mesh
- neko .nmsh file (converted from gmsh-re2-nmsh)
- case file

# In-situ data processing

For this example, we also show how to perform in-situ data processing, i.e., 
manipulating data sent from the Neko simulation while it is running.

In this particular case, three different fields are streamed and interpolated
on a 2D structured grid for visualization using `matplotlib`. The files related to insitu execution are:
- `cylinder_insitu.case` includes the `data_streamer` simulation component, 
  as well as a constant `strm_update_freq` constant which dictates how 
  frequently data should be sent. Not that other relevant simcomps use this
  parameter such that they are only computed when required.
- `insitu_task.py` is the python script that performs the in-situ
  data processing.

There are some dependencies for this to work:
- Neko must be installed with `Adios2` support (`--with-adios2=DIR`) such 
  that data can be transferred to python.
  @note pySEMTools provides installation wrappers for `adios2` and `mpi4py` 
  as third party libraries. 
- PySEMTools must be installed, as the python script relies heavily on it.
- `numpy`, `matplotlib`

A successful execution of this example will generate a set of `png`
images (total size ~38MB) in the `output_directory`, currently set to 
`"results_insitu"`. Below is a visualization of the snapshots compiled with
`ffmpeg`. 

![Animation of the generated in-situ snapshots](cylinder_insitu.gif)

## Executing neko in MPMD
The insitu run can be performed by executing neko and python in Multiple Program Multiple Data (mpmd) mode.

On personal computers with mpi, one can execute, for example:
```bash
mpirun -n 1 neko cylinder_insitu.case : -n 4 python3 insitu_task.py
```
which will run the neko simulation on 1 rank and the python processing script 
on 4 ranks.

On supercomputers, for the particular case of slurm, one can run in mpmd mode by executing a comand such as:
```bash
srun --multi-prog mpmd.conf
```
`mpmd.conf`, in its simplest form, is a file that indicates what application each rank executes. It can be
generated, for example, as:

```bash
cat > mpmd.conf << 'EOF'
0 ./select_gpu ./neko cylinder_insitu.case
1 python3 insitu_task.py
2 python3 insitu_task.py
3 python3 insitu_task.py
4 python3 insitu_task.py
EOF
```

This example has been succesfully executed on supercomputers such as lumi, leveraging the typically idling CPUs on a GPU compute node.
