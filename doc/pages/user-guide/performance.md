# Performance guidelines {#performance}

\tableofcontents

This is a short best practices guideline on how to achieve good
performance with Neko across various architectures. The guide covers
various configuration options, as well as advice on how to set up and
run simulations to achieve the best possible performance.

## Installation

When building Neko, it is important to enable full optimisation for
all configured backends. Configuring with `FCFLAGS` and `CFLAGS` to
`-O3` will provide an optimal build for CPUs and SX-Aurora backends,
while for accelerators, a user needs to set `CUDA_CFLAGS` or
`HIP_HIPCC_FLAG` depending on the backend type (OpenCL optimisation is
set in via the `CFLAGS` variable). Additional performance could also
be gained for both CPUs and GPUs by passing architecture-specific
configuration flags, either via `FCFLAGS` and `CFLAGS` (CPU, SX, and
OpenCL), `HIP_HIPCC_FLAGS` (HIP), or `CUDA_ARCH` (CUDA).

### Accelerator specific options

To avoid unnecessary data movement between host and device, it is
important to use a device-aware MPI implementation and configure Neko
with `--enable-device-mpi` for the CUDA and HIP backend. Furthermore,
to improve GPU utilisation, especially with few elements per device,
it could be beneficial to configure Neko to use OpenMP
`--enable-openmp`, and launch the simulation with two threads. This
would enable the task-parallel preconditioners inside Neko.

@note In the current release, OpenMP is only used to launch jobs
concurrently into different accelerator streams. Thus, no gains should
be expected from OpenMP when using a CPU backend.

## Simulation setup

This section explains how to set up and select options for a case to
achieve the best possible simulation performance.

### Load balancing

Load balancing is essential for good simulation performance. This can
be performed using graph partitioning of the mesh, either offline or
online (during a simulation with Neko). Both options require Neko to
be built with the optional dependency ParMetis (see @ref
installation).

The offline tool is called `prepart`, and is installed in the same
installation directory as `neko` and `makeneko`. To partition a mesh
into `nparts`, launch `prepart` as below:

```shell
./prepart <neko mesh> <nparts>
```

The output will be a new mesh, with an additional `_<nparts>` in the
filename to indicate that it has been balanced for `nparts`.

@note The offline tool can require substantial memory to partition
large meshes.

To enable load balancing (mesh partitioning) during a simulation, set
the `load_balancing` option to `true` (see \ref case-file). Runtime
load balancing will partition the mesh (in parallel) before the
time-stepping loop starts, and it will also save the balanced mesh
with an additional `_lb` in the filename.

@note Since mesh partitioning is non-deterministic, care has to be
taken when restarting a simulation such that one does not i) restart
with load balancing enabled and ii) use the original unbalanced
mesh. The suggested workflow is to run the first simulation with load
balancing enabled, and for all subsequent restarts, turn off load
balancing and use the balanced mesh with the additional `_lb` in the
filename.

 ### Parameters

Depending on the configured backend, some simulation parameters can
significantly impact performance, for example, the polynomial
order. For CPU backends, performance is less sensitive to the chosen
order and should be set depending on the accuracy and time-to-solution
needs of a given case. On accelerators, performance and device
utilisation are closely tied to the polynomial order. In general,
accelerators run best with the highest possible polynomial
order. However, using a very high order will force down the required
time-step size and, in turn, increase the time-to-solution for a given
case. A good rule of thumb is to use at least seventh-order
polynomials for accelerator backends.

Another important parameter to experiment with to improve performance
is the combination of linear solver and preconditioner, and their
associated tolerance criteria. The best combination is very case
dependent, thus it is best to experiment with the various types
provided (see #case-file).

## Running a simulation

When running a simulation, the only parameter a user has some control
over is the number of MPI ranks to use. Neko will always distribute a
mesh across all ranks in a job using a load-balanced linear
distribution. With too few elements per rank, communication costs will
start to dominate, reducing both scalability and performance. The
opposite, with too many elements per rank, will cause computational
cost (per rank) to increase at the cost of a reduced (overall)
parallel efficiency.

Neko has been optimised with strong scalability in mind, and has
demonstrated nearly ideal scaling and parallel efficiency with as few
as 10 elements per MPI rank on CPUs or 7000-10000 elements per MPI
rank on GPUs.

@note Given the significant computational capacity of GPUs, ideal
strong scalability, with as few elements per device as possible, is
not a guarantee of achieving the best time-to-solution. It is often
better to fill each GPU with as many elements as possible within the
device's memory. The indicative numbers given above (7000-10000
elements) should be seen as the capability of Neko to scale out a
problem across a large machine efficiently.

Finally, achieved performance depends on many factors, for example,
the interconnect; it is therefore advisable to invest time in
establishing a performance baseline the first time a new machine is
used.
