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

### CPU specific options

On CPU backends, the performance-critical path (math and solver
kernels, Krylov solvers and preconditioners, boundary conditions,
dealiasing, and the host gather-scatter) is threaded with OpenMP when
Neko is configured with `--enable-openmp`, making hybrid MPI+OpenMP a
supported execution mode: launch one rank per NUMA domain (e.g. one
per CMG on A64FX) and set `OMP_NUM_THREADS` to the number of cores in
it. Flat MPI (one rank per core, single-threaded) generally remains
the fastest configuration per node, but the hybrid mode reduces the
number of ranks -- and with it the memory footprint, initialisation
time, and communication metadata -- which pays off at very large
scale or under tight memory constraints.

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
with an additional `_lb_<nparts>` in the filename.

The load balancing is robust against restarts as it will check for the load
balanced mesh before performing the partitioning. If the load balanced mesh is
found, it will be used instead of the original mesh, and no partitioning will be
performed. It should be noted that a simulation restarted on a different number
of MPI ranks will trigger a new load balancing, and the balanced mesh will be
saved with the corresponding number of ranks in the filename.

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

### SEM operator auto-tuning {#performance-operator-autotuning}

On the CUDA and HIP backends the spectral element operators --- the
Helmholtz operator `Ax`, and `opgrad`, `dudxyz`, `cdtp`, `conv1`,
`convect_scalar` and `lambda2` --- each ship two kernel formulations,
and most of them a third --- with a fourth on Hopper, and a fifth for
the vector `Ax` --- and the best one depends on the polynomial order,
the element count and the hardware. Rather than fix a choice at build
time, Neko benchmarks them on the first call of each operator, using the
real element count of the running case, and caches the winner for the
rest of the run.

The two formulations are:

- the **1d** variant, which stages a whole element in shared memory and
  walks it with a flat thread block; and
- the **kstep** variant, which assigns one `(lx, lx)` thread plane per
  element and marches it in the third direction, holding the column in
  registers.

Five operators ship a third, which hands the tensor contractions of an
element to the matrix units of the device rather than to the ordinary
FMA pipes: on CUDA `Ax` (scalar and vector), `opgrad`, `dudxyz`, `conv1`
and `cdtp`; on HIP the same five less the vector `Ax`, which has no matrix
core kernel. `convect_scalar` and `lambda2` have no such variant on either
backend and keep tuning the two.

- the **dmma** variant on CUDA, which stages a cube in shared memory
  padded to `8^3` and issues each contraction on the fp64 tensor cores as
  an `8 x 64 x 8` GEMM built from `8 x 8 x 4` WMMA tiles. Its geometry
  parameter is the number of *warps per block* (2, 4 or 8), which stripe
  the eight tiles of a contraction among themselves.
- the **mfma** variant on HIP, which stages the element in LDS and issues
  each contraction on the matrix cores as a `D * U` GEMM with `M = lx`,
  `N = lx^2` and `K = lx`. Double precision uses the batched `4x4x4`
  tile, which is fully utilised in `M` for any `lx < 16`; single
  precision uses `16x16x4`, there being no `4x4x4` instruction for it.
  Its geometry parameter is the number of *wavefronts per block* (1, 2,
  4 or 8).

Padding the dmma cube to `8^3` keeps every tile full, so no lane ever has
to test an index, but it is exact only at `lx = 8`: at `lx = 4` just 64
of the 512 points are real and seven eighths of every contraction is
spent on zeros. Where `lx` divides 8 that waste is turned back into work
by packing `8/lx` elements along each axis --- `(8/lx)^3` per cube ---
and making the staged derivative matrix block diagonal, one `lx * lx`
copy per sub-cube. A contraction along any axis then couples only indices
within the same sub-cube, so the single padded contraction computes every
packed element independently and correctly. `lx = 8` packs one element
and is exactly what it was, `lx = 4` packs eight and `lx = 2` sixty-four,
none of it costing any extra shared memory. `lx = 3, 5, 6, 7` do not
divide 8 and keep one element and the old waste. The packing is aimed at
p-multigrid, which smooths at `lx = 4` and `lx = 2`, so low order `Ax` is
hot rather than incidental.

The four operators beyond `Ax` reuse that machinery unchanged on both
backends and differ only in what surrounds the contractions. `dudxyz`,
`opgrad` and `conv1` are the first half of `Ax` --- stage the field,
contract it three ways, and turn the result into the physical quantity
pointwise on the way out, with no second contraction --- so they stage
four cubes, 17.5 kB, exactly as the scalar `Ax` does. `cdtp` is the
other half: it forms three weighted fields first and then *accumulates*
three transposed contractions into a single output, the only one here
whose contractions follow its pointwise work rather than precede it.

The mfma variant has no padding to fill, but it meets the same imbalance
from the other side: a contraction offers only `ceil(lx^2/16)` column
groups of wavefront-parallel work --- one at `lx = 4`, four at `lx = 8`,
nine at `lx = 12` --- so a wavefront beyond that count has no matrix
core work left on the element. Rather than cap the sweep there, the
surplus wavefronts are given an element of their own: the block covers
`eb` elements, as many as the `ceil(lx^2/16)` groups leave wavefronts
for and rounded down to a power of two so that it divides `nwf`, and
`wpe = nwf / eb` wavefronts cooperate on each, the staging, geometry and
write-back passes parallelising over the whole block either way. At
`lx = 4` with eight wavefronts that is eight elements, one each, with
nothing idle; at `lx = 12` it is one element and eight cooperating
wavefronts.

On Hopper the dmma variant has two further forms, which differ from it
only in how an element reaches shared memory. The staging is nearly the
whole kernel --- at `lx = 8` an element is nine cubes of 4 kB, of which
the seven geometric factors alone are 78% of the read traffic --- and
plain dmma reads it with ordinary loads that nothing overlaps:

- the **dmma_tma** variant, which issues that traffic up front as bulk
  asynchronous copies on the Tensor Memory Accelerator, on two
  mbarriers: the cube of `u` on one, waited on immediately, and the
  seven geometric factor cubes on the other, waited on only at the
  pointwise step, so they are in flight underneath the first three
  contractions. The result cube goes back out as one bulk store rather
  than as scalar ones. Every operator that has a dmma variant has this
  one too, with its own set of cubes.
- the **dmma_tma_batch** variant, for the vector operator only, which
  stages all ten input cubes of an element --- the three components and
  the seven factors --- in a single batch of copies, each component
  keeping its own cube, in and then out, so that no load ever waits on a
  store. That is thirteen cubes, 54800 B, past the 48 kB a block gets
  without asking, so the allocation is dynamic and opted into once.

Neither buys the overlap for free. The seven extra cubes take the scalar
block from 17.5 kB to 45.5 kB, five resident blocks per multiprocessor
instead of eight, and the bet is that copies in flight replace the memory
level parallelism the lost warps were providing --- generated by one
thread per block rather than by thousands of loads. What that costs
differs per operator, because what has to stay resident does:

| operator              | copies batched | block   | blocks/SM |
| --------------------- | -------------- | ------- | --------- |
| `cdtp`                | 4              | 21.5 kB | 10        |
| `dudxyz`              | 5              | 33.5 kB | 6         |
| `Ax` scalar           | 8              | 45.5 kB | 5         |
| `opgrad`              | 10             | 53.5 kB | 4         |
| `Ax` vector (batched) | 10             | 53.5 kB | 4         |
| `conv1`               | 14             | 69.5 kB | 3         |

Anything past 48 kB is dynamic shared memory that the kernel opts into
once before its first launch, and the device is asked whether it will
grant it rather than assumed to. `cdtp` is the odd one out: it needs its
field *and* a factor before it can form anything, so a single wait would
leave the copies overlapping nothing, and its four are instead waited in
two stages, with half the factor traffic flying under the first
contraction.

On GH200 at `lx = 8` this measures out about 2% ahead of plain dmma for
the scalar `Ax`, while the batched vector form comes in some 6% ahead of
the component-at-a-time one, which only brings it level with the vector
dmma variant. For the gradient-type operators the margins are smaller
still: `dudxyz` picks dmma_tma, but by 0.8% over kstep, and `conv1` keeps
its 1d kernel with dmma_tma 0.7% behind it. Both of those are ties within
run-to-run drift. All of them are tuner candidates like every other one
and none is assumed to win.

@note Do not tune these by occupancy. Measured on GH200 at `lx = 8`, the
`dudxyz` and `conv1` TMA variants are flat to within 0.3% across the
whole 2, 4, 8 warp sweep --- for `conv1` that spans 6 to 24 warps per
multiprocessor, and its 6 warp configuration lands within 0.7% of a
1024 thread 1d kernel. The traffic is generated by one thread per block
regardless of how many warps exist, so once an operator is at its
bandwidth roof the block geometry stops mattering; the dmma variants
improve monotonically from 2 to 8 warps only because 8 is the first
candidate that fills the multiprocessor, and there is nothing beyond it
--- a contraction has eight tiles, so a sixteenth warp would have no
work.

@note The dmma variant only exists where the fp64 tensor cores do: a
double precision build, `2 <= lx <= 8` for the scalar `Ax`, `dudxyz`,
`opgrad`, `conv1` and `cdtp`, `4 <= lx <= 8` for the vector `Ax`, which
does not pack, and an sm_80
(A100) or sm_90 (H100 / GH200) device that the binary was actually
compiled for. Both the compile time guard and the runtime check
allow-list those two architectures rather than testing for "at least
sm_80", because consumer sm_86 / sm_89 assemble the instruction but have
no fp64 tensor hardware at all. A single precision build has no dmma
variant, and that is a hardware limit rather than an omission: NVIDIA
offers no fp32 tensor path, only TF32, whose 11 bit significand (unit
roundoff `4.9e-4`, against fp32's `6e-8`) is not usable for the operator
inside a Krylov solve.

@note The TMA variants need more than the tensor cores. The bulk copy
instructions are PTX 8.0, so a CUDA 12 or later toolkit is required, and
the architecture allow-list is sm_90 exactly --- sm_100 and later move
the bulk copy semantics and would have to be measured before being let
in. The order bound is `lx = 8` alone, because only there is a staged
cube the element's own contiguous chunk of global memory, which is what
lets a bulk copy address it with nothing but a base pointer and a length.
Every smaller order either strides into the padded cube or scatters
packed sub-cubes across it, which needs a tensor map re-encoded on the
host for every `Ax` call --- `u` and `w` are Krylov vectors and change on
each one --- and those are exactly the orders where dmma already measures
behind kstep. A bulk copy also requires 16 byte aligned addresses and
gives no diagnostic when it does not get them, so every pointer one will
touch is checked at tuning time rather than assumed: nine for the scalar
`Ax`, thirteen for the vector one, six for `dudxyz`, ten for `opgrad`,
fifteen for `conv1` and five for `cdtp`. The variants that stage past
48 kB additionally check that the device will grant a block the shared
memory they need.

@note The mfma variant covers `4 <= lx <= 12` in *both* precisions, on
gfx90a (MI250X) and gfx942 (MI300A / MI300X), for the scalar `Ax`,
`opgrad`, `dudxyz`, `conv1` and `cdtp` --- there is no vector mfma
kernel. `v_mfma_f32_16x16x4f32` is true IEEE fp32, so a single precision
build loses nothing in accuracy there; the upper order bound is the LDS
needed to keep the staged cubes resident, not the instruction.
Availability is confirmed by launching a probe kernel that reports
whether the device pass really was compiled for a matrix core
architecture, rather than by reading the device properties: a code
object selected from a fat binary built for something else would
otherwise turn the strategy into a silent no-op, launching, writing
nothing, and leaving stale values in the output.

@note Neither variant is a free win to expect. `Ax` is strongly bandwidth
bound, roughly 1.6 flop/byte at `lx = 8` against a ridge point near 9 on
GH200, so the matrix units can only pay by relieving shared memory and
issue pressure rather than by adding flops. The AMD side has already run
the experiment that removes precision from the argument: with a true fp32
matrix core, and therefore no accuracy compromise whatsoever, mfma still
measured 4.5% *behind* a plain 1d kernel at `lx = 8`. Whether either pays
on a given part is exactly what the tuner measures.

The vector (three component) Helmholtz operator runs its own search,
reported as a separate `Autotune Ax vector` section: the elements per
block of its kstep variant against, on CUDA, the warp counts of a vector
dmma variant and of the two TMA staged forms of it. It has no 1d
formulation, so `NEKO_AUTOTUNE=1D` selects kstep there, and no matrix
core variant on HIP, where the vector search is the elements per block
sweep alone. Blocking is not expected to pay much for the vector kernels
--- they sit at 254-255 registers, where a wider block changes threads
per block but not registers per thread, so it only saves the derivative
matrix loads --- but it is swept rather than assumed. Because the three
components share one set of geometric factors, the dmma variant reads
them once into registers and runs the components through the same staged
cubes rather than staging twelve of them. On GH200 at `lx = 8` that
loses narrowly to kstep --- holding the factors costs occupancy --- so
it is there for the low order end, where the register cost falls away
and the shared memory footprint does not.

Within each formulation the tuner also sweeps a geometry parameter. For
the kstep kernels this is the number of *elements per thread block*: a
single element is only half a warp at `lx = 4` and two warps at `lx = 8`,
so stacking several of them widens the block and lets the element
independent derivative matrices be loaded once per block instead of once
per element. For the 1d kernels it is the *chunk size*, which is both
the thread block size and the stride over the `lx^3` points of an
element; the historical value of 1024 leaves most of the block idle at
low order, where only 64 of 1024 threads have work to do at `lx = 4`.
Four sizes are measured (1024, 512, 256, 128), with any that would be
smaller than one derivative matrix falling back to 1024.

The tuner reports every candidate it measured, so the margins are
visible rather than implied. On a build or a device without the tensor
core variants that is the two formulations and their geometries:

```
---Autotune opgrad (lx: 4)----
  1D    ch=1024:     61.53 us/call
  1D    ch=512 :     32.10 us/call
  1D    ch=256 :     24.80 us/call
  1D    ch=128 :     21.40 us/call
  KSTEP eb=1   :     18.29 us/call
  KSTEP eb=2   :     15.40 us/call
  KSTEP eb=4   :     15.65 us/call
  Chose        : 2 (KSTEP, 2 elem/block)
```

The chosen line reports only the geometry that applies to the winning
formulation: an elements-per-block count when kstep wins, a chunk size
when the 1d variant does, and a warp or wavefront count when a matrix
core variant does, which on a part that has one adds lines of the
form

```
  DMMA  2w 8e  :     10.94 us/call
  DMMA  4w 8e  :     11.62 us/call
```

on CUDA, where the second field is the elements packed into one cube, and

```
  MFMA  4wf 1 e:    136.40 us/call
  MFMA  8wf 1 e:    136.40 us/call
```

on HIP, where it is the elements per block that the wavefront count
implies. The chosen line names the same pair, as
`Chose        : 3 (DMMA, 2 warps, 8 elem/blk)` or
`Chose        : 3 (MFMA, 4 wf, 1 elem/block)`.

At `lx = 8` on Hopper the TMA staged variants add lines of their own,
which carry no element count --- they stage a single element by
construction:

```
  TMA   2w     :     92.44 us/call
  TMA   4w     :     89.68 us/call
  TMAB  4w     :     89.61 us/call
```

`TMAB` is the batched vector form, so it appears in the
`Autotune Ax vector` section only, and the chosen line names them as
`Chose        : 4 (DMMA_TMA, 4 warps)` or
`Chose        : 5 (DMMA_TMA_BATCH, 4 warps)`.

A full section for one of the gradient-type operators, measured on GH200,
looks like this --- and is a fair illustration of how narrow the margins
between formulations can be once an operator is bandwidth bound:

```
----Autotune dudxyz (lx: 8)----
  1D    ch=1024:    233.25 us/call
  1D    ch=512 :    140.06 us/call
  1D    ch=256 :    147.39 us/call
  1D    ch=128 :    157.37 us/call
  KSTEP eb=1   :    112.70 us/call
  KSTEP eb=2   :    114.54 us/call
  KSTEP eb=4   :    123.90 us/call
  DMMA  2w 1e  :    119.13 us/call
  DMMA  4w 1e  :    115.52 us/call
  DMMA  8w 1e  :    113.11 us/call
  TMA   2w     :    112.18 us/call
  TMA   4w     :    111.82 us/call
  TMA   8w     :    112.06 us/call
  Chose        : 4 (DMMA_TMA, 4 warps)
```

@note The chunk sweep measures all four sizes rather than predicting a
best one. Thread utilisation alone suggests the smallest valid block,
but that is only right at the lowest orders: once the shared memory
footprint --- which grows as `lx^3` --- caps how many blocks fit on a
multiprocessor, a smaller block wastes thread slots instead of filling
them. Where the crossover falls depends on the shared memory per
multiprocessor, which differs across NVIDIA generations and again on
CDNA, so it is measured rather than assumed.

Sampling interleaves the variants over several rounds and keeps the best
round for each, rather than timing them one after another in a fixed
order --- on a machine whose clocks move during the sweep, a fixed order
biases the comparison by candidate position, which no amount of extra
iterations removes.

@note The measured ranking differs sharply between vendors, and by more
than the formulations alone. On an NVIDIA GH200 the kstep variant wins
and benefits from element blocking; on AMD gfx90a the same family
measures far worse than the 1d variant at `lx = 8`, and the blocked
kernels can overrun the register budget and spill to scratch. The
elements per block sweep is nevertheless enabled on both backends: a
verdict fixed at build time is one the tuner can never revisit, and
where the blocked variants do spill it simply rejects them. `NEKO_EB_TUNE`
turns the sweep off if the extra tuning time is not wanted.

The tuning behaviour is controlled by the environment variables
described in the @ref appendices_env-var reference: `NEKO_AUTOTUNE`
pins a formulation and skips the search entirely, `NEKO_EB_TUNE` and
`NEKO_MFMA_TUNE` enable or disable the elements per block and matrix
core sweeps, `NEKO_EB`, `NEKO_CHUNKS`, `NEKO_DMMA_NW`,
`NEKO_DMMA_TMA_NW` and `NEKO_MFMA_NWF` force a particular geometry when a
formulation is pinned, and
`NEKO_TUNE_ROUNDS` / `NEKO_TUNE_ITERS` control the sampling of both
sweeps. All
of them are useful mainly for A/B testing; the defaults are intended to
be right without intervention.

@note The search runs once per operator and polynomial order, at the
first call, and costs a few hundred kernel launches. On a case with a
very large element count this is visible in the startup time; reducing
`NEKO_TUNE_ITERS` shortens it, at the cost of noisier measurements.

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

### Gather-scatter communication backends

The dominant communication pattern in spectral element methods is the
gather-scatter (gs) operation that assembles partial element-local
contributions across element boundaries. Neko ships several
implementations of the off-process gs exchange, and the right choice
depends on the host/accelerator combination and the interconnect.

The active backend can be selected at runtime via the `NEKO_GS_COMM`
environment variable. If the variable is unset, a sensible default is
chosen based on the build configuration (device-aware MPI when device
MPI is available, and otherwise the fastest of the host backends the
build supports, picked by the autotuner described in
@ref performance-gs-autotuning). The supported values are:

| `NEKO_GS_COMM` | Backend | Requirement | Typical use |
|---|---|---|---|
| `MPI` | Host MPI (`gs_mpi`) | None | CPU runs, fallback for any system |
| `MPIGPU` | Device-aware MPI (`gs_device_mpi`) | `--enable-device-mpi`, GPU build | NVIDIA / AMD GPUs with GPU-aware MPI |
| `NCCL` | NCCL/RCCL (`gs_device_nccl`) | `--with-nccl=...` | Multi-GPU runs where NCCL outperforms MPI |
| `SHMEM` | NVSHMEM (`gs_device_shmem`) on GPU builds, OpenSHMEM (`gs_shmem`) on CPU builds | `--with-nvshmem=...` (GPU) or `--with-openshmem` (CPU, Cray) | NVIDIA GPUs with NVSHMEM-capable interconnect; CPU systems with a native OpenSHMEM (e.g. Cray OpenSHMEMX) |
| `CAF` | Coarray Fortran (`gs_caf`) | Coarray-capable Fortran compiler | Systems with a tuned coarray runtime |
| `NEIGHBOUR` | MPI neighbourhood collective (`gs_neighbour`) | None (MPI-3) | CPU runs where one collective per gs beats many point-to-point messages (e.g. Fugaku / Tofu) |
| `UTOFU` | Native Tofu RDMA (`gs_utofu`) | `--with-utofu` (Tofu-D, e.g. Fugaku) | CPU runs on Fugaku/Tofu; native RDMA, a faster successor to CAF |
| `MPIRMA` (or `RMA`) | MPI one-sided (`gs_mpi_rma`) | None (MPI-3) | CPU runs on systems with no OpenSHMEM, coarray or uTofu support, where one-sided is still worth trying; needs an MPI whose RMA progresses in hardware (see below) |
| `CRYSTAL` | Crystal router (`gs_crystal`) | None | CPU runs where the halo is many tiny messages: few elements per rank, coarse multigrid levels, or an MPI that handles many concurrent peers badly (see below) |
| `CRYSTALGPU` | Crystal router on the device (`gs_device_crystal`) | GPU-aware MPI, GPU build | The same, on GPU runs; the halo never leaves the device |

@note `MPIGPU` and `NCCL` require Neko to be built with GPU support
and the corresponding optional dependency. `SHMEM` picks the device
backend on GPU builds and the host backend on CPU builds, depending on
which of NVSHMEM / OpenSHMEM is present at configure time.

The backend can also be selected programmatically by passing the
`comm_bcknd` argument to `gs%init`, using the constants `GS_COMM_MPI`,
`GS_COMM_MPIGPU`, `GS_COMM_NCCL`, `GS_COMM_NVSHMEM`,
`GS_COMM_OPENSHMEM`, `GS_COMM_CAF`, `GS_COMM_NEIGHBOUR`, `GS_COMM_UTOFU`,
`GS_COMM_MPIRMA`, `GS_COMM_CRYSTAL` or `GS_COMM_CRYSTALGPU`
exposed by the `gather_scatter` module. The environment variable wins when
both are present.

#### The crystal router backends {#performance-gs-crystal}

`CRYSTAL` and `CRYSTALGPU` do not send one message per peer. They route
the halo to its destinations in recursive-bisection stages, keeping the
records already on the right side of the split and handing the rest to a
single partner, so one gather-scatter costs at most one message per stage
whatever the peer count is.

Because the schedule is fixed, the routing is worked out once at
initialisation: the partner ranks, the exact word counts and the index
lists that move the words which stay are all known before the first
exchange, so the runtime never scans destinations, negotiates sizes or
allocates. Stages with nothing to send or receive are dropped outright,
and since a partition's peers are close in rank index most of them are,
which leaves the number of stages set by the largest rank distance among
the peers rather than by @f$ \log_2 P @f$.

What it costs is store and forward. A word bound several stages away
crosses the network once per stage it survives and is copied locally each
time, and the stages are dependent, so only the first one overlaps the
local gather-scatter. This is a trade of bandwidth and latency for message
count: it pays where per-message overhead dominates and loses where the
halo is already large enough to be bandwidth bound. `CRYSTAL` is
benchmarked by the autotuning like any other host candidate, which is the
intended way to find out which regime a given run is in;
`NEKO_GS_TUNE=-CRYSTAL` drops it. `CRYSTALGPU` is not autotuned and has to
be asked for by name.

#### Runtime autotuning of the host backend {#performance-gs-autotuning}

Which host backend wins depends on the MPI implementation, the
interconnect and the halo of the particular decomposition, so the
choice is made by measurement rather than by a built-in rule. When
`NEKO_GS_COMM` is unset, no `comm_bcknd` argument is passed to
`gs%init` and the run has more than one rank, each `gs_t` instance
benchmarks the available host backends at initialisation and keeps the
fastest one. This mirrors the autotuning of the device MPI
synchronisation strategy (`NEKO_GS_STRTGY`) already done on
accelerator builds.

By default the candidates are every host backend the build supports
except `CAF`: `MPI`, `NEIGHBOUR`, `MPIRMA` and `CRYSTAL`, which need
nothing but MPI-3, plus `SHMEM` (OpenSHMEM) and `UTOFU` when the corresponding
support was built in -- a backend whose support is missing aborts in
its `init`, so it is left out of the list rather than tried. `SHMEM`
and `CAF` are additionally skipped when `NEKO_COMM` does not span every
process (`NEKO_COMM_ID`), since they address their peers by global PE /
image number; `UTOFU` and `MPIRMA` exchange their addresses over
`NEKO_COMM` itself and are kept in that case.

@note `MPIRMA` assumes the MPI implementation makes one-sided progress
without the target entering MPI. That holds for hardware-driven
components (Open MPI `osc/rdma` and `osc/ucx`, Cray MPICH over
libfabric) and not for `osc/pt2pt`, where its spin waits make it slow
rather than wrong -- it loses the benchmark and is dropped, at the cost
of the time spent measuring it. Use `NEKO_GS_TUNE=-MPIRMA` to skip it
outright on such a system.

##### Choosing the candidates

`NEKO_GS_TUNE` overrides that set. It takes a list of backend names,
spelled as for `NEKO_GS_COMM`, separated by commas or spaces and
matched case insensitively:

| `NEKO_GS_TUNE` | Candidates |
|---|---|
| unset | every supported host backend but `CAF` |
| `+CAF` | the default set, plus the coarray backend |
| `-MPIRMA` | the default set, without the MPI one-sided backend |
| `-CRYSTAL` | the default set, without the crystal router backend |
| `-SHMEM` | the default set, without OpenSHMEM |
| `+CAF,-SHMEM` | both deltas applied to the default set |
| `MPI,NEIGHBOUR` | exactly those two, whatever the build supports |
| `UTOFU` | uTofu alone -- nothing to compare, so it is simply used |

Names prefixed with `+`/`-` modify the default set, plain names
replace it, and mixing the two forms is an error, as is naming
something that is not a host gather-scatter backend. Backends selected
but unusable in this build or run are dropped from the comparison; a
backend that was named explicitly says so in the log:

```
 CAF          : unavailable
```

`CAF` is out of the default set because coarray support at configure
time says nothing about the job running more than one image: a
compiler may accept the syntax and still build every process as a
single-image program (gfortran's `-fcoarray=single`, or a coarray
runtime that was never linked in), where the backend has no peers to
put to. The case Neko can detect, `num_images() /= pe_size`, makes
`CAF` unusable no matter how it was asked for, and requesting the
backend outright with `NEKO_GS_COMM=CAF` then fails with an error
rather than exchanging nothing. When `CAF` is benchmarked, it runs in
the single signaling mode bound by `NEKO_GS_CAF_SIGNALING` (`sync`
unless set), unless that variable is set to `auto`, which tunes the
modes as well -- see @ref performance-caf-backend.

The gs schedule is computed once, by `gs_schedule`, and handed over
from one candidate backend to the next (`gs_comm_t%take_schedule`), so
tuning costs one extra backend setup plus a short burst of
gather-scatter operations (100 timed, 2 untimed, per candidate) on a
dummy vector, not a second pass over the connectivity. The handover is
split in two so that the outgoing backend is freed before the incoming
one is set up: backends holding scarce resources (uTofu VCQs,
symmetric memory, coarrays) never overlap. `GS_OP_MIN` is used for the
benchmark so that the working vector is left unchanged no matter how
many trials are run. The timings are averaged over all ranks with an
`MPI_Allreduce` before the winner is picked: setting up and driving a
backend is collective (the neighbourhood backend builds a graph
communicator, the one-sided backends allocate symmetric memory), so
every rank must arrive at the same decision. The log reports the
average time per gather-scatter operation for each candidate and the
backend that was kept:

```
 Comm         :         auto
 ...
 MPI          :  5.176E-05 s
 MPI neigh.   :  6.147E-05 s
 uTofu        :  1.732E-04 s
 Tuned comm   :          MPI
```

Tuning is skipped, keeping the host MPI backend, if any rank holds
zero dofs: such a rank skips the halo exchange entirely, which would
hang the other ranks in a collective-based backend. It is likewise
skipped when `NEKO_GS_TUNE` leaves fewer than two candidates, the
single candidate then simply being switched to. Since each `gs_t`
instance is tuned separately, a multigrid hierarchy tunes each level
on its own message sizes. Set `NEKO_GS_COMM` explicitly to pin a
backend and skip the benchmark, for example when comparing runs where
the tuning outcome itself should not vary.

@note Benchmarking `CAF` allocates the module-level receive coarray
shared by all `gs_caf_t` instances. That buffer is grown on demand and
retained for the lifetime of the program, so the memory stays
allocated even when the coarray backend loses the benchmark.

#### MPI neighbourhood collective backend {#performance-neighbour-backend}

The neighbourhood backend (`NEKO_GS_COMM=NEIGHBOUR`, the spelling
`NEIGHBOR` is also accepted) replaces the per-neighbour
`MPI_Isend`/`MPI_Irecv` pairs of the host `MPI` backend with a single
non-blocking neighbourhood collective, `MPI_Ineighbor_alltoallv`, over
a distributed-graph communicator. The graph is built once at
initialisation with `MPI_Dist_graph_create_adjacent` (sources are the
receive peers, destinations the send peers, `reorder = .false.`), so
the neighbour ordering matches the gs schedule and the concatenated
per-peer send / receive buffer layout of the host MPI backend is reused
unchanged. `nbsend` packs the send buffer and launches the collective
from the master thread; `nbwait` completes it and reduces each received
slab into `u` with the same serial-across-peers, parallel-within-peer
reduction as the host MPI backend.

The motivation is to collapse the N point-to-point calls of a halo
exchange into a single entry to the MPI runtime. On platforms where the
MPI library serialises concurrent calls internally (for example
Fujitsu MPI on Fugaku), this is one runtime crossing per gs op instead
of N, and it hands the striping of the transfers across the network
interfaces to the vendor runtime. The backend is host-only and requires
nothing beyond MPI-3, which every modern MPI implementation provides.

@note Whether the collective actually overlaps with the local gs work
depends on the MPI library making asynchronous progress; with the
default build of many implementations, progress happens at the
`MPI_Wait` in `nbwait`. Enable the library's asynchronous-progress
option to test the overlap benefit. Each `gs_t` instance owns its own
graph communicator, so multiple instances may coexist, but a given
instance must complete its `nbsend`/`nbwait` window before the next gs
op on the same instance.

#### NVSHMEM backend {#performance-nvshmem-backend}

The NVSHMEM backend (`NEKO_GS_COMM=SHMEM` on CUDA builds) keeps the
entire exchange on the GPU: a fused CUDA kernel
(`cuda_gs_pack_and_push`) packs the gathered shared dofs straight
from the device-resident solution buffer and issues
`nvshmemx_{float,double}_put_signal_nbi_block` to deliver each
neighbour's slice into the remote receive buffer together with a
per-sender completion signal. A `cuda_gs_unpack` kernel applies the
gather-scatter operation back into `u` on the device. The host never
sees the payload, so the backend is most effective on builds without
GPU-aware MPI or where NCCL/MPI go through host staging.

Each neighbour is driven on its own high-priority CUDA stream, with a
CUDA event recorded after the unpack so the caller's stream can wait
for completion -- this lets the receives complete out of order and
overlaps with the local gs work in `gs%bcknd`. Per-pair
`notifyDone`/`notifyReady` uint64 counters (allocated symmetrically
via `cudamalloc_nvshmem`) handle the double-buffered handshake so a
fast sender cannot overwrite a receive buffer the consumer has not
yet drained.

The backend is gated on `HAVE_CUDA` and `HAVE_NVSHMEM`. It is enabled
at configure time with `--with-nvshmem=DIR`, which requires a working
CUDA toolchain; non-CUDA builds and CUDA builds without NVSHMEM are
unaffected.

@note NVSHMEM requires a launcher and interconnect that support the
NVSHMEM transport (typically Slingshot, InfiniBand with GDR, or
NVLink).

#### OpenSHMEM backend {#performance-openshmem-backend}

The OpenSHMEM backend (`NEKO_GS_COMM=SHMEM` on CPU builds) uses
one-sided `shmem_putmem_signal_nbi` (OpenSHMEM 1.5) to deliver each
neighbour's contribution directly into a symmetric receive buffer
together with a per-sender completion signal. Receivers wait on
`shmem_signal_wait_until` per source rank and acknowledge the round
with `shmem_uint64_atomic_set`, which the sender consults before
overwriting the buffer on the next call. The handshake is therefore
fully per-pair and avoids any global barrier inside the gs op.

The backend is gated on the autoconf macro `HAVE_OPENSHMEM`. It is
currently wired up for Cray systems via `--with-openshmem` (selecting
either Cray OpenSHMEMX or Cray SHMEM, whichever is loaded in the
build environment); builds without a native OpenSHMEM library are
unaffected.

@note The OpenSHMEM backend allocates symmetric send and receive
buffers sized to the global maximum total send / receive count, plus
two symmetric `uint64` arrays of length `pe_size` per gs instance for
the data and ack signals. Per-rank slot indexing avoids the slot
exchange a neighbour-list scheme would require. Multiple `gs_t`
objects can coexist provided they are used strictly sequentially --
no overlapping nbsend / nbwait windows across instances.

#### Coarray Fortran backend {#performance-caf-backend}

The CAF backend (`NEKO_GS_COMM=CAF`) uses one-sided coarray puts
directly into a symmetric receive buffer, with a runtime-selectable
synchronisation strategy controlled by `NEKO_GS_CAF_SIGNALING`:

| `NEKO_GS_CAF_SIGNALING` | Standard | Mechanism | Notes |
|---|---|---|---|
| `sync` (default) | F2008 | `sync images` over the union of neighbour pairs, with a double-buffered receive coarray so only one rendezvous is needed per gs op | Most portable; works on every coarray-capable compiler |
| `atomic` | F2008 | Per-pair atomic counters via `atomic_define`/`atomic_ref` and a busy-wait spin | Avoids the image-set barrier; trade-off depends on the relative cost of pairwise atomics versus `sync images` on the target runtime |
| `event` | F2018 | F2018 events (`event post`/`event wait`) | Lowest theoretical overhead; requires a runtime that implements F2018 event semantics |
| `auto` | -- | Benchmarks the modes above that the build supports and binds the fastest | Only takes effect when the comm. backend is autotuned as well, with `CAF` among the candidates (`NEKO_GS_TUNE=+CAF`); see below |

The signaling mode is bound on the first gs initialisation and cannot
be changed thereafter. If the chosen mode requires a feature
unavailable in the build (e.g. `event` on a compiler without F2018
events), Neko aborts with a clear error.

With `NEKO_GS_CAF_SIGNALING=auto` the mode is chosen by measurement
instead, as part of the host backend autotuning described in
@ref performance-gs-autotuning -- which needs `CAF` among the tuning
candidates (`NEKO_GS_TUNE=+CAF`), as it is not one by default. When
that tuning reaches the `CAF` candidate, it times each available mode
on a freshly built backend -- the mode determines what per-instance
state `gs_caf_init` sets up, so
switching modes means rebuilding -- and the fastest one becomes the
`CAF` entry in the backend comparison:

```
 CAF sync     :  9.882E-05 s
 CAF atomic   :  7.431E-05 s
 Tuned CAF sig:       atomic
 MPI          :  5.176E-05 s
 MPI neigh.   :  6.147E-05 s
 CAF (atomic) :  7.431E-05 s
 Tuned comm   :          MPI
```

The per-mode lines are written while the `CAF` candidate is being
measured, so they precede the summary of the backend comparison.

The binding is program-wide -- every `gs_caf_t` instance reads the
same module-level mode -- so it is benchmarked only on the *first*
`gs_t` that tunes CAF and kept for the rest of the run. Re-binding it
later would swap the protocol under instances that are already live,
and in `atomic` mode desynchronise the round counters they share.
`auto` therefore has no effect when the comm. backend itself is not
autotuned (`NEKO_GS_COMM=CAF`, or a `comm_bcknd` argument), nor when
the autotuning runs without `CAF` among its candidates
(`NEKO_GS_TUNE` without `CAF`); the mode falls back to `sync` in those
cases.

@warning A signaling mode that the runtime does not really implement
fails by hanging, not by aborting, and the benchmark has no timeout.
`auto` is opt-in for that reason: `sync` remains the default, and
`atomic` / `event` are only exercised when explicitly asked for.

@note The CAF backend allocates a symmetric receive coarray sized to
the global maximum total receive count (doubled for the buffer
parity), shared by all gs instances. The buffer is grown on demand and
retained for the program lifetime. Multiple `gs_t` objects can coexist
provided they are used strictly sequentially -- no overlapping
nbsend/nbwait windows across instances.

#### uTofu backend {#performance-utofu-backend}

The uTofu backend (`NEKO_GS_COMM=UTOFU`) talks to the Tofu-D
interconnect directly through the native uTofu (`libtofucom`) API,
bypassing both MPI and the coarray runtime. It is the
performance-oriented successor to the CAF backend on Fugaku and other
Tofu systems, where coarray puts are roughly an order of magnitude
slower than MPI because every put is followed by a separate, serial
master-thread synchronisation step.

Each neighbour's contribution is packed into a per-instance send buffer
and delivered with a one-sided `utofu_put` straight into the remote
receive buffer. The slab tag is fused into the put's `edata` field, so
the put's arrival *is* the completion signal -- there is one network
operation per peer, with no separate atomic, credit message, or back
channel. This single-op-per-peer pattern is the main win over CAF.

Injection is non-blocking and thread-parallel. Packing each neighbour's
contribution into the send buffer is parallelised over the OpenMP team
in a single work-shared loop, and the `utofu_put` calls are as well:
one injection VCQ is created per OpenMP thread (with
`UTOFU_VCQ_FLAG_THREAD_SAFE`, which is what makes multithreaded uTofu
operation defined), dealt round-robin over the `NEKO_GS_UTOFU_NTNI`
TNIs (default `1`), and each thread fires the puts for its share of the
peers on its own VCQ. Since Fujitsu MPI does not provide a usable
`MPI_THREAD_MULTIPLE`, this is the only per-thread injection path
available on Fugaku; `NEKO_GS_UTOFU_MASTER_INJECT=1` restores serial
master-thread injection as an A/B reference. On the receive side each
slab is reduced into `u` as soon as its remote-put notice is observed
in the MRQ, mirroring the `Testsome`-style unpack-as-ready loop of the
MPI backend.

The send buffer is double-buffered (selected by the parity bit): each
put's completion is charged, via its `cbdata`, to a per-(VCQ, half)
counter, and the next `nbsend` drains only the half it is about to
reuse -- which was last used two rounds earlier, so the drain is a
non-blocking reap rather than a wait on the critical path. Each VCQ is
drained by exactly one thread per round, so the counters need no
atomics.

Correctness under skew is handled by a double-buffered receive region
selected by a parity bit (`edata = slab_tag*2 + parity`); a per-instance
stash holds notices that arrive from a neighbour running one round ahead
and replays them on the next round, so the unpack-as-ready path stays
correct without a global barrier or a back-pressure channel.

The backend also implements the fused vector (multi-component) halo
exchange used by `gs_op_r3`, so a three-component gather-scatter costs
one halo round instead of three. The vector path is a second, fully
independent instance of the same transport: its own registered
double-buffered slabs (sized for `GS_VEC_NC` components, each peer slab
holding `nc` consecutive component blocks), its own receive VCQ (and
thus private MRQ), and its own parity. The separation is what keeps the
skew handling sound -- a neighbour racing from a vector round into a
scalar round (or vice versa) lands its notice in the right queue. The
scalar per-peer layout is reused scaled by `nc`, with the destination
offsets recomputed each round from the receiver's advertised unscaled
layout.

The backend is gated on the autoconf macro `HAVE_UTOFU` and enabled at
configure time with `--with-utofu[=DIR]` (linking `-ltofucom`); builds
without uTofu are unaffected. Neighbour metadata (remote VCQ id,
receive STADD, the two parity offsets, and the slab tag) is exchanged
once at initialisation over `NEKO_COMM` via MPI; the data path
thereafter is pure uTofu.

@note The uTofu backend uses per-instance (not symmetric) send and
receive buffers, so it carries no global-maximum sizing constraint.
Each `gs_t` instance also gets its own receive VCQs (a set per path,
scalar and vector, spread over the TNIs with each receive peer bound
to one -- see `NEKO_GS_UTOFU_NRVCQ`), so an arrival notice for one
instance or path can never be consumed by another's poll -- multiple
`gs_t` objects can be used back to back safely. The injection VCQs
(one per OpenMP thread) are process-global and shared across
instances; send-completion accounting is per instance, per VCQ, and
per buffer-half.
