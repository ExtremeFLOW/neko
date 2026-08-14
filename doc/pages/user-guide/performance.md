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

@note `MPIGPU` and `NCCL` require Neko to be built with GPU support
and the corresponding optional dependency. `SHMEM` picks the device
backend on GPU builds and the host backend on CPU builds, depending on
which of NVSHMEM / OpenSHMEM is present at configure time.

The backend can also be selected programmatically by passing the
`comm_bcknd` argument to `gs%init`, using the constants `GS_COMM_MPI`,
`GS_COMM_MPIGPU`, `GS_COMM_NCCL`, `GS_COMM_NVSHMEM`,
`GS_COMM_OPENSHMEM`, `GS_COMM_CAF`, `GS_COMM_NEIGHBOUR` or `GS_COMM_UTOFU`
exposed by the `gather_scatter` module. The environment variable wins when
both are present.

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

The candidates are `MPI` and `NEIGHBOUR`, which need nothing but
MPI-3 and are therefore always benchmarked, plus `SHMEM`
(OpenSHMEM), `CAF` and `UTOFU` when the corresponding support was
built in -- a backend whose support is missing aborts in its `init`,
so it is left out of the list rather than tried. `SHMEM` and `CAF` are
additionally skipped when `NEKO_COMM` does not span every process
(`NEKO_COMM_ID`), since they address their peers by global PE / image
number; `UTOFU` exchanges its addresses over `NEKO_COMM` itself and is
kept in that case. `CAF` is benchmarked in the single signaling mode
bound by `NEKO_GS_CAF_SIGNALING` (`sync` unless set), unless that
variable is set to `auto`, which tunes the modes as well -- see
@ref performance-caf-backend.

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
hang the other ranks in a collective-based backend. Since each `gs_t`
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
| `auto` | -- | Benchmarks the modes above that the build supports and binds the fastest | Only takes effect when the comm. backend is autotuned as well; see below |

The signaling mode is bound on the first gs initialisation and cannot
be changed thereafter. If the chosen mode requires a feature
unavailable in the build (e.g. `event` on a compiler without F2018
events), Neko aborts with a clear error.

With `NEKO_GS_CAF_SIGNALING=auto` the mode is chosen by measurement
instead, as part of the host backend autotuning described in
@ref performance-gs-autotuning. When that tuning reaches the `CAF`
candidate, it times each available mode on a freshly built backend --
the mode determines what per-instance state `gs_caf_init` sets up, so
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
autotuned (`NEKO_GS_COMM=CAF`, or a `comm_bcknd` argument); the mode
falls back to `sync` in that case.

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
