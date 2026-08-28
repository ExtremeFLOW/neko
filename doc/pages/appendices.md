# Appendices {#appendices}

\tableofcontents

The appendices contain a few extra pages that are not directly related to the usage
of the code. But can be useful for users and developers alike.

- [Environmental variable reference](@ref appendices_env-var)
- \subpage governing-equations
- \subpage nmsh-format
- \subpage fld-format
- \subpage publications

## Environmental variable reference {#appendices_env-var}

| Name                     | Description                                                           | Default value |
| ------------------------ | --------------------------------------------------------------------- | ------------- |
| `NEKO_AUTOTUNE`          | Force SEM operator kernel formulation (``'1D'``,``'KSTEP'``,``'DMMA'``,``'DMMA_TMA'``,``'DMMA_TMA_BATCH'``,``'MFMA'``) | Unset |
| `NEKO_EB_TUNE`           | Sweep elements per block for the kstep kernels (boolean)              | 1             |
| `NEKO_EB`                | Elements per block candidate when `NEKO_AUTOTUNE=KSTEP` (0, 1 or 2)   | 0             |
| `NEKO_CHUNKS`            | Chunk size candidate when `NEKO_AUTOTUNE=1D` (0 to 3)                 | 0             |
| `NEKO_DMMA_NW`           | Warps per block candidate when `NEKO_AUTOTUNE=DMMA` (0, 1 or 2)       | 0             |
| `NEKO_DMMA_TMA_NW`       | Warps per block candidate when `NEKO_AUTOTUNE=DMMA_TMA` or `DMMA_TMA_BATCH` (0, 1 or 2) | 0 |
| `NEKO_MFMA_NWF`          | Wavefronts per block candidate when `NEKO_AUTOTUNE=MFMA` (0 to 3)     | 0             |
| `NEKO_MFMA_TUNE`         | Sweep the matrix core variant of the HIP Helmholtz operator (boolean) | 1             |
| `NEKO_TUNE_ROUNDS`       | Interleaved sampling rounds used by the operator auto-tuner           | 3             |
| `NEKO_TUNE_ITERS`        | Kernel launches timed per candidate per round                        | 100           |
| `NEKO_LOG_FILE`          | Log file name, uses `stdout` if not set.                              | Unset         |
| `NEKO_LOG_TAB_SIZE`      | Number of spaces added for each level of indentation in the log file. | 1             |
| `NEKO_LOG_LEVEL`         | Log verbosity level (integer > 0, default: 1)                         | Unset         |
| `NEKO_GS_STRTGY`         | Gather-scatter device MPI sync. strategy (0 < integer < 5 )           | Unset         |
| `NEKO_GS_COMM`           | Gather-scatter communication backend                                  | Unset         |
| `NEKO_GS_TUNE`           | Comm. backends the gather-scatter autotuning benchmarks (list)        | Unset (all but `CAF`) |
| `NEKO_GS_CAF_SIGNALING`  | Coarray Fortran gather-scatter signaling mode                         | Unset         |
| `NEKO_GS_RMA_FLUSH_ALL`  | Batch the MPI RMA gather-scatter payload flush (boolean)              | 1             |
| `NEKO_COMM_ID`           | Communicator id for this process (non-negative integer)               | 0             |
| `NEKO_HIP_ZEROCOPY`      | Zero-copy host/device mapping on unified memory (HIP), 1 enables      | 0             |
| `NEKO_METAL_ZEROCOPY`    | Zero-copy host/device mapping on unified memory (Metal), 0 disables   | 1             |
| `NEKO_MPI_THREAD_LEVEL`  | Requested MPI (and SHMEM) thread support level                        | Unset         |
| `NEKO_DEPRECATION_ERROR` | Whether to treat deprecated features as errors (boolean)              | Unset         |

### SEM operator auto-tuning details

On the CUDA and HIP backends each spectral element operator benchmarks
its available kernel formulations on first call and caches the winner
for the rest of the run; see @ref performance-operator-autotuning for
what is measured and why the defaults differ between vendors.

`NEKO_AUTOTUNE` pins a formulation and skips the search:

- `NEKO_AUTOTUNE=1D`    : always use the 1d variant, with the chunk size
  given by `NEKO_CHUNKS`.
- `NEKO_AUTOTUNE=KSTEP` : always use the kstep variant, with the
  elements per block given by `NEKO_EB`.
- `NEKO_AUTOTUNE=DMMA`  : always use the fp64 tensor core variant of the
  CUDA Helmholtz operator, scalar and vector, with the warps per block
  given by `NEKO_DMMA_NW`. Double precision only, `2 <= lx <= 8` for the
  scalar operator and `4 <= lx <= 8` for the vector one, on an sm_80 or
  sm_90 device. The other operators have no such variant and
  keep tuning as usual.
- `NEKO_AUTOTUNE=DMMA_TMA` : always use the TMA staged form of that
  variant, scalar and vector, with the warps per block given by
  `NEKO_DMMA_TMA_NW`. Double precision and `lx = 8` only, on an sm_90
  device, and needs a CUDA 12 or later toolkit.
- `NEKO_AUTOTUNE=DMMA_TMA_BATCH` : the batched form of the TMA variant,
  which stages all ten input cubes of an element at once. Same scope and
  same `NEKO_DMMA_TMA_NW` knob, but for the vector operator only --- the
  scalar operator has no such variant and reports the value as invalid.

- `NEKO_AUTOTUNE=MFMA`  : the HIP counterpart of `DMMA` --- always use the
  matrix core variant of the HIP Helmholtz operator, with the wavefronts per
  block given by `NEKO_MFMA_NWF`. Either precision, `4 <= lx <= 12`, on a
  gfx90a or gfx942 device, and for the scalar operator only.

The vector Helmholtz operator has no 1d formulation, so `1D` and `KSTEP`
both pin it to its kstep variant, and it has no matrix core variant on HIP.

Any other value is reported as an error and the search runs as usual.
`DMMA` on a build or a device without fp64 tensor cores, or `MFMA` on one
without matrix cores, or either at a polynomial order the variant does not
cover, is reported the same way. So are `DMMA_TMA` and `DMMA_TMA_BATCH`
on anything but an sm_90 device with a CUDA 12 toolkit, at any order but
`lx = 8`, or when the field pointers the operator is handed are not 16
byte aligned, which is checked per tune rather than assumed;
`DMMA_TMA_BATCH` additionally needs a device that will grant a block the
54800 B it stages into.

`NEKO_EB_TUNE` controls whether the elements per block dimension is
swept at all, and defaults to enabled on both backends. It was once off
on HIP, where the blocked kernels can exceed the register budget and
spill to scratch, but fixing that verdict at build time also stopped the
tuner from ever re-testing it. Where the blocked variants do spill the
tuner rejects them, so the cost of sweeping is tuning time rather than
run time. Set it to `0` to skip the sweep.

`NEKO_MFMA_TUNE` is the equivalent switch for the matrix core variant of
the HIP Helmholtz operator. It defaults to enabled wherever the hardware
and the polynomial order allow the variant at all; setting it to `0`
keeps the strategy out of the sweep while leaving
`NEKO_AUTOTUNE=MFMA` able to pin it explicitly.

`NEKO_CHUNKS` selects among the instantiated chunk sizes when the 1d
variant is pinned with `NEKO_AUTOTUNE=1D`. Candidates `0` to `3` are
1024, 512, 256 and 128 threads; candidate `0` is the historical value,
so `NEKO_AUTOTUNE=1D` on its own is the A/B baseline for the chunk
sweep. A candidate smaller than one derivative matrix (`lx*lx`) is
invalid and falls back to 1024.

`NEKO_EB` selects among the instantiated blocking candidates when a
formulation is pinned with `NEKO_AUTOTUNE=KSTEP`; candidate `0` is one
element per block, i.e. the unblocked geometry, which makes
`NEKO_AUTOTUNE=KSTEP` on its own the A/B baseline for the blocking.
Values outside the valid range fall back to `0`.

`NEKO_DMMA_NW` selects among the instantiated warps per block candidates
when the tensor core variant is pinned with `NEKO_AUTOTUNE=DMMA`;
candidates `0`, `1` and `2` are 2, 4 and 8 warps. Values outside the
valid range fall back to `0`. `NEKO_DMMA_TMA_NW` is the same selector for
the two TMA staged variants, with the same three candidates, and is read
by both `DMMA_TMA` and `DMMA_TMA_BATCH`. `NEKO_MFMA_NWF` is the HIP
equivalent, with candidates `0` to `3` selecting 1, 2, 4 and 8 wavefronts
per block.

On HIP that count also fixes the elements per block, and the two cannot
be set independently: `min(nwf, ceil(lx*lx/16))` wavefronts cooperate on
one element --- that being all the wavefront-parallel work a contraction
offers --- and the block covers `nwf` divided by that many elements. At
`lx = 8` the top candidate is eight wavefronts on one element, at `lx = 4`
it is eight elements with one wavefront each.

`NEKO_TUNE_ROUNDS` and `NEKO_TUNE_ITERS` control the sampling, and
apply to *both* sweeps --- they are not specific to the elements per
block dimension. Each round times every candidate once, chunk sizes
included, and the best round per candidate is kept, so raising the
round count guards against transients and clock drift during the sweep,
whereas raising the iteration count only reduces per-sample noise.

### Logging level details

A number of logging levels are supported.

- `NEKO_LOG_LEVEL=0`   : Quiet mode, minimal logging during execution.
- `NEKO_LOG_LEVEL=1`   : Default information mode, adding step informations.
- `NEKO_LOG_LEVEL=2`   : Verbose mode, logging extra details.
- `NEKO_LOG_LEVEL=5`   : Deprecated features will be logged if used.
- `NEKO_LOG_LEVEL=10`  : Debug mode.

### Gather-scatter communication backend details

A number of gather-scatter backends are supported.

- `NEKO_GS_COMM=MPI`    : Host based MPI
- `NEKO_GS_COMM=MPIGPU` : Device based MPI
- `NEKO_GS_COMM=NCCL`   : NCCL/RCCL using its point-to-point interface
- `NEKO_GS_COMM=SHMEM`  : NVSHMEM based on GPU builds (NVIDIA GPUs);
  OpenSHMEM based on CPU builds (requires a native OpenSHMEM library,
  e.g. Cray OpenSHMEMX, enabled at configure time with `--with-openshmem`)
- `NEKO_GS_COMM=CAF`    : Coarray Fortran (requires a coarray-capable compiler)
- `NEKO_GS_COMM=NEIGHBOUR` : Host MPI using an `MPI_Ineighbor_alltoallv`
  neighbourhood collective (`NEIGHBOR` is also accepted; MPI-3, host only)
- `NEKO_GS_COMM=UTOFU`  : Native Tofu interconnect (uTofu) one-sided RDMA
  (Fugaku and other Tofu-D systems; requires building with `--with-utofu`)
- `NEKO_GS_COMM=MPIRMA` : Host MPI one-sided puts into a passive-target
  window (`RMA` is also accepted; MPI-3, host only). Assumes an MPI whose
  RMA progresses without the target entering MPI, as the hardware-driven
  one-sided components do (Open MPI `osc/rdma` and `osc/ucx`, Cray MPICH
  over libfabric); `NEKO_GS_TUNE=-MPIRMA` drops it from the autotuning on
  systems where that does not hold

When `NEKO_GS_COMM` is unset and the build has no device-aware MPI,
the host backends are benchmarked at initialisation and the fastest
one is kept (see @ref performance-gs-autotuning). Which ones take part
is set by `NEKO_GS_TUNE`, a list of backend names spelled as for
`NEKO_GS_COMM` (comma or space separated, case insensitive). Unset, it
means every host backend the build supports except `CAF`, which a
compiler can accept at configure time while still giving the job a
single image the coarray backend cannot use.

- `NEKO_GS_TUNE=+CAF` adds the coarray backend to that default set,
  `NEKO_GS_TUNE=-SHMEM` removes a backend from it, and the two can be
  combined (`NEKO_GS_TUNE=+CAF,-SHMEM`).
- Plain names replace the set outright: `NEKO_GS_TUNE=MPI,NEIGHBOUR`
  benchmarks those two and nothing else. A single name pins that
  backend without benchmarking anything.
- The two forms cannot be mixed, and a name that is not a host
  gather-scatter backend is an error.

Backends that are selected but cannot run in this build or run are
dropped; those named explicitly are reported in the log as
`unavailable`. Set `NEKO_GS_COMM` to skip the benchmark and pin a
backend outright; that path ignores `NEKO_GS_TUNE`.

### MPI RMA backend details

`NEKO_GS_RMA_FLUSH_ALL` selects how the MPI one-sided backend
(`NEKO_GS_COMM=MPIRMA`) forces remote completion of the payload puts
before announcing them. MPI has no put-with-notify, so the data and the
signal are separate, unordered operations and the data must be flushed
in between.

- Unset or non-zero (the default): one `MPI_Win_flush_all` completes
  every put, then all the signals are issued.
- `NEKO_GS_RMA_FLUSH_ALL=0`: each peer gets its own `MPI_Win_flush` and
  is announced as soon as its put lands, which releases early receivers
  sooner but costs one flush call per peer.

The per-peer form is the better structure wherever a flush really is
scoped cheaply to the named target. On Cray MPICH over Slingshot it was
not -- the extra calls outweighed the pipelining they recovered -- so
the batched form is the default. It is worth re-measuring on an
unfamiliar MPI; the setting is read once, at the first gather-scatter
initialisation, so a run uses one strategy throughout.

### MPI thread level details

`NEKO_MPI_THREAD_LEVEL` overrides the default thread support requested
from the MPI runtime (and, when built with `--with-openshmem`, from
the SHMEM runtime as well). When unset, the level is chosen
automatically: `MPI_THREAD_MULTIPLE` when running with more than one
OpenMP thread per rank, falling back to `MPI_THREAD_FUNNELED` on
host-only backends, and `MPI_Init` otherwise. Device backends require
`MPI_THREAD_MULTIPLE` under the automatic policy. Setting the variable
bypasses these heuristics; if the runtime cannot honor the requested
level, initialisation aborts.

- `NEKO_MPI_THREAD_LEVEL=single`     : `MPI_THREAD_SINGLE` (calls `MPI_Init`).
- `NEKO_MPI_THREAD_LEVEL=funneled`   : `MPI_THREAD_FUNNELED`.
- `NEKO_MPI_THREAD_LEVEL=serialized` : `MPI_THREAD_SERIALIZED`.
- `NEKO_MPI_THREAD_LEVEL=multiple`   : `MPI_THREAD_MULTIPLE`.

### Coarray Fortran signaling mode details

When the CAF gather-scatter backend is used, the per-pair
synchronisation strategy is selected by `NEKO_GS_CAF_SIGNALING`. The
mode is bound on the first gather-scatter initialisation and cannot
change thereafter. The default (when unset) is `sync`.

- `NEKO_GS_CAF_SIGNALING=sync`   : `sync images` over the union of
  neighbour pairs, with a double-buffered receive coarray (F2008).
- `NEKO_GS_CAF_SIGNALING=atomic` : Per-pair atomic counters via
  `atomic_define`/`atomic_ref` with a busy-wait spin (F2008).
- `NEKO_GS_CAF_SIGNALING=event`  : F2018 events (`event post` /
  `event wait`); requires a runtime that implements F2018 event
  semantics.
- `NEKO_GS_CAF_SIGNALING=auto`   : Benchmark the modes above that the
  build supports and bind the fastest one. Only takes effect when the
  comm. backend is autotuned too (i.e. `NEKO_GS_COMM` unset), and only
  on the first gather-scatter that tunes CAF, as the mode is shared by
  every instance; falls back to `sync` otherwise. See
  @ref performance-gs-autotuning.

### uTofu injection details

When `NEKO_GS_COMM=UTOFU`, one injection virtual control queue (VCQ) is
created per OpenMP thread, dealt round-robin over the Tofu network
interfaces (TNIs) selected by `NEKO_GS_UTOFU_NTNI` (default `1`). Every
VCQ is created with `UTOFU_VCQ_FLAG_THREAD_SAFE`, and both the packing
of the send buffer and the `utofu_put` calls are parallelised over the
OpenMP team: each thread fires the puts for its share of the peers on
its own VCQ. This is the uTofu substitute for the per-thread
`MPI_Isend` that Fujitsu MPI's missing `MPI_THREAD_MULTIPLE` support
rules out. The counts are bound on the first gather-scatter
initialisation and cannot change thereafter; the TNI count is silently
clamped to the number of one-sided TNIs the hardware exposes, and the
VCQ count to whatever the TNIs' VCQ budget allows (threads beyond the
granted VCQ count do not inject).

- `NEKO_GS_UTOFU_NTNI=1` (default) : All injection VCQs sit on one TNI.
- `NEKO_GS_UTOFU_NTNI=N`           : Deal the per-thread VCQs across `N`
  TNIs (benefits large messages and high neighbour counts).
- `NEKO_GS_UTOFU_NVCQ=N`           : Cap the injection pool at `N` VCQs.
  The VCQ budget (roughly 8 per TNI per rank on Fugaku) is shared with
  the per-instance receive VCQs, so on tight budgets the pool must
  leave room.
- `NEKO_GS_UTOFU_NRVCQ=N`          : Receive VCQs per instance path
  (default: one per TNI; `1` restores a single receive queue). Each
  receive peer is bound to one of them, so an instance's incoming halo
  traffic is processed by several TNIs instead of funnelling through
  one.
- `NEKO_GS_UTOFU_MASTER_INJECT=1`  : Serialise all injection onto the
  master thread (still spread over every VCQ), as an A/B reference for
  the thread-parallel default.

The per-instance receive VCQs are placed round-robin over *all*
one-sided TNIs (not only the `NEKO_GS_UTOFU_NTNI` injection set),
sweeping past TNIs whose VCQ budget is exhausted; incoming puts are
routed by the advertised VCQ id, so this costs nothing and balances
receive processing across the interfaces. When the budget runs dry an
instance is granted fewer receive VCQs (at least one).

The startup log prints the granted counts, e.g.
`uTofu inj.   : 12 VCQs (12 threads) over 4 TNIs (6 requested)`.

The puts also request Tofu cache injection by default, which writes each
arriving slab straight into the receiver's cache. This helps when the
slab is consumed immediately but can pollute cache otherwise, so it can
be disabled for A/B testing:

- `NEKO_GS_UTOFU_CACHE_INJECT=1` (default, or unset) : cache injection on.
- `NEKO_GS_UTOFU_CACHE_INJECT=0`                     : cache injection off.

The fused vector (multi-component) exchange can be disabled with
`NEKO_GS_UTOFU_VEC=0`, in which case multi-component gather-scatter
falls back to independent scalar rounds (a validation/bisection aid).
