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
| `NEKO_AUTOTUNE`          | Force SEM operator kernel formulation (``'1D'``,``'KSTEP'``)          | Unset         |
| `NEKO_EB_TUNE`           | Sweep elements per block for the kstep kernels (boolean)              | CUDA: 1, HIP: 0 |
| `NEKO_EB`                | Elements per block candidate when `NEKO_AUTOTUNE=KSTEP` (0, 1 or 2)   | 0             |
| `NEKO_CHUNKS`            | Chunk size candidate when `NEKO_AUTOTUNE=1D` (0 to 3)                 | 0             |
| `NEKO_TUNE_ROUNDS`       | Interleaved sampling rounds used by the operator auto-tuner           | 3             |
| `NEKO_TUNE_ITERS`        | Kernel launches timed per candidate per round                        | 100           |
| `NEKO_LOG_FILE`          | Log file name, uses `stdout` if not set.                              | Unset         |
| `NEKO_LOG_TAB_SIZE`      | Number of spaces added for each level of indentation in the log file. | 1             |
| `NEKO_LOG_LEVEL`         | Log verbosity level (integer > 0, default: 1)                         | Unset         |
| `NEKO_GS_STRTGY`         | Gather-scatter device MPI sync. strategy (0 < integer < 5 )           | Unset         |
| `NEKO_GS_COMM`           | Gather-scatter communication backend                                  | Unset         |
| `NEKO_GS_CAF_SIGNALING`  | Coarray Fortran gather-scatter signaling mode                         | Unset         |
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

Any other value is reported as an error and the search runs as usual.

`NEKO_EB_TUNE` controls whether the elements per block dimension is
swept at all. It defaults to enabled on CUDA, where blocking measures a
clear gain, and disabled on HIP, where the blocked kernels exceed the
register budget and spill to scratch. Setting it to `0` or `1`
overrides the per-backend default.

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

When `NEKO_GS_COMM` is unset and the build has no device-aware MPI,
the host backends are benchmarked at initialisation and the fastest
one is kept (see @ref performance-gs-autotuning). `MPI` and
`NEIGHBOUR` are always benchmarked; `SHMEM`, `CAF` and `UTOFU` join
them if the build has the corresponding support. Set `NEKO_GS_COMM`
explicitly to skip the benchmark and pin a backend.

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
