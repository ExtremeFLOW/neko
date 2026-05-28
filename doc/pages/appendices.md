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
| `NEKO_AUTOTUNE`          | Force Ax auto-tuning strategy (``'1D'``,``'KSTEP'``)                  | Unset         |
| `NEKO_LOG_FILE`          | Log file name, uses `stdout` if not set.                              | Unset         |
| `NEKO_LOG_TAB_SIZE`      | Number of spaces added for each level of indentation in the log file. | 1             |
| `NEKO_LOG_LEVEL`         | Log verbosity level (integer > 0, default: 1)                         | Unset         |
| `NEKO_GS_STRTGY`         | Gather-scatter device MPI sync. strategy (0 < integer < 5 )           | Unset         |
| `NEKO_GS_COMM`           | Gather-scatter communication backend                                  | Unset         |
| `NEKO_GS_CAF_SIGNALING`  | Coarray Fortran gather-scatter signaling mode                         | Unset         |
| `NEKO_COMM_ID`           | Communicator id for this process (non-negative integer)               | 0             |
| `NEKO_MPI_THREAD_LEVEL`  | Requested MPI (and SHMEM) thread support level                        | Unset         |
| `NEKO_DEPRECATION_ERROR` | Whether to treat deprecated features as errors (boolean)              | Unset         |

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

When `NEKO_GS_COMM=CAF`, the per-pair synchronisation strategy is
selected by `NEKO_GS_CAF_SIGNALING`. The mode is bound on the first
gather-scatter initialisation and cannot change thereafter. The
default (when unset) is `sync`.

- `NEKO_GS_CAF_SIGNALING=sync`   : `sync images` over the union of
  neighbour pairs, with a double-buffered receive coarray (F2008).
- `NEKO_GS_CAF_SIGNALING=atomic` : Per-pair atomic counters via
  `atomic_define`/`atomic_ref` with a busy-wait spin (F2008).
- `NEKO_GS_CAF_SIGNALING=event`  : F2018 events (`event post` /
  `event wait`); requires a runtime that implements F2018 event
  semantics.
