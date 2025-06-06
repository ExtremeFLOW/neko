# Appendices {#appendices}

\tableofcontents

The appendices contain a few extra pages that are not directly related to the usage
of the code. But can be useful for users and developers alike.

- [Environmental variable reference](@ref appendices_env-var)
- \subpage governing-equations
- \subpage publications

## Environmental variable reference {#appendices_env-var}

| Name                | Description                                                           | Default value |
| ------------------- | --------------------------------------------------------------------- | ------------- |
| `NEKO_AUTOTUNE`     | Force Ax auto-tuning strategy (``'1D'``,``'KSTEP'``)                  | Unset         |
| `NEKO_LOG_FILE`     | Log file name, uses `stdout` if not set.                              | Unset         |
| `NEKO_LOG_TAB_SIZE` | Number of spaces added for each level of indentation in the log file. | 1             |
| `NEKO_LOG_LEVEL`    | Log verbosity level (integer > 0, default: 1)                         | Unset         |
| `NEKO_GS_STRTGY`    | Gather-scatter device MPI sync. strategy (0 < integer < 5 )           | Unset         |
| `NEKO_GS_COMM`      | Gather-scatter communication backend                                  | Unset         |
| `NEKO_COMM_ID`      | Communicator id for this process (non-negative integer)               | 0             |

### Logging level details

A number of logging levels are supported.

- `NEKO_LOG_LEVEL=0`   : Quiet mode, minimal logging during execution.
- `NEKO_LOG_LEVEL=1`   : Default information mode, adding step informations.
- `NEKO_LOG_LEVEL=2`   : Verbose mode, logging extra details.
- `NEKO_LOG_LEVEL=10`  : Debug mode.

### Gather-scatter communication backend details

A number of gather-scatter backends are supported.

- `NEKO_GS_COMM=MPI`    : Host based MPI
- `NEKO_GS_COMM=MPIGPU` : Device based MPI
- `NEKO_GS_COMM=NCCL`   : NCCL/RCCL using its point-to-point interface
- `NEKO_GS_COMM=SHMEM`  : NVSHMEM based (for NVIDIA GPUs)