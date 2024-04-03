# Appendix {#appendix}

## Environmental variable reference {#env-var}

Name                    | Description                                                   | Default value
----                    | -----------                                                   | -------------
`NEKO_AUTOTUNE`         | Force Ax auto-tuning strategy (``'1D'``,``'KSTEP'``)          | Unset
`NEKO_LOG_LEVEL`        | Log verbosity level (integer > 0, default: 1)                 | Unset
`NEKO_GS_STRTGY`        | Gather-scatter device MPI sync. strategy (0 < integer < 5 )   | Unset

### Logging level details

A number of logging levels are supported.

- `NEKO_LOG_LEVEL=0`   : Quiet mode, minimal logging during execution.
- `NEKO_LOG_LEVEL=1`   : Default information mode, adding step informations.
- `NEKO_LOG_LEVEL=2`   : Verbose mode, logging extra details.
- `NEKO_LOG_LEVEL=10`  : Debug mode.

## Profiling {#profiling}

Neko can be profiled in a couple of different ways,

- Using the built in profiler (enabled with `--enable-profiling`).
- Using external profiling tools.
  - NVIDIA Nsight Systems.
  - AMD ROCm Profiler.
  - Craypat.

### Using the built in profiler

The built in profiler can be enabled by configuring Neko with
`--enable-profiling`. This will enable the profiler and add the necessary flags
to the build system. The profiler will then use `MPI_Wtime` to measure the time
spent in a predefined set of regions of the code. The CPU time spent in these
regions will be accumulated across all MPI ranks and written to a file called
`profile.log` in the working directory, along with an accumulated number of
executions for each region.
