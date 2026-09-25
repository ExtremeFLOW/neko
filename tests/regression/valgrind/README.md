# Valgrind regression test

This test is intentionally separated from the ReFrame case-file tests in
`tests/reframe`. Its purpose is to validate leak regressions for a minimal TGV
run without coupling to the example-validation suite.

## Requirements

- `makeneko` available on `PATH`
- `valgrind` available on `PATH`
- `python3` available on `PATH`
- a runtime environment that already exposes all required shared libraries for
  the generated `neko` executable

The script does not assume any particular installation layout for external
dependencies. Library paths and any other runtime environment setup must be
provided by the user or calling test harness.

## Run

From the repository root:

```bash
bash tests/regression/valgrind/run_valgrind_check.sh
```

Optional threshold overrides:

```bash
NEKO_VALGRIND_MAX_DEFINITE=0
NEKO_VALGRIND_MAX_INDIRECT=0
NEKO_VALGRIND_MAX_POSSIBLE=0
NEKO_VALGRIND_MAX_REACHABLE=1300000
```

Backend-specific defaults used by the runner when overrides are not provided:

- `CPU`: `definite=0`, `indirect=0`, `possible=0`, `reachable=1200000`
- `CPU (libxsmm)`: `definite=0`, `indirect=0`, `possible=0`, `reachable=1200000`
- `SX-Aurora`: `definite=0`, `indirect=0`, `possible=0`, `reachable=1200000`
- `Accelerator (CUDA)`: `definite=0`, `indirect=0`, `possible=3000`, `reachable=9000000`
- `Accelerator (HIP)`: `definite=0`, `indirect=0`, `possible=3000`, `reachable=9000000`
- `Accelerator (OpenCL)`: `definite=0`, `indirect=0`, `possible=3000`, `reachable=9000000`
- `UNKNOWN` backend: `definite=0`, `indirect=0`, `possible=0`, `reachable=1200000`

You can override limits globally (`NEKO_VALGRIND_MAX_*`) or per backend using:

- `*_CPU`
- `*_LIBXSMM`
- `*_SX`
- `*_CUDA`
- `*_HIP`
- `*_OPENCL`
- `*_UNKNOWN`

Example:

```bash
NEKO_VALGRIND_MAX_REACHABLE_CUDA=9500000 \
NEKO_VALGRIND_MAX_REACHABLE_HIP=11000000 \
bash tests/regression/valgrind/run_valgrind_check.sh
```

Optional working directory override:

```bash
bash tests/regression/valgrind/run_valgrind_check.sh /tmp/neko-valgrind-regression
```
