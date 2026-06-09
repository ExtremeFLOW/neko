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

Optional working directory override:

```bash
bash tests/regression/valgrind/run_valgrind_check.sh /tmp/neko-valgrind-regression
```
