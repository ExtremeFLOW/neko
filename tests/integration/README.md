# Integration tests

Small cases, which tests that overall functionality have not been broken.
Ideally each test should not be longer than a few seconds.
The tests are run with `pytest`, you can install it using pip.

## Running

You need Python and the  `pytest` and `json5` Python packages to run the tests.

If you run Neko on the CPU  with `mpirun`, and `neko` and `makeneko` are in the
`$PATH`, you just execute
```
pytest
```
to run the tests

Otherwise there are three things you have to provide.

1. The backend to pytest using the `--backend` option. Can be one of `cpu`,
   `cuda`, `hip`, `opencl`. Defaults to `cpu`. This can be used to configure
   tests based on the active backend. 
2. A shell script that wraps the command used to execute neko. It should
   accept 3 arguments
   - Number of ranks.
   - The path to the neko executable.
   - The case file.
   There such launchers already provided. 
   - `default_cpu_launcher.sh` is used by default and uses `mpirun`.
   - `default_cuda_launcher.sh` may work with CUDA on some machines, but is more
     provided for inspiration.
3. The location of `neko` and `makeneko` via `NEKO_EXEC` and `MAKENEKO_EXEC`
   environmental variables.
4. Additionally, you can also provide `--max_nprocs` to override the number of
   ranks the tests are run on. Useful if you, for example, have 1 GPU and want
   to run tests, for which the number of ranks is set to something bigger.

So, running `pytest` is equivalent to
```
pytest --launcher-script=./default_cpu_launcher.sh --backend=cpu
```

## Adding new tests

- Take a look at `conftest.py` and `testlib.py` to get a sense of what fixtures 
  and convenience routines are available. 
- The `tests/test_demo/test_suite.py` has a simple test that runs Neko.
- LLMs are very good at `pytest`. If you upload the above-mentioned files to
  an LLM, it will help you a lot to write whatever test you want.
- You can get meshes and case files from the `examples` folder. Pick on one that
  is close to what you need, and manipulate it in Python. There is also a small
  cylinder mesh sitting in the `meshes`, which is recommended as the default
  one.

## Conventions.

- Each test suite should be in its own directory inside `tests`, where you can
  also add necessary stuff such as reference log files.
- Each test suite should be in a file test_\*.py, this way it is automatically
  discovered by `pytest`.
- Each test in the suite should be in a function which starts with `test_`.
- There is a fixture provided to generate a log file, that will be inside the 
  `logs` directory. This directory will be packaged by the CI for download and
  inspection.
