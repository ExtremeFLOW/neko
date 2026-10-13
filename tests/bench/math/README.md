# math_ops benchmark

Compares four ways of expressing the same work on Neko data:

| path          | call                                | notes                          |
| ------------- | ----------------------------------- | ------------------------------ |
| `math`        | `add2(a, b, n)`                     | the raw routine, called direct |
| `field_math`  | `field_add2(fa, fb, n)`             | `field_t`, per-DOF CFD data    |
| `vector_math` | `vector_add2(va, vb, n)`            | `vector_t`                     |
| `matrix_math` | `matrix_add2(ma, mb, n)`            | `matrix_t`, `nrows=n, ncols=1` |

Ops covered: `add2` (in-place elementwise), `col2` (in-place elementwise
product), `glsc3` (MPI-reduced triple inner product).

On a CPU build all four end in the same `math.f90` routine, so any measured
difference is the wrapper's dispatch cost. On a device build the wrappers
dispatch to `device_math` instead; the driver handles both with no source
change.

## Building and running

Standalone, like `tests/bench/{ax,gs,nekbone}` -- not part of `make check`.
It builds against an *installed* Neko via pkg-config:

```sh
export PKG_CONFIG_PATH=<neko-prefix>/lib/pkgconfig:$PKG_CONFIG_PATH
make
mpirun -np 1 ./mathbench ../nekbone/data/512.nmsh 2000
```

Arguments are `<mesh> [niter]`. **`niter = 0` runs verification only** --
cheap, and the right mode for the multi-rank sweep below.

If the binary fails with `libjsonfortran.so...: cannot open shared object
file`, set `LD_LIBRARY_PATH` to the json-fortran/HDF5/ParMETIS lib dirs that
`neko.pc` refers to.

Set `OMP_NUM_THREADS=1`. `math.f90`'s routines carry `!$omp parallel do`, and
thread-team startup otherwise confounds a per-call cost of a few nanoseconds.

## Verification, including multiple ranks

Correctness runs before any timing and aborts the run on mismatch. Two things
are checked:

1. All four paths agree, elementwise for `add2`/`col2` and as a scalar for
   `glsc3`.
2. `glsc3`'s MPI-reduced value is **rank-count independent**. Run the same
   mesh at several rank counts and compare:

   ```sh
   for np in 1 2 4 8; do
     mpirun --oversubscribe -np $np ./mathbench ../nekbone/data/512.nmsh 0
   done
   ```

   The printed `glsc3 =` value must agree across all of them to ~13
   significant digits (exact agreement is not expected -- floating-point
   reduction is not associative, so a different decomposition legitimately
   changes the last bits).

   This is the check a single-rank run cannot make. A wrapper reducing
   rank-local data, or using the wrong communicator, produces a
   rank-count-dependent answer and is caught here and nowhere else.

Values are filled from the dofmap's **global coordinates**, not from a
rank-local index. SEM partitioning distributes elements, each carrying its
full `lx**3` block, so a rank-local index names different elements at
different rank counts -- which would make `glsc3` rank-count-dependent for a
reason that has nothing to do with the code under test. Coordinates are
decomposition-independent, so the reduction genuinely is too.

## Reading the output

Each line reports `min`, `mean`, `stddev` and `Mdofs/s/pe`. **Use `min`.**
Scheduler preemption, page faults and frequency changes can only make an
iteration slower, so the fastest observed call is the cleanest estimate of the
true cost, and it is the only statistic that makes a nanosecond-scale dispatch
difference visible at all. `mean`/`stddev` are printed alongside so that a
noisy run is recognisable as noisy rather than silently reported as clean.
`Mdofs/s/pe` uses the ReFrame workrate formula already used for Neko
(`tests/reframe/checks.py`): `1e-3 * dofs / time / pes`.

Timing uses `-np 1` by default. For the in-place ops the operand is restored
from a reference between iterations, and that restore sits **outside** the
timed window, so the reported time is the operation alone.

## Inspecting inlining

The interesting inlining -- `math.f90`'s routine into the wrapper -- is a
cross-translation-unit decision, so it does **not** happen when Neko is built
normally. It requires Neko to be built with LTO:

```sh
./configure --prefix=<prefix> ... FCFLAGS="-g -O2 -flto -ffat-lto-objects"
make && make install
```

`-ffat-lto-objects` matters: `libneko` is a static archive with no shared
variant, and the fat objects are what let a consumer that does not itself pass
`-flto` still benefit.

The remarks appear when **this driver links**, not while Neko builds, because
that link is where the LTO plugin runs over the archive:

```sh
make FFLAGS="-g -O2 -fopt-info-inline-optimized $(pkg-config --cflags neko)" \
  2>&1 | tee inline.log
grep -oE "Inlined (add2|col2|glsc3)/[0-9]+ into (field|vector|matrix)_[a-z0-9]+/[0-9]+" inline.log
```

Do not use `-fopt-info-inline-optimized=<file>` under a parallel `make`;
concurrent compiles clobber the same path. The bare form writes to stderr.

### What was observed (gfortran, `-O2 -flto -ffat-lto-objects`)

- `add2` and `col2` **are** inlined into all three wrappers.
- `glsc3` is **not** inlined into any of them -- its body carries an
  `MPI_Allreduce` and an extended-precision accumulator, so the inliner
  declines it on size. The `%size()` type-bound call *is* still inlined into
  the `glsc3` wrappers.
- Throughput was unchanged by the inlining, within noise. These ops are
  memory-bandwidth-bound, so removing a call boundary does not move a number
  that is set by how fast memory can be streamed.

## What this benchmark can and cannot resolve

The wrappers cost a fixed few nanoseconds per call. The op itself costs
~1.3 us at `n = 4096` and ~690 us at `n = 884736`. Wrapper overhead is
therefore 0.001-0.1% of the work and sits below the measurement floor at every
size Neko actually uses -- measured spread across the four paths is 1-3%, which
is noise, not signal. Making the benchmark more careful will not change this;
the overhead is genuinely negligible rather than merely hard to see.

The sweep runs down to `lx = 2` precisely because that is the only end where a
per-call cost is a measurable fraction of the work. Even there it is not
resolvable.
