# Neko Reframe test suite

Requires Reframe to run (`spack install reframe` or `pip install reframe-hpc`).

```
reframe -C settings.py -c checks.py -r --performance-report
```

Substitute `-r` for `-l` to see what tests would run.

The device backend is chosen for partitions with a device with `type: gpu`.
Otherwise the cpu backend is used.

## Options

Run only a specific test:

```
# 2 iterations of hemi
-n MiniHemi

# Small TGV (512 elements)
-n Tgv8

# Medium TGV (32768 elements)
-n Tgv32

# math/device_math throughput and wrapper dispatch cost
-n MathOpsVerify -n MathOpsPerf
```

> **Attach references with `set_reference()`, never `reference.setdefault()`.**
> `reference` is a `ScopedDict`: assigning a whole `{scope: {var: ref}}`
> mapping flattens it to `scope:var` keys, but mutating the plain dict that
> `setdefault` returns just nests a dict under the scope key, where
> `check_performance` -- which looks up `scope:var` -- never finds it. The
> variable then has no bounds and passes unconditionally, while still
> appearing in the report. `Tgv8`/`Tgv32` attached `enstrophy_error` that way,
> so their accuracy gate was inert until this was fixed.

> **A performance variable that fails to evaluate is silently dropped.**
> ReFrame catches any exception from a perf function, logs `skipping
> evaluation of performance variable` at warning level and carries on
> (`check_performance` in `reframe/core/pipeline.py`). A regex that stops
> matching therefore does not fail the test -- the metric just vanishes from
> the report and its gate quietly stops existing. `MathOpsPerf` defends
> against this by asserting the expected `BENCH`/`RATIO` record counts in its
> sanity function; do the same for anything new.

> **A file that fails to load still exits 0.** ReFrame skips a check file it
> cannot import, prints `Found 0 check(s)` and returns success, so a single
> incompatibility turns the whole suite into a green no-op. This happened:
> `reframe-hpc` is installed unpinned in CI, and from 4.8 onwards
> `checks.py` failed to load because `OutOfSourceAutotools` redefined
> `configuredir`, which `ConfigureBasedBuildSystem` had gained. If you change
> `checks.py`, run it with `-l` first and check the count is non-zero;
> `check_reframe.yml` now does this automatically.

Control which configurations are used with environment variables:

```
# Real precision
NEKO_REAL=sp,dp

# Fluid scheme
NEKO_SCHEME=pnpn
```

## math_ops: tracking `math`/`device_math` and the wrapper libraries

`MathOpsVerify` and `MathOpsPerf` drive `tests/bench/math_ops`, comparing raw
`math` (or `device_math`, on a device build) against the `field_math`,
`vector_math` and `matrix_math` wrappers over `add2`, `col2` and `glsc3`.

**What is gated, and where:**

| variable | what it catches | gated on |
| --- | --- | --- |
| `overhead_mean_<op>_<wrapper>` | a wrapper starting to cost more than the direct call | everywhere |
| `overhead_spread_<op>_<wrapper>` | a regression confined to one problem size | everywhere |
| `glsc3_lx<N>` | a changed reduction, or one that is not rank-count invariant | everywhere, `dp` only |
| `workrate_{math,device_math}_<op>_lx<N>` | absolute throughput drift | nowhere yet -- see below |

The split is deliberate. The overhead variables are **ratios** of one timing
to another in the same run, so machine speed, compiler and precision all
cancel; they mean the same thing on a dedicated node and on a contended
2-core CI runner, and can gate every PR. Absolute throughput does not have
that property, which is why it is only recorded here.

### Filling in absolute throughput references

`workrate_*` is logged but has no reference, because a meaningful one has to
be measured on the machine it will gate. Follow what `Tgv8`/`Tgv32` do for
`total_runtime`: run `MathOpsPerf` a few times on a dedicated node, then add
to `MathOpsPerf.set_workrate_perf`

```python
if self.current_partition.fullname == 'dt:cpu':
    ref['workrate_math_add2_lx12'] = (<observed>, -0.10, None, 'Mdofs/s/pe')
```

Only the largest `lx` is worth gating -- the small sizes are dominated by
per-call effects and are noisier.

### Regenerating the pinned values

**`glsc3_lx<N>`** are values *of the pinned mesh at `dp`*. Regenerate after
changing the mesh, the fill, or the `lx` sweep:

```sh
cd tests/bench/math_ops
PKG_CONFIG_PATH=<neko-prefix>/lib/pkgconfig make
mpirun -np 1 ./mathbench ../nekbone/data/512.nmsh 0 | grep GLSC3
```

Then confirm the same values come back at `-np 2/4/8` before committing them:
that invariance is the property the reference encodes, so pinning numbers
that do not have it would bake in a bug. Observed spread across `-np 1/2/4/8`
was <= 2.4e-14 relative, against a 1e-12 tolerance.

**`overhead_mean_max` / `overhead_spread_cap`** are noise floors, so
recalibrate them only if the gate starts flaking, and use the worst case over
several runs rather than a single one:

```sh
for i in $(seq 1 8); do
  mpirun -np 1 ./mathbench ../nekbone/data/512.nmsh 200 | grep '^RATIO'
done
```

Group by `(op, path)`, take mean and stddev of `value` across the six `lx`
records within each run, then take the worst over runs. Calibrate against a
**contended** machine, not an idle one -- pinning to two cores against two
competing busy loops moved the worst case from 1.033/0.062 to 1.117/0.218.
The current caps (1.30 and 0.35) sit above that on purpose: the regressions
this can actually resolve are gross ones (an added copy roughly doubles the
ratio), while a lost inline is a few nanoseconds against an 8 us call and is
invisible at any threshold, so a loose gate that never flakes beats a tight
one that cries wolf.

### Checking the gate still fires

A gate never observed failing is not known to work. Add a deliberate cost to
a wrapper -- e.g. a redundant `copy` in `field_add2` -- rerun `MathOpsPerf`,
and confirm `overhead_mean_add2_field_math` fails before reverting.

## How to add a new system

The `settings.py` file defines how to run the suite on different systems.
Reference: https://reframe-hpc.readthedocs.io/en/stable/config_reference.html.

CPU partitions require processor information which can be automatically
generated using `reframe --detect-host-topology`. GPU partitions require a
device with `type: gpu`.

Also update `valid_systems` and potentially `valid_prog_environs` in `checks.py`.
