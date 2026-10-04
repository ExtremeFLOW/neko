# Dealiased advection check

`dealias_check.py` exercises the dealiased advection (`adv_dealias`) in every
combination of its options on two small cases shipped with Neko, the affine
TGV box (`examples/tgv/512.nmsh`) and the curved 2D cylinder
(`examples/2d_cylinder/2d_cylinder.nmsh`), and compares the end fields:

* `dealias_metrics`: `exact` or `interpolated`
* `dealias_store_metrics`: `true` or `false`
* `dealias_chunk_elements`: `-1` (all at once), `100` (forced chunking),
  `0` (automatic)

Optionally a reference `neko` built from `develop` is run on the same cases.
Fields are written in double precision with tight solver tolerances so that
the comparison is meaningful at rounding level. Only the Python standard
library is needed.

```
python3 contrib/dealias_check/dealias_check.py --neko $PWD/src/neko \
    --ref /path/to/develop/src/neko --launcher "mpirun -n {np}" --np 1
```

Expected outcome (measured on CPU, one rank):

* runs with the same `dealias_metrics` are bitwise identical regardless of
  storage and chunking;
* `exact` and `interpolated` agree to rounding on the TGV box (about 5e-14
  relative); on the cylinder, whose elements are only mildly curved, they
  differ by about 2e-12 in velocity and 1e-7 in pressure at order 7, which
  is the level the solvers amplify rounding to;
* the reference binary agrees with the `interpolated` runs to rounding
  (about 4e-15 on the TGV box, 3e-12 and 1e-8 on the cylinder).

With several MPI ranks the cylinder case shows run-to-run noise of the same
size in develop itself (summation order in the gather-scatter), so use one
rank, which is also all a GPU test needs: every operation in the dealiased
advection is element local and the parallel parts of the solver are
untouched.

The script pins the gather-scatter communication backend to `MPI` through
`NEKO_GS_COMM` unless that variable is already set: the autotuned choice
varies with the machine load, and a different backend sums in a different
order, which the solvers amplify to about 1e-7 in the pressure on these
tiny cases and would hide the differences the script is meant to show.

The per-step times printed are the mean `Fluid step time` over all but the
first step and are indicative only on these tiny cases. Use `--quick` to
skip the forced-chunking runs and `--steps N` to lengthen the runs for
timing.
