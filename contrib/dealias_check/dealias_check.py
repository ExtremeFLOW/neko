#!/usr/bin/env python3
"""Consistency and timing check for the dealiased advection.

Runs two small cases (the affine TGV box and the curved 2D cylinder) through
the dealiased advection in every combination of

  * dealias_metrics        : exact | interpolated
  * dealias_store_metrics  : true | false
  * dealias_chunk_elements : -1 (all at once) | 100 (forced chunking) | 0 (auto)

with double precision output and tight solver tolerances, and compares the
end fields of every run against the run that is algorithmically closest to
the previous implementation (interpolated metrics, stored, all at once), and
optionally against a reference `neko` binary built from develop.

Expected outcome (measured on CPU, one rank)
  * runs sharing the same `dealias_metrics` are bitwise identical whatever
    the storage or chunking: these only reorganise the same arithmetic;
  * `exact` vs `interpolated` agree to rounding on the TGV box (affine
    elements, ~5e-14 relative). On the cylinder, whose elements are only
    mildly curved, they differ by ~2e-12 in velocity and ~1e-7 in pressure
    at order 7, i.e. at the level the solvers amplify rounding to;
  * the reference binary agrees with the `interpolated` runs to rounding
    (~4e-15 on the TGV box, ~3e-12 / 1e-8 on the cylinder);
  * with several ranks the cylinder case carries run-to-run noise of the
    same ~1e-12 / 1e-7 size in develop as well (summation order in the
    gather-scatter), so use one rank for bitwise comparisons;
  * in single precision (--sp) only the bitwise identity within a metrics
    mode is a sharp test; differences between modes or against the
    reference sit at the ~1e-6 (velocity) to ~1e-4 (pressure) level the
    solvers amplify single-precision rounding to.

Only the standard library is needed. Typical use:

  python3 contrib/dealias_check/dealias_check.py --neko $PWD/src/neko \
      --ref /path/to/develop/src/neko --launcher "mpirun -n {np}" --np 1

See `--help` for the options.
"""
import argparse
import array
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

METRICS = ["exact", "interpolated"]
STORE = [True, False]
CHUNKS = [-1, 100, 0]
BASELINE = ("interpolated", True, -1)


# Solver tolerances: tight enough for rounding-level comparisons in double
# precision; in single precision an absolute 1e-10 is unreachable, the
# solvers would run to their iteration limit, and both timings and
# field differences would be meaningless, so --sp relaxes them.
TOL = {"dp": (1e-11, 1e-10), "sp": (1e-6, 1e-5)}


def tgv_case(steps, dt, order, tol=TOL["dp"]):
    return {
        "version": 1.0,
        "case": {
            "mesh_file": "512.nmsh",
            "output_boundary": False,
            "output_checkpoints": False,
            "output_at_end": True,
            "output_precision": "double",
            "time": {"end_time": steps * dt, "timestep": dt},
            "numerics": {"time_order": 3, "polynomial_order": order,
                         "dealias": True},
            "fluid": {
                "scheme": "pnpn",
                "Re": 1600.0,
                "initial_condition": {
                    "type": "expression",
                    "value": ["sin(x)*cos(y)*cos(z)",
                              "-cos(x)*sin(y)*cos(z)", "0"]},
                "velocity_solver": {
                    "type": "cg", "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tol[0], "max_iterations": 800},
                "pressure_solver": {
                    "type": "gmres", "preconditioner": {"type": "hsmg"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tol[1], "max_iterations": 800},
                "output_control": "never",
            },
        },
    }


def cyl_case(steps, dt, order, tol=TOL["dp"]):
    return {
        "version": 1.0,
        "case": {
            "mesh_file": "2d_cylinder.nmsh",
            "output_boundary": False,
            "output_checkpoints": False,
            "output_at_end": True,
            "output_precision": "double",
            "time": {"end_time": steps * dt, "timestep": dt},
            "numerics": {"time_order": 3, "polynomial_order": order,
                         "dealias": True},
            "fluid": {
                "scheme": "pnpn",
                "Re": 160.0,
                "initial_condition": {"type": "uniform",
                                      "value": [1.0, 0.0, 0.0]},
                "velocity_solver": {
                    "type": "cg", "preconditioner": {"type": "jacobi"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tol[0], "max_iterations": 800},
                "pressure_solver": {
                    "type": "gmres", "preconditioner": {"type": "hsmg"},
                    "projection_space_size": 0,
                    "absolute_tolerance": tol[1], "max_iterations": 800},
                "boundary_conditions": [
                    {"type": "no_slip", "zone_indices": [2]},
                    {"type": "velocity_value", "zone_indices": [3],
                     "value": [1.0, 0.0, 0.0]},
                    {"type": "outflow", "zone_indices": [4]}],
                "output_control": "never",
            },
        },
    }


CASES = {
    "tgv": (tgv_case, REPO / "examples/tgv/512.nmsh", 5, 1e-2),
    "cyl": (cyl_case, REPO / "examples/2d_cylinder/2d_cylinder.nmsh", 20,
            1e-3),
}


def read_fld(fn):
    """Read a Neko fld file written by a single-file output (any rank
    count) and return (header, idx, {'uvw':..., 'p':...}) as flat arrays
    ordered by global element index."""
    with open(fn, "rb") as f:
        hdr = f.read(132).decode()
        parts = hdr.split()
        wdsz, lx, ly, lz, nelv = (int(parts[1]), int(parts[2]),
                                  int(parts[3]), int(parts[4]), int(parts[5]))
        rdcode = parts[11]
        f.read(4)  # endian test pattern
        idx = array.array("i")
        idx.frombytes(f.read(4 * nelv))
        code = "d" if wdsz == 8 else "f"
        lxyz = lx * ly * lz
        out = {}

        def take(n):
            a = array.array(code)
            a.frombytes(f.read(wdsz * n))
            return a

        if "X" in rdcode:
            out["xyz"] = take(3 * lxyz * nelv)
        if "U" in rdcode:
            out["uvw"] = take(3 * lxyz * nelv)
        if "P" in rdcode:
            out["p"] = take(lxyz * nelv)
    # Reorder by global element id so that runs with different rank counts
    # (hence element orders) are comparable
    order = sorted(range(nelv), key=lambda e: idx[e])
    for k, per in (("uvw", 3 * lxyz), ("p", lxyz)):
        if k in out:
            a = out[k]
            b = array.array(code, bytes(len(a) * a.itemsize))
            for new_e, old_e in enumerate(order):
                b[new_e * per:(new_e + 1) * per] = a[old_e * per:(old_e + 1) * per]
            out[k] = b
    return hdr, out


def compare(a, b):
    res = {}
    for k in ("uvw", "p"):
        d = max(abs(x - y) for x, y in zip(a[k], b[k]))
        ref = max(abs(x) for x in a[k]) or 1.0
        res[k] = d / ref
    return res


def parse_log(log):
    txt = Path(log).read_text(errors="replace")
    m = re.findall(r"Fluid step time \(s\):\s+([0-9.E+-]+)", txt)
    steps = [float(x) for x in m]
    # skip the first step (warm-up, autotuning of kernels)
    per_step = sum(steps[1:]) / len(steps[1:]) if len(steps) > 1 else float("nan")
    tot = re.findall(r"Total elapsed time \(s\):\s+([0-9.E+-]+)", txt)
    total = float(tot[-1]) if tot else float("nan")
    deal = re.findall(r"Dealiasing\s*:\s*(.*)", txt)
    deal = "; ".join(d.strip() for d in deal) if deal else "(no dealiasing line)"
    err = "ERROR" in txt or "Error" in txt
    single = "single precision" in txt
    return per_step, total, deal, err, single


def run_one(neko, launcher, np_, workdir, name, case, mesh, env):
    d = Path(workdir) / name
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.f0*"):
        f.unlink()
    for f in d.glob("*.nek5000"):
        f.unlink()
    link = d / mesh.name
    if not link.exists():
        os.symlink(mesh, link)
    (d / "run.case").write_text(json.dumps(case, indent=2))
    cmd = launcher.format(np=np_).split() + [str(neko), "run.case"]
    t0 = time.time()
    with open(d / "log.txt", "w") as log:
        rc = subprocess.call(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT,
                             env=env)
    wall = time.time() - t0
    per_step, total, deal, err, single = parse_log(d / "log.txt")
    fld = d / "field0.f00000"
    ok = (rc == 0) and fld.exists() and not err
    return {"name": name, "dir": d, "rc": rc, "ok": ok, "wall": wall,
            "per_step": per_step, "total": total, "deal": deal,
            "single": single, "fld": fld if ok else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--neko", required=True, help="neko binary of this branch")
    ap.add_argument("--ref", help="reference neko binary (develop), optional")
    ap.add_argument("--launcher", default="mpirun -n {np}",
                    help="launcher prefix, {np} is replaced by --np")
    ap.add_argument("--np", type=int, default=1)
    ap.add_argument("--cases", default="tgv,cyl")
    ap.add_argument("--order", type=int, default=7)
    ap.add_argument("--steps", type=int, default=0,
                    help="override the number of time steps (0 = per case)")
    ap.add_argument("--workdir", default="dealias_check_out")
    ap.add_argument("--quick", action="store_true",
                    help="only exact/interpolated x stored/rebuilt, chunk -1 and 0")
    ap.add_argument("--sp", action="store_true",
                    help="binaries built in single precision: relax the solver "
                         "tolerances accordingly")
    args = ap.parse_args()

    env = dict(os.environ)
    # Pin the gather-scatter communication backend unless the caller chose
    # one: the autotuned choice differs from run to run with the machine
    # load, and a different backend sums in a different order, which the
    # solvers amplify to ~1e-7 in the pressure on these tiny cases. That
    # noise would hide the differences this script is meant to show.
    env.setdefault("NEKO_GS_COMM", "MPI")
    print(f"gather-scatter backend pinned to NEKO_GS_COMM={env['NEKO_GS_COMM']}")
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    chunks = [-1, 0] if args.quick else CHUNKS
    tol = TOL["sp" if args.sp else "dp"]
    summary = []
    warned_sp = False

    for cname in args.cases.split(","):
        make_case, mesh, steps, dt = CASES[cname]
        if args.steps > 0:
            steps = args.steps
        if not mesh.exists():
            print(f"mesh {mesh} not found, skipping {cname}")
            continue
        runs = {}
        print(f"\n=== case {cname}: {steps} steps of {dt}, order {args.order}")
        for metrics in METRICS:
            for store in STORE:
                for chunk in chunks:
                    case = make_case(steps, dt, args.order, tol)
                    num = case["case"]["numerics"]
                    num["dealias_metrics"] = metrics
                    num["dealias_store_metrics"] = store
                    num["dealias_chunk_elements"] = chunk
                    name = f"{cname}_{metrics}_{'store' if store else 'rebuild'}_chunk{chunk}"
                    r = run_one(args.neko, args.launcher, args.np, workdir,
                                name, case, mesh, env)
                    runs[(metrics, store, chunk)] = r
                    print(f"  {name:48s} rc={r['rc']} step={r['per_step']:.4f}s"
                          f"  {r['deal']}")
                    if not args.sp and not warned_sp and r["single"]:
                        warned_sp = True
                        print("  WARNING: this neko is built in single precision;"
                              " rerun with --sp, otherwise the tolerances are"
                              " unreachable and timings and differences are"
                              " meaningless")
        ref = None
        if args.ref:
            case = make_case(steps, dt, args.order, tol)
            ref = run_one(args.ref, args.launcher, args.np, workdir,
                          f"{cname}_reference", case, mesh, env)
            print(f"  {ref['name']:48s} rc={ref['rc']} step={ref['per_step']:.4f}s")

        base = runs.get(BASELINE)
        if base is None or not base["ok"]:
            # fall back to any successful run with interpolated metrics
            base = next((r for k, r in runs.items() if r["ok"] and k[0] == "interpolated"), None)
        print(f"\n  {'run':48s} {'uvw vs base':>12s} {'p vs base':>12s}"
              f" {'uvw vs ref':>12s} {'p vs ref':>12s} {'s/step':>8s}")
        base_f = read_fld(base["fld"])[1] if base and base["ok"] else None
        ref_f = read_fld(ref["fld"])[1] if ref and ref["ok"] else None
        for key, r in runs.items():
            if not r["ok"]:
                line = f"  {r['name']:48s} FAILED (rc={r['rc']}), see {r['dir']/'log.txt'}"
                print(line)
                summary.append(line)
                continue
            f = read_fld(r["fld"])[1]
            cb = compare(base_f, f) if base_f else {"uvw": float("nan"), "p": float("nan")}
            cr = compare(ref_f, f) if ref_f else {"uvw": float("nan"), "p": float("nan")}
            line = (f"  {r['name']:48s} {cb['uvw']:12.3e} {cb['p']:12.3e}"
                    f" {cr['uvw']:12.3e} {cr['p']:12.3e} {r['per_step']:8.4f}")
            print(line)
            summary.append(line)
        if ref and ref["ok"]:
            line = f"  {ref['name']:48s} {'-':>12s} {'-':>12s} {'-':>12s} {'-':>12s} {ref['per_step']:8.4f}"
            print(line)
            summary.append(line)

    (workdir / "summary.txt").write_text("\n".join(summary) + "\n")
    print(f"\nsummary written to {workdir/'summary.txt'}")
    print("baseline for 'vs base' is interpolated metrics, stored, all elements"
          " at once (closest to the previous implementation)")


if __name__ == "__main__":
    main()
