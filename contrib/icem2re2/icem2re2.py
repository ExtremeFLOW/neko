#!/usr/bin/env python3
"""icem2re2 — one-shot ICEM/Fluent .msh -> Neko .re2 converter.

Pipeline: .msh --(mshconvert)--> .rea --(pymech)--> .re2

Usage:
    icem2re2.py input.msh output.re2 --bcs bcs.json \\
        [--periodic periodic.json] [--keep-rea]

bcs.json maps non-periodic Fluent zone ids to Nek boundary-condition letters,
e.g.:
    {"13": "v", "14": "W", "15": "W", "16": "o"}

periodic.json (optional) declares translational-periodic pairs by zone id.
mshconvert handles the periodic linkage natively; do NOT list these zones in
bcs.json:
    [{"zones": [13, 14], "displacement": [0, 0, 0.05]}]
"""

import argparse
import json
import os
import shutil
import sys
import tempfile

# Make sibling modules importable regardless of how the script is launched.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def load_bcs(path):
    with open(path, 'r') as f:
        raw = json.load(f)
    bcs = {}
    for k, v in raw.items():
        # Accept decimal ("13") or hex ("0xd") zone ids — int(_, 0) auto-detects
        # the base from the prefix, so users can paste either form.
        try:
            zone = int(k, 0) if isinstance(k, str) else int(k)
        except (TypeError, ValueError):
            sys.exit(f"Error: zone id '{k}' in {path} is not an integer "
                     "(use decimal like \"13\" or hex like \"0xd\").")
        if not isinstance(v, str):
            sys.exit(f"Error: BC for zone {k} in {path} must be a string, got {type(v).__name__}.")
        v_clean = v.strip()
        if len(v_clean) != 1:
            sys.exit(f"Error: BC for zone {k} in {path} must be a single character, "
                     f"got {v!r}.")
        bcs[zone] = v_clean
    return bcs


def load_periodic(path):
    """Load periodic.json declaring translational-periodic zone pairs.

    JSON shape: a list of objects, each
        {"zones": [zone_a, zone_b], "displacement": [dx, dy, dz]}

    Zone ids accept decimal or hex (same parsing as bcs.json). Returns a
    dict shaped for mshconvert.convert(periodic_dx=...):
        {(zone_a, zone_b): (dx, dy, dz), ...}
    """
    with open(path, 'r') as f:
        raw = json.load(f)
    if not isinstance(raw, list):
        sys.exit(f"Error: {path} must be a JSON list of pair objects.")
    periodic_dx = {}
    for idx, item in enumerate(raw):
        if not isinstance(item, dict):
            sys.exit(f"Error: {path}[{idx}] must be an object with "
                     "'zones' and 'displacement'.")
        zones = item.get("zones")
        disp = item.get("displacement")
        if not isinstance(zones, list) or len(zones) != 2:
            sys.exit(f"Error: {path}[{idx}].zones must be a list of two zone ids.")
        try:
            zone_pair = tuple(int(z, 0) if isinstance(z, str) else int(z)
                              for z in zones)
        except (TypeError, ValueError):
            sys.exit(f"Error: {path}[{idx}].zones must be integers "
                     "(decimal or 0x-prefixed hex).")
        if not isinstance(disp, list) or not all(isinstance(x, (int, float)) for x in disp):
            sys.exit(f"Error: {path}[{idx}].displacement must be a numeric list.")
        periodic_dx[zone_pair] = tuple(float(x) for x in disp)
    return periodic_dx


def main():
    parser = argparse.ArgumentParser(
        prog='icem2re2',
        description='Convert an ICEM/Fluent .msh mesh to a Neko .re2 mesh.')
    parser.add_argument('input_msh', help='Input ANSYS Fluent .msh file (ICEM export).')
    parser.add_argument('output_re2', help='Output Neko .re2 file.')
    parser.add_argument('--bcs', required=True,
                        help='JSON file mapping zone ids to Nek BC letters (e.g. {"13":"v"}).')
    parser.add_argument('--periodic',
                        help='JSON file declaring translational-periodic zone pairs.')
    parser.add_argument('--keep-rea', action='store_true',
                        help='Keep the intermediate .rea file next to the output.')
    args = parser.parse_args()

    if not os.path.isfile(args.input_msh):
        sys.exit(f"Error: input mesh not found: {args.input_msh}")
    if not os.path.isfile(args.bcs):
        sys.exit(f"Error: BC file not found: {args.bcs}")
    if args.periodic and not os.path.isfile(args.periodic):
        sys.exit(f"Error: periodic file not found: {args.periodic}")

    bcs = load_bcs(args.bcs)
    periodic_dx = load_periodic(args.periodic) if args.periodic else {}

    try:
        from mshconvert import convert
    except ImportError as e:
        sys.exit(f"Error importing mshconvert: {e}\n"
                 "Required Python deps: numpy, scipy.")
    try:
        from rea2re2 import rea_to_re2
    except ImportError as e:
        sys.exit(f"Error importing rea2re2: {e}\n"
                 "Required Python deps: pymech.")

    # mshconvert.convert() writes '<basename>.rea' in the current working
    # directory next to the input file. Run it inside a temp dir so we can
    # control where the intermediate ends up and clean up on failure.
    workdir = tempfile.mkdtemp(prefix='icem2re2_')
    msh_abs = os.path.abspath(args.input_msh)
    re2_abs = os.path.abspath(args.output_re2)
    msh_basename = os.path.basename(args.input_msh)
    msh_stem = os.path.splitext(msh_basename)[0]
    msh_link = os.path.join(workdir, msh_basename)
    rea_tmp = os.path.join(workdir, msh_stem + '.rea')

    prev_cwd = os.getcwd()
    try:
        shutil.copy2(msh_abs, msh_link)
        os.chdir(workdir)
        print(f"[1/2] Converting {args.input_msh} -> {msh_stem}.rea")
        if periodic_dx:
            print(f"      Periodic pairs (handled by mshconvert): {periodic_dx}")
        convert(msh_basename, bcs=bcs, periodic_dx=periodic_dx)
        if not os.path.isfile(rea_tmp):
            sys.exit(f"Error: mshconvert did not produce {rea_tmp}.")

        print(f"[2/2] Converting {msh_stem}.rea -> {args.output_re2}")
        rea_to_re2(rea_tmp, re2_abs)

        if args.keep_rea:
            rea_final = os.path.splitext(re2_abs)[0] + '.rea'
            shutil.copy2(rea_tmp, rea_final)
            print(f"Intermediate .rea kept at: {rea_final}")
    finally:
        os.chdir(prev_cwd)
        shutil.rmtree(workdir, ignore_errors=True)

    print(f"\nConversion complete: {args.output_re2}")


if __name__ == '__main__':
    main()
