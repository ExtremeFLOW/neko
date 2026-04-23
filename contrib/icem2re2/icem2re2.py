#!/usr/bin/env python3
"""icem2re2 — one-shot ICEM/Fluent .msh -> Neko .re2 converter.

Pipeline: .msh --(mshconvert)--> .rea --(pymech)--> .re2

Usage:
    icem2re2.py input.msh output.re2 --bcs bcs.json [--keep-rea]

bcs.json maps Fluent zone ids to Nek boundary-condition letters, e.g.:
    {
        "13": "v",
        "14": "W",
        "15": "W",
        "16": "o"
    }
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
        try:
            bcs[int(k)] = v
        except ValueError:
            sys.exit(f"Error: zone id '{k}' in {path} is not an integer.")
    return bcs


def main():
    parser = argparse.ArgumentParser(
        prog='icem2re2',
        description='Convert an ICEM/Fluent .msh mesh to a Neko .re2 mesh.')
    parser.add_argument('input_msh', help='Input ANSYS Fluent .msh file (ICEM export).')
    parser.add_argument('output_re2', help='Output Neko .re2 file.')
    parser.add_argument('--bcs', required=True,
                        help='JSON file mapping zone ids to Nek BC letters (e.g. {"13":"v"}).')
    parser.add_argument('--keep-rea', action='store_true',
                        help='Keep the intermediate .rea file next to the output.')
    args = parser.parse_args()

    if not os.path.isfile(args.input_msh):
        sys.exit(f"Error: input mesh not found: {args.input_msh}")
    if not os.path.isfile(args.bcs):
        sys.exit(f"Error: BC file not found: {args.bcs}")

    bcs = load_bcs(args.bcs)

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
        convert(msh_basename, bcs=bcs)
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
