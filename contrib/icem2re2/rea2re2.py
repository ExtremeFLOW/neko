#!/usr/bin/env python3
"""Convert a Nek5000 .rea mesh (as produced by mshconvert) to binary .re2 via pymech."""

import sys
import time
import pymech
from pymech.core import HexaData


def parse_rea(filename):
    print(f"Reading {filename}...")
    t0 = time.time()

    with open(filename, 'r') as f:
        lines = f.readlines()

    nlines = len(lines)
    print(f"  File loaded: {nlines} lines in {time.time()-t0:.1f}s")

    idx = 0

    while idx < nlines:
        if 'DIMENSIONAL' in lines[idx]:
            ndim = int(lines[idx].split()[0])
            break
        idx += 1
    idx += 1
    print(f"  ndim = {ndim}")

    while idx < nlines:
        if 'MESH DATA' in lines[idx]:
            break
        idx += 1
    idx += 1

    parts = lines[idx].split()
    nel = abs(int(float(parts[0])))
    nelv = abs(int(float(parts[2])))
    idx += 1
    print(f"  nel = {nel}, nelv = {nelv}")

    if int(float(parts[0])) < 0:
        print("ERROR: Mesh data not in this .rea (negative NELT).")
        sys.exit(1)

    nfaces = 2 * ndim
    lr1 = [2, 2, 2 if ndim == 3 else 1]
    var = [ndim, 0, 0, 0, 0]

    print(f"  Allocating HexaData for {nel} elements...")
    data = HexaData(ndim, nel, lr1, var, nbc=1)
    data.wdsz = 8
    data.endian = 'little'

    # mshconvert .rea format is face-grouped:
    #   Line 1: x1 x2 x3 x4  (bottom face x)   Line 4: x5 x6 x7 x8  (top face x)
    #   Line 2: y1 y2 y3 y4                    Line 5: y5 y6 y7 y8
    #   Line 3: z1 z2 z3 z4                    Line 6: z5 z6 z7 z8
    # The 4 values per face map to pos[dim, iz, (iy,ix)] as
    #   val[0] -> pos[dim, iz, 0, 0], val[1] -> pos[dim, iz, 0, 1],
    #   val[2] -> pos[dim, iz, 1, 1], val[3] -> pos[dim, iz, 1, 0].

    print(f"  Reading {nel} element coordinates...")
    t2 = time.time()
    for iel in range(nel):
        if iel % 500000 == 0 and iel > 0:
            pct = 100.0 * iel / nel
            elapsed = time.time() - t2
            rate = iel / elapsed
            eta = (nel - iel) / rate / 60.0
            print(f"    ... {iel}/{nel} ({pct:.1f}%), ETA: {eta:.1f} min")

        elem = data.elem[iel]
        idx += 1  # skip element header

        for iz in range(ndim - 1):
            for idim in range(ndim):
                vals = lines[idx].split()
                idx += 1
                elem.pos[idim, iz, 0, 0] = float(vals[0])
                elem.pos[idim, iz, 0, 1] = float(vals[1])
                elem.pos[idim, iz, 1, 1] = float(vals[2])
                elem.pos[idim, iz, 1, 0] = float(vals[3])

    print(f"  Done reading coordinates in {time.time()-t2:.1f}s")

    while idx < nlines:
        if 'CURVED SIDE DATA' in lines[idx]:
            break
        idx += 1
    idx += 1

    ncurve = int(lines[idx].split()[0])
    idx += 1
    print(f"  Curved sides: {ncurve}")

    for _ in range(ncurve):
        parts = lines[idx].split()
        iedge = int(parts[0])
        ielem = int(parts[1]) - 1
        cvals = [float(parts[i]) for i in range(2, 7)]
        ccurve = parts[7] if len(parts) > 7 else ' '
        if 0 <= ielem < nel:
            data.elem[ielem].curv[iedge - 1] = cvals
            data.elem[ielem].ccurv[iedge - 1] = ccurve
        idx += 1

    while idx < nlines:
        if 'FLUID' in lines[idx] and 'BOUNDARY' in lines[idx]:
            break
        idx += 1
    idx += 1

    total_bc = nel * nfaces
    print(f"  Reading {total_bc} boundary conditions...")
    t4 = time.time()

    for ibc in range(total_bc):
        if ibc % 5000000 == 0 and ibc > 0:
            pct = 100.0 * ibc / total_bc
            elapsed = time.time() - t4
            rate = ibc / elapsed
            eta = (total_bc - ibc) / rate / 60.0
            print(f"    ... {ibc}/{total_bc} ({pct:.1f}%), ETA: {eta:.1f} min")

        line = lines[idx]
        idx += 1

        iel = ibc // nfaces
        iface = ibc % nfaces

        stripped = line.strip()
        bc = stripped[0] if stripped else ' '
        rest = stripped[1:].strip() if len(stripped) > 1 else ''
        parts = rest.split()

        data.elem[iel].bcs[0, iface][0] = bc
        data.elem[iel].bcs[0, iface][1] = iel + 1
        data.elem[iel].bcs[0, iface][2] = iface + 1

        if bc in ('E', 'e') and len(parts) >= 3:
            data.elem[iel].bcs[0, iface][3] = float(parts[1])
            data.elem[iel].bcs[0, iface][4] = float(parts[2])
            data.elem[iel].bcs[0, iface][5] = 0.0
            data.elem[iel].bcs[0, iface][6] = 0.0
            data.elem[iel].bcs[0, iface][7] = 0.0
        else:
            for ip in range(5):
                try:
                    data.elem[iel].bcs[0, iface][3 + ip] = float(parts[ip + 1]) if len(parts) > ip + 1 else 0.0
                except (ValueError, IndexError):
                    data.elem[iel].bcs[0, iface][3 + ip] = 0.0

    print(f"  Done reading BCs in {time.time()-t4:.1f}s")
    return data


def rea_to_re2(rea_file, re2_file):
    data = parse_rea(rea_file)
    print(f"\nWriting {re2_file} via pymech...")
    t0 = time.time()
    pymech.writere2(re2_file, data)
    print(f"Done! Output: {re2_file} ({time.time()-t0:.1f}s)")


def main():
    if len(sys.argv) != 3:
        print("Usage: rea2re2.py input.rea output.re2")
        sys.exit(1)
    rea_to_re2(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
