#!/usr/bin/env python3
"""train_final.msh (mm, hex8) -> hex27 curved -> metres -> gmsh2nek input.

Order matters: curve_to_cad snaps wall nodes onto the analytic wall(), which
is in mm -- so curving runs BEFORE scaling. Only wall groups (3,4,5) are
snapped; inlet/outlet are planar and stay put.
"""
import ogrid_cvd as og

og.curve_to_cad("train_final.msh", "train_o2.msh")
print("curved: train_o2.msh (hex27, walls on the analytic geometry, mm)")


def scale_msh(fin, fout, s=0.001):
    with open(fin) as f:
        lines = f.readlines()
    out, i = [], 0
    while i < len(lines):
        out.append(lines[i])
        if lines[i].startswith("$Nodes"):
            n = int(lines[i + 1]); out.append(lines[i + 1])
            for j in range(i + 2, i + 2 + n):
                p = lines[j].split()
                out.append("%s %.10g %.10g %.10g\n"
                           % (p[0], float(p[1])*s, float(p[2])*s, float(p[3])*s))
            i += 2 + n
            continue
        i += 1
    with open(fout, "w") as f:
        f.writelines(out)


scale_msh("train_o2.msh", "train_o2_m.msh")
print("scaled: train_o2_m.msh (metres)")
