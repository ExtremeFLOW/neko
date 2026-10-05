"""slab_table.py: E_r(20), E_v, worst violation over EVERY frame (README definition), svv off vs on;
plus psi statistics for psi_xi10 (the 'Measured, not just argued' table). 2026-10-05.
Needs the SVV-off runs in logs/svv_off_2026-10-02/, which is gitignored and local to the
workstation that ran them; run from examples/advecting_slab_1d/."""
import glob, os, sys, json, logging
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import norms as nm
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import Field
from pysemtools.datatypes.coef import Coef
logging.disable(logging.INFO); comm = MPI.COMM_WORLD
A = "logs/svv_off_2026-10-02/"
def meas(d, eps, stats=False):
    f = sorted(glob.glob(d + "/field0.f?????"))
    msh = Mesh(comm, data=preadnek(f[0], comm)); coef = Coef(msh, comm)
    wv, nv = 0.0, 0
    gmin, gmax, bm = np.inf, 0.0, []
    epin, epout, ephi, wrong = 0.0, 0.0, 0.0, 0
    for ff in f:
        o = Field(comm, data=preadnek(ff, comm)); p = o.fields["temp"][0]
        wv = max(wv, p.max() - 1, -p.min(), 0.0); nv = max(nv, nm.n_unbounded(p))
        if stats:
            ps = o.fields["scal"][0]
            gx = coef.dudxyz(ps, coef.drdx, coef.dsdx, coef.dtdx)
            g = np.abs(gx); b = p*(1 - p) > 1e-4
            gmin = min(gmin, g[b].min()); gmax = max(gmax, g[b].max()); bm.append(g[b].mean())
            pe = nm.slab_distance(msh.x, o.t)
            epin = max(epin, np.abs(ps - pe)[b].max()); epout = max(epout, np.abs(ps - pe)[~b].max())
            ephi = max(ephi, np.abs(p - nm.slab_phase(msh.x, eps, o.t)).max())
            wrong = max(wrong, int((np.sign(gx[b]) != np.sign(np.gradient(pe, axis=-1)[b] if False else -np.sign(msh.x[b] - (0.5 + o.t) % 1.0 + 0*gx[b]))).sum()) if False else 0)
    p0 = Field(comm, data=preadnek(f[0], comm)).fields["temp"][0]
    o = Field(comm, data=preadnek(f[-1], comm)); phi = o.fields["temp"][0]; pe = nm.slab_phase(msh.x, eps, t=o.t)
    r = dict(t=o.t, E_r=nm.E_r(phi, pe, coef, msh), E_v=nm.E_v(phi, p0, pe, coef, msh), wv=wv, nv=nv)
    if stats:
        r.update(gmin=gmin, gmax=gmax, bmin=min(bm), bmax=max(bm), epin=epin, epout=epout, ephi=ephi)
    return r
print("%-14s %5s %5s %4s | %9s %9s | %10s %10s | %6s %6s | %10s" % ("case", "xi", "gamma", "N", "E_r off", "E_r on", "viol off", "viol on", "n off", "n on", "E_v on"))
for c in sys.argv[1:]:
    cc = json.load(open(c + ".case"))["case"]; eps = cc["cdi"]["epsilon"]; g = cc["cdi"]["gamma"]; N = cc["numerics"]["polynomial_order"]
    a = meas(A + cc["output_directory"], eps); b = meas(cc["output_directory"], eps, stats=(c == "psi_xi10"))
    print("%-14s %5.2f %5.2f %4d | %9.5f %9.5f | %10.2e %10.2e | %6d %6d | %10.2e" % (c, eps*N/0.1, g, N, a["E_r"], b["E_r"], a["wv"], b["wv"], a["nv"], b["nv"], b["E_v"]), flush=True)
    if "gmin" in b:
        print("   psi_xi10 svv on: band |dpsi/dx| %.4f - %.3f, band mean %.3f - %.3f, max|psi-psi_ex| band %.5f out %.5f, max|phi-phi_ex| %.5f"
              % (b["gmin"], b["gmax"], b["bmin"], b["bmax"], b["epin"], b["epout"], b["ephi"]), flush=True)
