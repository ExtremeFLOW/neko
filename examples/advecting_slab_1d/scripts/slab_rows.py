"""slab_rows.py: every-frame E_r(20), E_v, worst violation, violating nodes, and psi statistics for
named output dirs (eps 0.01). 2026-10-05."""
import glob, os, sys, logging
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import norms as nm
from mpi4py import MPI
from pysemtools.io.ppymech.neksuite import preadnek
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import Field
from pysemtools.datatypes.coef import Coef
logging.disable(logging.INFO); comm = MPI.COMM_WORLD
tag, d = sys.argv[1], sys.argv[2]; eps = 0.01
f = sorted(glob.glob(d + "/field0.f?????"))
msh = Mesh(comm, data=preadnek(f[0], comm)); coef = Coef(msh, comm)
wv, nv, gmin, gmax, ein = 0.0, 0, np.inf, 0.0, 0.0
for ff in f:
    o = Field(comm, data=preadnek(ff, comm)); p, ps = o.fields["temp"][0], o.fields["scal"][0]
    wv = max(wv, p.max() - 1, -p.min(), 0.0); nv = max(nv, nm.n_unbounded(p))
    g = np.abs(coef.dudxyz(ps, coef.drdx, coef.dsdx, coef.dtdx)); b = p*(1 - p) > 1e-4
    gmin = min(gmin, g[b].min()); gmax = max(gmax, g[b].max())
    ein = max(ein, np.abs(ps - nm.slab_distance(msh.x, o.t))[b].max())
p0 = Field(comm, data=preadnek(f[0], comm)).fields["temp"][0]
o = Field(comm, data=preadnek(f[-1], comm)); phi = o.fields["temp"][0]; pe = nm.slab_phase(msh.x, eps, t=o.t)
print("%-28s t=%.2f E_r=%.5f E_v=%.2e viol=%.2e nviol=%d band|dpsi/dx| %.3f-%.3f max|psi-psi_ex|band %.5f"
      % (tag, o.t, nm.E_r(phi, pe, coef, msh), nm.E_v(phi, p0, pe, coef, msh), wv, nv, gmin, gmax, ein), flush=True)
