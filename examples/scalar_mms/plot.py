import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import pysemtools
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.field import Field, FieldRegistry
from pysemtools.io.ppymech.neksuite import pynekread
from pysemtools.interpolation.probes import Probes
from mpi4py import MPI

DATA_PATH = "."

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#%%
mesh = Mesh(comm, create_connectivity=True)
fld = FieldRegistry(comm)

pynekread(join(DATA_PATH, "field0.f00000"), comm, data_dtype=np.single, msh=mesh)
pynekread(join(DATA_PATH, "field0.f00001"), comm, data_dtype=np.single, fld=fld)

#%%
n = 500
x = np.linspace(0, 2 * np.pi, n)
y = np.zeros(n)
z = np.ones(n) * 0.5

if rank == 0:
    xyz = np.column_stack((x[:, None], y[:, None], z[:, None]))
else:
    xyz = None

#%%

probes = Probes(comm, probes = xyz, msh = mesh, 
                point_interpolator_type='multiple_point_legendre_numpy',
                max_pts=1000, find_points_comm_pattern='point_to_point')

probes.interpolate_from_field_list(0, [fld.registry['t']], comm, 
                                   write_data=False)

temp = probes.interpolated_fields[:, 1]
#%%
ic = np.sin(x)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(x, ic, '--k', label="ICs")
ax1.plot(x, temp, label=r"Final solution")

ax1.set_xlabel(r"$x$")
ax1.legend()
ax1.set_xlim(0, np.pi * 2)



ax2.plot(x, ic - temp, '--k', label="Absolute error")
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"absolute error")

plt.tight_layout()

print(f"Absolute error, Linf norm: {np.max(ic - temp)}")
plt.savefig(join(DATA_PATH, "results.png"), dpi=300)
