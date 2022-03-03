# Case File {#case-file}

The case file defines all the parameters of a simulation. It uses the plain-text
Fortran Namelist format and usually has a `.case` extension. A case file contains
two sections: `NEKO_CASE` defines the simulation case, and `NEKO_PARAMETERS` defines
the solver parameters.

## Example

This example is available in the source repository as file `examples/hemi/hemi.case`:

~~~
&NEKO_CASE
mesh_file= 'hemi.nmsh'
fluid_scheme='plan4'
lx = 6
source_term = 'noforce'
initial_condition = 'uniform'
/
&NEKO_PARAMETERS
dt = 1d-3
T_end = 2.0
nsamples = 100
uinf= 1.0,0.0,0.0
output_bdry = .true.
rho = 1
Re = 1400
abstol_vel = 1d-9
abstol_prs = 1d-9
ksp_vel = 'cg'
ksp_prs = 'gmres'
pc_vel = 'jacobi'
pc_prs = 'hsmg'
/
~~~

## Syntax

A section is started by a `&SECTION` line and ends with a `/` line. Inbetween
is a list of `name = value` assignments. The names are not case sensistive. The
values can be of type string (``'text'``), integer (`123`), floating point number (`1.2d3` \f$=1.2 \cdot 10^3\f$), or boolean (`.true.` or `.false.`).
Comments can be added by prepending a `!`.

## `NEKO_CASE` reference

Name                    | Description
----                    | -----------
`mesh_file`             | Path to the mesh file (`.nmsh` extension)
`fluid_scheme`          | Solver scheme (``'plan4'`` or ``'pnpn'``)
`lx`                    | Number of quadrature points points per element and direction, i.e. polynomial degree + 1
`source_term`           | Source term \f$ f \f$: default ``'noforce'`` for \f$ f=0 \f$, or ``'user'`` for user defined function
`initial_condition`     | Initial condition: ``'uniform'``, ``'blasius'`` or ``'user'``

## `NEKO_PARAMETERS` reference

Name                    | Description                                                   | Default value
----                    | -----------                                                   | -------------
`nsamples`              | Number of output samples                                      | `0`
`output_bdry`           | Output boundary markings                                      | `.false.`
`output_part`           | Output partitions                                             | `.false.`
`output_chkp`           | Output checkpoints                                            | `.false.`
`dt`                    | Time-step size                                                | -
`T_end`                 | Final time                                                    | -
`rho`                   | Density \f$ \rho \f$                                          | `1d0`
`mu`                    | Dynamic viscosity \f$ \mu \f$                                 | `1d0`
`Re`                    | Reynolds number                                               | `1d0`
`uinf`                  | Free-stream velocity \f$ u_\infty \f$                         | `0.0, 0.0, 0.0`
`abstol_vel`            | Tolerance for velocity solver                                 | `1d-9`
`abstol_prs`            | Tolerance for pressure solver                                 | `1d-9`
`ksp_vel`               | Krylov solver for velocity (``'cg'``, ``'gmres'``)            | ``'cg'``
`ksp_prs`               | Krylov solver for pressure (``'cg'``, ``'gmres'``)            | ``'gmres'``
`pc_vel`                | Preconditioner for velocity solver (`jacobi`, `hsmg`)         | ``'jacobi'``
`pc_prs`                | Preconditioner for pressure solver (`jacobi`, `hsmg`)                 | ``'hsmg'``
`fluid_inflow`          | Fluid inflow condition (``'default'``, ``'blasius'``, ``'user'``)     | ``'default'``
`vol_flow_dir`          | Direction of forced volume flow (none=0, x=1, y=2, z=3)               | `0`
`avflow`                | Use averaged flow for forced volume flow                              | `.true.`
`loadb`                 | Use load-balancing                                                    | `.false.`
`flow_rate`             | Volume flow speed                                                     | `0d0`
`proj_dim`              | Projection space for pressure solution                                | `20`
`time_order`            | Order of the time stepping                                            | `3`
`jlimit`                | Job limit in ``'HH:MM:SS'``                                           | -
`restart_file`          | Checkpoint filename                                                   | -
`stats_begin`           | Start time for statistics                                                                                     | `0d0`
`stats_mean_flow`       | Enable mean flow statistics                                                                                   | `.false.`
`output_mean_flow`      | Output mean flow field                                                                                        | `.false.`
`stats_mean_sqr_flow`   | Enable squared flow statistics                                                                                | `.false.`
`output_mean_sqr_flow`  | Output mean squared flow field                                                                                | `.false.`
`output_dir`            | Output directory                                                                                              | Current working directory
`dealias`               | Enable dealiasing                                                                                             | `.true.`
`lxd`                   | Size of dealiased space                                                                                       | `0`
`delta`                 | Boundary layer thickness \f$ \delta \f$                                                                       | `1d0`
`blasius_approx`        | Type of approximate Blasius profile (``'linear'``, ``'quadratic'``, ``'cubic'``, ``'quartic'``, ``'sin'``)    | ``'sin'``
`bc_labels`             | Type of boundary condition for each label                      | -