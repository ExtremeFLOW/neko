# The Guermond and Shen manufactured solution for the Stokes equations

## Demonstrated concepts
- User source term for the fluid.
- Outputing to a CSV in the user file.
- Working with the `neko_field_registry`.
- Initializing and freeing variables in the user file.
- Reading custom JSON keywords from the case file.
- Basic field arithmetics, global max and min reduction.
- The field_writer simcomp.

## The method of manufactured solutions

The MMS is a verification technique based on the following. Consider some
equation $L(s) = 0$. We come up with some desired analytical solution $\hat
s(x,t)$, and  then analytically compute $\hat f = L(\hat s)$. If we now set
$\hat f$ as the source term: $L(s) = \hat f$, we expect the solution to follow
$\hat s(x, t)$.
Initial conditions are obtained as $\hat s(x, 0)$ and boundary conditions should
be selected to be compatible with the behaviour of $\hat s$.

## The case

This case is known as the Guermond and Shen test for Stokes equations, although
it has been also expanded to the Navier-Stokes, see DOI:
`10.1090/S0025-5718-03-01621-1`.

The simulation domain is a square of size [0, 2]. In Neko, we add a redundant
third dimension with periodic conditions and a single element. No-slip
conditions are applied on all sides of the square. The generated solution is

$$
u = \pi \sin(t) \sin(2\pi y) \sin^2(\pi x) \\
v = -\pi \sin(t) \sin(2\pi x) \sin^2(\pi y) \\
p = \sin(t) \sin(\pi y) \cos(\pi x) \\
$$

Viscosity and density are assumed to be equal to 1.

This case can be used to analyze the error in the solution coming from different
contributions. For example, by increasing the mesh resolution one can relatively
quickly eliminate the spatial discretization error and look at the errors due to
time-stepping and pressure-velocity coupling. See, e.g. DOI
`10.1016/j.jcp.2018.06.062`. Thus, it is a useful tool when prototyping and
verifying the accuracy of solvers.


## Case setup
The mesh can generated with `genmeshbox`, e.g.
`genmeshbox 0 2 0 2 0 0.05 16 16 1  .false. .false. .true`
where the resolution in x and y can be changed as desired.

The user file is configured to compute the errors with respect to the true
solution and write the $L_\infty$ errors to a .csv file. Moreover, the
`field_writer` simcomp is used to output the errors as fields at the end of
the run. The frequency of the output can of course be adjusted.

Advection is turned off in the `fluid` so that the Stokes equations are solved.

Otherwise the case is flexible to the needs of the user and numerical and I/O
parameters can be tweaked as desired.
