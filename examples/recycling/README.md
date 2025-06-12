# Turbulent channel case at Re_tau=180 with 1600 elements and N=5 with recycling boundary conditions

Simulation of a turbulent channel in a domain of size
$(8\pi\delta,2\delta,4/3\pi\delta)$  where $\delta$ is the height of half then
channel. The bulk Reynolds number is $U_b=1$ which is enforced through recycling
boundary conditions. $Re_b=\delta U_b/\nu=2800$. The mesh is deformed and
refined closer to the wall in recycling.f90 (the original box.nmsh is of size
(8, 2, 1.5) and centred in (4, 0, 0.75)).

The channel starts with a slightly perturbed initial condition which becomes
turbulent. The transient for the pressure is quite long however. The intent of
this case is to illustrate how the global interpolation is used to set a
recycling boundary condition. This can be seen clearly as the value at x=0 and
at x=10 are the same at all times, with a slight scaling to ensure a bulk
velocity of 1.

This case is similar to the turbulent_channel case, but modified to illustrate
the use of global_interpolation in the recycling.f90 file. The global
interpolation automatically finds the value at some given xyz coordinates
regardless if it local to the calling process or distributed somewhere else
entirely in the domain. It find s the correct owning rank and element and
executes the necessary interpolation.

It can as such be used as a black box when one simply want to know the value of
the field u(x,y,z), but don't know how it is distributed among the different
ranks. x,y,z does not need to be on any collocation point either, but the value
is spectrally interpolated onto the point. 
