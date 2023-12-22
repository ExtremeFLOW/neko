Pipeflow at Re_b=5300 and Re_tau = 1800.
To perform a DNS of the case lx=8 should be used.
We have two case files, one called pipe.case using "flow_rate_force" to drive the flow,  and one called pipe_source.case that instead uses a constant source term.
The source term values is based on computing the shear stress as 2 * (2 * Re_tau/Re_B)^2, which assumes that the bulk velocity is 1. 



We use quite forgiving tolerances, and this case should be possible to run efficiently on consumer GPUs in single precision, specify '--enable-real=sp' when configuring.
