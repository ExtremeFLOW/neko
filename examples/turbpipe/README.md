Pipeflow at Re_b=5300 and Re_tau = 1800.
To perform a DNS of the case lx=8 should be used.
We have two case files, one called pipe.case using "vol_flow" to drive the flow,  and one called pipe_source.case that instead uses a source term. Please note that there is currently no support for "vol_flow" on GPUs when using the 'plan4' solver. 
We use quite forgiving tolerances, and this case should be possible to run effieinctly on consumer GPUs in single precision, specify '--enable-real=sp' when configuring.
