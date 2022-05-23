Pipeflow at Re_b=5300 and Re_tau = 1800.
To perform a DNS of the case lx=8 should be used.
We have two case files, one called pipe.case using "vol_flow" to drive the flow, this is currently not supported for GPUs.
We also have one called pipe_device.case that instead uses a source term, this can be run efficiently on GPUs.
We use quite forgiving tolerances, and this case should be possible to run effieinctly on consumer GPUs in single precision, specify '--enable-real=sp' when configuring.
