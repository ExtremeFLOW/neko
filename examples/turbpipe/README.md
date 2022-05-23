Pipeflow at Re=5300 and Re_tau = 1800.
To perform a DNS of the case lx=8 should be used.
We have two case files, one called pipe.case using "vol_flow" to drive the flow.
We also have one called pipe_device.case that instead uses a source term.
We use quite forgiving tolerances, and this case should be possible to run effieinctly consumer GPUs.
