include $(top_srcdir)/hip.mk
include $(top_srcdir)/cuda.mk

lib_LIBRARIES = libneko.a
libneko_a_SOURCES = \
	num_types.f90\
        structs.f90\
	parameters.f90\
	nmsh.f90\
	re2.f90\
	map.f90\
	log.f90\
	comm.f90\
        curve.f90\
	mpi_types.f90 \
	datadist.f90 \
	distdata.f90 \
	utils.f90\
	mxm_wrapper.F90\
	mxm_std.f90\
	speclib.f90 \
	math.f90\
	mathops.f90\
	fast3d.f90\
	tensor.f90\
	space.f90\
	dofmap.f90\
	coef.f90\
	gs_bcknd.f90\
	gs_cpu.f90\
	gs_hip.f90\
	gs_sx.f90\
	gs_ops.f90\
	gather_scatter.f90\
	entity.f90\
	point.f90\
	element.f90\
	ax.f90\
	precon.f90\
	krylov.f90\
	cg.f90\
	cg_hip.f90\
        ksp_bicgstab.f90\
	quad.f90\
	hex.f90\
	tet.f90\
	tri.f90\
	htable.f90\
	uset.f90\
	stack.f90\
	tuple.f90\
	zone.f90\
	mesh.f90\
	mesh_field.f90\
	tet_mesh.f90\
	field.f90\
	rea.f90\
	sx_cdtp.f90\
	sx_conv1.f90\
	sx_dudxyz.f90\
	sx_opgrad.f90\
	operators.f90\
	opr_cpu.f90\
	opr_sx.f90\
	opr_xsmm.F90\
	opr_device.F90\
	generic_file.f90\
	map_file.f90\
	re2_file.f90\
	rea_file.f90\
	fld_file.f90\
	vtk_file.f90\
	nmsh_file.f90\
	file.f90\
	output.f90\
	fluid_output.f90\
	sampler.f90\
	bc.f90\
	dirichlet.f90\
	wall.f90\
	inflow.f90\
	usr_inflow.f90\
	facet_normal.f90\
	symmetry.f90\
	gmres.f90\
	pc_jacobi.f90\
	pc_identity.f90\
	pc_jacobi_sx.f90\
	source.f90\
	abbdf.f90\
	neko_config.f90\
	case.f90\
	user_intf.f90\
	fluid_method.f90\
	fluid_plan4.f90\
	fluid_plan1.f90\
	fluid_schemes.f90\
	simulation.f90\
	ax_helm.f90\
	ax_helm_sx.f90\
	ax_helm_xsmm.F90\
	projection.f90\
	parmetis.F90\
	redist.f90\
	pc_hsmg.f90\
	fdm.f90\
	schwarz.f90\
	signal.f90\
	sighdl.c\
	jobctrl.f90\
	cuda_intf.F90\
	hip_intf.F90\
	device.F90\
	device_math.F90\
	device_dirichlet.F90\
	device_inflow.F90\
	device_wall.F90\
	device_symmetry.F90\
	device_facet_normal.F90\
	neko.f90

if ENABLE_PARMETIS
libneko_a_SOURCES += parmetis_wrapper.c
endif

if ENABLE_HIP
libneko_a_SOURCES += \
	hip/hip_wrapper.hip\
	hip/math.hip\
	hip/opr_dudxyz.hip\
	hip/gs_kernels.hip\
	hip/dirichlet_kernel.hip\
	hip/inflow_kernel.hip\
	hip/no_slip_wall_kernel.hip\
	hip/symmetry_kernel.hip\
	hip/facet_normal_kernel.hip
endif

if ENABLE_CUDA
libneko_a_SOURCES += cuda/cuda_wrapper.hip
endif

bin_PROGRAMS = neko
neko_SOURCES = $(libneko_a_SOURCES) driver.f90

pkginclude_HEADERS =  $(wildcard *.mod)

CLEANFILES =*.mod

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(libneko_a_SOURCES) driver.f90  > .depends
endif

include	.depends
