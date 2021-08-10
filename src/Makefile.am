lib_LIBRARIES = libneko.a
libneko_a_SOURCES = \
	config/num_types.f90\
	structs.f90\
	parameters.f90\
	io/format/nmsh.f90\
	io/format/re2.f90\
	io/format/map.f90\
	log.f90\
	comm.f90\
	mesh/curve.f90\
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
	gs_sx.f90\
	gs_ops.f90\
	gather_scatter.f90\
	mesh/entity.f90\
	mesh/point.f90\
	mesh/element.f90\
	ax.f90\
	precon.f90\
	krylov.f90\
	cg.f90\
	ksp_bicgstab.f90\
	mesh/quad.f90\
	mesh/hex.f90\
	mesh/tet.f90\
	mesh/tri.f90\
	htable.f90\
	uset.f90\
	stack.f90\
	tuple.f90\
	mesh/zone.f90\
	mesh/mesh.f90\
	mesh_field.f90\
	mesh/tet_mesh.f90\
	field.f90\
	io/format/rea.f90\
	sx_cdtp.f90\
	sx_conv1.f90\
	sx_dudxyz.f90\
	sx_opgrad.f90\
	operators.f90\
	opr_cpu.f90\
	opr_sx.f90\
	opr_xsmm.F90\
	checkpoint.f90\
	io/generic_file.f90\
	io/map_file.f90\
	io/re2_file.f90\
	io/rea_file.f90\
	io/fld_file.f90\
	io/vtk_file.f90\
	io/nmsh_file.f90\
	io/chkp_file.f90\
	io/file.f90\
	io/output.f90\
	io/fluid_output.f90\
	io/chkp_output.f90\
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
	config/neko_config.f90\
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
	neko.f90

if ENABLE_PARMETIS
libneko_a_SOURCES += parmetis_wrapper.c
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
