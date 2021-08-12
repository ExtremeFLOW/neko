include $(top_srcdir)/hip.mk
include $(top_srcdir)/cuda.mk

lib_LIBRARIES = libneko.a
libneko_a_SOURCES = \
	config/num_types.f90\
	common/structs.f90\
	common/parameters.f90\
	io/format/nmsh.f90\
	io/format/re2.f90\
	io/format/map.f90\
	common/log.f90\
	comm/comm.f90\
	mesh/curve.f90\
	comm/mpi_types.f90\
	common/datadist.f90\
	common/distdata.f90\
	common/utils.f90\
	math/mxm_wrapper.F90\
	math/bcknd/cpu/mxm_std.f90\
	sem/speclib.f90 \
	math/math.f90\
	math/mathops.f90\
	math/fast3d.f90\
	math/tensor.f90\
	sem/space.f90\
	sem/dofmap.f90\
	sem/coef.f90\
	gs/gs_bcknd.f90\
	gs/bcknd/cpu/gs_cpu.f90\
	gs/bcknd/sx/gs_sx.f90\
	gs/bcknd/device/gs_device.F90\
	gs/gs_ops.f90\
	gs/gather_scatter.f90\
	mesh/entity.f90\
	mesh/point.f90\
	mesh/element.f90\
	math/ax.f90\
	krylov/precon.f90\
	krylov/krylov.f90\
	krylov/cg.f90\
	krylov/bicgstab.f90\
	krylov/bcknd/device/cg_device.f90\
	mesh/quad.f90\
	mesh/hex.f90\
	mesh/tet.f90\
	mesh/tri.f90\
	adt/htable.f90\
	adt/uset.f90\
	adt/stack.f90\
	adt/tuple.f90\
	mesh/zone.f90\
	mesh/mesh.f90\
	field/mesh_field.f90\
	mesh/tet_mesh.f90\
	field/field.f90\
	io/format/rea.f90\
	math/bcknd/sx/sx_cdtp.f90\
	math/bcknd/sx/sx_conv1.f90\
	math/bcknd/sx/sx_dudxyz.f90\
	math/bcknd/sx/sx_opgrad.f90\
	math/operators.f90\
	math/bcknd/cpu/opr_cpu.f90\
	math/bcknd/sx/opr_sx.f90\
	math/bcknd/xsmm/opr_xsmm.F90\
	math/bcknd/device/opr_device.F90\
	common/checkpoint.f90\
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
	common/sampler.f90\
	bc/bc.f90\
	bc/dirichlet.f90\
	bc/wall.f90\
	bc/inflow.f90\
	bc/usr_inflow.f90\
	bc/facet_normal.f90\
	bc/symmetry.f90\
	krylov/gmres.f90\
	krylov/bcknd/cpu/pc_jacobi.f90\
	krylov/pc_identity.f90\
	krylov/bcknd/sx/pc_jacobi_sx.f90\
	common/source.f90\
	common/abbdf.f90\
	config/neko_config.f90\
	case.f90\
	common/user_intf.f90\
	fluid/fluid_method.f90\
	fluid/fluid_plan4.f90\
	fluid/fluid_plan1.f90\
	fluid/fluid_schemes.f90\
	simulation.f90\
	math/bcknd/cpu/ax_helm.f90\
	math/bcknd/sx/ax_helm_sx.f90\
	math/bcknd/xsmm/ax_helm_xsmm.F90\
	common/projection.f90\
	comm/parmetis.F90\
	comm/redist.f90\
	krylov/pc_hsmg.f90\
	math/fdm.f90\
	math/schwarz.f90\
	common/signal.f90\
	common/sighdl.c\
	common/jobctrl.f90\
	device/cuda_intf.F90\
	device/hip_intf.F90\
	device/device.F90\
	math/bcknd/device/device_math.F90\
	bc/bcknd/device/device_dirichlet.F90\
	bc/bcknd/device/device_inflow.F90\
	bc/bcknd/device/device_wall.F90\
	bc/bcknd/device/device_symmetry.F90\
	bc/bcknd/device/device_facet_normal.F90\
	neko.f90

if ENABLE_PARMETIS
libneko_a_SOURCES += comm/parmetis_wrapper.c
endif

if ENABLE_HIP
libneko_a_SOURCES += \
	device/hip/hip_wrapper.hip\
	math/bcknd/device/hip/math.hip\
	math/bcknd/device/hip/opr_dudxyz.hip\
	gs/bcknd/device/hip/gs_kernels.hip\
	bc/bcknd/device/hip/dirichlet_kernel.hip\
	bc/bcknd/device/hip/inflow_kernel.hip\
	bc/bcknd/device/hip/no_slip_wall_kernel.hip\
	bc/bcknd/device/hip/symmetry_kernel.hip\
	bc/bcknd/device/hip/facet_normal_kernel.hip
endif

if ENABLE_CUDA
libneko_a_SOURCES += device/cuda/cuda_wrapper.cu
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
