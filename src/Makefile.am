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
	comm/comm.F90\
	mesh/curve.f90\
	comm/mpi_types.f90\
	common/datadist.f90\
	common/distdata.f90\
	common/utils.f90\
	math/mxm_wrapper.F90\
	sem/speclib.f90 \
	math/math.f90\
	math/mathops.f90\
	math/fast3d.f90\
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
	mesh/tet_mesh.f90\
	field/field.f90\
	field/mesh_field.f90\
	field/mean_field.f90\
	field/mean_sqr_field.f90\
	io/format/rea.f90\
	math/bcknd/cpu/cdtp.f90\
	math/bcknd/cpu/conv1.f90\
	math/bcknd/cpu/dudxyz.f90\
	math/bcknd/cpu/opgrad.f90\
	math/bcknd/sx/sx_cdtp.f90\
	math/bcknd/sx/sx_conv1.f90\
	math/bcknd/sx/sx_dudxyz.f90\
	math/bcknd/sx/sx_opgrad.f90\
	math/operators.f90\
	math/bcknd/cpu/opr_cpu.f90\
	math/bcknd/sx/opr_sx.f90\
	math/bcknd/xsmm/opr_xsmm.F90\
	math/bcknd/device/opr_device.F90\
	math/tensor.f90\
	math/bcknd/cpu/tensor_cpu.f90\
	math/bcknd/sx/tensor_sx.f90\
	math/bcknd/xsmm/tensor_xsmm.F90\
	math/fdm.f90\
	math/bcknd/cpu/fdm_cpu.f90\
	math/bcknd/sx/fdm_sx.f90\
	math/bcknd/xsmm/fdm_xsmm.f90\
	math/schwarz.f90\
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
	io/mean_flow_output.f90\
	io/mean_sqr_flow_output.f90\
	common/sampler.f90\
	bc/bc.f90\
	bc/dirichlet.f90\
	bc/wall.f90\
	bc/inflow.f90\
	bc/usr_inflow.f90\
	bc/facet_normal.f90\
	bc/symmetry.f90\
	krylov/precon.f90\
	krylov/krylov.f90\
	krylov/pc_identity.f90\
	krylov/precon_fctry.f90\
	krylov/krylov_fctry.f90\
	krylov/bcknd/cpu/cg.f90\
	krylov/bcknd/cpu/cacg.f90\
	krylov/bcknd/cpu/pipecg.f90\
	krylov/bcknd/cpu/bicgstab.f90\
	krylov/bcknd/cpu/gmres.f90\
	krylov/bcknd/cpu/pc_jacobi.f90\
	krylov/bcknd/sx/cg_sx.f90\
	krylov/bcknd/sx/pipecg_sx.f90\
	krylov/bcknd/sx/gmres_sx.f90\
	krylov/bcknd/sx/pc_jacobi_sx.f90\
	krylov/bcknd/device/cg_device.f90\
	krylov/bcknd/device/pc_jacobi_device.F90\
	krylov/bcknd/device/pc_identity_device.f90\
	common/source.f90\
	common/abbdf.f90\
	common/stats_quant.f90\
	common/statistics.f90\
	config/neko_config.f90\
	case.f90\
	common/user_intf.f90\
	fluid/fluid_method.f90\
	fluid/fluid_plan4.f90\
	fluid/fluid_plan1.f90\
	fluid/fluid_fctry.f90\
	fluid/mean_flow.f90\
	fluid/mean_sqr_flow.f90\
	simulation.f90\
	math/ax_helm_fctry.f90\
	math/bcknd/cpu/ax_helm.f90\
	math/bcknd/sx/ax_helm_sx.f90\
	math/bcknd/xsmm/ax_helm_xsmm.F90\
	common/projection.f90\
	comm/parmetis.F90\
	comm/redist.f90\
	krylov/pc_hsmg.f90\
	common/signal.f90\
	common/sighdl.c\
	common/jobctrl.f90\
	device/cuda_intf.F90\
	device/hip_intf.F90\
	device/opencl_intf.F90\
	device/device.F90\
	math/bcknd/device/device_math.F90\
	math/bcknd/device/device_mathops.F90\
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
	math/bcknd/device/hip/math.hip\
	math/bcknd/device/hip/mathops.hip\
	math/bcknd/device/hip/opr_dudxyz.hip\
	math/bcknd/device/hip/opr_cdtp.hip\
	math/bcknd/device/hip/opr_conv1.hip\
	math/bcknd/device/hip/opr_opgrad.hip\
	krylov/bcknd/device/hip/pc_jacobi.hip\
	gs/bcknd/device/hip/gs.hip\
	bc/bcknd/device/hip/dirichlet.hip\
	bc/bcknd/device/hip/inflow.hip\
	bc/bcknd/device/hip/no_slip_wall.hip\
	bc/bcknd/device/hip/symmetry.hip\
	bc/bcknd/device/hip/facet_normal.hip
endif

if ENABLE_CUDA
libneko_a_SOURCES += \
	math/bcknd/device/cuda/math.cu\
	math/bcknd/device/cuda/mathops.cu\
	math/bcknd/device/cuda/opr_dudxyz.cu\
	math/bcknd/device/cuda/opr_cdtp.cu\
	math/bcknd/device/cuda/opr_conv1.cu\
	math/bcknd/device/cuda/opr_opgrad.cu\
	gs/bcknd/device/cuda/gs.cu\
	bc/bcknd/device/cuda/dirichlet.cu\
	bc/bcknd/device/cuda/inflow.cu\
	bc/bcknd/device/cuda/no_slip_wall.cu\
	bc/bcknd/device/cuda/symmetry.cu\
	bc/bcknd/device/cuda/facet_normal.cu
endif

if ENABLE_OPENCL
libneko_a_SOURCES += \
	math/bcknd/device/opencl/math.c

%.cl.h: %.cl
	./scripts/cltostring.sh $<
endif

bin_PROGRAMS = neko
neko_SOURCES = $(libneko_a_SOURCES) driver.f90

pkginclude_HEADERS = $(wildcard *.mod) device/device_config.h

CLEANFILES = *.mod

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(libneko_a_SOURCES) driver.f90  > .depends
endif

# Note we need to list any headers etc here for the dist. script
EXTRA_DIST = \
	bc/bcknd/device/cuda/dirichlet_kernel.h\
	bc/bcknd/device/cuda/inflow_kernel.h\
	bc/bcknd/device/cuda/no_slip_wall_kernel.h\
	bc/bcknd/device/cuda/symmetry_kernel.h\
	bc/bcknd/device/cuda/facet_normal_kernel.h\
	gs/bcknd/device/cuda/gs_kernels.h\
	math/bcknd/device/cuda/math_kernel.h\
	math/bcknd/device/cuda/mathops_kernel.h\
	math/bcknd/device/cuda/dudxyz_kernel.h\
	math/bcknd/device/cuda/cdtp_kernel.h\
	math/bcknd/device/cuda/conv1_kernel.h\
	math/bcknd/device/cuda/opgrad_kernel.h\
	bc/bcknd/device/hip/dirichlet_kernel.h\
	bc/bcknd/device/hip/inflow_kernel.h\
	bc/bcknd/device/hip/no_slip_wall_kernel.h\
	bc/bcknd/device/hip/symmetry_kernel.h\
	bc/bcknd/device/hip/facet_normal_kernel.h\
	gs/bcknd/device/hip/gs_kernels.h\
	math/bcknd/device/hip/math_kernel.h\
	math/bcknd/device/hip/mathops_kernel.h\
	math/bcknd/device/hip/dudxyz_kernel.h\
	math/bcknd/device/hip/cdtp_kernel.h\
	math/bcknd/device/hip/conv1_kernel.h\
	math/bcknd/device/hip/opgrad_kernel.h\
	krylov/bcknd/device/hip/pc_jacobi_kernel.h

include	.depends
