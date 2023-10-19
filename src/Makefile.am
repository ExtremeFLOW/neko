
lib_LIBRARIES = libneko.a
neko_fortran_SOURCES = \
	config/num_types.f90\
	common/structs.f90\
	io/format/nmsh.f90\
	io/format/re2.f90\
	io/format/map.f90\
	io/format/stl.f90\
	common/log.f90\
	comm/comm.F90\
	mesh/curve.f90\
	comm/mpi_types.f90\
	common/datadist.f90\
	common/distdata.f90\
	common/utils.f90\
	common/json_utils.f90\
	common/system.f90\
	math/mxm_wrapper.F90\
	sem/speclib.f90 \
	sem/local_interpolation.f90 \
	math/math.f90\
	math/mathops.f90\
	math/fast3d.f90\
	sem/space.f90\
	sem/dofmap.f90\
	sem/coef.f90\
	sem/cpr.f90\
	sem/spectral_error_indicator.f90\
	common/time_interpolator.f90\
	sem/interpolation.f90\
	sem/point_interpolator.f90\
	sem/bcknd/device/device_coef.F90\
	gs/gs_bcknd.f90\
	gs/bcknd/cpu/gs_cpu.f90\
	gs/bcknd/sx/gs_sx.f90\
	gs/bcknd/device/gs_device.F90\
	gs/gs_ops.f90\
	gs/gs_comm.f90\
	gs/gs_mpi.f90\
	gs/bcknd/device/gs_device_mpi.F90\
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
	mesh/facet_zone.f90\
	mesh/point_zone.f90\
	mesh/point_zones/sphere_point_zone.f90\
	mesh/point_zones/box_point_zone.f90\
	mesh/point_zone_factory.f90\
	mesh/point_zone_registry.f90\
	mesh/mesh.f90\
	mesh/octree.f90\
	mesh/tet_mesh.f90\
	mesh/tri_mesh.f90\
	field/field_registry.f90\
	field/scratch_registry.f90\
	field/field.f90\
	field/field_list.f90\
	field/field_series.f90\
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
	math/bcknd/sx/sx_cfl.f90\
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
	math/vector.f90\
	math/matrix.f90\
	common/checkpoint.f90\
	io/generic_file.f90\
	io/map_file.f90\
	io/re2_file.f90\
	io/rea_file.f90\
	io/fld_file.f90\
	io/fld_file_data.f90\
	io/vtk_file.f90\
	io/stl_file.f90\
	io/nmsh_file.f90\
	io/chkp_file.f90\
	io/csv_file.f90\
	io/file.f90\
	io/output.f90\
	io/fluid_output.f90\
	io/chkp_output.f90\
	io/mean_flow_output.f90\
	io/fluid_stats_output.f90\
	io/mean_sqr_flow_output.f90\
	io/data_streamer.F90\
	common/sampler.f90\
	common/global_interpolation.F90\
	common/profiler.F90\
	common/craypat.F90\
	bc/bc.f90\
	bc/dirichlet.f90\
	bc/dong_outflow.f90\
	bc/wall.f90\
	bc/inflow.f90\
	bc/usr_inflow.f90\
	bc/usr_scalar.f90\
	bc/facet_normal.f90\
	bc/symmetry.f90\
	bc/non_normal.f90\
	bc/blasius.f90\
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
	krylov/bcknd/device/pipecg_device.F90\
	krylov/bcknd/device/gmres_device.F90\
	krylov/bcknd/device/pc_jacobi_device.F90\
	krylov/bcknd/device/pc_identity_device.f90\
	time_schemes/time_scheme.f90\
	time_schemes/bdf_time_scheme.f90\
	time_schemes/ext_time_scheme.f90\
	time_schemes/ab_time_scheme.f90\
	time_schemes/time_scheme_controller.f90\
	common/time_based_controller.f90\
	common/stats_quant.f90\
	common/statistics.f90\
	common/rhs_maker.f90\
	common/rhs_maker_fctry.f90\
	simulation_components/probes.F90\
	common/bcknd/cpu/rhs_maker_cpu.f90\
	common/bcknd/sx/rhs_maker_sx.f90\
	common/bcknd/device/rhs_maker_device.F90\
	common/material_properties.f90\
	config/neko_config.f90\
	case.f90\
	common/user_intf.f90\
	fluid/fluid_scheme.f90\
	fluid/fluid_aux.f90\
	fluid/fluid_pnpn.f90\
	fluid/fluid_fctry.f90\
	fluid/fluid_volflow.f90\
	fluid/pnpn_res.f90\
	fluid/pnpn_res_fctry.f90\
	fluid/mean_flow.f90\
	fluid/fluid_stats.f90\
	fluid/mean_sqr_flow.f90\
	fluid/flow_profile.f90\
	fluid/flow_ic.f90\
	fluid/advection.f90\
	fluid/bcknd/cpu/pnpn_res_cpu.f90\
	fluid/bcknd/sx/pnpn_res_sx.f90\
	fluid/bcknd/device/pnpn_res_device.F90\
	fluid/fluid_user_source_term.f90\
	fluid/fluid_source_term.f90\
	simulation.f90\
	math/ax_helm_fctry.f90\
	math/bcknd/cpu/ax_helm.f90\
	math/bcknd/sx/ax_helm_sx.f90\
	math/bcknd/xsmm/ax_helm_xsmm.F90\
	math/bcknd/device/ax_helm_device.F90\
	common/projection.f90\
	common/bcknd/device/device_projection.F90\
	comm/parmetis.F90\
	comm/redist.f90\
	krylov/pc_hsmg.f90\
	common/signal.f90\
	common/jobctrl.f90\
	device/cuda_intf.F90\
	device/hip_intf.F90\
	device/cuda/nvtx.F90\
	device/hip/roctx.F90\
	device/opencl_intf.F90\
	device/opencl/prgm_lib.F90\
	device/dummy_device.F90\
	device/device.F90\
	math/bcknd/device/device_math.F90\
	math/bcknd/device/device_schwarz.F90\
	math/bcknd/device/tensor_device.F90\
	math/bcknd/device/fdm_device.F90\
	math/bcknd/device/device_mathops.F90\
	bc/bcknd/device/device_dirichlet.F90\
	bc/bcknd/device/device_inflow.F90\
	bc/bcknd/device/device_wall.F90\
	bc/bcknd/device/device_symmetry.F90\
	bc/bcknd/device/device_facet_normal.F90\
	bc/bcknd/device/device_inhom_dirichlet.F90\
	bc/bcknd/device/device_dong_outflow.F90\
	scalar/bcknd/device/scalar_residual_device.F90\
	scalar/scalar_scheme.f90\
	scalar/scalar_pnpn.f90\
	scalar/scalar_aux.f90\
	scalar/scalar_residual.f90\
	scalar/scalar_residual_fctry.f90\
	scalar/bcknd/cpu/scalar_residual_cpu.f90\
	scalar/bcknd/sx/scalar_residual_sx.f90\
	scalar/source_scalar.f90\
	simulation_components/simulation_component.f90\
	simulation_components/simulation_component_global.f90\
	simulation_components/vorticity.f90\
	simulation_components/lambda2.f90\
	simulation_components/simulation_component_factory.f90\
	source_terms/source_term.f90\
	source_terms/const_source_term.f90\
	source_terms/bcknd/cpu/const_source_term_cpu.f90\
	source_terms/bcknd/device/const_source_term_device.f90\
	source_terms/source_term_factory.f90\
	neko.f90

neko_c_SOURCES = comm/comm_wrapper.c\
		 common/sighdl.c\
		 common/cpuid.c\
                 gs/bcknd/device/device_mpi.c\
                 math/bcknd/device/device_mpi_reduce.c

libneko_a_SOURCES = $(neko_fortran_SOURCES)\
                    $(neko_c_SOURCES)

if ENABLE_PARMETIS
libneko_a_SOURCES += comm/parmetis_wrapper.c
endif

if ENABLE_HIP
libneko_a_SOURCES += \
	device/hip/check.hip\
	math/bcknd/device/hip/math.hip\
	math/bcknd/device/hip/schwarz.hip\
	math/bcknd/device/hip/tensor.hip\
	math/bcknd/device/hip/fdm.hip\
	math/bcknd/device/hip/mathops.hip\
	math/bcknd/device/hip/opr_dudxyz.hip\
	math/bcknd/device/hip/opr_cdtp.hip\
	math/bcknd/device/hip/opr_conv1.hip\
	math/bcknd/device/hip/opr_opgrad.hip\
	math/bcknd/device/hip/opr_lambda2.hip\
	math/bcknd/device/hip/opr_cfl.hip\
	math/bcknd/device/hip/ax_helm.hip\
	krylov/bcknd/device/hip/pc_jacobi.hip\
	krylov/bcknd/device/hip/pipecg_aux.hip\
	krylov/bcknd/device/hip/gmres_aux.hip\
	gs/bcknd/device/hip/gs.hip\
	bc/bcknd/device/hip/dirichlet.hip\
	bc/bcknd/device/hip/inflow.hip\
	bc/bcknd/device/hip/no_slip_wall.hip\
	bc/bcknd/device/hip/symmetry.hip\
	bc/bcknd/device/hip/facet_normal.hip\
	bc/bcknd/device/hip/inhom_dirichlet.hip\
	bc/bcknd/device/hip/dong_outflow.hip\
	common/bcknd/device/hip/rhs_maker.hip\
	common/bcknd/device/hip/projection.hip\
	fluid/bcknd/device/hip/pnpn_res.hip\
	scalar/bcknd/device/hip/scalar_residual.hip\
	sem/bcknd/device/hip/coef.hip
endif

if ENABLE_CUDA
libneko_a_SOURCES += \
	device/cuda/check.cu\
	math/bcknd/device/cuda/math.cu\
	math/bcknd/device/cuda/schwarz.cu\
	math/bcknd/device/cuda/tensor.cu\
	math/bcknd/device/cuda/fdm.cu\
	math/bcknd/device/cuda/mathops.cu\
	math/bcknd/device/cuda/opr_dudxyz.cu\
	math/bcknd/device/cuda/opr_cdtp.cu\
	math/bcknd/device/cuda/opr_conv1.cu\
	math/bcknd/device/cuda/opr_opgrad.cu\
	math/bcknd/device/cuda/opr_lambda2.cu\
	math/bcknd/device/cuda/opr_cfl.cu\
	math/bcknd/device/cuda/ax_helm.cu\
	krylov/bcknd/device/cuda/pc_jacobi.cu\
	krylov/bcknd/device/cuda/pipecg_aux.cu\
	krylov/bcknd/device/cuda/gmres_aux.cu\
	gs/bcknd/device/cuda/gs.cu\
	bc/bcknd/device/cuda/dirichlet.cu\
	bc/bcknd/device/cuda/inflow.cu\
	bc/bcknd/device/cuda/no_slip_wall.cu\
	bc/bcknd/device/cuda/symmetry.cu\
	bc/bcknd/device/cuda/facet_normal.cu\
	bc/bcknd/device/cuda/inhom_dirichlet.cu\
	bc/bcknd/device/cuda/dong_outflow.cu\
	common/bcknd/device/cuda/rhs_maker.cu\
	common/bcknd/device/cuda/projection.cu\
	fluid/bcknd/device/cuda/pnpn_res.cu\
	scalar/bcknd/device/cuda/scalar_residual.cu\
	sem/bcknd/device/cuda/coef.cu
endif

if ENABLE_OPENCL
libneko_a_SOURCES += \
	device/opencl/check.c\
	device/opencl/jit.c\
	math/bcknd/device/opencl/math.c\
	math/bcknd/device/opencl/schwarz.c\
	math/bcknd/device/opencl/tensor.c\
	math/bcknd/device/opencl/fdm.c\
	math/bcknd/device/opencl/mathops.c\
	math/bcknd/device/opencl/opr_dudxyz.c\
	math/bcknd/device/opencl/opr_cdtp.c\
	math/bcknd/device/opencl/opr_conv1.c\
	math/bcknd/device/opencl/opr_opgrad.c\
	math/bcknd/device/opencl/opr_lambda2.c\
	math/bcknd/device/opencl/opr_cfl.c\
	math/bcknd/device/opencl/ax_helm.c\
	krylov/bcknd/device/opencl/pc_jacobi.c\
	gs/bcknd/device/opencl/gs.c\
	bc/bcknd/device/opencl/dirichlet.c\
	bc/bcknd/device/opencl/inflow.c\
	bc/bcknd/device/opencl/no_slip_wall.c\
	bc/bcknd/device/opencl/symmetry.c\
	bc/bcknd/device/opencl/facet_normal.c\
	bc/bcknd/device/opencl/inhom_dirichlet.c\
	bc/bcknd/device/opencl/dong_outflow.c\
	common/bcknd/device/opencl/rhs_maker.c\
	fluid/bcknd/device/opencl/pnpn_res.c\
	scalar/bcknd/device/opencl/scalar_residual.c\
	sem/bcknd/device/opencl/coef.c

%.cl.h: %.cl
	bash ./scripts/cltostring.sh $<
endif

if ENABLE_ADIOS2
#
# Fortran code related to ADIOS2 can be added directly
# to libneko_fortran_SOURCES (as above)
#
# To avoid linker issues, any C++ code related to ADIOS2
# must be added to libneko_a_LIBADD and neko_LDADD,
# with a manual build rule as illustrated below
#
#io/adios2_example.o: io/adios2_example.cpp
#	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $^
#
#libneko_a_LIBADD = \
#	io/adios2_example.o
#
#neko_LDADD = \
#	io/adios2_example.o
#
#
io/nek_adios2.o: io/nek_adios2.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $^

libneko_a_LIBADD = \
	io/nek_adios2.o

neko_LDADD = \
	io/nek_adios2.o

endif


bin_PROGRAMS = neko
neko_SOURCES = $(libneko_a_SOURCES) driver.f90

pkginclude_HEADERS = $(wildcard *.mod)\
		     device/device_config.h\
		     common/neko_log.h

if ENABLE_HIP
pkginclude_HEADERS += device/hip/check.h
endif

if ENABLE_CUDA
pkginclude_HEADERS += device/cuda/check.h
endif

if ENABLE_OPENCL
pkginclude_HEADERS += device/opencl/check.h
endif


CLEANFILES = *.mod

if ENABLE_OPENCL
CLEANFILES += math/bcknd/device/opencl/*.cl.h\
	      gs/bcknd/device/opencl/*.cl.h\
	      bc/bcknd/device/opencl/*.cl.h\
	      krylov/bcknd/device/opencl/*.cl.h\
	      fluid/bcknd/device/opencl/*.cl.h\
	      scalar/bcknd/device/opencl/*.cl.h\
	      sem/bcknd/device/opencl/*.cl.h
endif

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(neko_fortran_SOURCES) driver.f90  > .depends
endif

# Note we need to list any headers etc here for the dist. script
EXTRA_DIST = \
	bc/bcknd/device/cuda/dirichlet_kernel.h\
	bc/bcknd/device/cuda/inflow_kernel.h\
	bc/bcknd/device/cuda/no_slip_wall_kernel.h\
	bc/bcknd/device/cuda/symmetry_kernel.h\
	bc/bcknd/device/cuda/facet_normal_kernel.h\
	bc/bcknd/device/cuda/inhom_dirichlet_kernel.h\
	bc/bcknd/device/cuda/dong_outflow_kernel.h\
	gs/bcknd/device/cuda/gs_kernels.h\
	math/bcknd/device/cuda/math_kernel.h\
	math/bcknd/device/cuda/tensor_kernel.h\
	math/bcknd/device/cuda/fdm_kernel.h\
	math/bcknd/device/cuda/mathops_kernel.h\
	math/bcknd/device/cuda/schwarz_kernel.h\
	math/bcknd/device/cuda/dudxyz_kernel.h\
	math/bcknd/device/cuda/cdtp_kernel.h\
	math/bcknd/device/cuda/conv1_kernel.h\
	math/bcknd/device/cuda/opgrad_kernel.h\
	math/bcknd/device/cuda/lambda2_kernel.h\
	math/bcknd/device/cuda/cfl_kernel.h\
	math/bcknd/device/cuda/ax_helm_kernel.h\
	common/bcknd/device/cuda/sumab_kernel.h\
	common/bcknd/device/cuda/makeext_kernel.h\
	common/bcknd/device/cuda/makebdf_kernel.h\
	common/bcknd/device/cuda/projection_kernel.h\
	fluid/bcknd/device/cuda/prs_res_kernel.h\
	fluid/bcknd/device/cuda/vel_res_update_kernel.h\
	scalar/bcknd/device/cuda/scalar_residual_update_kernel.h\
	sem/bcknd/device/cuda/coef_kernel.h\
	bc/bcknd/device/hip/dirichlet_kernel.h\
	bc/bcknd/device/hip/inflow_kernel.h\
	bc/bcknd/device/hip/no_slip_wall_kernel.h\
	bc/bcknd/device/hip/symmetry_kernel.h\
	bc/bcknd/device/hip/facet_normal_kernel.h\
	bc/bcknd/device/hip/inhom_dirichlet_kernel.h\
	bc/bcknd/device/hip/dong_outflow_kernel.h\
	gs/bcknd/device/hip/gs_kernels.h\
	math/bcknd/device/hip/math_kernel.h\
	math/bcknd/device/hip/tensor_kernel.h\
	math/bcknd/device/hip/schwarz_kernel.h\
	math/bcknd/device/hip/fdm_kernel.h\
	math/bcknd/device/hip/mathops_kernel.h\
	math/bcknd/device/hip/dudxyz_kernel.h\
	math/bcknd/device/hip/cdtp_kernel.h\
	math/bcknd/device/hip/conv1_kernel.h\
	math/bcknd/device/hip/opgrad_kernel.h\
	math/bcknd/device/hip/lambda2_kernel.h\
	math/bcknd/device/hip/cfl_kernel.h\
	math/bcknd/device/hip/ax_helm_kernel.h\
	common/bcknd/device/hip/sumab_kernel.h\
	common/bcknd/device/hip/makeext_kernel.h\
	common/bcknd/device/hip/makebdf_kernel.h\
	common/bcknd/device/hip/projection_kernel.h\
	fluid/bcknd/device/hip/prs_res_kernel.h\
	fluid/bcknd/device/hip/vel_res_update_kernel.h\
	scalar/bcknd/device/hip/scalar_residual_update_kernel.h\
	krylov/bcknd/device/cuda/pipecg_kernel.h\
	krylov/bcknd/device/hip/pipecg_kernel.h\
	krylov/bcknd/device/cuda/gmres_kernel.h\
	krylov/bcknd/device/hip/gmres_kernel.h\
	sem/bcknd/device/hip/coef_kernel.h\
	bc/bcknd/device/opencl/dirichlet_kernel.cl\
	bc/bcknd/device/opencl/inflow_kernel.cl\
	bc/bcknd/device/opencl/no_slip_wall_kernel.cl\
	bc/bcknd/device/opencl/symmetry_kernel.cl\
	bc/bcknd/device/opencl/facet_normal_kernel.cl\
	bc/bcknd/device/opencl/inhom_dirichlet_kernel.cl\
	bc/bcknd/device/opencl/dong_outflow_kernel.cl\
	math/bcknd/device/opencl/math_kernel.cl\
	math/bcknd/device/opencl/schwarz_kernel.cl\
	math/bcknd/device/opencl/tensor_kernel.cl\
	math/bcknd/device/opencl/fdm_kernel.cl\
	math/bcknd/device/opencl/mathops_kernel.cl\
	math/bcknd/device/opencl/dudxyz_kernel.cl\
	math/bcknd/device/opencl/ax_helm_kernel.cl\
	math/bcknd/device/opencl/cdtp_kernel.cl\
	math/bcknd/device/opencl/opgrad_kernel.cl\
	math/bcknd/device/opencl/conv1_kernel.cl\
	math/bcknd/device/opencl/lambda2_kernel.cl\
	math/bcknd/device/opencl/cfl_kernel.cl\
	common/bcknd/device/opencl/rhs_maker_kernel.cl\
	fluid/bcknd/device/opencl/pnpn_res_kernel.cl\
	scalar/bcknd/device/opencl/scalar_residual_kernel.cl\
	gs/bcknd/device/opencl/gs_kernels.cl\
	krylov/bcknd/device/opencl/jacobi_kernel.cl\
	sem/bcknd/device/opencl/coef_kernel.cl\
	math/bcknd/device/device_mpi_reduce.h\
	math/bcknd/device/device_mpi_op.h\
	device/opencl/jit.h\
	device/opencl/prgm_lib.h\
	device/opencl/check.h\
	device/hip/check.h\
	device/cuda/check.h\
	comm/comm.h\
	common/neko_log.h

# Fortran deps. (updated with makedepf90 via target 'depend')
include	.depends

# Device deps. (manually updated)
include .depends_device

# CUDA Build rules
.cu.o:
	$(NVCC) @CUDA_ARCH@ @CUDA_CFLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<

# HIP Build rules
.hip.o:
	$(HIPCC) @HIP_HIPCC_FLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<
