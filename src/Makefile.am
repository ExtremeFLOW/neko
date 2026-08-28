lib_LTLIBRARIES = libneko.la
libneko_la_LDFLAGS = -release ${PACKAGE_VERSION}
neko_fortran_SOURCES = \
	ale/ale_rigid_kinematics.f90\
	ale/ale_manager.f90\
	ale/bcknd/cpu/ale_routines_cpu.f90\
	ale/bcknd/device/ale_routines_device.F90\
	api/neko_api.f90\
	api/neko_api_user.f90\
	config/num_types.f90\
	common/structs.f90\
	io/format/nmsh.f90\
	io/format/re2.f90\
	io/format/map.f90\
	io/format/stl.f90\
	io/format/vtk.f90\
	io/buffer/buffer.F90\
	io/buffer/buffer_1d.F90\
	io/buffer/buffer_4d.F90\
	io/buffer/buffer_4d_npar.F90\
	common/log.f90\
	comm/comm.F90\
	comm/pfunit_comm_utils.F90\
	comm/crystal_router.f90\
	mesh/curve.f90\
	comm/mpi_types.f90\
	comm/shmem.F90\
	common/datadist.f90\
	common/distdata.f90\
	common/utils.f90\
	common/time_state.f90\
	common/json_utils.f90\
	common/case_file_utils.f90\
	common/mask.f90\
	common/expression.f90\
	common/system.f90\
	common/runtime_statistics.f90\
	common/user_access_singleton.f90\
	data_types/host_array.f90\
	data_types/device_array.f90\
	math/mxm_wrapper.F90\
	sem/speclib.f90 \
	qoi/drag_torque.f90\
	sem/local_interpolation.f90 \
	math/math.f90\
	math/mathops.f90\
	math/fast3d.f90\
	sem/space.f90\
	sem/map_1d.f90\
	sem/map_2d.f90\
	sem/dofmap.f90\
	sem/coef.f90\
	sem/cpr.f90\
	common/time_interpolator.f90\
	sem/interpolation.f90\
	sem/point_interpolator.f90\
	sem/bcknd/device/device_coef.F90\
	sem/bcknd/device/device_local_interpolation.F90\
	gs/gs_bcknd.f90\
	gs/bcknd/cpu/gs_cpu.f90\
	gs/bcknd/sx/gs_sx.f90\
	gs/bcknd/device/gs_device.F90\
	gs/gs_ops.f90\
	gs/gs_comm.f90\
	gs/gs_mpi.f90\
	gs/gs_mpi_rma.f90\
	gs/gs_shmem.F90\
	gs/gs_caf.F90\
	gs/gs_neighbour.f90\
	gs/gs_utofu.F90\
	gs/bcknd/device/gs_device_mpi.F90\
	gs/bcknd/device/gs_device_nccl.F90\
	gs/bcknd/device/gs_device_shmem.F90\
	gs/gather_scatter.f90\
	gs/gs_tune.f90\
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
	mesh/point_zones/cylinder_point_zone.f90\
	mesh/point_zones/combine_point_zone.f90\
	mesh/point_zone_fctry.f90\
	mesh/point_zone_registry.f90\
	mesh/mesh.f90\
	mesh/octree.f90\
	mesh/search_tree/aabb.f90\
	mesh/aabb_tree.f90\
	mesh/tet_mesh.f90\
	mesh/tri_mesh.f90\
	mesh/intersection_detector.f90\
	field/field.f90\
	field/field_array.f90\
	field/field_list.f90\
	field/field_series.f90\
	field/field_series_list.f90\
	field/mesh_field.f90\
	field/mean_field.f90\
	io/format/rea.f90\
	math/bcknd/cpu/cpu_cdtp.f90\
	math/bcknd/cpu/cpu_conv1.f90\
	math/bcknd/cpu/cpu_dudxyz.f90\
	math/bcknd/cpu/cpu_opgrad.f90\
	math/bcknd/cpu/convect_scalar.f90\
	math/bcknd/cpu/set_convect_rst.f90\
	math/bcknd/sx/sx_cdtp.f90\
	math/bcknd/sx/sx_conv1.f90\
	math/bcknd/sx/sx_dudxyz.f90\
	math/bcknd/sx/sx_opgrad.f90\
	math/bcknd/sx/sx_cfl.f90\
	math/bcknd/sx/sx_lambda2.f90\
	math/bcknd/sx/sx_convect_scalar.f90\
	math/bcknd/sx/sx_set_convect_rst.f90\
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
	common/particles.f90\
	math/vector_series.f90\
	math/vector_list.f90\
	math/vector_math.f90\
	registries/registry.f90\
	registries/registry_entry.f90\
	registries/scratch_registry.f90\
	math/matrix.f90\
	math/matrix_math.f90\
	math/signed_distance.f90\
	multigrid/tree_amg.f90\
	multigrid/tree_amg_utils.f90\
	multigrid/tree_amg_aggregate.f90\
	multigrid/tree_amg_multigrid.f90\
	multigrid/tree_amg_smoother.f90\
	multigrid/bcknd/device/device_tree_amg_smoother.F90\
	multigrid/phmg.f90\
	common/checkpoint.f90\
	io/generic_file.f90\
	io/map_file.f90\
	io/re2_file.f90\
	io/rea_file.f90\
	io/bp_file.F90\
	io/fld_file.f90\
	io/import_field_utils.f90\
	io/fld_file_data.f90\
	io/vtk_file.f90\
	io/stl_file.f90\
	io/nmsh_file.f90\
	io/chkp_file.f90\
	io/csv_file.f90\
	io/hdf5_session.F90\
	io/hdf5_file.F90\
	io/vtkhdf_file.F90\
	io/file.f90\
	io/output.f90\
	io/fluid_output.f90\
	io/field_output.f90\
	io/fld_file_output.f90\
	io/chkp_output.f90\
	io/fluid_stats_output.f90\
	io/fluid_sgs_stats_output.f90\
	io/scalar_stats_output.f90\
	io/scalar_sgs_stats_output.f90\
	io/spatial_average_output.f90\
	io/data_streamer.F90\
	io/output_controller.f90\
	global_interpolation/global_interpolation.f90\
	global_interpolation/global_interpolation_comm.f90\
	global_interpolation/aabb_pe_finder.f90\
	global_interpolation/aabb_el_finder.f90\
	global_interpolation/cartesian_el_finder.f90\
	global_interpolation/cartesian_pe_finder.f90\
	global_interpolation/el_finder.f90\
	global_interpolation/pe_finder.f90\
	global_interpolation/legendre_rst_finder.f90\
	io/mean_field_output.f90\
	common/profiler.F90\
	common/craypat.F90\
	bc/scalar_bc_projector.f90\
	bc/vector_bc_projector.f90\
	bc/bc.f90\
	bc/mixed_bc.f90\
	bc/bc_list.f90\
	bc/dirichlet.f90\
	bc/neumann.f90\
	bc/no_slip.f90\
	bc/shear_stress.f90\
	bc/dong_outflow.f90\
	bc/zero_dirichlet.f90\
	bc/inflow.f90\
	bc/field_dirichlet.f90\
	bc/field_neumann.f90\
	bc/field_dirichlet_vector.f90\
	bc/expression_dirichlet.f90\
	bc/expression_dirichlet_vector.f90\
	bc/overset_interface.f90\
	bc/overset_interface_vector.f90\
	bc/facet_normal.f90\
	bc/symmetry.f90\
	bc/symmetry_aligned.f90\
	bc/non_normal_aligned.f90\
	bc/non_normal.f90\
	bc/blasius.f90\
	bc/wall_model_bc.f90\
	krylov/precon.f90\
	krylov/krylov.f90\
	krylov/pc_identity.f90\
	krylov/precon_fctry.f90\
	krylov/krylov_fctry.f90\
	krylov/bcknd/cpu/cg.f90\
	krylov/bcknd/cpu/cacg.f90\
	krylov/bcknd/cpu/cheby.f90\
	krylov/bcknd/cpu/pipecg.f90\
	krylov/bcknd/cpu/bicgstab.f90\
	krylov/bcknd/cpu/gmres.f90\
	krylov/bcknd/cpu/pc_jacobi.f90\
	krylov/bcknd/cpu/cg_coupled.f90\
	krylov/bcknd/sx/cg_sx.f90\
	krylov/bcknd/sx/pipecg_sx.f90\
	krylov/bcknd/sx/gmres_sx.f90\
	krylov/bcknd/sx/pc_jacobi_sx.f90\
	krylov/bcknd/device/cg_device.f90\
	krylov/bcknd/device/cg_cpld_device.f90\
	krylov/bcknd/device/cheby_device.F90\
	krylov/bcknd/device/pipecg_device.F90\
	krylov/bcknd/device/fusedcg_device.F90\
	krylov/bcknd/device/fusedcg_cpld_device.F90\
	krylov/bcknd/device/gmres_device.F90\
	krylov/bcknd/device/pc_jacobi_device.F90\
	krylov/bcknd/device/pc_identity_device.f90\
	time_schemes/time_scheme.f90\
	time_schemes/bdf_time_scheme.f90\
	time_schemes/ext_time_scheme.f90\
	time_schemes/iextm_time_scheme.f90\
	time_schemes/ab_time_scheme.f90\
	time_schemes/time_scheme_controller.f90\
	time_schemes/runge_kutta_scheme.f90\
	common/boundary_data.f90\
	common/time_based_controller.f90\
	common/stats_quant.f90\
	common/statistics.f90\
	common/rhs_maker.f90\
	common/rhs_maker_fctry.f90\
	simulation_components/probes.F90\
	simulation_components/boundary_operation.f90\
	simulation_components/boundary_flux.f90\
	simulation_components/boundary_data_writer_simcomp.f90\
	simulation_components/data_streamer_simcomp.F90\
	simulation_components/field_subsampler.f90\
	simulation_components/force_torque.f90\
	simulation_components/fluid_stats_simcomp.f90\
	simulation_components/fluid_sgs_stats_simcomp.f90\
	simulation_components/scalar_stats_simcomp.f90\
	simulation_components/scalar_sgs_stats_simcomp.f90\
	simulation_components/spatial_average.f90\
	simulation_components/user_stats.f90\
	simulation_components/field_writer.f90\
	simulation_components/wall_shear_stress_simcomp.f90\
	common/bcknd/cpu/rhs_maker_cpu.f90\
	common/bcknd/sx/rhs_maker_sx.f90\
	common/bcknd/device/rhs_maker_device.F90\
	config/neko_config.f90\
	case.f90\
	common/user_intf.f90\
	fluid/stress_formulation/pnpn_res_stress_fctry.f90\
	fluid/stress_formulation/bcknd/cpu/pnpn_res_stress_cpu.f90\
	fluid/stress_formulation/bcknd/device/pnpn_res_stress_device.F90\
	fluid/fluid_scheme_base.f90\
	fluid/fluid_scheme_incompressible.f90\
	fluid/fluid_scheme_compressible.f90\
	fluid/regularization.f90\
	fluid/regularization_fctry.f90\
	fluid/regularization/entropy_viscosity.f90\
	fluid/regularization/bcknd/cpu/entropy_viscosity_cpu.f90\
	fluid/regularization/bcknd/device/entropy_viscosity_device.F90\
	fluid/fluid_scheme_compressible_ns.f90\
	fluid/fluid_aux.f90\
	fluid/fluid_pnpn.f90\
	fluid/fluid_base_fctry.f90\
	fluid/fluid_pnpn_bc_fctry.f90\
	fluid/compressible_ns_bc_fctry.f90\
	fluid/fluid_volflow.f90\
	fluid/pnpn_res.f90\
	fluid/pnpn_res_fctry.f90\
	fluid/compressible_res.f90\
	fluid/compressible_res_fctry.f90\
	fluid/fluid_stats.f90\
	fluid/flow_profile.f90\
	fluid/flow_ic.f90\
	fluid/advection.f90\
	fluid/advection_fctry.f90\
	fluid/bcknd/advection/adv_dealias.f90\
	fluid/bcknd/advection/adv_no_dealias.f90\
	fluid/bcknd/advection/adv_oifs.f90\
	fluid/bcknd/advection/adv_dummy.f90\
	fluid/bcknd/cpu/pnpn_res_cpu.f90\
	fluid/bcknd/sx/pnpn_res_sx.f90\
	fluid/bcknd/device/pnpn_res_device.F90\
	fluid/bcknd/cpu/compressible_res_cpu.f90\
	fluid/bcknd/device/compressible_res_device.F90\
	fluid/bcknd/cpu/compressible_ops_cpu.f90\
	fluid/bcknd/device/compressible_ops_device.F90\
	fluid/fluid_source_term.f90\
	common/time_step_controller.f90\
	simulation.f90\
	math/ax_helm_fctry.f90\
	math/ax_helm.f90\
	math/ax_helm_full.f90\
	math/bcknd/cpu/ax_helm_cpu.f90\
	math/bcknd/cpu/ax_helm_vector_cpu.f90\
	math/bcknd/cpu/ax_helm_full_cpu.f90\
	math/bcknd/sx/ax_helm_sx.f90\
	math/bcknd/xsmm/ax_helm_xsmm.F90\
	math/bcknd/device/ax_helm_device.F90\
	math/bcknd/device/ax_helm_full_device.F90\
	common/projection.f90\
	common/projection_vel.f90\
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
	device/metal_intf.F90\
	device/device.F90\
	math/field_math.f90\
	math/bcknd/device/device_math.F90\
	math/bcknd/device/hip/hip_math.f90\
	math/bcknd/device/cuda/cuda_math.f90\
	math/bcknd/device/opencl/opencl_math.f90\
	math/bcknd/device/metal/metal_math.f90\
	math/bcknd/device/device_schwarz.F90\
	math/bcknd/device/tensor_device.F90\
	math/bcknd/device/fdm_device.F90\
	math/bcknd/device/device_mathops.F90\
	bc/bcknd/device/device_dirichlet.F90\
	bc/bcknd/device/device_inflow.F90\
	bc/bcknd/device/device_zero_dirichlet.F90\
	bc/bcknd/device/device_symmetry_aligned.F90\
	bc/bcknd/device/device_facet_normal.F90\
	bc/bcknd/device/device_neumann.F90\
	bc/bcknd/device/device_inhom_dirichlet.F90\
	bc/bcknd/device/device_dong_outflow.F90\
	bc/bcknd/device/device_constrain_mixed_bc.F90\
	bc/bcknd/device/device_coupled_vector_bc_projector.F90\
	scalar/bcknd/device/scalar_residual_device.F90\
	scalar/scalar_scheme.f90\
	scalar/scalar_scheme_fctry.f90\
	scalar/scalars.f90\
	scalar/scalar_pnpn.f90\
	scalar/scalar_pnpn_bc_fctry.f90\
	scalar/scalar_aux.f90\
	scalar/scalar_residual.f90\
	scalar/scalar_residual_fctry.f90\
	scalar/scalar_ic.f90\
	scalar/scalar_source_term.f90\
	scalar/scalar_stats.f90\
	scalar/bcknd/cpu/scalar_residual_cpu.f90\
	scalar/bcknd/sx/scalar_residual_sx.f90\
	simulation_components/simulation_component.f90\
	simulation_components/simcomp_executor.f90\
	simulation_components/spectral_error.f90\
	simulation_components/lambda2.f90\
	lpt/bcknd/cpu/lpt_periodic_bc_cpu.f90\
	lpt/bcknd/cpu/lpt_wall_collision_cpu.f90\
	lpt/bcknd/device/lpt_periodic_bc_device.F90\
	lpt/bcknd/device/lpt_wall_collision_device.F90\
	lpt/lpt_periodic_bc.f90\
	lpt/lpt_output.f90\
	lpt/lpt.f90\
	simulation_components/lpt_simcomp.f90\
	lpt/lpt_migration.f90\
	lpt/lpt_wall_collision.f90\
	simulation_components/curl_simcomp.f90\
	simulation_components/divergence_simcomp.f90\
	simulation_components/weak_gradient_simcomp.f90\
	simulation_components/gradient_simcomp.f90\
	simulation_components/derivative_simcomp.f90\
	simulation_components/les_simcomp.f90\
	simulation_components/simulation_component_fctry.f90\
	source_terms/source_term.f90\
	source_terms/user_source_term.f90\
	source_terms/translation_source_term.f90\
	source_terms/coriolis_source_term.f90\
	source_terms/bcknd/cpu/coriolis_source_term_cpu.f90\
	source_terms/bcknd/device/coriolis_source_term_device.f90\
	source_terms/centrifugal_source_term.f90\
	source_terms/bcknd/cpu/centrifugal_source_term_cpu.f90\
	source_terms/bcknd/device/centrifugal_source_term_device.f90\
	source_terms/source_term_handler.f90\
	source_terms/const_source_term.f90\
	source_terms/field_source_term.f90\
	source_terms/bcknd/cpu/const_source_term_cpu.f90\
	source_terms/bcknd/device/const_source_term_device.f90\
	source_terms/boussinesq_source_term.f90\
	source_terms/bcknd/cpu/boussinesq_source_term_cpu.f90\
	source_terms/bcknd/device/boussinesq_source_term_device.f90\
	source_terms/hpfrt_source_term.f90\
	source_terms/source_term_fctry.f90\
	source_terms/sponge_source_term.f90\
	source_terms/brinkman_source_term.f90\
	filter/mappings.f90\
	filter/filter.f90\
	filter/filter_fctry.f90\
	filter/PDE_filter.f90\
	filter/elementwise_filter.f90\
	filter/bcknd/cpu/mappings_cpu.f90\
	filter/bcknd/device/mappings_device.F90\
	filter/bcknd/device/cuda/cuda_mappings.f90\
	filter/bcknd/device/hip/hip_mappings.f90\
	filter/bcknd/device/opencl/opencl_mappings.f90\
	filter/bcknd/device/metal/metal_mappings.f90\
	source_terms/gradient_jump_penalty.f90\
	source_terms/bcknd/device/device_gradient_jump_penalty.F90\
	les/les_model.f90\
	les/les_model_fctry.f90\
	les/fluid_sgs_stats.f90\
	les/scalar_sgs_stats.f90\
	les/vreman.f90\
	les/bcknd/cpu/vreman_cpu.f90\
	les/bcknd/device/vreman_device.f90\
	les/bcknd/device/device_vreman_nut.F90\
	les/smagorinsky.f90\
	les/bcknd/cpu/smagorinsky_cpu.f90\
	les/bcknd/device/smagorinsky_device.f90\
	les/bcknd/device/device_smagorinsky_nut.F90\
	les/dynamic_smagorinsky.f90\
	les/bcknd/cpu/dynamic_smagorinsky_cpu.f90\
	les/bcknd/device/dynamic_smagorinsky_device.f90\
	les/bcknd/device/device_dynamic_smagorinsky_nut.F90\
	les/sigma.f90\
	les/bcknd/cpu/sigma_cpu.f90\
	les/bcknd/device/sigma_device.f90\
	les/bcknd/device/device_sigma_nut.F90\
	les/wale.f90\
	les/bcknd/cpu/wale_cpu.f90\
	les/bcknd/device/wale_device.f90\
	les/bcknd/device/device_wale_nut.F90\
	les/deardorff.f90\
	les/bcknd/cpu/deardorff_cpu.f90\
	les/bcknd/device/deardorff_device.f90\
	les/bcknd/device/device_deardorff_nut.F90\
	wall_models/wall_sampler.f90\
	wall_models/wall_gll_sampler.f90\
	wall_models/wall_distance_sampler.f90\
	wall_models/wall_sampler_fctry.f90\
	wall_models/wall_model.f90\
	wall_models/bcknd/device/wall_model_device.F90\
	wall_models/spalding.f90\
	wall_models/cai_sagaut_model_ii.f90\
	wall_models/bcknd/cpu/spalding_cpu.f90\
	wall_models/bcknd/cpu/cai_sagaut_model_ii_cpu.f90\
	wall_models/bcknd/device/spalding_device.F90\
	wall_models/bcknd/device/cai_sagaut_model_ii_device.F90\
	wall_models/rough_log_law.f90\
	wall_models/bcknd/cpu/rough_log_law_cpu.f90\
	wall_models/bcknd/device/rough_log_law_device.F90\
	wall_models/most.f90\
	wall_models/bcknd/cpu/most_cpu.f90\
	wall_models/bcknd/device/most_device.F90\
	wall_models/richardson.f90\
	wall_models/bcknd/cpu/richardson_cpu.f90\
	wall_models/bcknd/device/richardson_device.F90\
	wall_models/wall_model_fctry.f90\
	neko.f90

neko_c_SOURCES = \
	comm/comm_wrapper.c\
	comm/comm_nccl.c\
	common/sighdl.c\
	common/cpuid.c\
	gs/bcknd/device/device_mpi.c\
	gs/bcknd/device/device_nccl.c\
	math/bcknd/device/device_mpi_reduce.c\
	math/bcknd/device/device_nccl_reduce.c

libneko_la_SOURCES = \
	$(neko_fortran_SOURCES)\
	$(neko_c_SOURCES)

libneko_la_LIBADD =

if ENABLE_PARMETIS
libneko_la_SOURCES += comm/parmetis_wrapper.c
endif

if ENABLE_UTOFU
libneko_la_SOURCES += gs/gs_utofu_aux.c
endif

if ENABLE_HIP
libneko_la_SOURCES += \
    ale/bcknd/device/hip/ale_kinematics.hip\
	device/hip/check.hip\
	device/hip/buffer.hip\
	device/hip/unified.hip\
	math/bcknd/device/hip/math.hip\
	math/bcknd/device/hip/schwarz.hip\
	math/bcknd/device/hip/tensor.hip\
	math/bcknd/device/hip/fdm.hip\
	math/bcknd/device/hip/mathops.hip\
	math/bcknd/device/hip/opr_dudxyz.hip\
	math/bcknd/device/hip/opr_cdtp.hip\
	math/bcknd/device/hip/opr_conv1.hip\
	math/bcknd/device/hip/opr_convect_scalar.hip\
	math/bcknd/device/hip/set_convect_rst.hip\
	math/bcknd/device/hip/opr_opgrad.hip\
	math/bcknd/device/hip/opr_lambda2.hip\
	math/bcknd/device/hip/opr_cfl.hip\
	math/bcknd/device/hip/opr_rotate_cyc.hip\
	math/bcknd/device/hip/ax_helm.hip\
	math/bcknd/device/hip/ax_helm_full.hip\
	krylov/bcknd/device/hip/pc_jacobi.hip\
	krylov/bcknd/device/hip/pipecg_aux.hip\
	krylov/bcknd/device/hip/fusedcg_aux.hip\
	krylov/bcknd/device/hip/fusedcg_cpld_aux.hip\
	krylov/bcknd/device/hip/gmres_aux.hip\
	gs/bcknd/device/hip/gs.hip\
	bc/bcknd/device/hip/dirichlet.hip\
	bc/bcknd/device/hip/inflow.hip\
	bc/bcknd/device/hip/zero_dirichlet.hip\
	bc/bcknd/device/hip/symmetry_aligned.hip\
	bc/bcknd/device/hip/constrain_mixed_bc.hip\
	bc/bcknd/device/hip/coupled_vector_bc_projector.hip\
	bc/bcknd/device/hip/facet_normal.hip\
	bc/bcknd/device/hip/neumann.hip\
	bc/bcknd/device/hip/inhom_dirichlet.hip\
	bc/bcknd/device/hip/dong_outflow.hip\
	common/bcknd/device/hip/rhs_maker.hip\
	common/bcknd/device/hip/projection.hip\
	fluid/bcknd/device/hip/pnpn_res.hip\
	fluid/bcknd/device/hip/compressible_res.hip\
	fluid/bcknd/device/hip/compressible_ops_compute_max_wave_speed.hip\
	fluid/bcknd/device/hip/compressible_ops_compute_entropy.hip\
	fluid/bcknd/device/hip/compressible_ops_update.hip\
	fluid/regularization/bcknd/device/hip/entropy_viscosity.hip\
	fluid/stress_formulation/bcknd/device/hip/pnpn_stress_res.hip\
	scalar/bcknd/device/hip/scalar_residual.hip\
	sem/bcknd/device/hip/coef.hip\
	sem/bcknd/device/hip/local_interpolation.hip\
	filter/bcknd/device/hip/mappings.hip\
	lpt/bcknd/device/hip/lpt_periodic_bc.hip\
	lpt/bcknd/device/hip/lpt_wall_collision.hip\
	source_terms/bcknd/device/hip/gradient_jump_penalty.hip\
	les/bcknd/device/hip/sigma_nut.hip\
	les/bcknd/device/hip/vreman_nut.hip\
	les/bcknd/device/hip/smagorinsky_nut.hip\
	les/bcknd/device/hip/dynamic_smagorinsky_nut.hip\
	les/bcknd/device/hip/wale_nut.hip\
	les/bcknd/device/hip/deardorff_nut.hip\
	wall_models/bcknd/device/hip/wall_model.hip\
	wall_models/bcknd/device/hip/spalding.hip\
	wall_models/bcknd/device/hip/cai_sagaut_model_ii.hip\
	wall_models/bcknd/device/hip/rough_log_law.hip\
	wall_models/bcknd/device/hip/most.hip\
	wall_models/bcknd/device/hip/richardson.hip\
	multigrid/bcknd/device/hip/amg_cheby.hip\
	krylov/bcknd/device/hip/cheby_aux.hip
endif

if ENABLE_CUDA
libneko_la_SOURCES += \
    ale/bcknd/device/cuda/ale_kinematics.cu\
	device/cuda/check.cu\
	device/cuda/buffer.cu\
	device/cuda/nvtx_wrapper.cu\
	math/bcknd/device/cuda/math.cu\
	math/bcknd/device/cuda/schwarz.cu\
	math/bcknd/device/cuda/tensor.cu\
	math/bcknd/device/cuda/fdm.cu\
	math/bcknd/device/cuda/mathops.cu\
	math/bcknd/device/cuda/opr_dudxyz.cu\
	math/bcknd/device/cuda/opr_cdtp.cu\
	math/bcknd/device/cuda/opr_conv1.cu\
	math/bcknd/device/cuda/opr_convect_scalar.cu\
	math/bcknd/device/cuda/set_convect_rst.cu\
	math/bcknd/device/cuda/opr_opgrad.cu\
	math/bcknd/device/cuda/opr_lambda2.cu\
	math/bcknd/device/cuda/opr_cfl.cu\
	math/bcknd/device/cuda/opr_rotate_cyc.cu\
	math/bcknd/device/cuda/ax_helm.cu\
	math/bcknd/device/cuda/ax_helm_full.cu\
	krylov/bcknd/device/cuda/pc_jacobi.cu\
	krylov/bcknd/device/cuda/pipecg_aux.cu\
	krylov/bcknd/device/cuda/fusedcg_aux.cu\
	krylov/bcknd/device/cuda/fusedcg_cpld_aux.cu\
	krylov/bcknd/device/cuda/gmres_aux.cu\
	gs/bcknd/device/cuda/gs.cu\
	bc/bcknd/device/cuda/dirichlet.cu\
	bc/bcknd/device/cuda/inflow.cu\
	bc/bcknd/device/cuda/zero_dirichlet.cu\
	bc/bcknd/device/cuda/symmetry_aligned.cu\
	bc/bcknd/device/cuda/facet_normal.cu\
	bc/bcknd/device/cuda/neumann.cu\
	bc/bcknd/device/cuda/inhom_dirichlet.cu\
	bc/bcknd/device/cuda/dong_outflow.cu\
	bc/bcknd/device/cuda/constrain_mixed_bc.cu\
	bc/bcknd/device/cuda/coupled_vector_bc_projector.cu\
	common/bcknd/device/cuda/rhs_maker.cu\
	common/bcknd/device/cuda/projection.cu\
	fluid/bcknd/device/cuda/pnpn_res.cu\
	fluid/bcknd/device/cuda/compressible_res.cu\
	fluid/bcknd/device/cuda/compressible_ops_compute_max_wave_speed.cu\
	fluid/bcknd/device/cuda/compressible_ops_compute_entropy.cu\
	fluid/bcknd/device/cuda/compressible_ops_update.cu\
	fluid/regularization/bcknd/device/cuda/entropy_viscosity.cu\
	fluid/stress_formulation/bcknd/device/cuda/pnpn_stress_res.cu\
	scalar/bcknd/device/cuda/scalar_residual.cu\
	sem/bcknd/device/cuda/coef.cu\
	sem/bcknd/device/cuda/local_interpolation.cu\
	filter/bcknd/device/cuda/mappings.cu\
	lpt/bcknd/device/cuda/lpt_periodic_bc.cu\
	lpt/bcknd/device/cuda/lpt_wall_collision.cu\
	source_terms/bcknd/device/cuda/gradient_jump_penalty.cu\
	les/bcknd/device/cuda/sigma_nut.cu\
	les/bcknd/device/cuda/vreman_nut.cu\
	les/bcknd/device/cuda/smagorinsky_nut.cu\
	les/bcknd/device/cuda/dynamic_smagorinsky_nut.cu\
	les/bcknd/device/cuda/wale_nut.cu\
	les/bcknd/device/cuda/deardorff_nut.cu\
	comm/bcknd/device/cuda/comm_nvshmem.cu\
	wall_models/bcknd/device/cuda/wall_model.cu\
	wall_models/bcknd/device/cuda/spalding.cu\
	wall_models/bcknd/device/cuda/cai_sagaut_model_ii.cu\
	wall_models/bcknd/device/cuda/rough_log_law.cu\
	wall_models/bcknd/device/cuda/most.cu\
	wall_models/bcknd/device/cuda/richardson.cu\
	multigrid/bcknd/device/cuda/amg_cheby.cu\
	krylov/bcknd/device/cuda/cheby_aux.cu

BUILT_SOURCES = @top_builddir@/libtool_cuda

if ENABLE_NVSHMEM
#
# Device-side NVSHMEM calls (gs_nvshmem_kernels.h) require relocatable
# device code and a device link against libnvshmem_device. The device
# code in gs_nvshmem.cu is self-contained, so the device link can be
# closed here, keeping nvcc out of the final link. Built under a
# different object name (and added to LIBADD, as with Metal) to not
# clash with the suffix rule and .depends_device
#
gs/bcknd/device/cuda/gs_nvshmem_rdc.lo: gs/bcknd/device/cuda/gs_nvshmem.cu gs/bcknd/device/cuda/gs_kernels.h gs/bcknd/device/cuda/gs_nvshmem_kernels.h @top_builddir@/libtool_cuda
	$(SHELL) @top_builddir@/libtool_cuda --tag=CC --mode=compile $(NVCC) -rdc=true @CUDA_ARCH@ @CUDA_CFLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<

#
# Wrap the device linked object with a libtool object
# (pic only; libtool falls back to the pic object for static archives)
#
gs/bcknd/device/cuda/gs_nvshmem_dlink.lo: gs/bcknd/device/cuda/gs_nvshmem_rdc.lo
	$(MKDIR_P) gs/bcknd/device/cuda/.libs
	obj=gs/bcknd/device/cuda/.libs/gs_nvshmem_rdc.o; \
	test -f $$obj || obj=gs/bcknd/device/cuda/gs_nvshmem_rdc.o; \
	$(NVCC) @CUDA_ARCH@ -dlink -Xcompiler -fPIC \
		-o gs/bcknd/device/cuda/.libs/gs_nvshmem_dlink.o \
		$$obj @NVSHMEM_LDFLAGS@ -lnvshmem_device
	printf "%s\n" \
		"# gs_nvshmem_dlink.lo - a libtool object file" \
		"# Generated by libtool (hand-written NVSHMEM device link)" \
		"pic_object='.libs/gs_nvshmem_dlink.o'" \
		"non_pic_object=none" > $@

libneko_la_LIBADD += \
	gs/bcknd/device/cuda/gs_nvshmem_rdc.lo\
	gs/bcknd/device/cuda/gs_nvshmem_dlink.lo
else
libneko_la_SOURCES += \
	gs/bcknd/device/cuda/gs_nvshmem.cu
endif

endif

if ENABLE_OPENCL
libneko_la_SOURCES += \
	device/opencl/check.c\
	device/opencl/jit.c\
	device/opencl/buffer.c\
	math/bcknd/device/opencl/math.c\
	math/bcknd/device/opencl/schwarz.c\
	math/bcknd/device/opencl/tensor.c\
	math/bcknd/device/opencl/fdm.c\
	math/bcknd/device/opencl/mathops.c\
	math/bcknd/device/opencl/opr_dudxyz.c\
	math/bcknd/device/opencl/opr_cdtp.c\
	math/bcknd/device/opencl/opr_conv1.c\
	math/bcknd/device/opencl/opr_convect_scalar.c\
	math/bcknd/device/opencl/set_convect_rst.c\
	math/bcknd/device/opencl/opr_opgrad.c\
	math/bcknd/device/opencl/opr_lambda2.c\
	math/bcknd/device/opencl/opr_cfl.c\
	math/bcknd/device/opencl/opr_rotate_cyc.c\
	math/bcknd/device/opencl/ax_helm.c\
	math/bcknd/device/opencl/ax_helm_full.c\
	krylov/bcknd/device/opencl/pc_jacobi.c\
	gs/bcknd/device/opencl/gs.c\
	bc/bcknd/device/opencl/dirichlet.c\
	bc/bcknd/device/opencl/inflow.c\
	bc/bcknd/device/opencl/zero_dirichlet.c\
	bc/bcknd/device/opencl/symmetry_aligned.c\
	bc/bcknd/device/opencl/constrain_mixed_bc.c\
	bc/bcknd/device/opencl/coupled_vector_bc_projector.c\
	bc/bcknd/device/opencl/facet_normal.c\
	bc/bcknd/device/opencl/neumann.c\
	bc/bcknd/device/opencl/inhom_dirichlet.c\
	bc/bcknd/device/opencl/dong_outflow.c\
	wall_models/bcknd/device/opencl/cai_sagaut_model_ii.c\
	common/bcknd/device/opencl/rhs_maker.c\
	fluid/bcknd/device/opencl/pnpn_res.c\
	fluid/stress_formulation/bcknd/device/opencl/pnpn_stress_res.c\
	fluid/bcknd/device/opencl/compressible_res.c\
	fluid/bcknd/device/opencl/compressible_ops_compute_max_wave_speed.c\
	fluid/bcknd/device/opencl/compressible_ops_compute_entropy.c\
	fluid/bcknd/device/opencl/compressible_ops_update.c\
	fluid/regularization/bcknd/device/opencl/entropy_viscosity.c\
	scalar/bcknd/device/opencl/scalar_residual.c\
	sem/bcknd/device/opencl/coef.c\
	sem/bcknd/device/opencl/local_interpolation.c\
	filter/bcknd/device/opencl/mappings.c

%.cl.h: %.cl
	bash ./scripts/cltostring.sh $<
endif

if ENABLE_METAL
#
# Objective-C sources can't be added to libneko_la_SOURCES (automake
# would require AC_PROG_OBJC and libtool lacks an OBJC tag), instead
# wrap them with libtool via the .m.lo rule and add them to LIBADD
#
#
# The embedded shader library is generated at build time and
# requires the Metal shader compiler (enforced by configure),
# hence nodist (it would not be buildable on non-macOS dist hosts)
#
nodist_libneko_la_SOURCES = device/metal/neko_metallib.c

METAL_OBJS = \
	device/metal/metal.lo\
	math/bcknd/device/metal/math.lo\
	math/bcknd/device/metal/ax_helm.lo\
	math/bcknd/device/metal/ax_helm_full.lo\
	math/bcknd/device/metal/tensor.lo\
	math/bcknd/device/metal/opr_cfl.lo\
	math/bcknd/device/metal/opr_dudxyz.lo\
	math/bcknd/device/metal/opr_opgrad.lo\
	math/bcknd/device/metal/opr_cdtp.lo\
	math/bcknd/device/metal/opr_conv1.lo\
	math/bcknd/device/metal/set_convect_rst.lo\
	math/bcknd/device/metal/opr_convect_scalar.lo\
	math/bcknd/device/metal/opr_lambda2.lo\
	math/bcknd/device/metal/opr_rotate_cyc.lo\
	math/bcknd/device/metal/schwarz.lo\
	math/bcknd/device/metal/fdm.lo\
	multigrid/bcknd/device/metal/amg_cheby.lo\
	krylov/bcknd/device/metal/cheby.lo\
	krylov/bcknd/device/metal/pc_jacobi.lo\
	scalar/bcknd/device/metal/scalar_residual.lo\
	math/bcknd/device/metal/mathops.lo\
	common/bcknd/device/metal/rhs_maker.lo\
	fluid/bcknd/device/metal/compressible_res.lo\
	fluid/bcknd/device/metal/pnpn_res.lo\
	fluid/stress_formulation/bcknd/device/metal/pnpn_stress_res.lo\
	fluid/bcknd/device/metal/compressible_ops_compute_max_wave_speed.lo\
	fluid/bcknd/device/metal/compressible_ops_compute_entropy.lo\
	fluid/bcknd/device/metal/compressible_ops_update.lo\
	fluid/regularization/bcknd/device/metal/entropy_viscosity.lo\
	filter/bcknd/device/metal/mappings.lo\
	gs/bcknd/device/metal/gs.lo\
	sem/bcknd/device/metal/coef.lo\
	sem/bcknd/device/metal/local_interpolation.lo\
	bc/bcknd/device/metal/dirichlet.lo\
	bc/bcknd/device/metal/inflow.lo\
	bc/bcknd/device/metal/zero_dirichlet.lo\
	bc/bcknd/device/metal/symmetry.lo\
	bc/bcknd/device/metal/facet_normal.lo\
	bc/bcknd/device/metal/neumann.lo\
	bc/bcknd/device/metal/inhom_dirichlet.lo\
	bc/bcknd/device/metal/dong_outflow.lo\
	bc/bcknd/device/metal/coupled_vector_bc_projector.lo

libneko_la_LIBADD += $(METAL_OBJS)

dist_noinst_DATA = \
	device/metal/metal.m\
	math/bcknd/device/metal/math.m\
	math/bcknd/device/metal/ax_helm.m\
	math/bcknd/device/metal/ax_helm_full.m\
	math/bcknd/device/metal/tensor.m\
	math/bcknd/device/metal/opr_cfl.m\
	math/bcknd/device/metal/opr_dudxyz.m\
	math/bcknd/device/metal/opr_opgrad.m\
	math/bcknd/device/metal/opr_cdtp.m\
	math/bcknd/device/metal/opr_conv1.m\
	math/bcknd/device/metal/set_convect_rst.m\
	math/bcknd/device/metal/opr_convect_scalar.m\
	math/bcknd/device/metal/opr_lambda2.m\
	math/bcknd/device/metal/opr_rotate_cyc.m\
	math/bcknd/device/metal/schwarz.m\
	math/bcknd/device/metal/fdm.m\
	multigrid/bcknd/device/metal/amg_cheby.m\
	krylov/bcknd/device/metal/pc_jacobi.m\
	krylov/bcknd/device/metal/cheby.m\
	scalar/bcknd/device/metal/scalar_residual.m\
	math/bcknd/device/metal/mathops.m\
	common/bcknd/device/metal/rhs_maker.m\
	fluid/bcknd/device/metal/compressible_res.m\
	fluid/bcknd/device/metal/pnpn_res.m\
	fluid/stress_formulation/bcknd/device/metal/pnpn_stress_res.m\
	fluid/bcknd/device/metal/compressible_ops_compute_max_wave_speed.m\
	fluid/bcknd/device/metal/compressible_ops_compute_entropy.m\
	fluid/bcknd/device/metal/compressible_ops_update.m\
	fluid/regularization/bcknd/device/metal/entropy_viscosity.m\
	filter/bcknd/device/metal/mappings.m\
	gs/bcknd/device/metal/gs.m\
	sem/bcknd/device/metal/coef.m\
	sem/bcknd/device/metal/local_interpolation.m\
	bc/bcknd/device/metal/dirichlet.m\
	bc/bcknd/device/metal/inflow.m\
	bc/bcknd/device/metal/zero_dirichlet.m\
	bc/bcknd/device/metal/symmetry.m\
	bc/bcknd/device/metal/facet_normal.m\
	bc/bcknd/device/metal/neumann.m\
	bc/bcknd/device/metal/inhom_dirichlet.m\
	bc/bcknd/device/metal/dong_outflow.m\
	bc/bcknd/device/metal/coupled_vector_bc_projector.m

METAL_KERNELS = \
	math/bcknd/device/metal/math_kernel.metal\
	math/bcknd/device/metal/ax_helm_kernel.metal\
	math/bcknd/device/metal/ax_helm_full_kernel.metal\
	math/bcknd/device/metal/tensor_kernel.metal\
	math/bcknd/device/metal/cfl_kernel.metal\
	math/bcknd/device/metal/dudxyz_kernel.metal\
	math/bcknd/device/metal/opgrad_kernel.metal\
	math/bcknd/device/metal/cdtp_kernel.metal\
	math/bcknd/device/metal/conv1_kernel.metal\
	math/bcknd/device/metal/set_convect_rst_kernel.metal\
	math/bcknd/device/metal/convect_scalar_kernel.metal\
	math/bcknd/device/metal/lambda2_kernel.metal\
	math/bcknd/device/metal/rotate_kernel.metal\
	math/bcknd/device/metal/schwarz_kernel.metal\
	math/bcknd/device/metal/fdm_kernel.metal\
	multigrid/bcknd/device/metal/amg_cheby_kernel.metal\
	krylov/bcknd/device/metal/pc_jacobi_kernel.metal\
	krylov/bcknd/device/metal/cheby_kernel.metal\
	scalar/bcknd/device/metal/scalar_residual_update_kernel.metal\
	math/bcknd/device/metal/mathops_kernel.metal\
	common/bcknd/device/metal/rhs_maker_kernel.metal\
	fluid/bcknd/device/metal/compressible_res_kernel.metal\
	fluid/bcknd/device/metal/pnpn_res_kernel.metal\
	fluid/stress_formulation/bcknd/device/metal/prs_stress_res_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_compute_max_wave_speed_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_compute_entropy_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_update_kernel.metal\
	fluid/regularization/bcknd/device/metal/entropy_viscosity_kernel.metal\
	filter/bcknd/device/metal/mapping_kernels.metal\
	gs/bcknd/device/metal/gs_kernels.metal\
	sem/bcknd/device/metal/coef_kernel.metal\
	sem/bcknd/device/metal/local_interpolation_kernel.metal\
	bc/bcknd/device/metal/dirichlet_kernel.metal\
	bc/bcknd/device/metal/inflow_kernel.metal\
	bc/bcknd/device/metal/zero_dirichlet_kernel.metal\
	bc/bcknd/device/metal/symmetry_kernel.metal\
	bc/bcknd/device/metal/facet_normal_kernel.metal\
	bc/bcknd/device/metal/neumann_kernel.metal\
	bc/bcknd/device/metal/inhom_dirichlet_kernel.metal\
	bc/bcknd/device/metal/dong_outflow_kernel.metal\
	bc/bcknd/device/metal/coupled_vector_bc_projector_kernel.metal

device/metal/neko_kernels.metallib: $(METAL_KERNELS)
	$(MKDIR_P) device/metal
	@XCRUN@ -sdk macosx metal -O2 -o $@ $^

device/metal/neko_metallib.c: device/metal/neko_kernels.metallib
	cd device/metal && xxd -i neko_kernels.metallib | \
		$(SED) -e 's/^unsigned/const unsigned/' > neko_metallib.c

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
# With the previous build system, the correct method was:
#
#io/nek_adios2.o: io/nek_adios2.cpp
#	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $^

#libneko_la_LIBADD = \
#	io/nek_adios2.o
#
# Now, use libtools to wrap everything

io/nek_adios2.lo: io/nek_adios2.cpp
	$(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

libneko_la_LIBADD += \
        io/nek_adios2.lo

#neko_LDADD = \
#	io/nek_adios2.o

endif


bin_PROGRAMS = neko
neko_SOURCES = driver.f90
neko_LDADD= libneko.la

pkginclude_HEADERS = \
	$(wildcard *.mod)

CLEANFILES = *.mod\
	     *.smod

if ENABLE_CUDA
CLEANFILES += @top_builddir@/libtool_cuda
endif

if ENABLE_NVSHMEM
CLEANFILES += \
	gs/bcknd/device/cuda/gs_nvshmem_rdc.lo\
	gs/bcknd/device/cuda/gs_nvshmem_rdc.o\
	gs/bcknd/device/cuda/.libs/gs_nvshmem_rdc.o\
	gs/bcknd/device/cuda/gs_nvshmem_dlink.lo\
	gs/bcknd/device/cuda/.libs/gs_nvshmem_dlink.o
endif

if ENABLE_METAL
CLEANFILES += $(METAL_OBJS)\
	device/metal/neko_kernels.metallib\
	device/metal/neko_metallib.c
endif

if ENABLE_OPENCL
CLEANFILES += \
	math/bcknd/device/opencl/*.cl.h\
	gs/bcknd/device/opencl/*.cl.h\
	bc/bcknd/device/opencl/*.cl.h\
	krylov/bcknd/device/opencl/*.cl.h\
	fluid/bcknd/device/opencl/*.cl.h\
	fluid/stress_formulation/bcknd/device/opencl/*.cl.h\
	fluid/regularization/bcknd/device/opencl/*.cl.h\
	scalar/bcknd/device/opencl/*.cl.h\
	sem/bcknd/device/opencl/*.cl.h\
	source_terms/bcknd/device/opencl/*.cl.h\
	wall_models/bcknd/device/opencl/*.cl.h\
	filter/bcknd/device/opencl/*.cl.h\
	common/bcknd/device/opencl/*.cl.h
endif

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(neko_fortran_SOURCES) | $(SED) 's/\.o/\.lo/g'  > .depends
endif

if ENABLE_PYTHON
SUBDIRS = python
endif

# Note we need to list any headers etc here for the dist. script
EXTRA_DIST = \
    ale/bcknd/device/cuda/ale_kinematics_kernel.h\
	bc/bcknd/device/cuda/dirichlet_kernel.h\
	bc/bcknd/device/cuda/inflow_kernel.h\
	bc/bcknd/device/cuda/zero_dirichlet_kernel.h\
	bc/bcknd/device/cuda/symmetry_aligned_kernel.h\
	bc/bcknd/device/cuda/facet_normal_kernel.h\
	bc/bcknd/device/cuda/neumann_kernel.h\
	bc/bcknd/device/cuda/bc_utils.h\
	bc/bcknd/device/cuda/inhom_dirichlet_kernel.h\
	bc/bcknd/device/cuda/dong_outflow_kernel.h\
	bc/bcknd/device/cuda/constrain_mixed_bc_kernel.h\
	bc/bcknd/device/cuda/coupled_vector_bc_projector_kernel.h\
	gs/bcknd/device/cuda/gs_kernels.h\
	gs/bcknd/device/cuda/gs_nvshmem_kernels.h\
	math/bcknd/device/cuda/math_kernel.h\
	math/bcknd/device/cuda/tensor_kernel.h\
	math/bcknd/device/cuda/fdm_kernel.h\
	math/bcknd/device/cuda/mathops_kernel.h\
	math/bcknd/device/cuda/schwarz_kernel.h\
	math/bcknd/device/cuda/dudxyz_kernel.h\
	math/bcknd/device/cuda/cdtp_kernel.h\
	math/bcknd/device/cuda/conv1_kernel.h\
	math/bcknd/device/cuda/convect_scalar_kernel.h\
	math/bcknd/device/cuda/set_convect_rst_kernel.h\
	math/bcknd/device/cuda/opgrad_kernel.h\
	math/bcknd/device/cuda/lambda2_kernel.h\
	math/bcknd/device/cuda/cfl_kernel.h\
	math/bcknd/device/cuda/rotate_kernel.h\
	math/bcknd/device/cuda/ax_helm_kernel.h\
	math/bcknd/device/cuda/ax_helm_full_kernel.h\
	math/bcknd/device/cuda/dmma_kernel.h\
	math/bcknd/device/cuda/dmma_tma_kernel.h\
	math/bcknd/device/cuda/elem_block.h\
	math/bcknd/device/cuda/elem_block_tune.h\
	common/bcknd/device/cuda/sumab_kernel.h\
	common/bcknd/device/cuda/makeext_kernel.h\
	common/bcknd/device/cuda/makebdf_kernel.h\
	common/bcknd/device/cuda/makeoifs_kernel.h\
	common/bcknd/device/cuda/projection_kernel.h\
	fluid/bcknd/device/cuda/prs_res_kernel.h\
	fluid/bcknd/device/cuda/vel_res_update_kernel.h\
	fluid/bcknd/device/cuda/compressible_res_kernel.h\
	fluid/bcknd/device/cuda/compressible_ops_kernel.h\
	fluid/regularization/bcknd/device/cuda/entropy_viscosity_kernel.h\
	fluid/stress_formulation/bcknd/device/cuda/prs_stress_res_kernel.h\
	scalar/bcknd/device/cuda/scalar_residual_update_kernel.h\
	sem/bcknd/device/cuda/coef_kernel.h\
	lpt/bcknd/device/cuda/lpt_periodic_bc_kernel.h\
	lpt/bcknd/device/cuda/lpt_wall_collision_kernel.h\
	les/bcknd/device/cuda/sigma_nut_kernel.h\
	les/bcknd/device/cuda/vreman_nut_kernel.h\
	les/bcknd/device/cuda/smagorinsky_nut_kernel.h\
	les/bcknd/device/cuda/dynamic_smagorinsky_nut_kernel.h\
	les/bcknd/device/cuda/wale_nut_kernel.h\
	les/bcknd/device/cuda/deardorff_nut_kernel.h\
	wall_models/bcknd/device/cuda/wall_model_kernel.h\
	wall_models/bcknd/device/cuda/spalding_kernel.h\
	wall_models/bcknd/device/cuda/cai_sagaut_model_ii_kernel.h\
	wall_models/bcknd/device/cuda/rough_log_law_kernel.h\
    ale/bcknd/device/hip/ale_kinematics_kernel.h\
	wall_models/bcknd/device/cuda/most_kernel.h\
	wall_models/bcknd/device/cuda/richardson_kernel.h\
	bc/bcknd/device/hip/dirichlet_kernel.h\
	bc/bcknd/device/hip/inflow_kernel.h\
	bc/bcknd/device/hip/zero_dirichlet_kernel.h\
	bc/bcknd/device/hip/symmetry_aligned_kernel.h\
	bc/bcknd/device/hip/constrain_mixed_bc_kernel.h\
	bc/bcknd/device/hip/coupled_vector_bc_projector_kernel.h\
	bc/bcknd/device/hip/facet_normal_kernel.h\
	bc/bcknd/device/hip/neumann_kernel.h\
	bc/bcknd/device/hip/bc_utils.h\
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
	math/bcknd/device/hip/convect_scalar_kernel.h\
	math/bcknd/device/hip/set_convect_rst_kernel.h\
	math/bcknd/device/hip/opgrad_kernel.h\
	math/bcknd/device/hip/lambda2_kernel.h\
	math/bcknd/device/hip/cfl_kernel.h\
	math/bcknd/device/hip/rotate_kernel.h\
	math/bcknd/device/hip/ax_helm_kernel.h\
	math/bcknd/device/hip/ax_helm_full_kernel.h\
	math/bcknd/device/hip/mfma_kernel.h\
	math/bcknd/device/hip/elem_block.h\
	math/bcknd/device/hip/wave.h\
	math/bcknd/device/hip/elem_block_tune.h\
	common/bcknd/device/hip/sumab_kernel.h\
	common/bcknd/device/hip/makeext_kernel.h\
	common/bcknd/device/hip/makebdf_kernel.h\
	common/bcknd/device/hip/makeoifs_kernel.h\
	common/bcknd/device/hip/projection_kernel.h\
	fluid/bcknd/device/hip/prs_res_kernel.h\
	fluid/bcknd/device/hip/vel_res_update_kernel.h\
	fluid/bcknd/device/hip/compressible_res_kernel.h\
	fluid/bcknd/device/hip/compressible_ops_kernel.h\
	fluid/regularization/bcknd/device/hip/entropy_viscosity_kernel.h\
	fluid/stress_formulation/bcknd/device/hip/prs_stress_res_kernel.h\
	scalar/bcknd/device/hip/scalar_residual_update_kernel.h\
	krylov/bcknd/device/cuda/pipecg_kernel.h\
	krylov/bcknd/device/hip/pipecg_kernel.h\
	krylov/bcknd/device/cuda/fusedcg_kernel.h\
	krylov/bcknd/device/hip/fusedcg_kernel.h\
	krylov/bcknd/device/cuda/fusedcg_cpld_kernel.h\
	krylov/bcknd/device/hip/fusedcg_cpld_kernel.h\
	krylov/bcknd/device/cuda/gmres_kernel.h\
	krylov/bcknd/device/hip/gmres_kernel.h\
	sem/bcknd/device/hip/coef_kernel.h\
	lpt/bcknd/device/hip/lpt_periodic_bc_kernel.h\
	lpt/bcknd/device/hip/lpt_wall_collision_kernel.h\
	les/bcknd/device/hip/sigma_nut_kernel.h\
	les/bcknd/device/hip/vreman_nut_kernel.h\
	les/bcknd/device/hip/smagorinsky_nut_kernel.h\
	les/bcknd/device/hip/dynamic_smagorinsky_nut_kernel.h\
	les/bcknd/device/hip/wale_nut_kernel.h\
	les/bcknd/device/hip/deardorff_nut_kernel.h\
	wall_models/bcknd/device/hip/wall_model_kernel.h\
	wall_models/bcknd/device/hip/spalding_kernel.h\
	wall_models/bcknd/device/hip/cai_sagaut_model_ii_kernel.h\
	wall_models/bcknd/device/hip/rough_log_law_kernel.h\
	wall_models/bcknd/device/hip/most_kernel.h\
	wall_models/bcknd/device/hip/richardson_kernel.h\
	bc/bcknd/device/opencl/dirichlet_kernel.cl\
	bc/bcknd/device/opencl/inflow_kernel.cl\
	bc/bcknd/device/opencl/zero_dirichlet_kernel.cl\
	bc/bcknd/device/opencl/symmetry_aligned_kernel.cl\
	bc/bcknd/device/opencl/constrain_mixed_bc_kernel.cl\
	bc/bcknd/device/opencl/coupled_vector_bc_projector_kernel.cl\
	bc/bcknd/device/opencl/facet_normal_kernel.cl\
	bc/bcknd/device/opencl/neumann_kernel.cl\
	bc/bcknd/device/opencl/bc_utils.h\
	bc/bcknd/device/opencl/inhom_dirichlet_kernel.cl\
	bc/bcknd/device/opencl/dong_outflow_kernel.cl\
	wall_models/bcknd/device/opencl/cai_sagaut_model_ii_kernel.cl\
	math/bcknd/device/opencl/math_kernel.cl\
	math/bcknd/device/opencl/schwarz_kernel.cl\
	math/bcknd/device/opencl/tensor_kernel.cl\
	math/bcknd/device/opencl/fdm_kernel.cl\
	math/bcknd/device/opencl/mathops_kernel.cl\
	math/bcknd/device/opencl/dudxyz_kernel.cl\
	math/bcknd/device/opencl/ax_helm_kernel.cl\
	math/bcknd/device/opencl/ax_helm_full_kernel.cl\
	math/bcknd/device/opencl/cdtp_kernel.cl\
	math/bcknd/device/opencl/opgrad_kernel.cl\
	math/bcknd/device/opencl/conv1_kernel.cl\
	math/bcknd/device/opencl/convect_scalar_kernel.cl\
	math/bcknd/device/opencl/set_convect_rst_kernel.cl\
	math/bcknd/device/opencl/lambda2_kernel.cl\
	math/bcknd/device/opencl/cfl_kernel.cl\
	math/bcknd/device/opencl/rotate_kernel.cl\
	common/bcknd/device/opencl/rhs_maker_kernel.cl\
	fluid/bcknd/device/opencl/pnpn_res_kernel.cl\
	fluid/bcknd/device/opencl/compressible_res_kernel.cl\
	fluid/stress_formulation/bcknd/device/opencl/prs_stress_res_kernel.cl\
	fluid/bcknd/device/opencl/compressible_ops_compute_max_wave_speed_kernel.cl\
	fluid/bcknd/device/opencl/compressible_ops_compute_entropy_kernel.cl\
	fluid/bcknd/device/opencl/compressible_ops_update_kernel.cl\
	fluid/regularization/bcknd/device/opencl/entropy_viscosity_kernel.cl\
	scalar/bcknd/device/opencl/scalar_residual_kernel.cl\
	gs/bcknd/device/opencl/gs_kernels.cl\
	krylov/bcknd/device/opencl/jacobi_kernel.cl\
	sem/bcknd/device/opencl/coef_kernel.cl\
	sem/bcknd/device/opencl/local_interpolation_kernel.cl\
	math/bcknd/device/device_mpi_reduce.h\
	math/bcknd/device/device_mpi_op.h\
	math/bcknd/device/device_nccl_reduce.h\
	math/bcknd/device/device_nccl_op.h\
	filter/bcknd/device/cuda/mapping_kernels.h\
	filter/bcknd/device/hip/mapping_kernels.h\
	filter/bcknd/device/opencl/mapping_kernels.cl\
	source_terms/bcknd/device/cuda/gradient_jump_penalty_kernel.h\
	source_terms/bcknd/device/hip/gradient_jump_penalty_kernel.h\
	math/bcknd/device/metal/math_kernel.metal\
	math/bcknd/device/metal/ax_helm_kernel.metal\
	math/bcknd/device/metal/ax_helm_full_kernel.metal\
	math/bcknd/device/metal/tensor_kernel.metal\
	math/bcknd/device/metal/cfl_kernel.metal\
	math/bcknd/device/metal/dudxyz_kernel.metal\
	math/bcknd/device/metal/opgrad_kernel.metal\
	math/bcknd/device/metal/cdtp_kernel.metal\
	math/bcknd/device/metal/conv1_kernel.metal\
	math/bcknd/device/metal/set_convect_rst_kernel.metal\
	math/bcknd/device/metal/convect_scalar_kernel.metal\
	math/bcknd/device/metal/lambda2_kernel.metal\
	math/bcknd/device/metal/rotate_kernel.metal\
	math/bcknd/device/metal/schwarz_kernel.metal\
	math/bcknd/device/metal/fdm_kernel.metal\
	multigrid/bcknd/device/metal/amg_cheby_kernel.metal\
	krylov/bcknd/device/metal/cheby_kernel.metal\
	krylov/bcknd/device/metal/pc_jacobi_kernel.metal\
	scalar/bcknd/device/metal/scalar_residual_update_kernel.metal\
	math/bcknd/device/metal/mathops_kernel.metal\
	common/bcknd/device/metal/rhs_maker_kernel.metal\
	fluid/bcknd/device/metal/compressible_res_kernel.metal\
	fluid/bcknd/device/metal/pnpn_res_kernel.metal\
	fluid/stress_formulation/bcknd/device/metal/prs_stress_res_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_compute_max_wave_speed_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_compute_entropy_kernel.metal\
	fluid/bcknd/device/metal/compressible_ops_update_kernel.metal\
	fluid/regularization/bcknd/device/metal/entropy_viscosity_kernel.metal\
	filter/bcknd/device/metal/mapping_kernels.metal\
	gs/bcknd/device/metal/gs_kernels.metal\
	sem/bcknd/device/metal/coef_kernel.metal\
	sem/bcknd/device/metal/local_interpolation_kernel.metal\
	bc/bcknd/device/metal/dirichlet_kernel.metal\
	bc/bcknd/device/metal/inflow_kernel.metal\
	bc/bcknd/device/metal/zero_dirichlet_kernel.metal\
	bc/bcknd/device/metal/symmetry_kernel.metal\
	bc/bcknd/device/metal/facet_normal_kernel.metal\
	bc/bcknd/device/metal/neumann_kernel.metal\
	bc/bcknd/device/metal/inhom_dirichlet_kernel.metal\
	bc/bcknd/device/metal/dong_outflow_kernel.metal\
	bc/bcknd/device/metal/coupled_vector_bc_projector_kernel.metal\
	bc/bcknd/device/metal/bc_utils.h\
	bc/bcknd/device/metal/bc_shader_utils.h\
	device/metal/check.h\
	device/device_config.h\
	device/opencl/jit.h\
	device/opencl/prgm_lib.h\
	device/opencl/check.h\
	device/hip/check.h\
	device/cuda/check.h\
	device/cuda/buffer.h\
	device/hip/buffer.h\
	device/hip/unified.h\
	device/opencl/buffer.h\
	comm/comm_nccl.h\
	comm/comm.h\
	common/neko_log.h


# Expose the listed headers in the package
nobase_pkginclude_HEADERS = ${EXTRA_DIST}

# Expose API header
include_HEADERS = api/neko.h

# Fortran deps. (updated with makedepf90 via target 'depend')
include	.depends

# Device deps. (manually updated)
include .depends_device

# TurboNeko deps. (manually updated)
include .depends_turboneko

# CUDA Build rules
.cu.lo:
	$(SHELL) @top_builddir@/libtool_cuda --tag=CC --mode=compile $(NVCC) @CUDA_ARCH@ @CUDA_CFLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<

@top_builddir@/libtool_cuda: @top_builddir@/libtool
	$(SED) 's/pic_flag=.*/pic_flag="-Xcompiler -fPIC"/g' @top_builddir@/libtool > $@

# HIP Build rules
.hip.lo:
	$(LIBTOOL) --tag=CC --mode=compile $(HIPCC) @HIP_HIPCC_FLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<

# Metal (Objective-C) Build rules
.m.lo:
	$(MKDIR_P) $(@D)
	$(LIBTOOL) --tag=CC --mode=compile $(OBJC) $(DEFS) $(CPPFLAGS) $(OBJCFLAGS) @METAL_OBJCFLAGS@ -I@top_builddir@/src -I@top_srcdir@/src -o $@ -c $<
