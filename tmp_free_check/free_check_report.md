# Free-check report

Scanned files: 515

## src/adt/htable.f90
- Type `h_tuple_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `htable_iter_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/bc/bc.f90
- Type `bc_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `bc_alloc_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/case.f90
- Type `case_t`: free missing deallocate for: output_directory
## src/common/checkpoint.f90
- Type `chkp_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/common/mask.f90
- Type `mask_t`: free missing nullify for: mask_d
## src/common/structs.f90
- Type `array_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/field/field.f90
- Type `field_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/field/field_series.f90
- Type `field_series_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/fluid/fluid_scheme_base.f90
- Type `fluid_scheme_base_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/fluid/fluid_scheme_compressible.f90
- Type `fluid_scheme_compressible_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/fluid/fluid_scheme_incompressible.f90
- Type `fluid_scheme_incompressible_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/global_interpolation/global_interpolation.f90
- Type `global_interpolation_t`: free missing deallocate for: xyz, rst, pe_owner, el_owner0, el_owner0_local, rst_local, xyz_local
## src/global_interpolation/global_interpolation_comm.f90
- Type `glb_intrp_comm_mpi_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `glb_intrp_comm_t`: free missing deallocate for: send_dof, recv_dof, send_pe, recv_pe
## src/gs/bcknd/device/gs_device_mpi.F90
- Type `gs_device_mpi_t`: free missing deallocate for: event
## src/gs/bcknd/device/gs_device_nccl.F90
- Type `gs_device_nccl_t`: free missing deallocate for: event
## src/gs/bcknd/device/gs_device_shmem.F90
- Type `gs_device_shmem_buf_t`: free missing deallocate for: remote_offset
- Type `gs_device_shmem_t`: free missing deallocate for: event, notifyDone, notifyReady
## src/gs/gs_comm.f90
- Type `gs_comm_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/gs/gs_mpi.f90
- Type `gs_comm_mpi_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/buffer/buffer_1d.F90
- Type `buffer_1d_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/buffer/buffer_4d.F90
- Type `buffer_4d_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/buffer/buffer_4d_npar.F90
- Type `buffer_4d_npar_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/chkp_file.f90
- Type `chkp_file_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/chkp_output.f90
- Type `chkp_output_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/data_streamer.F90
- Type `data_streamer_t`: free missing deallocate for: lglel
## src/io/fluid_stats_output.f90
- Type `fluid_stats_output_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/format/rea.f90
- Type `rea_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/mean_field_output.f90
- Type `mean_field_output_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/output.f90
- Type `output_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/io/scalar_stats_output.f90
- Type `scalar_stats_output_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/krylov/bcknd/cpu/cheby.f90
- Type `cheby_t`: free missing nullify for: schwarz
## src/krylov/bcknd/device/cheby_device.F90
- Type `cheby_device_t`: free missing nullify for: schwarz
## src/krylov/bcknd/device/fusedcg_cpld_device.F90
- Type `fusedcg_cpld_device_t`: free missing deallocate for: p1_d, p2_d, p3_d
## src/krylov/bcknd/device/fusedcg_device.F90
- Type `fusedcg_device_t`: free missing deallocate for: p_d
## src/krylov/bcknd/device/gmres_device.F90
- Type `gmres_device_t`: free missing deallocate for: z_d, h_d, v_d
## src/krylov/bcknd/device/pipecg_device.F90
- Type `pipecg_device_t`: free missing deallocate for: u_d
## src/krylov/krylov.f90
- Type `ksp_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/krylov/pc_hsmg.f90
- Type `multigrid_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `hsmg_t`: free missing deallocate for: amg_solver, pc_crs
- Type `hsmg_t`: free missing nullify for: msh
## src/math/matrix.f90
- Type `matrix_t`: free missing nullify for: x_d
- Type `matrix_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/math/vector.f90
- Type `vector_ptr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/mesh/aabb_tree.f90
- Type `aabb_tree_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/mesh/mesh.f90
- Type `mesh_element_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/mesh/octree.f90
- Type `octree_t`: free missing nullify for: root
## src/mesh/point.f90
- Type `point_ptr`: MISSING type-bound `free` (has alloc/pointer components)
## src/mesh/point_zone.f90
- Type `point_zone_wrapper_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `point_zone_pointer_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/multigrid/phmg.f90
- Type `phmg_lvl_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `phmg_hrchy_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `phmg_t`: free missing deallocate for: ax, intrp
- Type `phmg_t`: free missing nullify for: msh
## src/multigrid/tree_amg_multigrid.f90
- Type `tamg_wrk_t`: MISSING type-bound `free` (has alloc/pointer components)
- Type `tamg_solver_t`: free missing deallocate for: jsmoo
## src/scalar/scalar_scheme.f90
- Type `scalar_scheme_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/sem/cpr.f90
- Type `cpr_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/sem/map_2d.f90
- Type `map_2d_t`: MISSING type-bound `free` (has alloc/pointer components)
## src/simulation_components/simulation_component.f90
- Type `simulation_component_wrapper_t`: MISSING type-bound `free` (has alloc/pointer components)
