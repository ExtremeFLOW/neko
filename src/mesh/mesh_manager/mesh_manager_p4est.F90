! Copyright (c) 2019-2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Main interface for data exchange and manipulation between p4est and neko
module mesh_manager_p4est
  use mpi_f08
  use num_types, only : i4, i8, rp, dp
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use profiler, only : profiler_start_region, profiler_end_region
  use stack, only : stack_i4_t
  use hex, only : hex_t
  use mesh, only : mesh_t
  use manager_mesh, only : manager_mesh_t
  use manager_geom_p4est, only : manager_geom_node_ind_p4est_t, &
       manager_geom_p4est_t
  use manager_conn_p4est, only : manager_conn_obj_p4est_t, manager_conn_p4est_t
  use manager_mesh_p4est, only : manager_mesh_p4est_t
  use mesh_manager, only : mesh_manager_t
  use mesh_manager_transfer_p4est, only : mesh_manager_transfer_p4est_t

  implicit none

  private

  ! Neko is 3D only
#define N_DIM 3
#if N_DIM == 2
#undef P4_TO_P8
#else
#define P4_TO_P8
#endif
  ! refinement flag definition
#include "mesh/mesh_manager/amr.h"

  !> p4est mesh manager
  type, public, extends(mesh_manager_t) :: mesh_manager_p4est_t
     !> Tree file
     character(len=:), allocatable :: tree_file
     !> Non-periodic connectivity file
     character(len=:), allocatable :: cnn_np_file
     !> Is domain periodic
     logical :: is_periodic
     !> Max refinement level
     integer :: ref_level_max
     !> Log level for p4est
     integer :: log_level
   contains
     !> Start p4est
     procedure, pass(this) :: start => p4est_start
     !> Stop p4est
     procedure, pass(this) :: stop => p4est_stop
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => p4est_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => p4est_free
     !> Import mesh data into current type
     procedure, pass(this) :: import => p4est_import
     !> Apply data from nmsh file to mesh manager structures
     procedure, pass(this) :: mesh_file_apply => p4est_mesh_file_apply
     !> Perform refinement/coarsening on the mesh manager side
     procedure, pass(this) :: refine => p4est_refine
     !> Construct neko mesh type based on mesh manager data
     procedure, pass(this) :: mesh_construct => p4est_mesh_construct
#ifdef HAVE_P4EST
     !> The costrucructor from type components.
     procedure, pass(this) :: init_from_component => &
          p4est_init_from_components
#endif
  end type mesh_manager_p4est_t

#ifdef HAVE_P4EST
  ! set of interfaces to call p4est functions
  interface
     subroutine wp4est_init(fmpicomm, catch_signals, print_backtrace, &
          log_threshold) bind(c, name = 'wp4est_init')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: fmpicomm, catch_signals, &
            print_backtrace, log_threshold
     end subroutine wp4est_init

     subroutine wp4est_finalize() bind(c, name = 'wp4est_finalize')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_finalize

     subroutine wp4est_is_initialized(is_init) &
          bind(c, name = 'wp4est_is_initialized')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_init
     end subroutine wp4est_is_initialized

#ifdef P4_TO_P8
     subroutine wp4est_cnn_new(num_vertices, num_trees, num_edges, num_corners,&
          vertices, tree_to_vertex, tree_to_tree, tree_to_face, tree_to_edge,&
          ett_offset, edge_to_tree, edge_to_edge, tree_to_corner, ctt_offset,&
          corner_to_tree, corner_to_corner) bind(c, name = 'wp4est_cnn_new')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: num_vertices, num_trees, num_edges, num_corners
       type(c_ptr), value :: vertices, tree_to_vertex, tree_to_tree, &
            tree_to_face, tree_to_edge, ett_offset, edge_to_tree, edge_to_edge,&
            tree_to_corner, ctt_offset, corner_to_tree, corner_to_corner
     end subroutine wp4est_cnn_new
#else
     subroutine wp4est_cnn_new(num_vertices, num_trees, num_corners,&
          vertices, tree_to_vertex, tree_to_tree, tree_to_face, tree_to_corner,&
          ctt_offset,  corner_to_tree, corner_to_corner) &
          bind(c, name = 'wp4est_cnn_new')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: num_vertices, num_trees, num_corners
       type(c_ptr), value :: vertices, tree_to_vertex, tree_to_tree, &
            tree_to_face, tree_to_corner, ctt_offset, corner_to_tree,&
            corner_to_corner
     end subroutine wp4est_cnn_new
#endif

     subroutine wp4est_cnn_brick(nx, ny, nz) bind(c, name = 'wp4est_cnn_brick')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: nx, ny, nz
     end subroutine wp4est_cnn_brick

     subroutine wp4est_cnn_unit_cube() bind(c, name = 'wp4est_cnn_unit_cube')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_unit_cube

     subroutine wp4est_cnn_unit_cube_periodic() &
          bind(c, name = 'wp4est_cnn_unit_cube_periodic')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_unit_cube_periodic

     subroutine wp4est_cnn_rot_cubes() bind(c, name = 'wp4est_cnn_rot_cubes')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_rot_cubes

     subroutine wp4est_cnn_del() bind(c, name = 'wp4est_cnn_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_del

     subroutine wp4est_cnn_valid(is_valid) bind(c, name = 'wp4est_cnn_valid')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_valid
     end subroutine wp4est_cnn_valid

     subroutine wp4est_cnn_np_valid(is_valid) &
          bind(c, name = 'wp4est_cnn_np_valid')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_valid
     end subroutine wp4est_cnn_np_valid

     subroutine wp4est_cnn_attr(enable_tree_attr) &
          bind(c, name = 'wp4est_cnn_attr')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: enable_tree_attr
     end subroutine wp4est_cnn_attr

     subroutine wp4est_cnn_complete() bind(c, name = 'wp4est_cnn_complete')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_complete

     subroutine wp4est_cnn_save(filename) &
          bind(c, name = 'wp4est_cnn_save')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_cnn_save

     subroutine wp4est_cnn_load(filename) &
          bind(c, name = 'wp4est_cnn_load')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_cnn_load

     subroutine wp4est_cnn_np_load(filename) &
          bind(c, name = 'wp4est_cnn_np_load')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_cnn_np_load

     subroutine wp4est_geom_del() bind(c, name = 'wp4est_geom_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_geom_del

     subroutine wp4est_tree_new() bind(c, name = 'wp4est_tree_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_new

     subroutine wp4est_tree_del() bind(c, name = 'wp4est_tree_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_del

     subroutine wp4est_tree_valid(is_valid) bind(c, name = 'wp4est_tree_valid')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_valid
     end subroutine wp4est_tree_valid

     subroutine wp4est_tree_save(filename) &
          bind(c, name = 'wp4est_tree_save')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_tree_save

     subroutine wp4est_tree_load(filename) &
          bind(c, name = 'wp4est_tree_load')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_tree_load

     subroutine wp4est_tree_cnn_swap() bind(c, name = 'wp4est_tree_cnn_swap')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_cnn_swap

     subroutine wp4est_tree_cnn_swap_back() &
          bind(c, name = 'wp4est_tree_cnn_swap_back')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_cnn_swap_back

     subroutine wp4est_ghost_new() bind(c, name = 'wp4est_ghost_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_ghost_new

     subroutine wp4est_ghost_del() bind(c, name = 'wp4est_ghost_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_ghost_del

     subroutine wp4est_mesh_new() bind(c, name = 'wp4est_mesh_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_mesh_new

     subroutine wp4est_mesh_del() bind(c, name = 'wp4est_mesh_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_mesh_del

     subroutine wp4est_nodes_new() bind(c, name = 'wp4est_nodes_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_nodes_new

     subroutine wp4est_nodes_del() bind(c, name = 'wp4est_nodes_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_nodes_del

     subroutine wp4est_lnodes_new(degree) bind(c, name = 'wp4est_lnodes_new')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: degree
     end subroutine wp4est_lnodes_new

     subroutine wp4est_lnodes_edge(npts) bind(c, name = 'wp4est_lnodes_edge')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: npts
     end subroutine wp4est_lnodes_edge

     subroutine wp4est_lnodes_del() bind(c, name = 'wp4est_lnodes_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_lnodes_del

     subroutine wp4est_part() bind(c, name = 'wp4est_part')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_part

     subroutine wp4est_elm_ini_dat(gidx, imsh, igrp, crv, bc) &
          bind(c, name = 'wp4est_elm_ini_dat')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: gidx, imsh, igrp, crv, bc
     end subroutine wp4est_elm_ini_dat

     subroutine wp4est_bc_check() bind(c, name = 'wp4est_bc_check')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_bc_check

     subroutine wp4est_msh_get_size(mdim, nelgt, nelgto, nelt, nmsh, ngrp, &
          maxl) bind(c, name = 'wp4est_msh_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: mdim, nmsh, ngrp, maxl
       integer(c_int32_t) :: nelt
       integer(c_int64_t) :: nelgt, nelgto
     end subroutine wp4est_msh_get_size

     subroutine wp4est_nds_get_size(nowin, nowsh, oowin, nin, nhf, nhe) &
          bind(c, name = 'wp4est_nds_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: nowin, nowsh, oowin, nin, nhf, nhe
     end subroutine wp4est_nds_get_size

     subroutine wp4est_nds_get_ind(nglid, nown, ncoord)&
          bind(c, name = 'wp4est_nds_get_ind')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: nglid, nown, ncoord
     end subroutine wp4est_nds_get_ind

     subroutine wp4est_nds_get_hfc(depend, ncoord)&
          bind(c, name = 'wp4est_nds_get_hfc')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: depend, ncoord
     end subroutine wp4est_nds_get_hfc

     subroutine wp4est_nds_get_hed(depend, ncoord)&
          bind(c, name = 'wp4est_nds_get_hed')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: depend, ncoord
     end subroutine wp4est_nds_get_hed

     subroutine wp4est_nds_get_vmap(vmap)&
          bind(c, name = 'wp4est_nds_get_vmap')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: vmap
     end subroutine wp4est_nds_get_vmap

     subroutine wp4est_elm_get_dat(gidx, level, igrp, crv, bc, falg) &
          bind(c, name = 'wp4est_elm_get_dat')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: gidx, level, igrp, crv, bc, falg
     end subroutine wp4est_elm_get_dat

     subroutine wp4est_elm_get_lnode(lnnum, lnown, lnoff, lnodes) &
          bind(c, name = 'wp4est_elm_get_lnode')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: lnnum, lnown
       integer(c_int64_t) :: lnoff
       type(c_ptr), value :: lnodes
     end subroutine wp4est_elm_get_lnode

     subroutine wp4est_sharers_get_size(nrank, nshare) &
          bind(c, name = 'wp4est_sharers_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: nrank, nshare
     end subroutine wp4est_sharers_get_size

     subroutine wp4est_sharers_get_ind(nglid, lrank, loff, lshare) &
          bind(c, name = 'wp4est_sharers_get_ind')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: nglid, lrank, loff, lshare
     end subroutine wp4est_sharers_get_ind

     subroutine wp4est_hang_get_info(hang_elm, hang_fsc, hang_edg) &
          bind(c, name = 'wp4est_hang_get_info')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: hang_elm, hang_fsc, hang_edg
     end subroutine wp4est_hang_get_info

     subroutine wp4est_fml_get_info(family, nelf) &
          bind(c, name = 'wp4est_fml_get_info')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: family
       integer(c_int) :: nelf
     end subroutine wp4est_fml_get_info

     subroutine wp4est_refine(max_level) &
          bind(c, name = 'wp4est_refine')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: max_level
     end subroutine wp4est_refine

     subroutine wp4est_coarsen() &
          bind(c, name = 'wp4est_coarsen')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_coarsen

     subroutine wp4est_balance() &
          bind(c, name = 'wp4est_balance')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_balance

     subroutine wp4est_tree_compare_copy(quad_data) &
          bind(c, name = 'wp4est_tree_compare_copy')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_compare_copy

     subroutine wp4est_tree_compare_del() &
          bind(c, name = 'wp4est_tree_compare_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_compare_del

     subroutine wp4est_tree_compare_check(check, quad_data) &
          bind(c, name = 'wp4est_tree_compare_check')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: check
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_compare_check

     subroutine wp4est_refm_put(ref_mark) &
          bind(c, name = 'wp4est_refm_put')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: ref_mark
     end subroutine wp4est_refm_put

     subroutine wp4est_egmap_put(el_gnum, el_lnum, el_nid) &
          bind(c, name = 'wp4est_egmap_put')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: el_gnum, el_lnum, el_nid
     end subroutine wp4est_egmap_put

     subroutine wp4est_msh_get_hst_size(map_nr, rfn_nr, crs_nr) &
          bind(c, name = 'wp4est_msh_get_hst_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: map_nr, rfn_nr, crs_nr
     end subroutine wp4est_msh_get_hst_size

     subroutine wp4est_msh_get_hst(elgl_mapg, elgl_map, elgl_rfng, elgl_rfn, &
          elgl_crsg, elgl_crs) &
          bind(c, name = 'wp4est_msh_get_hst')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: elgl_mapg, elgl_map, elgl_rfng, elgl_rfn, &
            elgl_crsg, elgl_crs
     end subroutine wp4est_msh_get_hst

     subroutine wp4est_vtk_write(filename) &
          bind(c, name = 'wp4est_vtk_write')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_vtk_write

  end interface
#endif

contains

#ifdef HAVE_P4EST
  !> Start p4est
  !! @param[out]  ierr  error flag
  subroutine p4est_start(this, json, ierr)
    class(mesh_manager_p4est_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(out) :: ierr
    integer :: catch_signals, print_backtrace, log_level
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_start, t_end

    t_start = MPI_WTIME()
    write(log_buf, '(a)') 'Starting p4est.'
    call neko_log%message(log_buf, NEKO_LOG_INFO)

    call wp4est_is_initialized(log_level)
    if (log_level == 0) then
       ! done by neko, so not needed here
       catch_signals = 0
       print_backtrace = 0

       ! p4est has internal logging system with the following levels
       ! DEFAULT   (-1)   Selects the SC default threshold.
       ! ALWAYS      0    Log absolutely everything.
       ! TRACE       1    Prefix file and line number.
       ! DEBUG       2    Any information on the internal state.
       ! VERBOSE     3    Information on conditions, decisions.
       ! INFO        4    Most relevant things a function is doing.
       ! STATISTICS  5    Important for consistency/performance.
       ! PRODUCTION  6    A few lines at most for a major api function.
       ! ESSENTIAL   7    Log a few lines max (version info) per program.
       ! ERROR       8    Log errors only.  This is suggested over SILENT.
       ! SILENT      9    Never log anything.  Instead suggesting ERROR.
       call json_get_or_default(json, 'log_level', log_level, 8)

       if ((log_level >= 0) .and. (log_level <= 9)) then
          this%log_level = log_level
       else
          call neko_warning('The p4est log level out of bounds [0..9] ' // &
               'with 0 beeing the most verbose. ' // &
               'Resetting to 8 (error only).')
          this%log_level = 8
       end if

       ! start p4est
       call wp4est_init(NEKO_COMM%mpi_val, catch_signals, print_backtrace, &
            this%log_level)
    else
       call neko_error('p4est already initialised; unknown MPI communicator')
    end if

    call wp4est_is_initialized(log_level)
    if (log_level == 0) then
       call neko_error('Failure starting p4est')
    else
       this%ifstarted = .true.
       ierr = 0
    end if

    t_end = MPI_WTIME()
    write(log_buf, '(A,F9.6)') 'Mesh manager starting time (s): ', &
         t_end - t_start
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

  end subroutine p4est_start

  !> Stop p4est
  !! @param[out]  ierr  error flag
  subroutine p4est_stop(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    ! stop p4est
    ! There is some memory issue here, so I comment it for now.
    ! As it is one of the last operations in the code I leave the investigation
    ! for the future.
    !call wp4est_finalize()
    this%ifstarted = .false.

  end subroutine p4est_stop

  !> The common constructor using a JSON object.
  !! @param json       The JSON object.
  !! @param type_name  Manager type name
  subroutine p4est_init_from_json(this, json)
    class(mesh_manager_p4est_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: tree_file, cnn_file
    integer :: ref_level_max

    ! Extract runtime parameters
    ! tree_file is mandatory
    call json_get_or_default(json, 'tree_file', tree_file, 'no tree')
    if (trim(tree_file) .eq. 'no tree') then
       call neko_error('The tree_file keyword could not be found in the .' // &
            'case file. Often caused by incorrectly formatted json.')
    end if

    ! check if there is any non-periodic connectivity file
    call json_get_or_default(json, 'connectivity_file', cnn_file, 'no cnn')

    ! p4est supports AMR; get maximum allowed refinement level
    ! Notice, p4est has internal max refinement level (P8EST_QMAXLEVEL
    ! defined in p8est.h) that cannot be exceeded. Moreover, negative
    ! value of ref_level_max will automatically set the restriction to
    ! P8EST_QMAXLEVEL. By default we set no refinement.
    call json_get_or_default(json, "ref_level_max", ref_level_max, 0)

    call p4est_init_from_components(this, tree_file, cnn_file, ref_level_max)

    if (allocated(tree_file)) deallocate(tree_file)
    if (allocated(cnn_file)) deallocate(cnn_file)

  end subroutine p4est_init_from_json

  !> Destructor.
  subroutine p4est_free(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(a)') 'Finalising p4est.'
    call neko_log%message(log_buf, NEKO_LOG_INFO)

    ! clean the memory in C side
    call wp4est_lnodes_del()
    call wp4est_nodes_del()
    call wp4est_mesh_del()
    call wp4est_ghost_del()
    call wp4est_tree_compare_del()
    call wp4est_tree_del()
    call wp4est_geom_del()
    call wp4est_cnn_del()

    call this%free_base()

    this%is_periodic = .false.
    this%ref_level_max = 0
    this%log_level = 0

    if (allocated(this%tree_file)) deallocate(this%tree_file)
    if (allocated(this%cnn_np_file)) deallocate(this%cnn_np_file)

  end subroutine p4est_free

  !> Import mesh data
  subroutine p4est_import(this)
    class(mesh_manager_p4est_t), intent(inout) :: this
    type(manager_mesh_p4est_t) :: mesh_new

    ! allocate types and import data from p4est
    call mesh_new%init()
    call mesh_new%free_data()
    call p4est_import_data(mesh_new, this%is_periodic)

    ! fill mesh information
    call this%mesh%init_type(mesh_new)
    call mesh_new%free()

    ! get new element distribution; get partitioning
    call this%transfer%elem_dist_construct(this%mesh)

  end subroutine p4est_import

  !> The constructor from type components.
  subroutine p4est_init_from_components(this, tree_file, cnn_file, &
       ref_level_max)
    class(mesh_manager_p4est_t), intent(inout) :: this
    character(len=*), intent(in) :: tree_file, cnn_file
    integer, intent(in) :: ref_level_max
    integer(i4) :: is_valid, ierr
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_start, t_end

    t_start = MPI_WTIME()
    ! allocate mesh types
    if (allocated(this%mesh))then
       call this%mesh%free()
       deallocate(this%mesh)
    end if
    allocate(manager_mesh_p4est_t::this%mesh)
    call this%mesh%init()
    call this%mesh%free_data()

    ! is p4est started?
    if (.not.this%ifstarted) call neko_error('p4est not started')

    this%tree_file = trim(tree_file)

    ! If domain is periodic special treatment of independent nodes is needed
    if (trim(cnn_file) .eq. 'no cnn') then
       this%is_periodic = .false.
    else
       this%cnn_np_file = trim(cnn_file)
       this%is_periodic = .true.
    end if

    ! AMR related stuff
    this%ref_level_max = ref_level_max

    ! for testing
!    call wp4est_cnn_rot_cubes()
!    call wp4est_cnn_save('rot_cubes.cnn')
!    call wp4est_tree_new()
!    call wp4est_tree_save('rot_cubes.tree')
!    call wp4est_cnn_brick(2, 2, 2)
!    call wp4est_cnn_unit_cube_periodic()
!    call wp4est_cnn_save('test.cnn')
!    call wp4est_cnn_load('unit_cube_periodic.cnn')
!    call wp4est_cnn_valid(is_valid)
!    write(*,*) 'TESTcnn', is_valid
!    call wp4est_cnn_np_load('unit_cube.cnn')
!    call wp4est_cnn_np_valid(is_valid)
!    write(*,*) 'TESTcnn_np', is_valid
!    call wp4est_tree_new()
!    call wp4est_tree_save('unit_cube_periodic.tree')
!    call wp4est_tree_load('512.tree')
!    call wp4est_vtk_write('test')
!    call wp4est_tree_del()
!    call wp4est_cnn_del()

    write(log_buf, '(a)') 'Reading p4est tree/connectivity data.'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)
    ! read the tree file
    call wp4est_tree_load(trim(this%tree_file))
    call wp4est_tree_valid(is_valid)
    if (is_valid == 0) call neko_error('Invalid p4est tree')
    call wp4est_cnn_valid(is_valid)
    if (is_valid == 0) call neko_error('Invalid p4est connectivity')
    ! for periodic domains load non-periodic connectivity
    if (this%is_periodic) then
       call wp4est_cnn_np_load(trim(this%cnn_np_file))
       call wp4est_cnn_np_valid(is_valid)
       if (is_valid == 0) call neko_error('Invalid non-periodic connectivity')
    end if
    ! perform partitioning on p4est side
    call wp4est_part()

!    call wp4est_cnn_save('test.cnn')
!    call wp4est_tree_save('test.tree')

    !call MPI_Barrier(NEKO_COMM, ierr)
    t_end = MPI_WTIME()
    write(log_buf, '(A,F9.6)') &
         'Mesh manager initialisation (including file reading) time (s): ', &
         t_end - t_start
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! initial importing of mesh data; it may be not complete
    call this%import()
    t_start = MPI_WTIME()
    write(log_buf, '(A,F9.6)') &
         'Mesh manager initial import time (s): ', &
         t_start - t_end
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

  end subroutine p4est_init_from_components

  !> Import data from p4est
  subroutine p4est_import_data(mesh_new, is_periodic)
    type(manager_mesh_p4est_t), intent(inout) :: mesh_new
    logical, intent(in) :: is_periodic
    character(len=LOG_SIZE) :: log_buf
    character(len=*), parameter :: frmt1="('mesh: element number =', i6,&
         &', max ref. lev. = ',i2)"
    character(len=*), parameter :: frmt2="('geometry: dim = ', i1, ';&
         & independent node number=', i6)"
    character(len=*), parameter :: frmt3="('connectivity: dim = ', i1, ';&
         & number of: vrt = ', i6, ', fcs = ', i6, ', edg = ', i6)"
    integer(i4) :: nvert, nface, nedge, ierr, gdim, nelt, nelv, ngrp, maxl, &
         maxg, ndep, lown, lshr, loff, lnum_in, lnum_fh, lnum_eh, nrank, nshare
    integer(i8) :: itmp8, gnelt, gnelto, goff, gnum
    integer(i8), dimension(3) :: itmp8lv
    integer(i8), allocatable, target, dimension(:) :: itmp8v1
    integer(i8), allocatable, target, dimension(:, :) :: itmp8v21
    integer(i4), allocatable, target, dimension(:) :: itmp4v1, itmp4v2, &
         itmp4v3, hngel
    integer(i4), allocatable, target, dimension(:, :) :: itmp4v21, itmp4v22, &
         itmp4v23, hngfc, hnged, vmap, fmap, emap, ealgn
    real(dp), allocatable, target, dimension(:, :) :: rtmpv1
!!$    integer(i4) :: il, jl

    call profiler_start_region("p4est import", 32)

    write(log_buf, '(a)') 'Importing p4est data'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! if mesh has periodic boundaries swap connectivity in tree to get
    ! geometrical information
    if (is_periodic) call wp4est_tree_cnn_swap()

    ! create p4est ghost zones and nodes
    call wp4est_ghost_new()
    call wp4est_nodes_new()

    ! get mesh size and distribution information
    call wp4est_msh_get_size(gdim, gnelt, gnelto, nelt, nelv, ngrp, maxl)

    ! object number
    nvert = 2**gdim
    nface = 2 * gdim
    nedge = 12 * (gdim - 2)

    ! get max refinement level across all ranks
    call MPI_Allreduce(maxl, maxg, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
    ! group information not used for now

    write(log_buf, frmt1) gnelt, maxg
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

!    if (nelt == 0) then
!    end if

    ! get geometry info
    select type (geom => mesh_new%geom)
    type is (manager_geom_p4est_t)

       ! get nodes and their coordinates
       call wp4est_nds_get_size(lown, lshr, loff, lnum_in, lnum_fh, lnum_eh)
       itmp8 = lown
       call MPI_Allreduce(itmp8, gnum, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, &
            ierr)
       write(log_buf, frmt2) gdim, gnum
       call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

!!$       write(*,*) 'TESTsize ', pe_rank, nelv, ngrp, nelt
!!$       write(*,*) 'TESTsizen ', pe_rank, lown, lshr, loff, lnum_in, &
!!$            lnum_fh, lnum_eh

       ! independent nodes
       allocate(itmp8v1(lnum_in), itmp4v1(lnum_in), rtmpv1(gdim, lnum_in))
       call wp4est_nds_get_ind(c_loc(itmp8v1), c_loc(itmp4v1),c_loc(rtmpv1))

!!$       do il = 1, lnum_in
!!$          write(*,*) 'TESTind ', pe_rank, il, itmp8v1(il), &
!!$               itmp4v1(il), rtmpv1(:,il)
!!$       end do

       select type (ind => geom%ind)
       type is(manager_geom_node_ind_p4est_t)
          call ind%init_data(lown, lshr, loff, lnum_in, gdim, &
               itmp8v1, itmp4v1, rtmpv1)

          if (lnum_fh > 0) then
             ! face hanging nodes
             ndep = 4
             allocate(itmp4v21(ndep, lnum_fh), rtmpv1(gdim, lnum_fh))
             call wp4est_nds_get_hfc(c_loc(itmp4v21), c_loc(rtmpv1))

             ! get global numbering of hanging nodes
             call p4est_hng_node_gnum_get(itmp8v1, itmp4v21, lnum_fh, ndep, &
                  ind%gidx, ind%ndown, ind%lnum)

             call geom%hng_fcs%init_data(lnum_fh, gdim, ndep, &
                  itmp8v1, itmp4v21, rtmpv1)
          end if

          if (lnum_eh > 0) then
             ! edge hanging nodes
             ndep = 2
             allocate(itmp4v21(ndep, lnum_eh), rtmpv1(gdim, lnum_eh))
             call wp4est_nds_get_hed(c_loc(itmp4v21), c_loc(rtmpv1))

             ! get global numbering of hanging nodes
             call p4est_hng_node_gnum_get(itmp8v1, itmp4v21, lnum_eh, ndep, &
                  ind%gidx, ind%ndown, ind%lnum)

             call geom%hng_edg%init_data(lnum_eh, gdim, ndep, &
                  itmp8v1, itmp4v21, rtmpv1)
          end if
       end select

       ! get element vertex mapping to nodes
       allocate(itmp4v21(nvert, nelt))
       call wp4est_nds_get_vmap(c_loc(itmp4v21))

!!$       do il = 1, nelt
!!$          do jl = 1, nvert
!!$             write(*,*) 'TESTvmap ', pe_rank, il, jl, itmp4v21(jl, il)
!!$          end do
!!$       end do

       ! geometry type saving already imported data
       call geom%init_data(gdim, nelt, itmp4v21, .true.)

    end select

    call wp4est_nodes_del()

    ! for periodic domain reconstruct ghosts and nodes to include periodicity
    ! in connectivity data
    if (is_periodic) then
       call wp4est_ghost_del()

       call wp4est_tree_cnn_swap_back()

       call wp4est_ghost_new()
    end if

    ! get connectivity info
    select type (conn => mesh_new%conn)
    type is (manager_conn_p4est_t)

       ! vertices
       call wp4est_lnodes_new(1)

       ! get hanging object info; based on lnode information
       allocate(hngel(nelt), hngfc(nface, nelt), hnged(nedge, nelt))
       call wp4est_hang_get_info(c_loc(hngel), c_loc(hngfc), c_loc(hnged))

       ! vertex connectivity
       allocate(vmap(nvert, nelt))
       call wp4est_elm_get_lnode(lnum_in, lown, goff, c_loc(vmap))
       call wp4est_sharers_get_size(nrank, nshare)
       allocate(itmp8v1(lnum_in), itmp4v1(nrank), itmp4v2(nrank+1), &
            itmp4v3(nshare))
       call wp4est_sharers_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
            c_loc(itmp4v2), c_loc(itmp4v3))
       itmp8 = lown
       call MPI_Allreduce(itmp8, gnum, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, &
            ierr)
       itmp8lv(1) = gnum ! for stamping log
       select type (vrt => conn%vrt)
       type is (manager_conn_obj_p4est_t)
          call vrt%init_data(lnum_in, lown, goff, gnum, nrank, nshare, &
               itmp8v1, itmp4v1, itmp4v3, itmp4v2)
       end select

!!$       write(*,*) 'TESTcvrt', lnum_in, lown, goff, gnum, nrank, nshare

       call wp4est_lnodes_del()

       ! faces
       call wp4est_lnodes_new(-1)

       ! face connectivity
       allocate(fmap(nface, nelt))
       call wp4est_elm_get_lnode(lnum_in, lown, goff, c_loc(fmap))
       call wp4est_sharers_get_size(nrank, nshare)
       allocate(itmp8v1(lnum_in), itmp4v1(nrank), itmp4v2(nrank+1), &
            itmp4v3(nshare))
       call wp4est_sharers_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
            c_loc(itmp4v2), c_loc(itmp4v3))
       itmp8 = lown
       call MPI_Allreduce(itmp8, gnum, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, &
            ierr)
       itmp8lv(2) = gnum ! for stamping log
       select type (fcs => conn%fcs)
       type is (manager_conn_obj_p4est_t)
          call fcs%init_data(lnum_in, lown, goff, gnum, &
               nrank, nshare, itmp8v1, itmp4v1, itmp4v3, itmp4v2)
       end select

!!$       write(*,*) 'TESTcfcs', lnum_in, lown, goff, gnum, nrank, nshare

       call wp4est_lnodes_del()

       ! edges
       call wp4est_lnodes_edge(1)

       ! edge connectivity
       allocate(emap(nedge, nelt))
       call wp4est_elm_get_lnode(lnum_in, lown, goff, c_loc(emap))
       call wp4est_sharers_get_size(nrank, nshare)
       allocate(itmp8v1(lnum_in), itmp4v1(nrank), itmp4v2(nrank+1), &
            itmp4v3(nshare))
       call wp4est_sharers_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
            c_loc(itmp4v2), c_loc(itmp4v3))
       itmp8 = lown
       call MPI_Allreduce(itmp8, gnum, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, &
            ierr)
       itmp8lv(3) = gnum ! for stamping log
       select type (edg => conn%edg)
       type is (manager_conn_obj_p4est_t)
          call edg%init_data(lnum_in, lown, goff, gnum, &
               nrank, nshare, itmp8v1, itmp4v1, itmp4v3, itmp4v2)
       end select

       write(log_buf, frmt3) gdim, itmp8lv
       call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

!!$       write(*,*) 'TESTcedg', lnum_in, lown, goff, gnum, nrank, nshare

       call wp4est_lnodes_del()

       ! get element data
       allocate(itmp8v1(nelt), itmp4v1(nelt), itmp4v2(nelt), &
            & itmp4v21(nface, nelt), itmp4v22(nface, nelt), &
            & itmp4v23(nface, nelt))
       call wp4est_elm_get_dat(c_loc(itmp8v1), c_loc(itmp4v1), &
            & c_loc(itmp4v2), c_loc(itmp4v21), c_loc(itmp4v22), &
            & c_loc(itmp4v23))

       ! get edge alignment
       call p4est_edge_alignment_get(ealgn, nelt, nedge)

!!$       write(*,*) 'TESTvrt'
!!$       do il = 1, nelt
!!$          do jl = 1, nvert
!!$             write(*,*) 'TESTvmap', il, jl, vmap(jl, il)
!!$          end do
!!$       end do
!!$       write(*,*) 'TESTfcs'
!!$       do il = 1, nelt
!!$          do jl = 1, nface
!!$             write(*,*) 'TESTfmap', il, jl, fmap(jl, il), itmp4v23(jl, il)
!!$          end do
!!$       end do
!!$       write(*,*) 'TESTedg'
!!$       do il = 1, nelt
!!$          do jl = 1, nedge
!!$             write(*,*) 'TESTemap', il, jl, emap(jl, il), ealgn(jl, il)
!!$          end do
!!$       end do

       ! element connectivity mappings saving already imported data
       call conn%init_data(gdim, nelt, vmap, fmap, itmp4v23, emap, ealgn, &
            hngel, hngfc, hnged, .true.)

    end select

    ! element family
    call p4est_family_get(itmp8v21, nelt, nvert, gnelto)

    ! import element general information saving already imported data
    call mesh_new%init_data(nelt, nelv, gnelt, gnelto, maxg, gdim, itmp8v1, &
            itmp4v1, itmp4v2, itmp4v21, itmp4v22, itmp8v21, .true.)

    ! do not destroy p4est ghost cells here, as they could be needed
!    call wp4est_ghost_del()

    call profiler_end_region("p4est import", 32)

  end subroutine p4est_import_data

  !> Get global numbering of hanging nodes
  !! @details p4est is a connectivity code, so it does not care about global
  !! numbering of hanging nodes just mapping them to independent nodes
  subroutine p4est_hng_node_gnum_get(gidx, map, lnum, ndep, ingidx, inown, &
       inlnum)
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), dimension(:, :), intent(in) :: map
    integer(i4), intent(in) :: lnum, ndep
    integer(i8), dimension(:), intent(in) :: ingidx
    integer(i4), dimension(:), intent(in) :: inown
    integer(i4), intent(in) :: inlnum

    allocate(gidx(lnum))

    ! This is not finished yet
    call neko_error('Hanging node numbering not done yet')

  end subroutine p4est_hng_node_gnum_get

  !> Extract edge orientation
  !! @detail Unfortunately p4est does not keep global edge orientation, as
  !! it distinguishes between face edge and true edge keeping relative
  !! orientation of the lat category (in pairs) only. To get it I use a tric
  !! related to the lnodes position.
  !! @param[inout]  ealgn    edge alignment
  !! @param[in]     nelt     number of elements
  !! @param[in]     nedge    number of edges
  subroutine p4est_edge_alignment_get(ealgn, nelt, nedge)
    integer(i4), allocatable, dimension(:, :), intent(inout) :: ealgn
    integer(i4), intent(in) :: nelt, nedge
    integer(i4) :: lnum, lown, il, jl
    integer(i8) :: goff
    integer(i4), allocatable, target, dimension(:, :, :) :: emap

    allocate(ealgn(nedge, nelt), emap(2, nedge, nelt))
    ! edge lnodes
    call wp4est_lnodes_edge(2)
    call wp4est_elm_get_lnode(lnum, lown, goff, c_loc(emap))

    ! get alignment
    do  il = 1, nelt
       do jl = 1, nedge
          if (emap(1, jl, il) < emap(2, jl, il)) then
             ealgn(jl, il) = 0
          else
             ealgn(jl, il) = 1
          end if
       end do
    end do

    call wp4est_lnodes_del()

    deallocate(emap)

  end subroutine p4est_edge_alignment_get

  !> Extract family information
  !! @details This routine provides information about sets of elements that
  !! share the same parent and could be destroyed togeter during coarsening
  !! step.
  !! @param[inout]  family   family flag
  !! @param[in]     nelt     number of elements
  !! @param[in]     nvert    number of vertices
  !! @param[in]     gnelto   global element offset
  subroutine p4est_family_get(family, nelt, nvert, gnelto)
    integer(i8), allocatable, dimension(:, :), intent(inout) :: family
    integer(i4), intent(in) :: nelt, nvert
    integer(i8), intent(in) :: gnelto
    integer(i4) :: ierr, il, jl
    integer(i8) :: itmp8
    integer(i4) :: nlfam ! local number of families
    integer(i8) :: nlfam_off ! global family nuimber offset
    integer(i8), allocatable, target, dimension(:) :: itmp4v1

    ! import family mark
    allocate(itmp4v1(nelt))
    call wp4est_fml_get_info(c_loc(itmp4v1), nlfam)

    ! test number of families; it must be multiply of number of vertices
    if ((nlfam > nelt).or.(mod(nlfam, nvert) /= 0)) &
         call neko_error('Invalid number of families.')

    ! get global family offset
    nlfam = nlfam/nvert
    itmp8 = nlfam
    call MPI_Scan(itmp8, nlfam_off, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
    nlfam_off = nlfam_off - itmp8

    ! mark all local families
    allocate(family(2, nelt))
    family(:, :) = 0
    do il=1, nlfam
       nlfam_off = nlfam_off + 1
       do jl = 1, nvert
          itmp8 = itmp4v1(jl + (il-1)*nvert) - gnelto
          ! one could test if 0<itmp8<=nelt
          family(1, int(itmp8,i4)) = nlfam_off
          family(2, int(itmp8,i4)) = nvert + 1 - jl
       end do
    end do

    deallocate(itmp4v1)

  end subroutine p4est_family_get

  !> Apply data read from mesh file to mesh manager structures
  subroutine p4est_mesh_file_apply(this)
    class(mesh_manager_p4est_t), intent(inout) :: this
    integer(i4) :: il, iel, ifc, ibc, nin, nhf, nhe
    integer(i8) :: itmp8
    integer(i4), target, allocatable, dimension(:) :: imsh

    select type(mesh => this%mesh)
    type is (manager_mesh_p4est_t)

       ! Fill all the information that may be not properly initialised
       ! neko supports just V-type elements
       allocate(imsh(mesh%nelt))
       imsh(:) = 0
       !> group flag is not used for now
       mesh%igrp(:) = 0
       ! face curvature data is not used for now; projection will be based on bc
       mesh%crv(:, :) = 0
       ! fill in boundary condition
       ! start with marking everything internal
       mesh%bc(:, :) = 0
       do il = 1, this%nmsh_mesh%nzone
          itmp8 = int(this%nmsh_mesh%zone(il)%e, i8)
          iel =  int(itmp8 - this%nmsh_mesh%offset_el, i4)
          if (itmp8 .ne. mesh%gidx(iel)) &
               call neko_error('Inconsistent global element number')
          select case(this%nmsh_mesh%zone(il)%type)
          case (5)
             ibc = -1
          case (7)
             ibc = this%nmsh_mesh%zone(il)%p_f
          end select
          mesh%bc(this%nmsh_mesh%zone(il)%f, iel) = ibc
       end do

       ! transfer data to p4est
       call wp4est_elm_ini_dat(c_loc(mesh%gidx), c_loc(imsh), &
            c_loc(mesh%igrp), c_loc(mesh%crv), c_loc(mesh%bc))

       ! check if applied boundary conditions are consistent with tree structure
       ! the ghost mesh was not destroyed at the end of data import, so no need
       ! to set it here
!       call wp4est_ghost_new()
       call wp4est_bc_check()
!       call wp4est_ghost_del()

       deallocate(imsh)
    end select

    ! correct vertex information in the mesh manager geometry
    !>
    !! @todo This is just a hack and should be done better in the future
    ! For now just local copy of vertices coordinates; it does not mater if the
    ! mesh file is fine, but in case of corrupted file could be important
    select type(geom => this%mesh%geom)
    type is (manager_geom_p4est_t)
       ! get limits for different node types
       nin = geom%ind%lnum
       nhf = nin + geom%hng_fcs%lnum
       nhe = nhf + geom%hng_edg%lnum
       do il = 1, geom%nel
          do ifc = 1, geom%nvrt
             ibc = geom%vmap(ifc, il)
             if (ibc .le. nin) then
                geom%ind%coord(:, ibc) = this%nmsh_mesh%hex(il)%v(ifc)%v_xyz(:)
             else if (ibc .le. nhf) then
                geom%hng_fcs%coord(:, ibc - nin) = &
                     this%nmsh_mesh%hex(il)%v(ifc)%v_xyz(:)
             else if (ibc .le. nhe) then
                geom%hng_edg%coord(:, ibc - nhf) = &
                     this%nmsh_mesh%hex(il)%v(ifc)%v_xyz(:)
             else
                call neko_error('Inconsistent vertex to node mapping')
             end if
          end do
       end do

    end select

  end subroutine p4est_mesh_file_apply

  !> Perform refinement/coarsening on the mesh manager side
  !! @param  ref_mark     refinement flag
  !! @param[out]  ifmod        mesh modification flag
  subroutine p4est_refine(this, ref_mark, ifmod)
    class(mesh_manager_p4est_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: ref_mark
    character(len=LOG_SIZE) :: log_buf
    logical, intent(out) :: ifmod
    integer(i8), target, allocatable, dimension(:) :: pel_gnum
    integer(i4), target, allocatable, dimension(:) :: pref_mark, pel_lnum, &
         pel_nid
    integer(i4), parameter :: p4est_compare = 0
    integer(i4) :: il, itmp
    ! element restructure data
    integer(i4) :: map_nr, rfn_nr, crs_nr
    integer(i8), target, allocatable, dimension(:) :: map_gidx
    integer(i8), target, allocatable, dimension(:, :) :: rfn_gidx
    integer(i8), target, allocatable, dimension(:, :, :) :: crs_gidx
    integer(i4), target, allocatable, dimension(:, :) :: map, rfn
    integer(i4), target, allocatable, dimension(:, :, :) :: crs

    call profiler_start_region("p4est refine", 31)

    write(log_buf, '(a)') 'p4est refinement/coarsening start'
    call neko_log%message(log_buf, NEKO_LOG_INFO)

    ifmod = .false.

    select type(transfer => this%transfer)
    type is (mesh_manager_transfer_p4est_t)
       ! check ref_mark size
       if (size(ref_mark) .ne. transfer%nelt_neko) &
            call neko_error('Inconsistent ref_mark array size')

       ! transfer refinement flag between neko and p4est distributions
       call transfer%ref_mark_transfer(ref_mark, pref_mark)

       ! PLACE FOR LOCAL CONSISTENCY CHECKS

       ! feed p4est with refinement data
       call wp4est_refm_put(c_loc(pref_mark))

       ! set element distribution info in p4est
       call transfer%neko_elem_dist_set(pel_gnum, pel_lnum, pel_nid)
       call wp4est_egmap_put(c_loc(pel_gnum), c_loc(pel_lnum), c_loc(pel_nid))

       ! destroy p4est ghost cells
       call wp4est_ghost_del()

       ! perform local refine/coarsen/balance on p4est side
       call wp4est_tree_compare_copy(p4est_compare)
       call wp4est_refine(this%ref_level_max)
       call wp4est_coarsen()
       call wp4est_balance()
       ! perform partitioning on p4est side
       call wp4est_part()
       call wp4est_tree_compare_check(itmp, p4est_compare)

       if (itmp == 0 ) then
          ! set refinement flag
          ifmod = .true.

          write(log_buf, '(a)') 'p4est refinement; mesh changed'
          call neko_log%message(log_buf, NEKO_LOG_INFO)

          ! import new data
          call this%import()

          ! get data to reconstruct fields
          ! data sizes
          call wp4est_msh_get_hst_size(map_nr, rfn_nr, crs_nr)
          ! number of children is equal to number of vertices
          itmp = 2**this%mesh%tdim
          allocate(map_gidx(this%mesh%nelt), map(2, this%mesh%nelt), &
               rfn_gidx(2, rfn_nr), rfn(3, rfn_nr), &
               crs_gidx(2, itmp, crs_nr), crs(2, itmp, crs_nr))
          call wp4est_msh_get_hst(c_loc(map_gidx), c_loc(map), &
               c_loc(rfn_gidx), c_loc(rfn), c_loc(crs_gidx), c_loc(crs))
          ! store data in the type
          call transfer%reconstruct_data_set(map_nr, map_gidx, map, rfn_nr, &
               rfn_gidx, rfn, crs_nr, crs_gidx, crs)
       else
          ! regenerate the ghost layer
          call wp4est_ghost_new()
          write(log_buf, '(a)') 'p4est refinement; mesh not changed'
          call neko_log%message(log_buf, NEKO_LOG_INFO)
       end if

       deallocate(pref_mark, pel_gnum, pel_lnum, pel_nid)
    class default
       call neko_error('Wrong data redistribution type')
    end select

    write(log_buf, '(a)') 'p4est refinement/coarsening end'
    call neko_log%message(log_buf, NEKO_LOG_INFO)

    call profiler_end_region("p4est refine", 31)

  end subroutine p4est_refine

  !> Construct neko mesh type based on mesh manager data
  !! @param[inout]   mesh     neko mesh type
  subroutine p4est_mesh_construct(this, mesh)
    class(mesh_manager_p4est_t), intent(in) :: this
    type(mesh_t), intent(inout) :: mesh
    integer :: il, itmp
    integer(i4) :: nidx
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_start, t_end

    t_start = MPI_WTIME()
    write(log_buf, '(a)') 'Constructing neko mesh type'
    call neko_log%message(log_buf, NEKO_LOG_INFO)

    ! What follows is just a hack, as I'm trying to minimise changes outside
    ! this file
    ! There are some differences in a concept, so not everything can be filled
    ! in properly

    ! Free neko mesh type
    call mesh%free()

    ! General mesh info
    ! There is no distinction between geometry and topology dimension
    mesh%gdim = this%mesh%tdim
    mesh%npts = this%mesh%geom%nvrt ! number of element vertices

    ! general element info
    mesh%nelv = this%mesh%nelt ! local number of elements
    mesh%glb_nelv = int(this%mesh%gnelt, i4) ! global number of elelments
    mesh%offset_el = int(this%mesh%gnelto, i4) ! global element offset;

    ! Geometrical context
    ! Local information only; global one will be related to connectivity
    select type (geom => this%mesh%geom)
    type is (manager_geom_p4est_t)
       ! local number of unique vertices
       mesh%mpts = geom%ind%lnum + geom%hng_fcs%lnum + geom%hng_edg%lnum

       ! fill the points
       ! order matters due to the element vertex mapping
       allocate(mesh%points(mesh%mpts))
       itmp = 0
       ! independent nodes
       if (geom%ind%lnum .gt. 0) then
          do il = 1, geom%ind%lnum
             itmp = itmp + 1
             nidx = int(geom%ind%gidx(il), i4)
             call mesh%points(itmp)%init(geom%ind%coord(:, il), nidx)
          end do
       end if
       ! face hanging nodes
       if (geom%hng_fcs%lnum .gt. 0) then
          do il = 1, geom%hng_fcs%lnum
             itmp = itmp + 1
             nidx = int(geom%hng_fcs%gidx(il), i4)
             call mesh%points(itmp)%init(geom%hng_fcs%coord(:, il), nidx)
          end do
       end if
       ! edge hanging nodes
       if (geom%hng_edg%lnum .gt. 0) then
          do il = 1, geom%hng_edg%lnum
             itmp = itmp + 1
             nidx = int(geom%hng_edg%gidx(il), i4)
             call mesh%points(itmp)%init(geom%hng_edg%coord(:, il), nidx)
          end do
       end if

       ! Global to local element mapping
       call mesh%htel%init(mesh%nelv, il)

       ! Fill in element array
       allocate(mesh%elements(mesh%nelv))
       do il = 1, mesh%nelv
          allocate(hex_t::mesh%elements(il)%e)
          nidx = int(this%mesh%gidx(il), i4)
          select type(ep => mesh%elements(il)%e)
          type is (hex_t)
             call ep%init(nidx, &
                  mesh%points(geom%vmap(1, il)), &
                  mesh%points(geom%vmap(2, il)), &
                  mesh%points(geom%vmap(3, il)), &
                  mesh%points(geom%vmap(4, il)), &
                  mesh%points(geom%vmap(5, il)), &
                  mesh%points(geom%vmap(6, il)), &
                  mesh%points(geom%vmap(7, il)), &
                  mesh%points(geom%vmap(8, il)))
          end select
          ! Global to local element mapping
          itmp = il
          call mesh%htel%set(nidx, itmp)
       end do
    end select

    ! Element deformation
    ! deformed element flag; It is set in mesh_generate_flags
    ! (dependent on vertex ordering)!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(mesh%dfrmd_el(mesh%nelv))
    
    ! for now set everything true
    mesh%dfrmd_el = .true.
    

    ! Connectivity context
    ! missing number of vertices in connectivity context
    mesh%mfcs = this%mesh%conn%fcs%lnum ! local number of unique faces
    mesh%meds = this%mesh%conn%edg%lnum ! local number of unique edges

    ! global number of unique vertices, faces and edges
    mesh%glb_mpts = int(this%mesh%conn%vrt%gnum, i4)
    mesh%glb_mfcs = int(this%mesh%conn%fcs%gnum,  i4)
    mesh%glb_meds = int(this%mesh%conn%edg%gnum, i4)
    ! this doesn't seem to be used in original version and
    ! I do not need it for p4est either
    !mesh%max_pts_id =

    ! mesh%htp - not used as vertices global id is already provided by p4est
    ! mesh%htf - not used as faces global id is already provided by p4est
    ! mesh%hte - not used as edge global id is already provided by p4est

    ! Fill in neighbour information
    select type (conn => this%mesh%conn)
    type is (manager_conn_p4est_t)
       ! USED IN gather_scatter.f90;
       ! MORE INVESTIGATION NEEDED
       ! get neighbour rank info; based on vertex sharing info
       call p4est_rank_neighbour_fill(mesh, conn)

       ! get vertex neighbours; not used outside mesh.f90
       call p4est_vertex_neighbour_fill(mesh, conn)

       ! USED IN gather_scatter.f90, tree_amg_multigrid.f90;
       ! MORE INVESTIGATION NEEDED
       ! get face neighbours
       call p4est_face_neighbour_fill(mesh, conn)
    end select

    write(*, *) 'TEST mesh construct', pe_rank, mesh%glb_mpts, mesh%glb_mfcs, &
         mesh%glb_meds, mesh%mfcs, mesh%meds, this%mesh%conn%vrt%lnum

    !call MPI_Barrier(NEKO_COMM, ierr)
    t_end = MPI_WTIME()
    write(log_buf, '(A,F9.6)') 'Mesh construction time (s): ', &
         t_end - t_start
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

  end subroutine p4est_mesh_construct

  !> Fill the mesh type with rank neighbour information
  !! @details Based on vertex connectivity info
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      conn    connectivity
  subroutine p4est_rank_neighbour_fill(mesh, conn)
    type(mesh_t), intent(inout) :: mesh
    type(manager_conn_p4est_t), intent(in) :: conn
    integer :: il, jl, src, dst
    type(stack_i4_t) :: neigh_order

    ! initialize neighbour ranks
    allocate(mesh%neigh(0:pe_size-1))
    mesh%neigh = .false.

    if (pe_size > 1) then
       select type (vrt => conn%vrt)
       type is (manager_conn_obj_p4est_t)
          ! get rank list based on vertex information
          do il= 1, vrt%nrank ! mpi rank loop
             ! not sure I have to exclude myself
             if (vrt%rank(il) /= pe_rank) mesh%neigh(vrt%rank(il)) = .true.
          end do
       end select

       ! Generate neighbour exchange order
       call neigh_order%init(pe_size)

       do il = 1, pe_size - 1
          src = modulo(pe_rank - il + pe_size, pe_size)
          dst = modulo(pe_rank + il, pe_size)
          if (mesh%neigh(src) .or. mesh%neigh(dst)) then
             jl = il ! adhere to standards...
             call neigh_order%push(jl)
          end if
       end do

       allocate(mesh%neigh_order(neigh_order%size()))
       select type(order => neigh_order%data)
       type is (integer)
          do il = 1, neigh_order%size()
             mesh%neigh_order(il) = order(il)
          end do
       end select
       call neigh_order%free()
    else
       allocate(mesh%neigh_order(1))
       mesh%neigh_order = 1
    end if

  end subroutine p4est_rank_neighbour_fill

  !> Fill the mesh type with vertex neighbour information
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      conn    connectivity
  subroutine p4est_vertex_neighbour_fill(mesh, conn)
    type(mesh_t), intent(inout) :: mesh
    type(manager_conn_p4est_t), intent(in) :: conn

    write(*, *) 'TEST vertex neigh', pe_rank
  end subroutine p4est_vertex_neighbour_fill

  !> Fill the mesh type with face neighbour information
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      conn    connectivity
  subroutine p4est_face_neighbour_fill(mesh, conn)
    type(mesh_t), intent(inout) :: mesh
    type(manager_conn_p4est_t), intent(in) :: conn

    write(*, *) 'TEST face neigh', pe_rank
  end subroutine p4est_face_neighbour_fill

#else

  !> Start p4est
  !! @param[out]  ierr  error flag
  subroutine p4est_start(this, json, ierr)
    class(mesh_manager_p4est_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(out) :: ierr

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_start

  !> Stop p4est
  !! @param[out]  ierr  error flag
  subroutine p4est_stop(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_stop

  !> The common constructor using a JSON object.
  !! @param json       The JSON object.
  !! @param type_name  Manager type name
  subroutine p4est_init_from_json(this, json)
    class(mesh_manager_p4est_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_init_from_json

  !> Destructor.
  subroutine p4est_free(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_free

  !> Import mesh data
  subroutine p4est_import(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_import

  !> Apply data read from mesh file to mesh manager structures
  subroutine p4est_mesh_file_apply(this)
    class(mesh_manager_p4est_t), intent(inout) :: this

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_mesh_file_apply

  !> Perform refinement/coarsening on the mesh manager side
  !! @param  ref_mark     refinement flag
  !! @param[out]  ifmod        mesh modification flag
  subroutine p4est_refine(this, ref_mark, ifmod)
    class(mesh_manager_p4est_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: ref_mark
    logical, intent(out) :: ifmod

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_refine

  !> Construct neko mesh type based on mesh manager data
  !! @param[inout]   mesh     neko mesh type
  subroutine p4est_mesh_construct(this, mesh)
    class(mesh_manager_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: mesh
  end subroutine p4est_mesh_construct
#endif

end module mesh_manager_p4est
