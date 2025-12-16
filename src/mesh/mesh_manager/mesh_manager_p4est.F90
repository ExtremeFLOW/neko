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
  use tuple, only : tuple_i4_t, tuple4_i4_t
  use stack, only : stack_i4_t, stack_i4t2_t
  use hex, only : hex_t
  use mesh_conn, only : mesh_conn_t
  use mesh, only : mesh_t, mesh_generate_flags, NEKO_MSH_MAX_ZLBLS
  use nmsh, only: nmsh_mesh_t
  use manager_mesh, only : manager_mesh_t
  use manager_geom_p4est, only : manager_geom_node_ind_p4est_t, &
       manager_geom_p4est_t
  use manager_conn_p4est, only : manager_conn_obj_p4est_t, manager_conn_p4est_t
  use manager_mesh_p4est, only : manager_mesh_p4est_t
  use mesh_manager, only : mesh_manager_t
  use mesh_manager_transfer_p4est, only : mesh_manager_transfer_p4est_t
  use, intrinsic :: iso_c_binding, only : c_null_char

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

     subroutine wp4est_nds_get_vmap(vmap, vcoord, tol) &
          bind(c, name = 'wp4est_nds_get_vmap')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: vmap, vcoord
       real(c_double) :: tol
     end subroutine wp4est_nds_get_vmap

     subroutine wp4est_nds_get_vcoord(vcoord) &
          bind(c, name = 'wp4est_nds_get_vcoord')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: vcoord
     end subroutine wp4est_nds_get_vcoord

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
    if (log_level .eq. 0) then
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

       if ((log_level .ge. 0) .and. (log_level .le. 9)) then
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
    if (log_level .eq. 0) then
       call neko_error('Failure starting p4est')
    else
       this%ifstarted = .true.
       ierr = 0
    end if

    ! this is AMR framework
    this%isamr = .true.

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
    character(len=:), allocatable :: tree_file
    integer :: ref_level_max

    ! Extract runtime parameters
    ! tree_file is mandatory
    call json_get_or_default(json, 'tree_file', tree_file, 'no tree')
    if (trim(tree_file) .eq. 'no tree') then
       call neko_error('The tree_file keyword could not be found in the .' // &
            'case file. Often caused by incorrectly formatted json.')
    end if

    ! p4est supports AMR; get maximum allowed refinement level
    ! Notice, p4est has internal max refinement level (P8EST_QMAXLEVEL
    ! defined in p8est.h) that cannot be exceeded. Moreover, negative
    ! value of ref_level_max will automatically set the restriction to
    ! P8EST_QMAXLEVEL. By default we set no refinement.
    call json_get_or_default(json, "ref_level_max", ref_level_max, 0)

    call p4est_init_from_components(this, trim(tree_file), ref_level_max)

    if (allocated(tree_file)) deallocate(tree_file)

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

    this%ref_level_max = 0
    this%log_level = 0

    if (allocated(this%tree_file)) deallocate(this%tree_file)

  end subroutine p4est_free

  !> Import mesh data
  subroutine p4est_import(this, ifcomplete)
    class(mesh_manager_p4est_t), intent(inout) :: this
    logical, intent(in) :: ifcomplete
    type(manager_mesh_p4est_t) :: mesh_new

    ! allocate types and import data from p4est
    call mesh_new%init()
    call mesh_new%free_data()
    call p4est_import_data(mesh_new, ifcomplete)

    ! fill mesh information
    call this%mesh%init_type(mesh_new)
    call mesh_new%free()

    ! get new element distribution; get partitioning
    call this%transfer%elem_dist_construct(this%mesh)

  end subroutine p4est_import

  !> The constructor from type components.
  subroutine p4est_init_from_components(this, tree_file, ref_level_max)
    class(mesh_manager_p4est_t), intent(inout) :: this
    character(len=*), intent(in) :: tree_file
    integer, intent(in) :: ref_level_max
    integer(i4) :: is_valid, ierr
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_start, t_end
    logical :: exist

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

    ! AMR related stuff
    this%ref_level_max = ref_level_max

    ! for testing
!    call wp4est_cnn_rot_cubes()
!    call wp4est_cnn_brick(2, 2, 2)
!    call wp4est_cnn_unit_cube_periodic()
!    call wp4est_tree_new()
!    call wp4est_vtk_write('test')
!    call wp4est_tree_del()
!    call wp4est_cnn_del()

    call neko_log%message('Reading p4est tree data from the file '// &
         trim(this%tree_file), NEKO_LOG_VERBOSE)
    inquire(file = trim(this%tree_file), exist = exist)
    if (.not. exist) call neko_error('Missing p4est tree file.')
    ! read the tree file
    call wp4est_tree_load(trim(this%tree_file)//c_null_char)
    call wp4est_tree_valid(is_valid)
    if (is_valid .eq. 0) call neko_error('Invalid p4est tree')
    call wp4est_cnn_valid(is_valid)
    if (is_valid .eq. 0) call neko_error('Invalid p4est connectivity')
    ! perform partitioning on p4est side
    call wp4est_part()

!    call wp4est_cnn_save('test.cnn'//c_null_char)
!    call wp4est_tree_save('test.tree'//c_null_char)
!    call wp4est_vtk_write('test'//c_null_char)

    !call MPI_Barrier(NEKO_COMM, ierr)
    t_end = MPI_WTIME()
    write(log_buf, '(A,F9.6)') &
         'Mesh manager initialisation (including file reading) time (s): ', &
         t_end - t_start
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

  end subroutine p4est_init_from_components

  !> Import data from p4est
  subroutine p4est_import_data(mesh_new, ifcomplete)
    type(manager_mesh_p4est_t), intent(inout) :: mesh_new
    logical, intent(in) :: ifcomplete
    character(len=LOG_SIZE) :: log_buf
    character(len=*), parameter :: frmt1="('mesh: element number =', i6,&
         &', max ref. lev. = ',i2)"
    character(len=*), parameter :: frmt2="('geometry: dim = ', i1, ';&
         & independent node number=', i6)"
    character(len=*), parameter :: frmt3="('connectivity: dim = ', i1, ';&
         & number of: vrt = ', i6, ', fcs = ', i6, ', edg = ', i6)"
    integer(i4) :: nvert, nface, nedge, ierr, gdim, nelt, nelv, ngrp, maxl, &
         maxg, ndep, lown, lshr, loff, lnum_in, lnum_fh, lnum_eh, nrank, &
         nshare, nindp, nfhngp, nehngp
    integer(i8) :: itmp8, gnelt, gnelto, goff, gnum
    integer(i8), dimension(3) :: itmp8lv
    integer(i8), allocatable, target, dimension(:) :: itmp8v1, itmp8v2
    integer(i8), allocatable, target, dimension(:, :) :: itmp8v21
    integer(i4), allocatable, target, dimension(:) :: itmp4v1, itmp4v2, &
         itmp4v3, hngei, ind_map, fhng_map, ehng_map
    logical, allocatable, target, dimension(:) :: hngel
    integer(i4), allocatable, target, dimension(:, :) :: itmp4v21, itmp4v22, &
         itmp4v23, hngfc, hnged, vmap, fmap, emap, ealgn
    real(dp), allocatable, target, dimension(:, :) :: rtmpv1, ind_coord, &
         fhng_coord, ehng_coord
    real(dp), allocatable, target, dimension(:, :, :) :: vcoord
    real(dp), parameter :: tol = 1.0D-10

    call profiler_start_region("p4est import", 32)

    write(log_buf, '(a)') 'Importing p4est data'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! create p4est ghost zones
    call wp4est_ghost_new()

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

!    if (nelt .eq. 0) then
!    end if

    ! get geometry info
    select type (geom => mesh_new%geom)
    type is (manager_geom_p4est_t)

       ! get complete global information
       if (ifcomplete) then
          ! create p4est nodes
          call wp4est_nodes_new()

          ! get nodes and their coordinates
          call wp4est_nds_get_size(lown, lshr, loff, lnum_in, lnum_fh, lnum_eh)
          itmp8 = lown
          call MPI_Allreduce(itmp8, gnum, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, &
               ierr)
          write(log_buf, frmt2) gdim, gnum
          call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

          ! get element vertex mapping to nodes and identify periodic ones
          allocate(itmp4v21(nvert, nelt), vcoord(gdim, nvert, nelt))
          vcoord(:, :, :) = 0.0
          call wp4est_nds_get_vmap(c_loc(itmp4v21), c_loc(vcoord), tol)

          ! Get periodic nodes and correct vertex mapping
          call p4est_periodic_get(itmp4v21, vcoord, nindp, ind_map, ind_coord, &
               nfhngp, fhng_map, fhng_coord, nehngp, ehng_map, ehng_coord, &
               nvert, nelt, lnum_in, lnum_fh, lnum_eh, gdim)
          deallocate(vcoord)

          ! geometry type saving already imported data
          call geom%init_data(gdim, nelt, itmp4v21, .true.)

          select type (ind => geom%ind)
          type is(manager_geom_node_ind_p4est_t)
             ! independent nodes
             allocate(itmp8v1(lnum_in), itmp4v1(lnum_in), &
                  rtmpv1(gdim, lnum_in), itmp8v2(nindp))
             call wp4est_nds_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
                  c_loc(rtmpv1))

             ! for now just a placeholder for global numbering of periodic nodes
             call ind%init_data(lown, lshr, loff, lnum_in, gdim, &
                  itmp8v1, itmp4v1, rtmpv1, nindp, itmp8v2, ind_coord, ind_map)

             if (lnum_fh > 0) then
                ! face hanging nodes
                ndep = 4
                allocate(itmp8v1(lnum_fh), itmp4v21(ndep, lnum_fh), &
                     rtmpv1(gdim, lnum_fh), itmp8v2(nfhngp))
                call wp4est_nds_get_hfc(c_loc(itmp4v21), c_loc(rtmpv1))

                ! for now just a placeholder for global numbering of nodes
                call geom%hng_fcs%init_data(lnum_fh, gdim, ndep, &
                     itmp8v1, itmp4v21, rtmpv1, nfhngp, itmp8v2, fhng_coord, &
                     fhng_map)
             end if

             if (lnum_eh > 0) then
                ! edge hanging nodes
                ndep = 2
                allocate(itmp8v1(lnum_eh), itmp4v21(ndep, lnum_eh), &
                     rtmpv1(gdim, lnum_eh), itmp8v2(nehngp))
                call wp4est_nds_get_hed(c_loc(itmp4v21), c_loc(rtmpv1))

                ! for now just a placeholder for global numbering of nodes
                call geom%hng_edg%init_data(lnum_eh, gdim, ndep, &
                     itmp8v1, itmp4v21, rtmpv1, nehngp, itmp8v2, ehng_coord, &
                     ehng_map)
             end if
          end select

          call wp4est_nodes_del()
       else ! ifcomplete
          ! just a simple set of information containing approximate (linear
          ! interpolation) coordinates of element vertices, but no global
          ! mapping
          ! get coordinates of element vertices
          allocate(vcoord(gdim, nvert, nelt))
          call wp4est_nds_get_vcoord(c_loc(vcoord))

          ! geometry type simple initialisation
          call geom%init_simple(gdim, nelt, vcoord)

          call neko_log%message('Extracted local geometry representation', &
               NEKO_LOG_VERBOSE)
       end if
    end select

    ! get connectivity info
    select type (conn => mesh_new%conn)
    type is (manager_conn_p4est_t)

       ! vertices
       call wp4est_lnodes_new(1)

       ! get hanging object info; based on lnode information
       allocate(hngei(nelt), hngel(nelt), hngfc(nface, nelt), &
            hnged(nedge, nelt))
       call wp4est_hang_get_info(c_loc(hngei), c_loc(hngfc), c_loc(hnged))
       hngel(:) = (hngei(:) == 1)
       deallocate(hngei)

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

       ! element connectivity mappings saving already imported data
       call conn%init_data(gdim, nelt, vmap, fmap, itmp4v23, emap, ealgn, &
            hngel, hngfc, hnged, .true.)

    end select

    ! element family
    call p4est_family_get(itmp8v21, nelt, nvert, gnelto)

    ! import element general information saving already imported data
    call mesh_new%init_data(nelt, nelv, gnelt, gnelto, maxg, gdim, itmp8v1, &
         itmp4v1, itmp4v2, itmp4v21, itmp4v22, itmp8v21, ifcomplete, .true.)

    ! Get global numbering of hanging and periodic nodes including reordering
    ! of faces with respect of multiplicity
    call p4est_hng_periodic_gnum_get(mesh_new)

    ! do not destroy p4est ghost cells here, as they could be needed
!    call wp4est_ghost_del()

    call profiler_end_region("p4est import", 32)

  end subroutine p4est_import_data

  !> Get periodic nodes and correct vertex mapping
  !! @details p4est is a connectivity code, so it does not distinguish between
  !! independent node and the node set on periodic boundary. This can result in
  !! a wrong element shape, so I have to identify them on my own
  subroutine p4est_periodic_get(vmap, vcoord, nind, ind_map, ind_coord, &
       nfhng, fhng_map, fhng_coord, nehng, ehng_map, ehng_coord, nvert, nelt, &
       lnum_in, lnum_fh, lnum_eh, gdim)
    integer(i4), allocatable, target, dimension(:, :), intent(inout) :: vmap
    real(dp),allocatable, target, dimension(:, :, :), intent(inout) :: vcoord
    integer(i4), allocatable, target, dimension(:), intent(inout) :: ind_map, &
         fhng_map, ehng_map
    real(dp), allocatable, target, dimension(:, :), intent(inout) :: &
         ind_coord, fhng_coord, ehng_coord
    integer(i4), intent(out) :: nind, nfhng, nehng
    integer(i4), intent(in) :: nvert, nelt, lnum_in, lnum_fh, lnum_eh, gdim
    integer :: il, jl, itmp, bndf, bnde, bndfp, bndep
    real(dp), allocatable, target, dimension(:, :) :: lind_coord, &
         lfhng_coord, lehng_coord
    integer(i4), allocatable, target, dimension(:) :: lind_map, lfhng_map, &
         lehng_map

    ! THIS IS NOT UPDATED YET!!! OLD METHOD WORKED WITH SINGLE PERIODICITY ONLY
    call neko_error('Periodic nodes not finalised yet')
    ! count periodic nodes
    ! independent
    nind = 0
    do il = 1, lnum_in
       if (ind_map(il) .gt. 0) then
          nind = nind + 1
          ind_map(il) = nind
       end if
    end do
    ! face hanging
    nfhng = 0
    do il = 1, lnum_fh
       if (fhng_map(il) .gt. 0) then
          nfhng = nfhng + 1
          fhng_map(il) = nfhng
       end if
    end do
    ! edge hanging
    nehng = 0
    do il = 1, lnum_eh
       if (ehng_map(il) .gt. 0) then
          nehng = nehng + 1
          ehng_map(il) = nehng
       end if
    end do

    ! shrink coordinate lists
    allocate(lind_coord(gdim, nind), lfhng_coord(gdim, nfhng), &
         lehng_coord(gdim, nehng))
    do il = 1, lnum_in
       if (ind_map(il) .gt. 0) then
          lind_coord(:, ind_map(il)) = ind_coord(:, il)
       end if
    end do
    call move_alloc(lind_coord, ind_coord)
    do il = 1, lnum_fh
       if (fhng_map(il) .gt. 0) then
          lfhng_coord(:, fhng_map(il)) = fhng_coord(:, il)
       end if
    end do
    call move_alloc(lfhng_coord, fhng_coord)
    do il = 1, lnum_eh
       if (ehng_map(il) .gt. 0) then
          lehng_coord(:, ehng_map(il)) = ehng_coord(:, il)
       end if
    end do
    call move_alloc(lehng_coord, ehng_coord)

    ! update vertex mapping
    bndf = lnum_in + lnum_fh
    bnde = lnum_in + lnum_fh + lnum_eh
    bndfp = lnum_in + nind + lnum_fh
    bndep = lnum_in + nind + lnum_fh + nfhng + lnum_eh
    do il = 1, nelt
       do jl = 1, nvert
          itmp = abs(vmap(jl, il))
          if (itmp .le. lnum_in) then
             ! independent node
             if (vmap(jl, il) .lt. 0) then
                ! periodic
                vmap(jl, il) = lnum_in + ind_map(itmp)
             else
                ! nothing to do
             end if
          else if (itmp .le. bndf) then
             ! face hanging node
             if (vmap(jl, il) .lt. 0) then
                ! periodic
                itmp = itmp - lnum_in
                vmap(jl, il) = bndfp + fhng_map(itmp)
             else
                ! update position including periodic nodes
                vmap(jl, il) = vmap(jl, il) + nind
             end if
          else if (itmp .le. bnde) then
             ! edge hanging node
             if (vmap(jl, il) .lt. 0) then
                ! periodic
                itmp = itmp - bndf
                vmap(jl, il) = bndep + ehng_map(itmp)
             else
                ! update position including periodic nodes
                vmap(jl, il) = vmap(jl, il) + nind + nfhng
             end if
          else
             call neko_error('Wrong vertex mapping index.')
          end if
       end do
    end do

    ! get periodic node to already existing node mapping
    allocate(lind_map(nind), lfhng_map(nfhng), lehng_map(nehng))
    do il = 1, lnum_in
       if (ind_map(il) .gt. 0) then
          lind_map(ind_map(il)) = il
       end if
    end do
    call move_alloc(lind_map, ind_map)
    do il = 1, lnum_fh
       if (fhng_map(il) .gt. 0) then
          lfhng_map(fhng_map(il)) = il
       end if
    end do
    call move_alloc(lfhng_map, fhng_map)
    do il = 1, lnum_eh
       if (ehng_map(il) .gt. 0) then
          lehng_map(ehng_map(il)) = il
       end if
    end do
    call move_alloc(lehng_map, ehng_map)

    ! missing global numbering of periodic nodes; will be done later

  end subroutine p4est_periodic_get

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
    if ((nlfam .gt. nelt).or.(mod(nlfam, nvert) .ne. 0)) &
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

  !> Get global numbering of hanging and periodic nodes
  !! @details p4est is a connectivity code, so it does not care about global
  !! numbering of hanging or periodic nodes just mapping them to independent
  !! nodes
  subroutine p4est_hng_periodic_gnum_get(mesh_new)
    type(manager_mesh_p4est_t), intent(inout) :: mesh_new

    ! This is not finished yet
    !call neko_error('Nodes global numbering not finalised yet')

    ! periodic node numbering
    ! hanging node numbering
    ! reordering of faces taking into account their multiplicity

  end subroutine p4est_hng_periodic_gnum_get

  !> Apply data read from mesh file to mesh manager structures
  subroutine p4est_mesh_file_apply(this)
    class(mesh_manager_p4est_t), intent(inout) :: this
    integer(i4) :: il, iel, ifc, ibc, nin, ninp, nhf, nhfp, nhe, nhep
    integer(i8) :: itmp8
    integer(i8), target, allocatable, dimension(:) :: gidx
    integer(i4), target, allocatable, dimension(:) :: imsh, igrp
    integer(i4), target, allocatable, dimension(:, :) :: crv, bc
    ! conversion from circular to symmetric notation
    integer(i4) , dimension(8), parameter :: cs_cnv = &
         [1, 2, 4, 3, 5, 6, 8, 7]

    select type(mesh => this%mesh)
    type is (manager_mesh_p4est_t)

       ! Fill all the information that may be not properly initialised
       ! neko supports just V-type elements
       allocate(gidx(mesh%nelt), imsh(mesh%nelt), igrp(mesh%nelt), &
            crv(mesh%nfcs, mesh%nelt), bc(mesh%nfcs, mesh%nelt))
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
       gidx(:) = mesh%gidx(:)
       igrp(:) = mesh%igrp(:)
       crv(:, :) = mesh%crv(:, :)
       bc(:, :) = mesh%bc(:, :)
       call wp4est_elm_ini_dat(c_loc(gidx), c_loc(imsh), c_loc(igrp), &
            c_loc(crv), c_loc(bc))

       ! check if applied boundary conditions are consistent with tree structure
       ! the ghost mesh was not destroyed at the end of data import, so no need
       ! to set it here
!       call wp4est_ghost_new()
       call wp4est_bc_check()
!       call wp4est_ghost_del()

       deallocate(gidx, imsh, igrp, crv, bc)
    end select

    ! correct vertex information in the mesh manager geometry
    !>
    !! @todo This is just a hack and should be done better in the future
    ! For now just local copy of vertices coordinates; it does not mater if the
    ! mesh file is fine, but in case of corrupted file could be important
    select type(geom => this%mesh%geom)
    type is (manager_geom_p4est_t)
       ! complete global geometry information
       if (geom%ifcomplete) then
          select type (ind => geom%ind)
          type is(manager_geom_node_ind_p4est_t)
             ! get limits for different node types
             nin = ind%lnum
             ninp = nin + ind%pnum
             nhf = ninp + geom%hng_fcs%lnum
             nhfp = nhf + geom%hng_fcs%pnum
             nhe = nhfp + geom%hng_edg%lnum
             nhep = nhe + geom%hng_edg%pnum
             do il = 1, geom%nel
                do ifc = 1, geom%nvrt
                   ibc = geom%vmap(ifc, il)
                   if (ibc .le. nin) then
                      ind%coord(:, ibc) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else if (ibc .le. ninp) then
                      ind%pcoord(:, ibc - nin) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else if (ibc .le. nhf) then
                      geom%hng_fcs%coord(:, ibc - ninp) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else if (ibc .le. nhfp) then
                      geom%hng_fcs%pcoord(:, ibc - nhf) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else if (ibc .le. nhe) then
                      geom%hng_edg%coord(:, ibc - nhfp) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else if (ibc .le. nhep) then
                      geom%hng_edg%pcoord(:, ibc - nhe) = &
                           this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
                   else
                      call neko_error('Inconsistent vertex to node mapping')
                   end if
                end do
             end do
          end select
       else ! ifcomplete
          ! only simplified info; local coordinates of element vertices
          do il = 1, geom%nel
             do ifc = 1, geom%nvrt
                geom%vcoord(:, ifc, il) = &
                     this%nmsh_mesh%hex(il)%v(cs_cnv(ifc))%v_xyz(:)
             end do
          end do
       end if
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

    write(log_buf, '(a)') 'p4est refinement/coarsening'
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

       if (itmp .eq. 0 ) then
          ! set refinement flag
          ifmod = .true.

          write(log_buf, '(a)') 'p4est refinement; mesh changed'
          call neko_log%message(log_buf, NEKO_LOG_INFO)

          ! import new data; just a simple geometry representation
          call this%import(.false.)

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


    call profiler_end_region("p4est refine", 31)

  end subroutine p4est_refine

  !> Construct neko mesh type based on mesh manager data
  !! @param[inout]   mesh     neko mesh type
  !! @param[in]      ifnmsh   use curvature information form nmsh file
  subroutine p4est_mesh_construct(this, mesh, ifnmsh)
    class(mesh_manager_p4est_t), intent(in) :: this
    type(mesh_t), intent(inout) :: mesh
    logical, intent(in) :: ifnmsh
    type(stack_i4t2_t), allocatable, dimension(:) :: obj
    type(tuple_i4_t), pointer, dimension(:) :: neighl
    integer :: il, jl, itmp, neighn
    character(len=LOG_SIZE) :: log_buf

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
    mesh%glb_nelv = int(this%mesh%gnelt, i4) ! global number of elements
    mesh%offset_el = int(this%mesh%gnelto, i4) ! global element offset;

    ! Geometrical context
    ! Local information only; global one will be related to connectivity
    select type (geom => this%mesh%geom)
    type is (manager_geom_p4est_t)
       ! get node information
       call p4est_node_fill(mesh, geom)

       ! get element information
       call p4est_element_fill(mesh, geom, this%mesh%gidx)
    end select

    ! Element deformation
    allocate(mesh%dfrmd_el(mesh%nelv))
    call mesh%curve%init(mesh%nelv)
    if (ifnmsh) then
       ! Get curvature information
       call p4est_curve_fill(mesh, this%nmsh_mesh)
    end if
    call mesh%curve%finalize() ! curved sides finalisation
    ! deformed element flag is not set 100% correct, as it neglects curve flags
    call mesh_generate_flags(mesh) ! deformation flag


    ! Connectivity context
    mesh%mpts = this%mesh%conn%vrt%lnum ! local number of unique vertices
    mesh%mfcs = this%mesh%conn%fcs%lnum ! local number of unique faces
    mesh%meds = this%mesh%conn%edg%lnum ! local number of unique edges

    ! global number of unique vertices, faces and edges
    mesh%glb_mpts = int(this%mesh%conn%vrt%gnum, i4)
    mesh%glb_mfcs = int(this%mesh%conn%fcs%gnum,  i4)
    mesh%glb_meds = int(this%mesh%conn%edg%gnum, i4)
    ! this doesn't seem to be used in original version and
    ! I do not need it for p4est either
    !mesh%max_pts_id

    ! mesh%htp - not used as vertices global id is already provided by p4est
    ! mesh%htf - not used as faces global id is already provided by p4est
    ! mesh%hte - not used as edge global id is already provided by p4est

    ! used internally to get connectivity only; not added
    ! mesh%facet_map

    ! Fill in neighbour and mesh distribution information
    select type (conn => this%mesh%conn)
    type is (manager_conn_p4est_t)
       ! used in: gather_scatter.f90;
       ! get neighbour rank info
       call p4est_rank_neighbour_fill(mesh, conn)

!!$       ! not used outside mesh.f90
!!$       ! get vertex neighbours
!!$       ! MORE INVESTIGATION NEEDED
!!$       select type (vrt => conn%vrt)
!!$       type is (manager_conn_obj_p4est_t)
!!$          call p4est_object_neighbour_fill(obj, mesh%nelv, &
!!$               this%mesh%gidx, vrt, conn%vmap, conn%nvrt)
!!$          ! extract global element index
!!$          if (allocated(obj)) then
!!$             allocate(mesh%point_neigh(size(obj)))
!!$             do il = 1, size(obj)
!!$                call mesh%point_neigh(il)%init()
!!$                neighn = obj(il)%size()
!!$                neighl => obj(il)%array()
!!$                do jl = 1, neighn
!!$                   call mesh%point_neigh(il)%push(neighl(jl)%x(2))
!!$                end do
!!$                call obj(il)%free()
!!$             end do
!!$             deallocate(obj)
!!$          end if
!!$       end select

       ! used in: gather_scatter.f90, tree_amg_multigrid.f90;
       ! MORE INVESTIGATION NEEDED
       ! get face neighbours
       select type (fcs => conn%fcs)
       type is (manager_conn_obj_p4est_t)
          call p4est_object_neighbour_fill(obj, mesh%nelv, &
               this%mesh%gidx, fcs, conn%fmap, conn%nfcs)
          ! extract global element index
          ! THIS WORKS FOR CONFORMING MESHES ONLY
          ! IN GENERAL THIS IS NOT A FAST AND SMOOTH CODE
          if (allocated(obj)) then
             allocate(mesh%facet_neigh(conn%nfcs, conn%nel))
             do il = 1, conn%nel
                do jl = 1, conn%nfcs
                   itmp = conn%fmap(jl, il)
                   neighn = obj(itmp)%size()
                   neighl => obj(itmp)%array()
                   select case(neighn)
                   case (1)
                      ! just one face
                      mesh%facet_neigh(jl, il) = 0
                   case (2)
                      if (neighl(1)%x(2) .eq. int(this%mesh%gidx(il), i4)) then
                         mesh%facet_neigh(jl, il) = neighl(2)%x(2)
                      else
                         mesh%facet_neigh(jl, il) = neighl(1)%x(2)
                      end if
                   case default
                      ! none or too many neighbours
                      call neko_error('This works for conformal mesh only.')
                   end select
                end do
             end do
             do il = 1, size(obj)
                call obj(il)%free()
             end do
             deallocate(obj)
          end if
       end select

!!$       ! used in dofmap.f90 through get_global and is_shared methods
!!$       ! MORE INVESTIGATION NEEDED
!!$       ! replaced by new mesh connectivity type
!!$       ! get mesh data distribution
!!$       call p4est_distdata_fill(mesh, conn)

       ! fill connectivity mapping information
       call p4est_conn_mapping_fill(mesh%conn, conn)
    end select

    ! Boundary conditions
    select type (meshmm => this%mesh)
    type is (manager_mesh_p4est_t)
       ! get boundary condition information
       call p4est_bc_fill(mesh, meshmm)
    end select

    ! final flags
    mesh%lconn = .true.
    ! not used
    !mesh%lnumr = .true.

  end subroutine p4est_mesh_construct

  !> Fill the geometrical nodes information
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      geom    geometry
  subroutine p4est_node_fill(mesh, geom)
    type(mesh_t), intent(inout) :: mesh
    type(manager_geom_p4est_t), intent(in) :: geom
    integer :: il, jl, itmp
    integer(i4) :: nidx

    ! complete global geometry information
    if (geom%ifcomplete) then
       select type (ind => geom%ind)
       type is(manager_geom_node_ind_p4est_t)
          ! local number of unique vertices
          mesh%gpts = ind%lnum +  ind%pnum + geom%hng_fcs%lnum + &
               geom%hng_fcs%pnum + geom%hng_edg%lnum+ geom%hng_edg%pnum

          ! fill the points
          ! order matters due to the element vertex mapping
          allocate(mesh%points(mesh%gpts))
          itmp = 0
          ! independent nodes
          if (ind%lnum .gt. 0) then
             do il = 1, ind%lnum
                itmp = itmp + 1
                nidx = int(ind%gidx(il), i4)
                call mesh%points(itmp)%init(ind%coord(:, il), nidx)
             end do
          end if
          ! independent periodic
          if (ind%pnum .gt. 0) then
             do il = 1, ind%pnum
                itmp = itmp + 1
                nidx = int(ind%pgidx(il), i4)
                call mesh%points(itmp)%init(ind%pcoord(:, il), nidx)
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
          ! face hanging periodic
          if (geom%hng_fcs%pnum .gt. 0) then
             do il = 1, geom%hng_fcs%pnum
                itmp = itmp + 1
                nidx = int(geom%hng_fcs%pgidx(il), i4)
                call mesh%points(itmp)%init(geom%hng_fcs%pcoord(:, il), nidx)
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
          ! edge hanging periodic
          if (geom%hng_edg%pnum .gt. 0) then
             do il = 1, geom%hng_edg%pnum
                itmp = itmp + 1
                nidx = int(geom%hng_edg%pgidx(il), i4)
                call mesh%points(itmp)%init(geom%hng_edg%pcoord(:, il), nidx)
             end do
          end if
       end select
    else ! ifcomplete
       ! only simplified info; local coordinates of element vertices

       ! THIS IS JUST A HACK, AS POINTS IN THE POINT LIST ARE NOT UNIQUE
       mesh%gpts = geom%nvrt * geom%nel
       ! fill the points
       ! order matters due to the way elements are initialised
       allocate(mesh%points(mesh%gpts))
       itmp = 0
       do il = 1, geom%nel
          do jl = 1, geom%nvrt
             itmp = itmp + 1
             ! No real global numbering
             call mesh%points(itmp)%init(geom%vcoord(:, jl, il), itmp)
          end do
       end do
    end if

  end subroutine p4est_node_fill

  !> Fill the geometrical element information
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      geom    geometry
  !! @param[in]      gidx    element global index
  subroutine p4est_element_fill(mesh, geom, gidx)
    type(mesh_t), intent(inout) :: mesh
    type(manager_geom_p4est_t), intent(in) :: geom
    integer(i8), dimension(:), intent(in) :: gidx
    integer :: il, itmp
    integer(i4) :: nidx

    ! Global to local element mapping
    call mesh%htel%init(mesh%nelv, il)

    ! Fill in element array
    allocate(mesh%elements(mesh%nelv))
    ! complete global geometry information
    if (geom%ifcomplete) then
       do il = 1, mesh%nelv
          allocate(hex_t::mesh%elements(il)%e)
          nidx = int(gidx(il), i4)
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
    else ! ifcomplete
       ! only simplified info; local coordinates of element vertices

       ! THIS IS JUST A HACK, AS POINTS IN THE POINT LIST ARE NOT UNIQUE
       do il = 1, mesh%nelv
          itmp = (il - 1) * geom%nvrt
          allocate(hex_t::mesh%elements(il)%e)
          nidx = int(gidx(il), i4)
          select type(ep => mesh%elements(il)%e)
          type is (hex_t)
             call ep%init(nidx, mesh%points(1 + itmp), mesh%points(2 + itmp), &
                  mesh%points(3 + itmp), mesh%points(4 + itmp), &
                  mesh%points(5 + itmp), mesh%points(6 + itmp), &
                  mesh%points(7 + itmp), mesh%points(8 + itmp))
          end select
          ! Global to local element mapping
          itmp = il
          call mesh%htel%set(nidx, itmp)
       end do
    end if

  end subroutine p4est_element_fill

  !> Fill the mesh type with curvature information
  !! @param[inout]   mesh      neko mesh type
  !! @param[in]      nmsh      nmsh mesh information
  subroutine p4est_curve_fill(mesh, nmsh)
    type(mesh_t), intent(inout) :: mesh
    type(nmsh_mesh_t), intent(in) :: nmsh
    integer :: il, element_offset4, itmp

    ! add curved sides
    if (nmsh%ncurve .gt. 0) then
       element_offset4 = int(nmsh%offset_el, i4)
       do il = 1, nmsh%ncurve
          itmp = nmsh%curve(il)%e - element_offset4
          call mesh%mark_curve_element(itmp, nmsh%curve(il)%curve_data, &
               nmsh%curve(il)%type)
       end do
    end if

  end subroutine p4est_curve_fill

  !> Fill the mesh type with rank neighbour information
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
       ! get rank list
       ! vertex information
       select type (vrt => conn%vrt)
       type is (manager_conn_obj_p4est_t)
          do il= 1, vrt%nrank ! mpi rank loop
             ! not sure I have to exclude myself
             if (vrt%rank(il) .ne. pe_rank) mesh%neigh(vrt%rank(il)) = .true.
          end do
       end select
       ! face information
       select type (fcs => conn%fcs)
       type is (manager_conn_obj_p4est_t)
          do il= 1, fcs%nrank ! mpi rank loop
             ! not sure I have to exclude myself
             if (fcs%rank(il) .ne. pe_rank) mesh%neigh(fcs%rank(il)) = .true.
          end do
       end select
       ! edge information
       select type (edg => conn%edg)
       type is (manager_conn_obj_p4est_t)
          do il= 1, edg%nrank ! mpi rank loop
             ! not sure I have to exclude myself
             if (edg%rank(il) .ne. pe_rank) mesh%neigh(edg%rank(il)) = .true.
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

  !> Fill the mesh type with vertex, face and edge neighbour information
  !! @param[inout]   obj     object (position, neighbour) list
  !! @param[in]      nel     local number of elements
  !! @param[in]      gidx    element global index
  !! @param[in]      objmm   mesh manager info about the object
  !! @param[in]      map     object mapping
  !! @param[in]      nmap    map array size
  subroutine p4est_object_neighbour_fill(obj, nel, gidx, objmm, map, nmap)
    type(stack_i4t2_t), allocatable, dimension(:), intent(inout) :: obj
    integer, intent(in) :: nel, nmap
    integer(i8), dimension(:), intent(in) :: gidx
    type(manager_conn_obj_p4est_t), intent(in) :: objmm
    integer(i4), dimension(:,:), intent(in) :: map

    integer(i4) :: il, jl, kl, ll, ml
    integer(i4) :: itmp, ierr
    type(tuple_i4_t) :: ttmp
    ! offest in send/receive buffers
    integer(i4), allocatable, dimension(:) :: cmoff, cmoffr
    integer(i4) :: neighn
    type(tuple_i4_t), pointer, dimension(:) :: neighl
    integer(i4), allocatable, dimension(:, :) :: rbuf, sbuf ! snd/rcv buffers
    type(MPI_Request), allocatable, dimension(:) :: request ! MPI request
    type(MPI_Status), allocatable, dimension(:) :: status ! MPI status

    ! Get list of neighbours (elements a given object belongs to)
    if (allocated(obj)) then
       do il = 1, size(obj)
          call obj(il)%free()
       end do
       deallocate(obj)
    end if
    ! local object neighbours
    allocate(obj(objmm%lnum))
    do il = 1, objmm%lnum
       call obj(il)%init()
    end do
    do il = 1, nel
       ! get global element id
       itmp = int(gidx(il), i4)
       do jl = 1, nmap
          ttmp%x = [jl, itmp]
          call obj(map(jl, il))%push(ttmp)
       end do
    end do

    ! SHOULD THE COMMUNICATION PATTERN BE CHANGED?
    ! non-local object neighbours
    ! count size of send/receive buffers
    allocate(cmoff(objmm%nrank), cmoffr(objmm%nrank))
    itmp = 1
    cmoff(itmp) = 1
    do il= 1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! mpi rank loop
          cmoff(itmp) = cmoff(itmp - 1)
          do jl = objmm%off(il), objmm%off(il + 1) -1
             ! sum number of local object neighbours for a given rank
             neighn = obj(objmm%share(jl))%size()
             cmoff(itmp) = cmoff(itmp) + 1 + neighn
          end do
       end if
    end do
    allocate(request(objmm%nrank - 1), status(objmm%nrank - 1))

    ! Exchange size information
    ! set non-blocking receive
    itmp = 0
    do il=1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! count messages
          call MPI_IRecv(cmoffr(itmp), 1, MPI_INTEGER, objmm%rank(il), &
               objmm%rank(il), NEKO_COMM, request(itmp), ierr)
       end if
    end do

    ! send buffer size
    itmp = 0
    do il=1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! count messages
          jl = cmoff(itmp + 1) - cmoff(itmp)
          call MPI_Send(jl, 1, MPI_INTEGER, objmm%rank(il), &
               pe_rank, NEKO_COMM, ierr)
       end if
    end do

    ! finalise communication
    call MPI_Waitall(objmm%nrank - 1, request, status, ierr)

    itmp = cmoffr(1)
    cmoffr(1) = 1
    do il=2, objmm%nrank ! mpi rank loop
       ierr = cmoffr(il)
       cmoffr(il) = cmoffr(il - 1) + itmp
       itmp = ierr
    end do
    allocate(rbuf(2, cmoffr(objmm%nrank)), sbuf(2, cmoff(objmm%nrank)))

    ! set non-blocking receive
    itmp = 0
    do il=1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! count messages
          jl = 2 * (cmoffr(itmp + 1) - cmoffr(itmp))
          call MPI_IRecv(rbuf(:, cmoffr(itmp) : cmoffr(itmp + 1) - 1), jl, &
               MPI_INTEGER, objmm%rank(il), jl, NEKO_COMM, request(itmp), ierr)
       end if
    end do

    ! redistribute object information
    itmp = 0
    do il = 1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! count messages
          ll = cmoff(itmp)
          do jl = objmm%off(il), objmm%off(il + 1) - 1
             ! extract local neighbours
             neighn = obj(objmm%share(jl))%size()
             neighl => obj(objmm%share(jl))%array()
             sbuf(1, ll) = neighn ! local number of neighbours
             ! object global id
             sbuf(2, ll) = int(objmm%gidx(objmm%share(jl)), i4)
             ll = ll + 1
             do kl = 1, neighn
                sbuf(:, ll) = neighl(kl)%x(:)
                ll = ll + 1
             end do
          end do
          ! sanity check
          if (ll .ne. cmoff(itmp + 1)) call neko_error('Inconsistent number of &
               &local object neighbours; receive.')
          jl = 2 * (cmoff(itmp+1) - cmoff(itmp))
          call MPI_Send(sbuf(:, cmoff(itmp) : cmoff(itmp + 1) - 1), jl, &
               MPI_INTEGER, objmm%rank(il), jl, NEKO_COMM, ierr)
       end if
    end do

    ! finalise communication
    call MPI_Waitall(objmm%nrank - 1, request, status, ierr)

    ! extract data
    itmp = 0
    do il = 1, objmm%nrank ! mpi rank loop
       if (objmm%rank(il) .ne. pe_rank) then
          itmp = itmp + 1 ! count messages
          ll = cmoffr(itmp)
          do jl = objmm%off(il), objmm%off(il + 1) - 1
             ! check object global id
             if (objmm%gidx(objmm%share(jl)) .ne. rbuf(2, ll)) &
                  call neko_error('Inconsistent global vertex number in &
                  &neighbour exchange.')
             neighn = rbuf(1, ll)
             ll = ll + 1
             do kl = 1, neighn
                ttmp%x(:) = -rbuf(:, ll)
                call obj(objmm%share(jl))%push(ttmp)
                ll = ll + 1
             end do
          end do
          ! sanity check
          if (ll .ne. cmoffr(itmp + 1)) call neko_error('Inconsistent number &
               &of local vertex neighbours; extract.')
       end if
    end do

    deallocate(cmoff, request, status, rbuf, sbuf)

  end subroutine p4est_object_neighbour_fill

  !> Fill the mesh type with mesh data distribution information
  !! @param[inout]   mesh    neko mesh type
  !! @param[in]      conn    connectivity
  subroutine p4est_distdata_fill(mesh, conn)
    type(mesh_t), intent(inout) :: mesh
    type(manager_conn_p4est_t), intent(in) :: conn
    integer :: il, jl, itmp, neighn
    type(stack_i4t2_t), allocatable, dimension(:) :: obj
    type(tuple_i4_t) :: ttmp
    type(tuple_i4_t), pointer, dimension(:) :: neighl

    ! initialise connectivity information
    call mesh%ddata%init()

    ! Add shared vertices
    select type (vrt => conn%vrt)
    type is (manager_conn_obj_p4est_t)
       ! Local id
       do il = 1, vrt%nrank
          if (vrt%rank(il) == pe_rank) then
             do jl = vrt%off(il), vrt%off(il + 1) - 1
                call mesh%ddata%set_shared_point(vrt%share(jl))
             end do
             exit
          end if
       end do
    end select

    ! Add shared faces
    select type (fcs => conn%fcs)
    type is (manager_conn_obj_p4est_t)
       ! Local id
       do il = 1, fcs%nrank
          if (fcs%rank(il) == pe_rank) then
             itmp = il
             do jl = fcs%off(il), fcs%off(il + 1) - 1
                call mesh%ddata%set_shared_facet(fcs%share(jl))
             end do
             exit
          end if
       end do
       ! Get faces that are shared
       ! This part is to some extent already done in p4est_object_neighbour_fill
       ! (local element number missing, so maybe could be taken from previous
       ! steps.
       allocate(obj(fcs%lnum))
       do il = 1, fcs%lnum
          call obj(il)%init()
       end do
       do il = 1, conn%nel
          do jl = 1, conn%nfcs
             ttmp%x = [il, jl]
             call obj(conn%fmap(jl, il))%push(ttmp)
          end do
       end do
       ! extract positions
       do il = fcs%off(itmp), fcs%off(itmp + 1) - 1
          jl = fcs%share(il)
          neighn = obj(jl)%size()
          neighl => obj(jl)%array()
          do jl = 1, neighn
             call mesh%ddata%set_shared_el_facet(neighl(jl)%x(1), &
                  neighl(jl)%x(2))
          end do
       end do
       ! Local to global id mapping
       allocate(mesh%ddata%local_to_global_facet(fcs%lnum))
       do il = 1, fcs%lnum
          ! global face id; notice type casting
          itmp = int(fcs%gidx(il),i4)
          call mesh%ddata%set_local_to_global_facet(il, itmp)
       end do

       do il = 1, size(obj)
          call obj(il)%free()
       end do
       deallocate(obj)
    end select

    ! Add shared edges
    select type (edg => conn%edg)
    type is (manager_conn_obj_p4est_t)
       ! Local id
       do il = 1, edg%nrank
          if (edg%rank(il) == pe_rank) then
             do jl = edg%off(il), edg%off(il + 1) - 1
                call mesh%ddata%set_shared_edge(edg%share(jl))
             end do
             exit
          end if
       end do
       ! Local to global id mapping
       allocate(mesh%ddata%local_to_global_edge(edg%lnum))
       do il = 1, edg%lnum
          ! global face id; notice type casting
          itmp = int(edg%gidx(il),i4)
          call mesh%ddata%set_local_to_global_edge(il, itmp)
       end do
    end select

    ! not used
    !mesh%ldist = .true.

  end subroutine p4est_distdata_fill

  !> Fill the mesh type with connectivity mapping information
  !! @param[inout]   conn    neko mesh type connectivity
  !! @param[in]      connmm  mesh manager  connectivity
  subroutine p4est_conn_mapping_fill(conn, connmm)
    type(mesh_conn_t), intent(inout) :: conn
    type(manager_conn_p4est_t), intent(in) :: connmm
    logical, allocatable, dimension(:) :: share
    integer :: il, jl

    call conn%init(connmm%tdim, connmm%nel, connmm%hngel)

    select type (vrt => connmm%vrt)
    type is (manager_conn_obj_p4est_t)
       allocate(share(vrt%lnum))
       share(:) = .false.
       do il = 1, vrt%nrank
          if (vrt%rank(il) == pe_rank) then
             do jl = vrt%off(il), vrt%off(il + 1) - 1
                share(vrt%share(jl)) = .true.
             end do
             exit
          end if
       end do
       call conn%vrt%init(vrt%lnum, vrt%gnum, connmm%nel, connmm%nvrt, &
            vrt%gidx, share, connmm%vmap)
       deallocate(share)
    end select

    select type (edg => connmm%edg)
    type is (manager_conn_obj_p4est_t)
       allocate(share(edg%lnum))
       share(:) = .false.
       do il = 1, edg%nrank
          if (edg%rank(il) == pe_rank) then
             do jl = edg%off(il), edg%off(il + 1) - 1
                share(edg%share(jl)) = .true.
             end do
             exit
          end if
       end do
       call conn%edg%init(edg%lnum, edg%gnum, connmm%nel, connmm%nedg, &
            edg%gidx, share, connmm%emap, connmm%ealgn, connmm%hnged)
       deallocate(share)
    end select

    select type (fcs => connmm%fcs)
    type is (manager_conn_obj_p4est_t)
       allocate(share(fcs%lnum))
       share(:) = .false.
       do il = 1, fcs%nrank
          if (fcs%rank(il) == pe_rank) then
             do jl = fcs%off(il), fcs%off(il + 1) - 1
                share(fcs%share(jl)) = .true.
             end do
             exit
          end if
       end do
       call conn%fcs%init(fcs%lnum, fcs%gnum, connmm%nel, connmm%nfcs, &
            fcs%gidx, share, connmm%fmap, connmm%falgn, connmm%hngfc)
       deallocate(share)
    end select

  end subroutine p4est_conn_mapping_fill

  !> Fill the mesh type with boundary condition information
  !! @param[inout]   mesh      neko mesh type
  !! @param[in]      meshmm    mesh manager mesh information
  subroutine p4est_bc_fill(mesh, meshmm)
    type(mesh_t), intent(inout) :: mesh
    type(manager_mesh_p4est_t), intent(in) :: meshmm
    integer :: il, jl, itmp
    integer(i4) :: p_f, p_e ! dummy variables for setting periodic bc
    integer(i4), dimension(4) :: pt_id ! dummy array for setting periodic bc

    ! Initialize element boundary condition
    allocate(mesh%facet_type(meshmm%nfcs, meshmm%nelt))
    mesh%facet_type(:, :) = 0

    allocate(mesh%labeled_zones(NEKO_MSH_MAX_ZLBLS))
    do il = 1, NEKO_MSH_MAX_ZLBLS
       call mesh%labeled_zones(il)%init(mesh%nelv)
    end do

    call mesh%periodic%init(mesh%nelv)

    ! Set BC
    do il = 1, meshmm%nelt
       do jl = 1, meshmm%nfcs
          select case(meshmm%bc(jl, il))
          case(-1)
             ! this is my periodic bc mark, which is inconsistent with neko (5).
             ! I don't think I have to do here much more than just mark a face,
             ! as periodicity is already taken into account in communicator
             ! IT WOULD BE SUFFICIENT TO HAVE THIS ZONE TO BE LIKE ANY OTHER ONE
             call mesh%mark_periodic_facet(jl, il, p_f, p_e, pt_id)
             ! There is no need to call mesh apply, as I do not stick to
             ! node numbering anyhow.
          case (1 : NEKO_MSH_MAX_ZLBLS) ! everything elese marked as labeled bc
             itmp = meshmm%bc(jl,il) ! once again problem with inout attribute
             call mesh%mark_labeled_facet(jl, il, itmp)
          end select
       end do
    end do

    ! Finalise boundary conditions
    call mesh%periodic%finalize()

    do il = 1, NEKO_MSH_MAX_ZLBLS
       call mesh%labeled_zones(il)%finalize()
    end do

  end subroutine p4est_bc_fill

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
  subroutine p4est_import(this, ifcomplete)
    class(mesh_manager_p4est_t), intent(inout) :: this
    logical, intent(in) :: ifcomplete

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
  !! @param[in]      ifnmsh   use curvature information form nmsh file
  subroutine p4est_mesh_construct(this, mesh, ifnmsh)
    class(mesh_manager_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: mesh
    logical, intent(in) :: ifnmsh

    call neko_error('p4est mesh manager must be compiled with p4est support.')

  end subroutine p4est_mesh_construct
#endif

end module mesh_manager_p4est
