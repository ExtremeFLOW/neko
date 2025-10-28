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
module p4est
  use mpi_f08

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

  ! connectivity parameter arrays
  ! face vertices
  integer, parameter, dimension(4,6) :: p4_vface = reshape(&
       (/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),&
       shape(p4_vface))

  ! edge vertices
  integer, parameter, dimension(2,12) :: p4_vedge  = reshape(&
       (/ 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7&
       , 4,8 /),shape(p4_vedge))

  ! edge related faces
  integer, parameter, dimension(2,12) :: p4_eface  = reshape(&
       (/ 3,5 , 4,5 , 3,6 , 4,6 , 1,5 , 2,5 , 1,6 , 2,6 , 1,3 , 2,3 , 1,4&
       , 2,4 /),shape(p4_eface))

  ! corner related faces
  integer, parameter, dimension(3,8) :: p4_cface = reshape(&
       (/ 1,3,5 , 2,3,5 , 1,4,5 , 2,4,5 , 1,3,6 , 2,3,6 , 1,4,6 , 2,4,6 /),&
       shape(p4_cface))

  ! corner related edges
  integer, parameter, dimension(3,8) :: p4_cedge = reshape(&
       (/ 1,5,9 , 1,6,10 , 2,5,11 , 2,6,12 , 3,7,9 , 3,8,10 , 4,7,11 ,&
        4,8,12 /),shape(p4_cedge))

  ! corner to face corner
  integer, parameter, dimension(6,8) :: p4_cfcrn = reshape(&
       (/ 1,-1, 1,-1, 1,-1 , -1, 1, 2,-1, 2,-1 ,  2,-1,-1, 1, 3,-1 &
       , -1, 2,-1, 2, 4,-1 ,  3,-1, 3,-1,-1, 1 , -1, 3, 4,-1,-1, 2 &
       ,  4,-1,-1, 3,-1, 3 , -1, 4,-1, 4,-1, 4 /),shape(p4_cfcrn))

  ! to calculate neighbour face corner
  integer, parameter, dimension(6,6) :: p4_rt =reshape( &
       (/ 1,2,2,1,1,2 , 3,1,1,2,2,1 , 3,1,1,2,2,1 , 1,3,3,1,1,2 &
       , 1,3,3,1,1,2 , 3,1,1,3,3,1 /),shape(p4_rt))
  integer, parameter, dimension(4,3) :: p4_qt = reshape(&
       (/ 2,3,6,7 , 1,4,5,8 , 1,5,4,8 /),shape(p4_qt))
  integer, parameter, dimension(4,8) :: p4_pt = reshape(&
       (/ 1,2,3,4 , 1,3,2,4 , 2,1,4,3 , 2,4,1,3 , 3,1,4,2 , 3,4,1,2 &
       , 4,2,3,1 , 4,3,2,1 /),shape(p4_pt))

  ! default log threshold - production; for more info see sc.h
  integer, parameter :: p4_lp_production = 6

  ! set of interfaces to call p4est functions
  interface
     subroutine wp4est_init(fmpicomm, catch_signals, print_backtrace, &
          log_threshold) bind(c, name = 'wp4est_init')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: fmpicomm, catch_signals, &
            print_backtrace, log_threshold
     end subroutine wp4est_init

     subroutine wp4est_finalize(log_priority) bind(c, name = 'wp4est_finalize')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: log_priority
     end subroutine wp4est_finalize

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

     subroutine wp4est_msh_get_size(mdim, nelgt, nelgto, nelt, nelv, maxl) &
          bind(c, name = 'wp4est_msh_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: mdim, nelv, maxl
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

     subroutine wp4est_elm_get_dat(gidx, level, igrp, crv, bc, coord, falg) &
          bind(c, name = 'wp4est_elm_get_dat')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: gidx, level, igrp, crv, bc, coord, falg
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

     subroutine wp4est_tree_copy(quad_data) &
          bind(c, name = 'wp4est_tree_copy')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_copy

     subroutine wp4est_tree_check(check, quad_data) &
          bind(c, name = 'wp4est_tree_check')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: check
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_check

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

     subroutine wp4est_msh_get_hst(map_nr, rfn_nr, crs_nr, &
          elgl_map, elgl_rfn, elgl_crs) &
          bind(c, name = 'wp4est_msh_get_hst')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: map_nr, rfn_nr, crs_nr
       type(c_ptr), value :: elgl_map, elgl_rfn, elgl_crs
     end subroutine wp4est_msh_get_hst

     subroutine wp4est_vtk_write(filename) &
          bind(c, name = 'wp4est_vtk_write')
       USE, INTRINSIC :: ISO_C_BINDING
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_vtk_write

  end interface

contains


end module p4est
