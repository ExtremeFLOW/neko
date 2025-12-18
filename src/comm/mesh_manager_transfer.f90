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
!> Abstract types for mesh manager data redistribution
module mesh_manager_transfer
  use mpi_f08
  use num_types, only : rp
  use comm, only : NEKO_COMM, pe_rank
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use manager_mesh, only : manager_mesh_t

  implicit none

  private

  !> Base abstract type for mesh manager data redistribution
  type, abstract, public :: mesh_manager_transfer_t
     !> Mesh partitioning flag
     logical :: ifpartition
   contains
     !> Initialise the base type.
     procedure, pass(this) :: init_base => mesh_manager_init_base
     !> Free the base type.
     procedure, pass(this) :: free_base => mesh_manager_free_base
     !> Destructor.
     procedure(mesh_manager_free), pass(this), deferred :: free
     !> Construct new element distribution
     procedure(mesh_manager_element_dist), pass(this), deferred :: &
          elem_dist_construct
     !> Get refinement/coarsening vector sizes and mappings
     procedure(mesh_manager_vector_map), pass(this), deferred :: &
          vector_map
     !> Free refinement/coarsening vector sizes and mappings
     procedure(mesh_manager_free), pass(this), deferred :: &
          vector_map_free
     !> Construct vectors for refinement/coarsening
     procedure(mesh_manager_vector_constr), pass(this), deferred :: &
          vector_constr
  end type mesh_manager_transfer_t

  abstract interface
     !> Destructor
     subroutine mesh_manager_free(this)
       import mesh_manager_transfer_t
       class(mesh_manager_transfer_t), intent(inout) :: this
     end subroutine mesh_manager_free

     !> Get element distribution
     !! @param[in]   mesh     mesh manager mesh
     subroutine mesh_manager_element_dist(this, mesh)
       import mesh_manager_transfer_t, manager_mesh_t
       class(mesh_manager_transfer_t), intent(inout) :: this
       class(manager_mesh_t), intent(in) :: mesh
     end subroutine mesh_manager_element_dist

     !> Get refinement/coarsening vectors sizes and mappings
     !! @param[out]    nold    old element number
     !! @param[out]    nnew    new element number
     !! @param[out]    nref    refinement mapping size
     !! @param[out]    ncrs    coarsening mapping size
     !! @param[inout]  rmap    refinement mapping (element and child position)
     !! @param[inout]  cmap    coarsening mapping (element position)
     subroutine mesh_manager_vector_map(this, nold, nnew, nref, ncrs, rmap, &
          cmap)
       import mesh_manager_transfer_t, rp
       class(mesh_manager_transfer_t), intent(inout) :: this
       integer, intent(out) :: nold, nnew, nref, ncrs
       integer, dimension(:, :), allocatable, intent(inout) :: rmap
       integer, dimension(:), allocatable, intent(inout) :: cmap
     end subroutine mesh_manager_vector_map

     !> Construct vectors for refinement/coarsening
     !! @param[in]   vin     original vector
     !! @param[out]  vout    output vector for refinement
     !! @param[out]  vcrs    vector with additional data for coarsening
     subroutine mesh_manager_vector_constr(this, vin, vout, vcrs)
       import mesh_manager_transfer_t, rp
       class(mesh_manager_transfer_t), intent(inout) :: this
       real(rp), dimension(:, :, :, :), intent(in) :: vin
       real(rp), dimension(:, :, :, :), intent(out) :: vout, vcrs
     end subroutine mesh_manager_vector_constr
  end interface

contains

  !> Initialise the base type.
  !! @param[in]   ifpartition     partitioning flag
  subroutine mesh_manager_init_base(this, ifpartition)
    class(mesh_manager_transfer_t), intent(inout) :: this
    logical, intent(in) :: ifpartition

    call this%free_base()

    this%ifpartition = ifpartition

  end subroutine mesh_manager_init_base

  !> Free the base type.
  subroutine mesh_manager_free_base(this)
    class(mesh_manager_transfer_t), intent(inout) :: this

    this%ifpartition = .false.

  end subroutine mesh_manager_free_base

end module mesh_manager_transfer
