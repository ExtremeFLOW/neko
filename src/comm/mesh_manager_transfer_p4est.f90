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
!> Types for mesh manager data redistribution
module mesh_manager_transfer_p4est
  use mpi_f08
  use num_types, only : i4, i8, rp, dp
  use comm, only : NEKO_COMM, pe_rank
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use profiler, only : profiler_start_region, profiler_end_region
  use mesh_manager_transfer, only : mesh_manager_transfer_t
  use manager_mesh, only : manager_mesh_t
  use manager_mesh_p4est, only : manager_mesh_p4est_t

  implicit none

  private

  !> p4est mesh manager data redistribution routines
  type, public, extends(mesh_manager_transfer_t) :: &
       mesh_manager_transfer_p4est_t
     !> Local element number on mesh manager side
     integer(i4) :: nelt_mm
     !> Element global id on mesh manager side
     integer(i8), allocatable, dimension(:) :: gidx_mm
     !> Local element number on neko side
     integer(i4) :: nelt_neko
     !> Element global id on neko side
     integer(i8), allocatable, dimension(:) :: gidx_neko
     !> Forward communication (mm->neko); local element number and MPI rank
     integer(i4), allocatable, dimension(:, :) :: fwd_cmm
     !> Backward communication (mm<-neko); local element number and MPI rank
     integer(i4), allocatable, dimension(:, :) :: bwd_cmm
   contains
     !> Destructor.
     procedure, pass(this) :: free => p4est_free
     !> Get element distribution
     procedure, pass(this) :: elem_dist_get => p4est_elem_dist_get
  end type mesh_manager_transfer_p4est_t

contains

  !> Destructor.
  subroutine p4est_free(this)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nelt_mm = 0
    this%nelt_neko = 0

    if (allocated(this%gidx_mm)) deallocate(this%gidx_mm)
    if (allocated(this%gidx_neko)) deallocate(this%gidx_neko)
    if (allocated(this%fwd_cmm)) deallocate(this%fwd_cmm)
    if (allocated(this%bwd_cmm)) deallocate(this%bwd_cmm)

  end subroutine p4est_free

  !> Get element distribution
  !! @param[in]   mesh     mesh manager mesh
  subroutine p4est_elem_dist_get(this, mesh)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    class(manager_mesh_t), intent(in) :: mesh
    integer :: il

    ! clean old element distribution
    this%nelt_mm = 0
    this%nelt_neko = 0
    if (allocated(this%gidx_mm)) deallocate(this%gidx_mm)
    if (allocated(this%gidx_neko)) deallocate(this%gidx_neko)
    if (allocated(this%fwd_cmm)) deallocate(this%fwd_cmm)
    if (allocated(this%bwd_cmm)) deallocate(this%bwd_cmm)

    ! get new one
    select type (mesh)
    type is (manager_mesh_p4est_t)
       ! is partitioning required
       if (this%ifpartition) then
          this%nelt_mm = mesh%nelt
          allocate(this%gidx_mm(this%nelt_mm), this%fwd_cmm(2, this%nelt_mm))
          this%gidx_mm(:) = mesh%gidx(:)

          call neko_error('Nothing done yet')
       else
          ! no partitioning; mesh manager and neko share element distribution
          this%nelt_mm = mesh%nelt
          this%nelt_neko = mesh%nelt
          allocate(this%gidx_mm(this%nelt_mm), this%fwd_cmm(2, this%nelt_mm), &
               this%gidx_neko(this%nelt_neko), this%bwd_cmm(2, this%nelt_neko))
          do il = 1, mesh%nelt
             this%fwd_cmm(1, il) = il
             this%fwd_cmm(2, il) = pe_rank
          end do
          this%gidx_mm(:) = mesh%gidx(:)
          this%gidx_neko(:) = this%gidx_mm(:)
          this%bwd_cmm(:, :) = this%fwd_cmm(:, :)
       end if
    class default
       call neko_error('Wrong mesh type')
    end select

  end subroutine p4est_elem_dist_get

end module mesh_manager_transfer_p4est
