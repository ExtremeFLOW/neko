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
     ! element reconstruction data
     !> Number of untouched elements
     integer(i4) :: same_nr
     !> Old global id of untouched elements
     integer(i8), allocatable, dimension(:) :: same_gidx
     !> Mapping of untouched elements
     ! 1 - old local id; 2 - old MPI rank
     integer(i4), allocatable, dimension(:, :) :: same
     !> Number of refined elements
     integer(i4) :: refine_nr
     !> Global id of refined elements
     ! 1 - current global id; 2 - old parent global id
     integer(i8), allocatable, dimension(:, :) :: refine_gidx
     !> Mapping of refined elements
     ! 1 - old parent local id; 2 - old parent MPI rank; 3 - child position
     integer(i4), allocatable, dimension(:, :) :: refine
     !> Number of coarsened elements
     integer(i4) :: coarsen_nr
     !> Global id of coarsened elements
     ! 1 - new global id; 2 - old child global id
     integer(i8), allocatable, dimension(:, :, :) :: coarsen_gidx
     !> Mapping of coarsened elements
     ! 1 - old child local id; 2 - old child MPI rank
     integer(i4), allocatable, dimension(:, :, :) :: coarsen
   contains
     !> Destructor.
     procedure, pass(this) :: free => p4est_free
     !> Get element distribution
     procedure, pass(this) :: elem_dist_construct => p4est_elem_dist_construct
     !> Backward communication for element refinement flag
     procedure, pass(this) :: ref_mark_transfer => p4est_ref_mark_transfer
     !> Provide neko element distribution to p4est
     procedure, pass(this) :: neko_elem_dist_set => p4est_neko_elem_dist_set
     !> Set element distribution for field reconstruction
     procedure, pass(this) :: reconstruct_data_set => p4est_reconstruct_data_set
  end type mesh_manager_transfer_p4est_t

contains

  !> Destructor.
  subroutine p4est_free(this)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nelt_mm = 0
    this%nelt_neko = 0
    this%same_nr = 0
    this%refine_nr = 0
    this%coarsen_nr = 0

    if (allocated(this%gidx_mm)) deallocate(this%gidx_mm)
    if (allocated(this%gidx_neko)) deallocate(this%gidx_neko)
    if (allocated(this%fwd_cmm)) deallocate(this%fwd_cmm)
    if (allocated(this%bwd_cmm)) deallocate(this%bwd_cmm)
    if (allocated(this%same_gidx)) deallocate(this%same_gidx)
    if (allocated(this%same)) deallocate(this%same)
    if (allocated(this%refine_gidx)) deallocate(this%refine_gidx)
    if (allocated(this%refine)) deallocate(this%refine)
    if (allocated(this%coarsen_gidx)) deallocate(this%coarsen_gidx)
    if (allocated(this%coarsen)) deallocate(this%coarsen)

  end subroutine p4est_free

  !> Get element distribution
  !! @param[in]   mesh     mesh manager mesh
  subroutine p4est_elem_dist_construct(this, mesh)
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

  end subroutine p4est_elem_dist_construct

  !> Backward communication for element refinement flag
  !! @parameter[in]   ref_mark    element refinement flag; neko distribution
  !! @parameter[out]  pref_mark   element refinement flag; p4est distribution
  subroutine p4est_ref_mark_transfer(this, ref_mark, pref_mark)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: ref_mark
    integer(i4), target, allocatable, dimension(:), intent(out) :: pref_mark

    ! check size on neko side
    if (size(ref_mark) .ne. this%nelt_neko) &
         call neko_error('Inconsistent ref_mark array size')

    ! p4est distribution
    allocate(pref_mark(this%nelt_mm))

    ! partitioned mesh
    if (this%ifpartition) then
       ! backward communication

       call neko_error('Nothing done yet')
    else
       ! the same distribution
       pref_mark(:) = ref_mark(:)
    end if

  end subroutine p4est_ref_mark_transfer

  !> Provide neko element distribution to p4est
  !! @parameter[out]  pel_gnum       element global number
  !! @parameter[out]  pel_lnum       element local number on neko side
  !! @parameter[out]  pel_nid        element owner id on neko side
  subroutine p4est_neko_elem_dist_set(this, pel_gnum, pel_lnum, pel_nid)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    integer(i8), target, allocatable, dimension(:), intent(out) :: pel_gnum
    integer(i4), target, allocatable, dimension(:), intent(out) :: pel_lnum, &
         pel_nid

    allocate(pel_gnum(this%nelt_mm), pel_lnum(this%nelt_mm), &
         pel_nid(this%nelt_mm))
    pel_gnum(:) = this%gidx_mm(:)
    pel_lnum(:) = this%fwd_cmm(1,:)
    pel_nid(:) = this%fwd_cmm(2,:)

  end subroutine p4est_neko_elem_dist_set

  !> Get element distribution for field reconstruction
  subroutine p4est_reconstruct_data_set(this, same_nr, same_gidx, same, &
       refine_nr, refine_gidx, refine, coarsen_nr, coarsen_gidx, coarsen)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: same_nr, refine_nr, coarsen_nr
    integer(i8), allocatable, dimension(:), intent(inout) :: same_gidx
    integer(i8), allocatable, dimension(:, :), intent(inout) :: refine_gidx
    integer(i8), allocatable, dimension(:, :, :), intent(inout) :: coarsen_gidx
    integer(i4), allocatable, dimension(:, :), intent(inout) :: same, refine
    integer(i4), allocatable, dimension(:, :, :), intent(inout) :: coarsen

    this%same_nr = same_nr
    this%refine_nr = refine_nr
    this%coarsen_nr = coarsen_nr

    call move_alloc(same_gidx, this%same_gidx)
    call move_alloc(same, this%same)
    call move_alloc(refine_gidx, this%refine_gidx)
    call move_alloc(refine, this%refine)
    call move_alloc(coarsen_gidx, this%coarsen_gidx)
    call move_alloc(coarsen, this%coarsen)

  end subroutine p4est_reconstruct_data_set

end module mesh_manager_transfer_p4est
