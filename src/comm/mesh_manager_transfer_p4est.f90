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
  use mpi_f08, only : MPI_Allreduce, MPI_Sendrecv, MPI_Irecv, MPI_Isend, &
       MPI_Wait, MPI_Get_count, MPI_IN_PLACE, MPI_LOGICAL, MPI_INTEGER, &
       MPI_INTEGER8, MPI_LOR, MPI_MAX, MPI_STATUS_IGNORE, MPI_Status, &
       MPI_Request
  use num_types, only : i4, i8, rp, dp
  use comm, only : NEKO_COMM, pe_rank, pe_size, MPI_REAL_PRECISION
  use logger, only : neko_log, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEBUG, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use profiler, only : profiler_start_region, profiler_end_region
  use math, only : swap, sort, sort_tuple
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
     !> Old local element number on neko side
     integer(i4) :: nelt_old
     !> Old element global id on neko side
     integer(i8), allocatable, dimension(:) :: gidx_old
     !> Number of untouched elements
     integer(i4) :: same_nr
     !> Old global id of untouched elements
     ! of size nelt_neko
     integer(i8), allocatable, dimension(:) :: same_gidx
     !> Mapping of untouched elements
     ! of size nelt_neko
     ! if same_gidx /= 0:
     ! 1 - old local id; 2 - old MPI rank
     ! if same_gidx == 0:
     ! 1 - corresponding position in refine;
     ! 2 - corresponding position in coarsen
     integer(i4), allocatable, dimension(:, :) :: same
     !> Number of refined elements
     integer(i4) :: refine_nr
     !> Global id of refined elements
     ! of size refine_nr
     ! 1 - current global id; 2 - old parent global id
     integer(i8), allocatable, dimension(:, :) :: refine_gidx
     !> Mapping of refined elements
     ! of size refine_nr
     ! 1 - old parent local id; 2 - old parent MPI rank; 3 - child position
     integer(i4), allocatable, dimension(:, :) :: refine
     !> Number of coarsened elements
     integer(i4) :: coarsen_nr
     !> Number of children
     integer(i4) :: nchildren
     !> Global id of coarsened elements
     ! of size coarsen_nr
     ! 1 - new global id; 2 - old child global id
     integer(i8), allocatable, dimension(:, :, :) :: coarsen_gidx
     !> Mapping of coarsened elements
     ! of size coarsen_nr
     ! 1 - old child local id; 2 - old child MPI rank
     integer(i4), allocatable, dimension(:, :, :) :: coarsen
     ! Communication part
     !> Communication flag for vector refinement/coarsening step
     logical :: ifcomm
     !> Mapping for filling output vector with same and refine elements
     ! of size nelt_neko
     ! positive - element from local vector
     ! negative - element fetched from other rank
     ! 0 - coarsening; do not fill
     integer(i4), allocatable, dimension(:) :: same_ref_fill_map
     !> Mapping for filling coarsening vector
     ! of size coarsen_nr
     ! positive - element from local vector
     ! negative - element fetched from other rank
     integer(i4), allocatable, dimension(:, :) :: crs_fill_map
     ! Receive/send part
     !> Number of MPI ranks; receive
     integer(i4) :: nrank_rcv
     !> MPI rank list; receive
     integer(i4), allocatable, dimension(:) :: rank_lst_rcv
     !> Ordered neighbour numbering and permutation to simplify communication
     integer(i4), allocatable, dimension(:) :: ngh_rcv, ind_rcv
     !> Element list (local element number); receive
     integer(i4), allocatable, dimension(:) :: elem_lst_rcv
     !> Offset in element list; receive
     integer(i4), allocatable, dimension(:) :: off_rcv
     !> Number of MPI ranks; send
     integer(i4) :: nrank_snd
     !> MPI rank list; send
     integer(i4), allocatable, dimension(:) :: rank_lst_snd
     !> Ordered neighbour numbering and permutation to simplify communication
     integer(i4), allocatable, dimension(:) :: ngh_snd, ind_snd
     !> Element list (local element number); send
     integer(i4), allocatable, dimension(:) :: elem_lst_snd
     !> Offset in element list; send
     integer(i4), allocatable, dimension(:) :: off_snd
     !> Combined neighbour list length
     integer(i4) :: nngh
     !> Combined neighbour list
     integer(i4), allocatable, dimension(:) :: ngh
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
     !> Free element distribution for field reconstruction
     procedure, pass(this) :: reconstruct_data_free => &
          p4est_reconstruct_data_free
     !> Get refinement/coarsening vector sizes and mappings
     procedure, pass(this) :: vector_map => p4est_vector_map
     !> Free refinement/coarsening vector mappings
     procedure, pass(this) :: vector_map_free => p4est_vector_map_free
     !> Construct vectors for refinement/coarsening
     procedure, pass(this) :: vector_constr => p4est_vector_constr
  end type mesh_manager_transfer_p4est_t

contains

  !> Destructor.
  subroutine p4est_free(this)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this

    call this%free_base()
    call this%vector_map_free()
    call this%reconstruct_data_free()

    this%nelt_mm = 0
    this%nelt_neko = 0
    this%nelt_old = 0

    if (allocated(this%gidx_mm)) deallocate(this%gidx_mm)
    if (allocated(this%gidx_neko)) deallocate(this%gidx_neko)
    if (allocated(this%fwd_cmm)) deallocate(this%fwd_cmm)
    if (allocated(this%bwd_cmm)) deallocate(this%bwd_cmm)
    if (allocated(this%gidx_old)) deallocate(this%gidx_old)

  end subroutine p4est_free

  !> Get element distribution
  !! @param[in]   mesh     mesh manager mesh
  subroutine p4est_elem_dist_construct(this, mesh)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    class(manager_mesh_t), intent(in) :: mesh
    integer :: il

    ! save old element distribution
    this%nelt_old = this%nelt_neko
    if (allocated(this%gidx_neko)) call move_alloc(this%gidx_neko, &
         this%gidx_old)

    ! clean old element distribution
    this%nelt_mm = 0
    this%nelt_neko = 0
    if (allocated(this%gidx_mm)) deallocate(this%gidx_mm)
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

          call neko_error('Nothing done yet; partitioning')
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
       ! ref_mark

       call neko_error('Nothing done yet; mark_transfer')
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
    integer(i8), allocatable, dimension(:), intent(in) :: same_gidx
    integer(i8), allocatable, dimension(:, :), intent(in) :: refine_gidx
    integer(i8), allocatable, dimension(:, :, :), intent(in) :: coarsen_gidx
    integer(i4), allocatable, dimension(:, :), intent(in) :: same, refine
    integer(i4), allocatable, dimension(:, :, :), intent(in) :: coarsen
    integer :: il, iref, icrs

    call this%reconstruct_data_free()

    this%same_nr = same_nr
    this%refine_nr = refine_nr
    this%coarsen_nr = coarsen_nr
    this%nchildren = size(coarsen, 2)

    ! same has never zero size
    allocate(this%same_gidx, source = same_gidx)
    allocate(this%same, source = same)
    if (this%refine_nr .gt. 0) then
       allocate(this%refine_gidx, source = refine_gidx)
       allocate(this%refine, source = refine)
    end if
    if (this%coarsen_nr .gt. 0) then
       allocate(this%coarsen_gidx, source = coarsen_gidx)
       allocate(this%coarsen, source = coarsen)
    end if

    ! check mapping consistency
    iref = 0
    icrs = 0
    do il = 1, this%nelt_neko
       if (this%same_gidx(il) .eq. 0) then
          if (this%same(1, il) .eq. 0) then
             icrs = icrs + 1
          else
             iref = iref + 1
          end if
       end if
    end do
    if (icrs .ne. this%coarsen_nr .or. iref .ne. this%refine_nr) &
         call neko_error('Inconsistent number of refined/coarsened elements')

  end subroutine p4est_reconstruct_data_set

  !> Free element distribution data for field reconstruction
  subroutine p4est_reconstruct_data_free(this)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this

    this%same_nr = 0
    this%refine_nr = 0
    this%coarsen_nr = 0
    this%nchildren = 0

    if (allocated(this%same_gidx)) deallocate(this%same_gidx)
    if (allocated(this%same)) deallocate(this%same)
    if (allocated(this%refine_gidx)) deallocate(this%refine_gidx)
    if (allocated(this%refine)) deallocate(this%refine)
    if (allocated(this%coarsen_gidx)) deallocate(this%coarsen_gidx)
    if (allocated(this%coarsen)) deallocate(this%coarsen)

  end subroutine p4est_reconstruct_data_free

  !> Get refinement/coarsening vectors sizes and mappings
  !! @param[out]    nold       old element number
  !! @param[out]    nnew       new element number
  !! @param[out]    nref       refinement mapping size
  !! @param[out]    ncrs       coarsening mapping size
  !! @param[inout]  rmap       refinement mapping (elem. and child position)
  !! @param[inout]  cmap       coarsening mapping (elem. and child position)
  !! @param[out]    nchildren  mesh manager children number
  !! @param[out]    ifchange   mesh change flag
  !! @param[out]    nrcv       receive buffer size
  !! @param[out]    nsnd       send buffer size
  subroutine p4est_vector_map(this, nold, nnew, nref, ncrs, rmap, cmap, &
       nchildren, ifchange, nrcv, nsnd)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    integer, intent(out) :: nold, nnew, nref, ncrs, nchildren, nrcv, nsnd
    integer, dimension(:, :), allocatable, intent(inout) :: rmap, cmap
    logical, intent(out) :: ifchange
    integer :: il, jl, kl, ml, ll, itmp1, itmp2, ierr, iref, icrs, cmmsame, &
         cmmref, cmmcrs
    integer, dimension(:, :), allocatable :: cmapl, cmmapl
    integer, dimension(:), allocatable :: bnd

    ! reset communication
    call this%vector_map_free()
    ifchange = .false.

    ! keep information about relative child position in the parent
    if (this%coarsen_nr .gt. 0) then
       allocate(cmapl(this%nchildren, this%coarsen_nr))
       do il = 1, this%coarsen_nr
          do jl = 1, this%nchildren
             cmapl(jl, il) = jl - 1
          end do
       end do
    end if

    ! partitioned mesh
    if (this%ifpartition) then
       ! forward communication
       ! same_gidx, same, refine_gidx, refine, coarsen_gidx, coarsen, cmapl

       
       call neko_error('Nothing done yet; vector_map')
       
    end if

    ! has the local mesh changed
    ! refinement/coarsening
    if (this%refine_nr .gt. 0 .or. this%coarsen_nr .gt. 0) ifchange = .true.
    ! same element shift
    if (.not. ifchange) then
       do il = 1, this%nelt_neko
          if (this%same(1, il) .ne. il .or. this%same(2, il) .ne. pe_rank) then
             ifchange = .true.
             exit
          end if
       end do
    end if

    ! Is communication with other ranks needed?
    do il = 1, this%nelt_neko
       if (this%same_gidx(il) .ne. 0 .and. &
            this%same(2, il) .ne. pe_rank) then
          this%ifcomm = .true.
          exit
       end if
    end do
    if (.not. this%ifcomm .and. this%refine_nr .gt. 0) then
       do il = 1, this%refine_nr
          if (this%refine(2, il) .ne. pe_rank) then
             this%ifcomm = .true.
             exit
          end if
       end do
    end if
    if (.not. this%ifcomm .and. this%coarsen_nr .gt. 0) then
       element : do il = 1, this%coarsen_nr
          do jl = 1, this%nchildren
             if (this%coarsen(2, jl, il) .ne. pe_rank) then
                this%ifcomm = .true.
                exit element
             end if
          end do
       end do element
    end if
    call MPI_Allreduce(MPI_IN_PLACE, this%ifcomm, 1, MPI_LOGICAL, MPI_LOR, &
         NEKO_COMM, ierr)

    if (allocated(rmap)) deallocate(rmap)
    if (allocated(cmap)) deallocate(cmap)

    nold = this%nelt_old
    nnew = this%nelt_neko
    nref = this%refine_nr
    ncrs = this%coarsen_nr
    nchildren = this%nchildren

    ! exchange information between processors
    if (this%ifcomm) call p4est_vector_map_comm(this, this%nchildren, cmmsame, &
         cmmref, cmmcrs, cmmapl, bnd)

    ! fill in map arrays
    if (this%refine_nr .gt. 0) then
       allocate(rmap(2, this%refine_nr))
       rmap(:, :) = 0
       iref = 0
       do il = 1, this%nelt_neko
          if (this%same_gidx(il) .eq. 0 .and. this%same(1, il) .ne. 0) then
             iref = iref + 1
             rmap(1, il) = il
             rmap(2, il) = this%refine(3, this%same(1, il))
          end if
       end do
    end if
    if (this%coarsen_nr .gt. 0) then
       allocate(cmap(1 + nchildren, this%coarsen_nr))
       cmap(:, :) = 0
       icrs = 0
       do il = 1, this%nelt_neko
          if (this%same_gidx(il) .eq. 0 .and. this%same(1, il) .eq. 0) then
             icrs = icrs + 1
             cmap(1, icrs) = il
          end if
       end do
    end if

    ! get element mapping
    allocate(this%same_ref_fill_map(this%nelt_neko))
    this%same_ref_fill_map(:) = 0
    if (this%coarsen_nr .gt. 0) then
       allocate(this%crs_fill_map(this%nchildren, this%coarsen_nr))
       this%crs_fill_map(:, :) = 0
    end if

    ! exchanged data
    if (this%ifcomm) then
       itmp1 = cmmsame + cmmref
       itmp2 = cmmsame + cmmref + cmmcrs
       if (itmp2 .gt. 0) then
          ! count exchanged elements
          kl = 0
          do il = 1, itmp2
             if (bnd(il) .ne. 0) kl = kl + 1
             if (cmmapl(3 ,il) .le. itmp1) then
                ! same and refine elements
                this%same_ref_fill_map(cmmapl(4, il)) = - kl
             else if (cmmapl(3 ,il) .le. itmp2) then
                ! coarsening elements
                this%crs_fill_map(cmmapl(5, il), cmmapl(4, il)) = - kl
             else
                call neko_error('Wrong storage position')
             end if
          end do
       end if
    end if

    ! local elements
    ! same and refined elements
    do il = 1, this%nelt_neko
       if (this%same_gidx(il) .ne. 0) then
          ! unchanged element
          if (this%same(2, il) .eq. pe_rank) then
             ! local element
             this%same_ref_fill_map(il) = this%same(1, il)
          else
             ! fetched data
             if (this%same_ref_fill_map(il) .ge. 0) &
                  call neko_error('Same: wrong storage position')
          end if
       else
          ! changed element
          if (this%same(1, il) .ne. 0) then
             ! refinement
             if (this%refine(2, this%same(1, il)) .eq. pe_rank) then
                ! local element
                this%same_ref_fill_map(il) = this%refine(1, this%same(1, il))
             else
                ! fetched data
                if (this%same_ref_fill_map(il) .ge. 0) &
                     call neko_error('Ref: wrong storage position')
             end if
          end if
       end if
    end do

    ! coarsening elements
    if (this%coarsen_nr .gt. 0) then
       do il = 1, this%coarsen_nr
          do jl = 1, this%nchildren
             if (this%coarsen(2, jl, il) .eq. pe_rank) then
                ! local element
                this%crs_fill_map(jl, il) = this%coarsen(1, jl, il)
                cmap(1 + jl, il) = cmapl(jl, il)
             else
                ! fetched data
                if (this%crs_fill_map(jl, il) .ge. 0) &
                     call neko_error('CRS: wrong storage position')
             end if
          end do
       end do
    end if

    ! receive/send buffers size and sanity check
    nrcv = 0
    nsnd = 0
    if (this%nrank_rcv .gt. 0) then
       nrcv = this%off_rcv(this%nrank_rcv + 1) - 1
       if (nrcv .eq. 0) &
            call neko_error('Inconsistent receive rank and offset arrays')
    end if
    if (this%nrank_snd .gt. 0) then
       nsnd = this%off_snd(this%nrank_snd + 1) - 1
       if (nsnd .eq. 0) &
            call neko_error('Inconsistent send rank and offset arrays')
    end if

    if (allocated(cmapl)) deallocate(cmapl)
    if (allocated(cmmapl)) deallocate(cmmapl)
    if (allocated(bnd)) deallocate(bnd)

  end subroutine p4est_vector_map

  !> Exchange data between MPI ranks
  !! @param[in]     nchildren mesh manager children number
  !! @param[out]    cmmsame   number of same elements for communication
  !! @param[out]    cmmref    number of refined elements for communication
  !! @param[out]    cmmcrs    number of coarsened elements for communication
  !! @param[inout]  cmmapl    element info
  !! @param[inout]  ind       rank/element section boundary
  subroutine p4est_vector_map_comm(this, nchildren, cmmsame, cmmref, cmmcrs, &
       cmmapl, ind)
    type(mesh_manager_transfer_p4est_t), intent(inout) :: this
    integer, intent(in) :: nchildren
    integer, intent(out) :: cmmsame, cmmref, cmmcrs
    integer, dimension(:, :), allocatable, intent(inout) :: cmmapl
    integer, dimension(:), allocatable, intent(inout) :: ind
    integer :: il, jl, kl, ll, ml, itmp, bmax, src, dst, ierr, n_recv
    integer(i8), dimension(:), allocatable :: cmmgidxl, cmmgidxls
    integer(i8), dimension(:, :), allocatable :: rbuf, sbuf
    integer, dimension(:), allocatable :: ngh_dst, ind_dst, vtmp
    integer, parameter :: lda = 5 ! tuple length
    integer, dimension(lda) :: aa ! tmp array fro sorting
    integer, parameter :: nkey = 2 ! number of keys
    integer, dimension(nkey) :: key
    type(MPI_Status) :: status
    logical :: ifconsistent

    if (allocated(cmmapl)) deallocate(cmmapl)
    if (allocated(ind)) deallocate(ind)
    cmmsame = 0
    cmmref = 0
    cmmcrs = 0

    ! build data exchange information
    ! count same elements requiring communication
    do il = 1, this%nelt_neko
       if (this%same_gidx(il) .ne. 0 .and. &
            this%same(2, il) .ne. pe_rank) then
          cmmsame = cmmsame + 1
       end if
    end do
    ! count refined elements requiring communication
    if (this%refine_nr .gt. 0) then
       do il = 1, this%refine_nr
          if (this%refine(2, il) .ne. pe_rank) then
             cmmref = cmmref + 1
          end if
       end do
    end if
    ! count coarsened elements requiring communication
    if (this%coarsen_nr .gt. 0) then
       do il = 1, this%coarsen_nr
          do jl = 1, nchildren
             if (this%coarsen(2, jl, il) .ne. pe_rank) then
                cmmcrs = cmmcrs + 1
             end if
          end do
       end do
    end if

    ! get element list for receive
    itmp = cmmsame + cmmref + cmmcrs
    if (itmp .gt. 0) then
       allocate (cmmapl(lda, itmp), cmmgidxl(itmp), ind(itmp))
       ! same
       if (cmmsame .gt. 0) then
          cmmsame = 0
          do il = 1, this%nelt_neko
             if (this%same_gidx(il) .ne. 0 .and. &
                  this%same(2, il) .ne. pe_rank) then
                cmmsame = cmmsame + 1
                cmmapl(1, cmmsame) = this%same(2, il) ! MPI rank
                cmmapl(2, cmmsame) = this%same(1, il) ! old loc. id
                cmmapl(3, cmmsame) = cmmsame ! storage position
                cmmapl(4, cmmsame) = il ! position in rmap array
                cmmapl(5, cmmsame) = 0
                cmmgidxl(cmmsame) = this%same_gidx(il) ! old gl. id
             end if
          end do
       end if
       ! refine
       if (cmmref .gt. 0) then
          cmmref = 0
          do il = 1, this%refine_nr
             if (this%refine(2, il) .ne. pe_rank) then
                cmmref = cmmref + 1
                itmp = cmmsame + cmmref
                cmmapl(1, itmp) = this%refine(2, il) ! MPI rank
                cmmapl(2, itmp) = this%refine(1, il) ! old loc. id
                cmmapl(3, itmp) = itmp ! storage position
                cmmapl(4, itmp) = il ! position in rmap array
                cmmapl(5, itmp) = 0
                cmmgidxl(itmp) = this%refine_gidx(2, il) ! old gl. id
             end if
          end do
       end if
       ! coarsen
       if (cmmcrs .gt. 0) then
          cmmcrs = 0
          do il = 1, this%coarsen_nr
             do jl = 1, nchildren
                if (this%coarsen(2, jl, il) .ne. pe_rank) then
                   cmmcrs = cmmcrs + 1
                   itmp = cmmsame + cmmref + cmmcrs
                   cmmapl(1, itmp) = this%coarsen(2, jl, il) ! MPI rank
                   cmmapl(2, itmp) = this%coarsen(1, jl, il) ! old loc. id
                   cmmapl(3, itmp) = itmp ! storage position
                   cmmapl(4, itmp) = il ! position in cmap array
                   cmmapl(5, itmp) = jl ! child position
                   cmmgidxl(itmp) = this%coarsen_gidx(2, jl, il) !old gl. id
                end if
             end do
          end do
       end if

       ! sort elements with respect of MPI rank and local element number
       itmp = cmmsame + cmmref + cmmcrs
       key(1) = 1
       key(2) = 2
       call sort_tuple(cmmapl, lda, itmp, key, nkey, ind, aa)
       ! sort old element global id
       call swap(cmmgidxl, ind, itmp)

       ! count various receive ranks and get receive mappings
       ind(:) = 0
       ind(1) = 1
       jl = cmmapl(1, 1) ! MPI rank
       kl = cmmapl(2, 1) ! local element number
       ll = 1 ! MPI rank count
       ml = 1 ! element count
       do il = 2, itmp
          if (jl .ne. cmmapl(1, il)) then
             ind(il) = 1 ! rank border
             jl = cmmapl(1, il)
             kl = cmmapl(2, il)
             ll = ll + 1
             ml = ml + 1
          else
             if (kl .ne. cmmapl(2, il)) then
                ind(il) = 2 ! element border
                kl = cmmapl(2, il)
                ml = ml + 1
             end if
          end if
       end do
       this%nrank_rcv = ll
       allocate(this%rank_lst_rcv(ll), this%off_rcv(ll + 1), &
            this%elem_lst_rcv(ml), cmmgidxls(ml))
       this%rank_lst_rcv(1) = cmmapl(1, 1)
       this%off_rcv(1) = 1
       this%elem_lst_rcv(1) = cmmapl(2, 1)
       cmmgidxls(1) = cmmgidxl(1)
       kl = 1 ! rank count
       ll = 1 ! element count
       do il = 2, itmp
          if (ind(il) .eq. 1) then ! rank boundary
             kl = kl + 1
             ll = ll + 1
             this%rank_lst_rcv(kl) = cmmapl(1, il)
             this%off_rcv(kl) = ll
             this%elem_lst_rcv(ll) = cmmapl(2, il)
             cmmgidxls(ll) = cmmgidxl(il)
          else if (ind(il) .eq. 2) then ! element boundary
             ll = ll + 1
             this%elem_lst_rcv(ll) = cmmapl(2, il)
             cmmgidxls(ll) = cmmgidxl(il)
          end if
       end do
       kl = kl + 1
       this%off_rcv(kl) = ll + 1
    else
       this%nrank_rcv = 0
    end if ! (itmp .gt. 0)

    ! exchange information with other ranks
    ! get max value of the buffer
    bmax = 0
    if (this%nrank_rcv .gt. 0) then
       do il = 1, this%nrank_rcv
          bmax = max(bmax, this%off_rcv(il + 1) - this%off_rcv(il))
       end do
    end if
    call MPI_Allreduce(MPI_IN_PLACE, bmax, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    if (bmax .eq. 0) call neko_error('p4est_vector_map_comm: inconsistent &
         &communication flag')
    allocate(rbuf(2, bmax), sbuf(2, bmax))
    ! take int account data  amount is doubled
    bmax = bmax * 2

    ! get list of destination neighbours
    if (this%nrank_rcv .gt. 0) then
       allocate(ngh_dst(this%nrank_rcv), ind_dst(this%nrank_rcv))
       do il = 1, this%nrank_rcv
          ngh_dst(il) = mod(this%rank_lst_rcv(il) - pe_rank + pe_size, pe_size)
       end do
       ! sort destination neighbours
       call sort(ngh_dst, ind_dst, this%nrank_rcv)
    end if

    ! count send requests
    this%nrank_snd = 0
    ! count destinations
    jl = 1
    do il = 1, pe_size - 1
       src = modulo(pe_rank - il + pe_size, pe_size)
       dst = modulo(pe_rank + il, pe_size)
       ! destination depends on receiving requests, take advantage of the
       ! fact destinations are ordered
       ll = 0
       if (this%nrank_rcv .gt. 0) then
          if (ngh_dst(jl) .eq. il) then
             ll = this%off_rcv(ind_dst(jl) + 1) - this%off_rcv(ind_dst(jl))
             do kl = 1, ll
                sbuf(1, kl) = & ! old element local index
                     this%elem_lst_rcv(this%off_rcv(ind_dst(jl)) + kl - 1)
                sbuf(2, kl) = & ! old element global index
                     cmmgidxls(this%off_rcv(ind_dst(jl)) + kl - 1)
             end do
             if (jl .lt. this%nrank_rcv) jl = jl + 1
          end if
       end if
       ll = ll * 2
       call MPI_Sendrecv(sbuf, ll, MPI_INTEGER8, dst, 0, &
            rbuf, bmax, MPI_INTEGER8, src, 0, NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)

       if (n_recv .gt. 0) then
          n_recv = n_recv / 2
          ! check if the received data is consistent with local one by
          ! comparing global element numbers
          ifconsistent = .true.
          do kl = 1, n_recv
             if (rbuf(2, kl) .ne. this%gidx_old(int(rbuf(1, kl), i4))) then
                ifconsistent = .false.
                exit
             end if
          end do
          if (.not. ifconsistent) &
               call neko_error('Inconsistent global element number; &
               &vector_map')
          ! update send processor data
          if (this%nrank_snd .eq. 0) then
             this%nrank_snd = 1
             allocate(this%rank_lst_snd(1))
             this%rank_lst_snd(1) = src
             allocate(this%off_snd(2))
             this%off_snd(1) = 1
             this%off_snd(2) = n_recv + 1
             allocate(this%elem_lst_snd(n_recv))
             do kl = 1, n_recv
                this%elem_lst_snd(kl) = int(rbuf(1, kl), i4)
             end do
          else
             this%nrank_snd = this%nrank_snd + 1
             allocate(vtmp(this%nrank_snd))
             do kl = 1, this%nrank_snd - 1
                vtmp(kl) = this%rank_lst_snd(kl)
             end do
             vtmp(this%nrank_snd) = src
             call move_alloc(vtmp, this%rank_lst_snd)
             allocate(vtmp(this%nrank_snd + 1))
             do kl = 1, this%nrank_snd
                vtmp(kl) = this%off_snd(kl)
             end do
             vtmp(this%nrank_snd + 1) = vtmp(this%nrank_snd) + n_recv
             call move_alloc(vtmp, this%off_snd)
             allocate(vtmp(this%off_snd(this%nrank_snd + 1) -1))
             do kl = 1, this%off_snd(this%nrank_snd) - 1
                vtmp(kl) = this%elem_lst_snd(kl)
             end do
             do kl = 1, n_recv
                vtmp(this%off_snd(this%nrank_snd) - 1 + kl) = &
                     int(rbuf(1, kl), i4)
             end do
             call move_alloc(vtmp, this%elem_lst_snd)
          end if
       end if

    end do ! processor number

    ! get ordered neighbour numbering to simplify communication call
    if (this%nrank_rcv .gt. 0) then
       allocate(this%ngh_rcv(this%nrank_rcv), this%ind_rcv(this%nrank_rcv))
       do il = 1, this%nrank_rcv
          this%ngh_rcv(il) = mod(pe_rank - this%rank_lst_rcv(il) + pe_size, &
               pe_size)
       end do
       ! sort receive neighbours
       call sort(this%ngh_rcv, this%ind_rcv, this%nrank_rcv)
    end if

    if (this%nrank_snd .gt. 0) then
       allocate(this%ngh_snd(this%nrank_snd), this%ind_snd(this%nrank_snd))
       do il = 1, this%nrank_snd
          this%ngh_snd(il) = mod(this%rank_lst_snd(il) - pe_rank + pe_size, &
               pe_size)
       end do
       ! sort send neighbours
       call sort(this%ngh_snd, this%ind_snd, this%nrank_snd)
    end if

    ! combine neighbour lists
    this%nngh = 0
    if (this%nrank_rcv .gt. 0 .and. this%nrank_snd .gt. 0) then
       ! take advantage of the fact both lists are sorted and positive
       if (allocated(vtmp)) deallocate(vtmp)
       allocate(vtmp(this%nrank_rcv + this%nrank_snd))
       il = 1
       jl = 1
       kl = 1
       vtmp(kl) = min(this%ngh_rcv(il), this%ngh_snd(jl))
       do
          if (this%ngh_rcv(il) .eq. this%ngh_snd(jl)) then
             if (vtmp(kl) .lt. this%ngh_rcv(il)) then
                kl = kl + 1
                vtmp(kl) = this%ngh_rcv(il)
             end if
             il = il + 1
          else if (this%ngh_rcv(il) .gt. this%ngh_snd(jl)) then
             if (vtmp(kl) .lt. this%ngh_snd(jl)) then
                kl = kl + 1
                vtmp(kl) = this%ngh_snd(jl)
             end if
             jl = jl + 1
          else
             if (vtmp(kl) .lt. this%ngh_rcv(il)) then
                kl = kl + 1
                vtmp(kl) = this%ngh_rcv(il)
             end if
             il = il + 1
          end if
          if (il .gt. this%nrank_rcv .or. jl .gt. this%nrank_snd) exit
       end do
       if (il .gt. this%nrank_rcv) then
          do il = jl, this%nrank_snd
             if (vtmp(kl) .lt. this%ngh_snd(il)) then
                kl = kl + 1
                vtmp(kl) = this%ngh_snd(il)
             end if
          end do
       else if (jl .gt. this%nrank_snd) then
          do jl = il, this%nrank_rcv
             if (vtmp(kl) .lt. this%ngh_rcv(jl)) then
                kl = kl + 1
                vtmp(kl) = this%ngh_rcv(jl)
             end if
          end do
       end if
       this%nngh = kl
       allocate(this%ngh(kl))
       do il = 1, kl
          this%ngh(il) = vtmp(il)
       end do
    else if (this%nrank_rcv .gt. 0) then
       this%nngh = this%nrank_rcv
       allocate(this%ngh(this%nrank_rcv))
       this%ngh(:) = this%ngh_rcv(:)
    else if (this%nrank_snd .gt. 0) then
       this%nngh = this%nrank_snd
       allocate(this%ngh(this%nrank_snd))
       this%ngh(:) = this%ngh_snd(:)
    end if

    if (allocated(cmmgidxl)) deallocate(cmmgidxl)
    if (allocated(cmmgidxls)) deallocate(cmmgidxls)
    if (allocated(rbuf)) deallocate(rbuf)
    if (allocated(sbuf)) deallocate(sbuf)
    if (allocated(ngh_dst)) deallocate(ngh_dst)
    if (allocated(ind_dst)) deallocate(ind_dst)
    if (allocated(vtmp)) deallocate(vtmp)

  end subroutine p4est_vector_map_comm

  !> Free refinement/coarsening vector mappings
  subroutine p4est_vector_map_free(this)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this

    this%ifcomm = .false.
    this%nrank_rcv = 0
    this%nrank_snd = 0
    this%nngh = 0

    if (allocated(this%same_ref_fill_map)) deallocate(this%same_ref_fill_map)
    if (allocated(this%crs_fill_map)) deallocate(this%crs_fill_map)
    if (allocated(this%rank_lst_rcv)) deallocate(this%rank_lst_rcv)
    if (allocated(this%ngh_rcv)) deallocate(this%ngh_rcv)
    if (allocated(this%ind_rcv)) deallocate(this%ind_rcv)
    if (allocated(this%elem_lst_rcv)) deallocate(this%elem_lst_rcv)
    if (allocated(this%off_rcv)) deallocate(this%off_rcv)
    if (allocated(this%rank_lst_snd)) deallocate(this%rank_lst_snd)
    if (allocated(this%ngh_snd)) deallocate(this%ngh_snd)
    if (allocated(this%ind_snd)) deallocate(this%ind_snd)
    if (allocated(this%elem_lst_snd)) deallocate(this%elem_lst_snd)
    if (allocated(this%off_snd)) deallocate(this%off_snd)
    if (allocated(this%ngh)) deallocate(this%ngh)

  end subroutine p4est_vector_map_free

  !> Construct vectors for refinement/coarsening
  !! @param[in]   vin       original vector
  !! @param[out]  vout      output vector for refinement
  !! @param[out]  vcrs      vector with additional data for coarsening
  !! @param[out]  buff_rcv  receive buffer
  !! @param[out]  buff_snd  send buffer
  !! @param[in]   elsize    element size
  subroutine p4est_vector_constr(this, vin, vout, vcrs, buff_rcv, buff_snd, &
       elsize)
    class(mesh_manager_transfer_p4est_t), intent(inout) :: this
    real(rp), dimension(:, :, :, :), intent(in) :: vin
    real(rp), dimension(:, :, :, :), intent(out) :: vout, buff_rcv, buff_snd
    real(rp), dimension(:, :, :, :, :), intent(out) :: vcrs
    integer, intent(in) :: elsize
    integer :: il, jl, kl, src, dst, start, dim, ierr
    logical :: ifwait_rcv, ifwait_snd
    type(MPI_Request) :: snd_req, rcv_req

    ! sanity check
    if (size(vin, 4) .ne. this%nelt_old .or. &
         size(vout, 4) .ne. this%nelt_neko) &
         call neko_error('Inconsistent input/output vector size')
    if (this%coarsen_nr .gt. 0) then
       if (size(vcrs, 5) .ne. this%coarsen_nr) &
            call neko_error('Inconsistent coarsening vector size')
    end if
    if (this%nrank_rcv .gt. 0) then
       if (size(buff_rcv, 4) .ne. this%off_rcv(this%nrank_rcv + 1) - 1) &
            call neko_error('Inconsistent receive buffer size')
    end if
    if (this%nrank_snd .gt. 0) then
       if (size(buff_snd, 4) .ne. this%off_snd(this%nrank_snd + 1) - 1) &
            call neko_error('Inconsistent send buffer size')
    end if

    if (this%ifcomm) then
       ! fill send buffer
       if (this%nrank_snd .gt. 0) then
          do il = 1, this%off_snd(this%nrank_snd + 1) - 1
             buff_snd(:, :, :, il) = vin(:, :, :, this%elem_lst_snd(il))
          end do
       end if

       ! exchange data with other MPI ranks
       ! count receives and sends
       jl = 1
       kl = 1
       do il = 1, this%nngh
          ! receive
          if (this%nrank_rcv .gt. 0) then
             if (this%ngh(il) .eq. this%ngh_rcv(jl)) then
                src = this%rank_lst_rcv(this%ind_rcv(jl))
                start = this%off_rcv(this%ind_rcv(jl))
                dim = elsize * (this%off_rcv(this%ind_rcv(jl) + 1) - start)
                call MPI_Irecv(buff_rcv(: , :, :, start), dim, &
                     MPI_REAL_PRECISION, src, 0, NEKO_COMM, rcv_req, ierr)
                ifwait_rcv = .true.
             else
                ifwait_rcv = .false.
             end if
          else
             ifwait_rcv = .false.
          end if

          ! send
          if (this%nrank_snd .gt. 0) then
             if (this%ngh(il) .eq. this%ngh_snd(kl)) then
                dst = this%rank_lst_snd(this%ind_snd(kl))
                start = this%off_snd(this%ind_snd(kl))
                dim = elsize * (this%off_snd(this%ind_snd(kl) + 1) - start)
                call MPI_Isend(buff_snd(: , :, :, start), dim, &
                     MPI_REAL_PRECISION, dst, 0, NEKO_COMM, snd_req, ierr)
                ifwait_snd = .true.
             else
                ifwait_snd = .false.
             end if
          else
             ifwait_snd = .false.
          end if

          ! receive
          if (ifwait_rcv) then
             call MPI_Wait(rcv_req, MPI_STATUS_IGNORE, ierr)
             if (jl .lt. this%nrank_rcv) jl = jl + 1
          end if

          ! send
          if (ifwait_snd) then
             call MPI_Wait(snd_req, MPI_STATUS_IGNORE, ierr)
             if (kl .lt. this%nrank_snd) kl = kl + 1
          end if
       end do
    end if

    ! Fill same and refined elements
    do il = 1, this%nelt_neko
       if (this%same_ref_fill_map(il) .gt. 0) then
          ! local array
          vout(:, :, :, il) = vin(:, :, :, this%same_ref_fill_map(il))
       else if (this%same_ref_fill_map(il) .lt. 0) then
          ! fetched data
          vout(:, :, :, il) = &
               buff_rcv(:, :, :, abs(this%same_ref_fill_map(il)))
       end if
    end do

    ! Fill coarsening data
    if (this%coarsen_nr .gt. 0) then
       do il = 1, this%coarsen_nr
          do jl = 1, this%nchildren
             if (this%crs_fill_map(jl, il) .gt. 0) then
                ! local array
                vcrs(:, :, :, jl, il) = vin(:, :, :, this%crs_fill_map(jl, il))
             else if (this%crs_fill_map(jl, il) .lt. 0) then
                ! fetched data
                vcrs(:, :, :, jl, il) = &
                     buff_rcv(:, :, :, abs(this%crs_fill_map(jl, il)))
             end if
          end do
       end do
    end if

  end subroutine p4est_vector_constr

end module mesh_manager_transfer_p4est
