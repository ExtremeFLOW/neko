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
!> Implementation of the mesh connectivity type for p4est mesh manager
module manager_conn_p4est
  use num_types, only : i4, i8, rp, dp
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use utils, only : neko_error
  use math, only : sort, sort_tuple
  use manager_conn, only : manager_conn_obj_t, manager_conn_t
  use mpi_f08, only : MPI_Allreduce, MPI_Irecv, MPI_Isend, MPI_Wait, &
       MPI_Get_count, MPI_IN_PLACE, MPI_STATUS_IGNORE, MPI_INTEGER, &
       MPI_INTEGER8, MPI_MAX, MPI_Status, MPI_Request, MPI_Barrier
  implicit none
  private

  !> Extended type for connectivity information regarding vertices, faces and
  !! edges.
  !! @details It adds sharing information for communication
  type, extends(manager_conn_obj_t), public :: manager_conn_obj_p4est_t
     !> Number of MPI ranks sharing objects
     integer(i4) :: nrank
     !> Length of the share array
     integer(i4) :: nshare
     !> List of ranks sharing objects
     integer(i4), allocatable, dimension(:) :: rank
     !> List of shared objects
     integer(i4), allocatable, dimension(:) :: share
     !> Mapping of the shared objects local number to the remote local number
     ! the values related to the local rank are set to -1
     integer(i4), allocatable, dimension(:) :: sharemap
     !> Offset in the share list
     integer(i4), allocatable, dimension(:) :: off
     !> Local mapping of object to element component
     ! inverse to map in manager_conn_t; related to the local numbering
     ! first dimension: 1- local element number, 2 - position in element
     integer(i4), allocatable, dimension(:,:) :: lmap
     !> Offset in the local mapping
     integer(i4), allocatable, dimension(:) :: lmapoff
     !> List of shared objects (is included in share under pe_rank as well)
     integer(i4), allocatable, dimension(:) :: glist
     !> Global mapping of object to element component
     ! This array does not include local mapping
     ! first dimension: 1 - MPI rank, 2- local element number,
     ! 3 - position in element
     integer(i4), allocatable, dimension(:,:) :: gmap
     !> Offset in the global mapping
     integer(i4), allocatable, dimension(:) :: gmapoff
   contains
     procedure, pass(this) :: init_data => manager_conn_obj_init_data_p4est
     procedure, pass(this) :: init_type => manager_conn_obj_init_type_p4est
     procedure, pass(this) :: free => manager_conn_obj_free_p4est
  end type manager_conn_obj_p4est_t

  !> Type for element connectivity information
  type, extends(manager_conn_t), public :: manager_conn_p4est_t
     ! Face orientation.
     ! In 2D case permutation array takes 2 values
     !          0 => no permutations
     !          1 => row permutations
     ! In 3D case there are 8 possible transformations for faces
     ! Symmetric vertex numbering  (0,1,2,3)
     !                       2--3
     !                       |  |
     !                       0--1
     ! They are numbered according to p4est by 0..7 and stored in falg
     ! { 0, 1, 2, 3 }      0 => identity
     ! { 0, 2, 1, 3 }      1 => T
     ! { 1, 0, 3, 2 }      2 => P_x
     ! { 1, 3, 0, 2 }      3 => P_x T
     ! { 2, 0, 3, 1 }      4 => P_y T
     ! { 2, 3, 0, 1 }      5 => P_y
     ! { 3, 1, 2, 0 }      6 => P_y P_x T
     ! { 3, 2, 1, 0 }      7 => P_y P_x
     ! where   T - transpose
     !               P_x - permutation in x (rows)
     !               P_y - permutation in y (columns)
     !> Face alignment
     integer(i4), allocatable, dimension(:,:) :: falgn

     ! Edge orientation is similar to 2D face
     !> Edge alignment
     integer(i4), allocatable, dimension(:,:) :: ealgn

     ! Topology information: hanging elements, faces, edges and vertices
     ! The element is marked as hanging if it contains at least one hanging
     ! object.
     ! Hanging face
     ! 2D
     !   = -1 if the face is not hanging,
     !   = 0 if the face is the first half,
     !   = 1 if the face is the second half.
     ! 3D
     !   = -1 if the face is not hanging,
     !   = the corner of the full face that it touches:
     ! Hanging edge
     !   = -1 if the edge is not hanging,
     !   =  0 if the edge is the first half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  1 if the edge is the second half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  2 if the edge is the first half of a full edge
     !        and is on the boundary of a full face,
     !   =  3 if the edge is the second half of a full edge
     !        and is on the boundary of a full face,
     !   =  4 if the edge is in the middle of a full face.
     !        See the diagram below for clarification.
     !  o...............o o...............o +---2---+.......o o.......+---3---+
     !  :               : :               : |       |       : :       |       |
     !  :               : :               : 3   2   4       : :       4   3   3
     !  :               : :               : |       |       : :       |       |
     !  +---4---+       : :       +---4---+ +---4---+       : :       +---4---+
     !  |       |       : :       |       | :               : :               :
     !  2   0   4       : :       4   1   2 :               : :               :
     !  |       |       : :       |       | :               : :               :
     !  +---2---+.......o o.......+---3---+ o...............o o...............o
     !
     !                     o                  +-------+
     !                     :                  |\       \
     !                     :                  1 \       \
     !                     :                  |  +-------+
     !                     +-------+          +  |       |
     !                     |\       \         :\ |       |
     !                     0 \       \        : \|       |
     !                     |  +-------+       :  +-------+
     !                     +  |       |       o
     !                      \ |       |
     !                       \|       |
     !                        +-------+
     !> Hanging element list; 1 if at least one hanging face or edge otherwise 0
     logical, allocatable, dimension(:) :: hngel
     !> Hanging face list; position of hanging face; otherwise -1
     integer(i4), allocatable, dimension(:,:) :: hngfc
     !> Hanging edge list;
     integer(i4), allocatable, dimension(:,:) :: hnged

   contains
     procedure, pass(this) :: init => manager_conn_init_p4est
     procedure, pass(this) :: init_data => manager_conn_init_data_p4est
     procedure, pass(this) :: init_type => manager_conn_init_type_p4est
     procedure, pass(this) :: free_data => manager_conn_free_data_p4est
     procedure, pass(this) :: free => manager_conn_free_p4est
  end type manager_conn_p4est_t

contains

  !> Initialise connectivity object type
  !! @param[in]    lnum    local number of objects
  !! @param[in]    lown    number of owned objects
  !! @param[in]    goff    global object offset
  !! @param[in]    gnum    global number of objects
  !! @param[in]    nrank   number of MPI ranks sharing objects and
  !! @param[in]    nshare  number of shared objects
  !! @param[inout] gidx    global object index
  !! @param[inout] rank    list of ranks sharing object
  !! @param[inout] share   list of shared objects
  !! @param[inout] off     offset in share list
  !! @param[in]    nobj    number of objects in the element
  !! @param[in]    nelt    number of elements
  !! @param[in]    map     element to object mapping
  subroutine manager_conn_obj_init_data_p4est(this, lnum, lown, goff, gnum, &
       nrank, nshare, gidx, rank, share, off, nobj, nelt, map)
    class(manager_conn_obj_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown, nrank, nshare, nobj, nelt
    integer(i8), intent(in) :: goff, gnum
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: rank, share, off
    integer(i4), dimension(nobj, nelt), intent(in) :: map

    ! sanity check
    if ((lnum .ne. size(gidx)) .or. (nrank .ne. size(rank)) .or. &
         (nshare .ne. size(share)) .or. (nrank + 1 .ne. size(off))) &
         call neko_error('Inconsistent array sizes; p4est%conn_obj')

    call this%free()
    call this%init_data_base(lnum, lown, goff, gnum, gidx)

    this%nrank = nrank
    this%nshare = nshare

    if (allocated(rank)) call move_alloc(rank, this%rank)
    if (allocated(share)) call move_alloc(share, this%share)
    if (allocated(off)) call move_alloc(off, this%off)

    ! Create global inverse mapping
    call manager_conn_obj_global_map(this, nobj, nelt, map)

  end subroutine manager_conn_obj_init_data_p4est

  !> Build global object to element mapping
  !! @param[inout] this    object type
  !! @param[in]    nobj    number of objects in the element
  !! @param[in]    nelt    number of elements
  !! @param[in]    map     element to object mapping
  subroutine manager_conn_obj_global_map(this, nobj, nelt, map)
    type(manager_conn_obj_p4est_t), intent(inout) :: this
    integer, intent(in) :: nobj, nelt
    integer, dimension(nobj, nelt), intent(in) :: map
    integer :: il, jl, kl, ll, ml, bmax, ierr, nngh, itmp, src, dst, start, &
         dim, n_recv, ncomm
    integer(i8), allocatable, dimension(:, :) :: rbuf_gidx, sbuf_gidx, tmp_gidx
    integer, allocatable, dimension(:, :) :: ngh_src, ngh_dst, mapl, rbuf, &
         sbuf, vtmp
    integer, dimension(:), allocatable :: ind_src, ind_dst, ngh, ind, llist
    integer, parameter :: lda1 = 2, lda2 = 3 , lda3 = 4! tuple length
    integer, dimension(lda1) :: aa1 ! tmp array for sorting
    integer(i8), dimension(lda1) :: aa1_gidx ! tmp array for sorting
    integer, dimension(lda2) :: aa2 ! tmp array for sorting
    integer, dimension(lda3) :: aa3 ! tmp array for sorting
    integer, parameter :: nkey = 3 ! max. number of keys
    integer, dimension(nkey) :: key
    type(MPI_Request) :: src_req, dst_req
    type(MPI_Status) :: status
    logical :: ifwait_src, ifwait_dst

    ! Get mapping of the shared objects
    ! get the buffer
    bmax = 0
    if (this%nrank .gt. 0) then
       do il = 1, this%nrank
          if (this%rank(il) .ne. pe_rank) &
               bmax = max(bmax, this%off(il + 1) - this%off(il))
       end do
    end if
    call MPI_Allreduce(MPI_IN_PLACE, bmax, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    if (this%nrank .gt. 0) then
       allocate(rbuf_gidx(2, bmax), sbuf_gidx(2, bmax), tmp_gidx(2, bmax))
       ! get neighbour index skipping own rank
       allocate(ngh_src(lda1, this%nrank -1), ind_src(this%nrank -1), &
            ngh_dst(lda1, this%nrank -1), ind_dst(this%nrank - 1))
       jl = 0
       do il = 1, this%nrank
          if (this%rank(il) .ne. pe_rank) then
             jl = jl + 1
             ngh_src(1, jl) = mod(this%rank(il) - pe_rank + pe_size, pe_size)
             ngh_src(2, jl) = il
             ngh_dst(1, jl) = mod(pe_rank - this%rank(il) + pe_size, pe_size)
             ngh_dst(2, jl) = il
          end if
       end do
       itmp = this%nrank -1
       key(1) = 1
       call sort_tuple(ngh_src, lda1, itmp, key, 1, ind_src, aa1)
       call sort_tuple(ngh_dst, lda1, itmp, key, 1, ind_dst, aa1)
       ! combine neighbour lists
       nngh = 0
       ! take advantage of the fact both lists are sorted and positive
       allocate(ind(2 * itmp))
       il = 1
       jl = 1
       kl = 1
       ind(kl) = min(ngh_src(1, il), ngh_dst(1, jl))
       do
          if (ngh_src(1, il) .eq. ngh_dst(1, jl)) then
             if (ind(kl) .lt. ngh_src(1, il)) then
                kl = kl + 1
                ind(kl) = ngh_src(1, il)
             end if
             il = il + 1
          else if (ngh_src(1, il) .gt. ngh_dst(1, jl)) then
             if (ind(kl) .lt. ngh_dst(1, jl)) then
                kl = kl + 1
                ind(kl) = ngh_dst(1, jl)
             end if
             jl = jl + 1
          else
             if (ind(kl) .lt. ngh_src(1, il)) then
                kl = kl + 1
                ind(kl) = ngh_src(1, il)
             end if
             il = il + 1
          end if
          if (il .gt. itmp .or. jl .gt. itmp) exit
       end do
       if (il .gt. itmp) then
          do il = jl, itmp
             if (ind(kl) .lt. ngh_dst(1, il)) then
                kl = kl + 1
                ind(kl) = ngh_dst(1, il)
             end if
          end do
       else if (jl .gt. itmp) then
          do jl = il, itmp
             if (ind(kl) .lt. ngh_src(1, jl)) then
                kl = kl + 1
                ind(kl) = ngh_src(1, jl)
             end if
          end do
       end if
       nngh = kl
       allocate(ngh(kl))
       do il = 1, kl
          ngh(il) = ind(il)
       end do
       deallocate(ind)

       allocate(this%sharemap(this%nshare))
       this%sharemap(:) = -1

       ! while counting sends/receives one cannot go beyond array size
       ncomm = this%nrank -1

       ! Exchange global and local object numbers with other ranks
       allocate(ind(bmax))
       jl = 1
       kl = 1
       do il = 1, nngh
          ! receive
          if (ngh(il) .eq. ngh_src(1, jl)) then
             src = this%rank(ngh_src(2, jl))
             start = this%off(ngh_src(2, jl))
             dim = 2 * (this%off(ngh_src(2, jl) + 1) - start)
             call MPI_Irecv(rbuf_gidx, dim, MPI_INTEGER8, src, 0, NEKO_COMM, &
                  src_req, ierr)
             ifwait_src = .true.
          else
             ifwait_src = .false.
          end if

          ! send
          if (ngh(il) .eq. ngh_dst(1, kl)) then
             dst = this%rank(ngh_dst(2, kl))
             start = this%off(ngh_dst(2, kl))
             dim = this%off(ngh_dst(2, kl) + 1) - start
             start = start - 1
             do ll = 1, dim
                sbuf_gidx(1, ll) = this%share(start + ll)
                sbuf_gidx(2, ll) = this%gidx(this%share(start + ll))
             end do
             dim = 2 * dim
             call MPI_Isend(sbuf_gidx, dim, MPI_INTEGER8, dst, 0, NEKO_COMM, &
                  dst_req, ierr)
             ifwait_dst = .true.
          else
             ifwait_dst = .false.
          end if

          ! receive
          if (ifwait_src) then
             call MPI_Wait(src_req, MPI_STATUS_IGNORE, ierr)
             ! find object mapping
             ! local objects
             start = this%off(ngh_src(2, jl))
             dim = this%off(ngh_src(2, jl) + 1) - start
             start = start - 1
             do ll = 1, dim
                tmp_gidx(1, ll) = start + ll ! position
                tmp_gidx(2, ll) = this%gidx(this%share(start + ll)) ! gidx
             end do
             ! sort local and remote objects with respect to global id
             key(1) = 2
             call sort_tuple(tmp_gidx, lda1, dim, key, 1, ind, aa1_gidx)
             call sort_tuple(rbuf_gidx, lda1, dim, key, 1, ind, aa1_gidx)
             ! compare object global id
             do ll = 1, dim
                if (tmp_gidx(2, ll) .eq. rbuf_gidx(2, ll)) then
                   this%sharemap(int(tmp_gidx(1, ll), i4)) = &
                        int(rbuf_gidx(1, ll) , i4)
                else
                   call neko_error('Inconsistent global id of remote object')
                end if
             end do
             if (jl .lt. ncomm) jl = jl + 1
          end if

          ! send
          if (ifwait_dst) then
             call MPI_Wait(dst_req, MPI_STATUS_IGNORE, ierr)
             if (kl .lt. ncomm) kl = kl + 1
          end if
       end do
       deallocate(ind, tmp_gidx)
    end if

    ! Get local inverse mapping
    itmp = nobj * nelt
    allocate(mapl(lda2 , itmp), ind(itmp))
    do il = 1, nelt
       do jl = 1, nobj
          itmp = jl + (il - 1) * nobj
          mapl(1, itmp) = map(jl, il) ! local object number
          mapl(2, itmp) = il ! local element number
          mapl(3, itmp) = jl ! position in the element
       end do
    end do
    ! sort tuple with respect to object and element number
    itmp = nobj * nelt
    key(1) = 1
    key(2) = 2
    call sort_tuple(mapl, lda2, itmp, key, 2, ind, aa2)

    ! extract information
    allocate(this%lmap(2, itmp), this%lmapoff(this%lnum + 1))
    this%lmapoff(1) = 1
    this%lmap(1, 1) = mapl(2, 1) ! local element number
    this%lmap(2, 1) = mapl(3, 1) ! position in the element
    ! count local objects
    jl = 1
    ! count positions in map array
    kl = 1
    do il = 2, itmp
       if (jl .eq. mapl(1, il)) then
          kl = kl + 1
          this%lmap(1, kl) = mapl(2, il) ! local element number
          this%lmap(2, kl) = mapl(3, il) ! position in the element
       else
          jl = jl + 1
          ! every node should be present
          if (jl .ne. mapl(1, il)) &
               call neko_error('Local object not mapped to element')
          kl = kl + 1
          this%lmapoff(jl) = kl
          this%lmap(1, kl) = mapl(2, il) ! local element number
          this%lmap(2, kl) = mapl(3, il) ! position in the element
       end if
       if (jl .gt. this%lnum) &
            call neko_error('Too many local object in object mapping')
    end do
    if (jl .ne. this%lnum) &
         call neko_error('Wrong local object number in object mapping')
    if (kl .ne. itmp) &
         call neko_error('Wrong position in object mapping')
    this%lmapoff(this%lnum + 1) = itmp + 1
    deallocate(mapl, ind)

    ! Get global object mapping
    ! Get the buffer
    bmax = 0
    if (this%nrank .gt. 0) then
       do il = 1, this%nrank
          if (this%rank(il) .ne. pe_rank) then
             itmp = 0
             do jl = this%off(il), this%off(il + 1) - 1
                itmp = itmp + this%lmapoff(this%share(jl) + 1) - &
                     this%lmapoff(this%share(jl))
             end do
             bmax = max(bmax, itmp)
          end if
       end do
    end if
    call MPI_Allreduce(MPI_IN_PLACE, bmax, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    if (this%nrank .gt. 0) then
       allocate(rbuf(lda3, bmax), sbuf(lda3, bmax))
       ! exchange mapping data
       jl = 1
       kl = 1
       do il = 1, nngh
          ! receive
          if (ngh(il) .eq. ngh_src(1, jl)) then
             src = this%rank(ngh_src(2, jl))
             dim = lda3 * bmax
             call MPI_Irecv(rbuf, dim, MPI_INTEGER, src, 0, NEKO_COMM, &
                  src_req, ierr)
             ifwait_src = .true.
          else
             ifwait_src = .false.
          end if

          ! send
          if (ngh(il) .eq. ngh_dst(1, kl)) then
             dst = this%rank(ngh_dst(2, kl))
             dim = 0
             do ll = this%off(ngh_dst(2, kl)), this%off(ngh_dst(2, kl) + 1) - 1
                itmp = this%sharemap(ll) ! object number on receiver side
                do ml = this%lmapoff(this%share(ll)), &
                     this%lmapoff(this%share(ll) + 1) - 1
                   dim = dim + 1
                   sbuf(1, dim) = itmp ! remote local object id
                   sbuf(2, dim) = pe_rank ! MPI rank
                   sbuf(3, dim) = this%lmap(1, ml) ! local element number
                   sbuf(4, dim) = this%lmap(2, ml) ! position in the element
                end do
             end do
             dim = dim * lda3
             call MPI_Isend(sbuf, dim, MPI_INTEGER, dst, 0, NEKO_COMM, &
                  dst_req, ierr)
             ifwait_dst = .true.
          else
             ifwait_dst = .false.
          end if

          ! receive
          if (ifwait_src) then
             call MPI_Wait(src_req, status, ierr)
             call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)
             ! array length
             dim = n_recv / lda3
             ! concatenate data data
             if (allocated(mapl)) then
                itmp = size(mapl, 2)
             else
                itmp = 0
             end if
             allocate(vtmp(lda3, itmp + dim))
             if (itmp .gt. 0) then
                vtmp(:, 1: itmp) = mapl(:, 1: itmp)
             end if
             vtmp(:, itmp + 1: itmp + dim) = rbuf(:, 1: dim)
             ! Following code due to the problem with Cray Fortran :
             ! Version 18.0.1 on Dardel
             !call move_alloc(vtmp, mapl)
             if (allocated(mapl)) deallocate(mapl)
             allocate(mapl, source = vtmp)
             deallocate(vtmp)

             if (jl .lt. ncomm) jl = jl + 1
          end if

          ! send
          if (ifwait_dst) then
             call MPI_Wait(dst_req, MPI_STATUS_IGNORE, ierr)
             if (kl .lt. ncomm) kl = kl + 1
          end if
       end do

       ! order objects
       if (allocated(mapl)) then
          itmp = size(mapl, 2)
          allocate(ind(itmp))
          key(1) = 1
          key(2) = 2
          key(3) = 3
          call sort_tuple(mapl, lda3, itmp, key, 3, ind, aa3)
          ! count objects and mark object sections boundaries
          ind(:) = 0
          jl = 1
          kl = mapl(1, 1)
          ind(1) = 1
          do il = 2, itmp
             if (kl .ne. mapl(1, il)) then
                jl = jl + 1
                kl = mapl(1, il)
                ind(il) = 1
             end if
          end do
          ! sanity check; number of shared objects
          do il = 1, this%nrank
             if (this%rank(il) .eq. pe_rank) then
                ll = il
                exit
             end if
          end do
          if (jl .ne. this%off(ll + 1) - this%off(ll)) &
               call neko_error('Inconsistent local number of shared objects.')

          allocate(this%glist(jl), this%gmap(3, itmp), this%gmapoff(jl + 1))
          ! extract mapping
          jl = 0
          do il = 1, itmp
             this%gmap(:, il) = mapl(2: 4, il)
             if (ind(il) .eq.1) then
                jl = jl + 1
                this%glist(jl) = mapl(1, il)
                this%gmapoff(jl) = il
             end if
          end do
          this%gmapoff(jl + 1) = itmp + 1

          ! sanity check; local id of shared objects
          allocate(llist(jl))
          start = this%off(ll) - 1
          do il = 1, jl
             llist(il) = this%share(start + il)
          end do
          call sort(llist, ind, jl)
          do il = 1, jl
             if (llist(il) .ne. this%glist(il)) &
                  call neko_error('Inconsistent lists of shared objects')
          end do

          deallocate(mapl, ind, llist)
       end if
    end if

    if (this%nrank .gt. 0) then
       deallocate(rbuf_gidx, sbuf_gidx, ngh_src, ind_src, ngh_dst, ind_dst, &
            ngh, rbuf, sbuf)
    end if

  end subroutine manager_conn_obj_global_map

  !> Initialise connectivity object type based on another connectivity type
  !! @param[inout] conn   connectivity object data
  subroutine manager_conn_obj_init_type_p4est(this, conn)
    class(manager_conn_obj_p4est_t), intent(inout) :: this
    class(manager_conn_obj_t), intent(inout) :: conn

    call this%free()
    call this%init_type_base(conn)

    select type (conn)
    type is(manager_conn_obj_p4est_t)
       this%nrank = conn%nrank
       this%nshare = conn%nshare

       if (allocated(conn%rank)) call move_alloc(conn%rank, this%rank)
       if (allocated(conn%share)) call move_alloc(conn%share, this%share)
       if (allocated(conn%sharemap)) &
            call move_alloc(conn%sharemap, this%sharemap)
       if (allocated(conn%off)) call move_alloc(conn%off, this%off)
       if (allocated(conn%lmap)) call move_alloc(conn%lmap, this%lmap)
       if (allocated(conn%lmapoff)) call move_alloc(conn%lmapoff, this%lmapoff)
       if (allocated(conn%glist)) call move_alloc(conn%glist, this%glist)
       if (allocated(conn%gmap)) call move_alloc(conn%gmap, this%gmap)
       if (allocated(conn%gmapoff)) call move_alloc(conn%gmapoff, this%gmapoff)
    end select

  end subroutine manager_conn_obj_init_type_p4est

  !> Free connectivity object type
  subroutine manager_conn_obj_free_p4est(this)
    class(manager_conn_obj_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nrank = 0
    this%nshare = 0

    if (allocated(this%rank)) deallocate(this%rank)
    if (allocated(this%share)) deallocate(this%share)
    if (allocated(this%sharemap)) deallocate(this%sharemap)
    if (allocated(this%off)) deallocate(this%off)
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%lmapoff)) deallocate(this%lmapoff)
    if (allocated(this%glist)) deallocate(this%glist)
    if (allocated(this%gmap)) deallocate(this%gmap)
    if (allocated(this%gmapoff)) deallocate(this%gmapoff)

  end subroutine manager_conn_obj_free_p4est

  !> Allocate types
  subroutine manager_conn_init_p4est(this)
    class(manager_conn_p4est_t), intent(inout) :: this

    if (allocated(this%vrt))then
       call this%vrt%free()
       deallocate(this%vrt)
    end if
    allocate(manager_conn_obj_p4est_t::this%vrt)

    if (allocated(this%fcs))then
       call this%fcs%free()
       deallocate(this%fcs)
    end if
    allocate(manager_conn_obj_p4est_t::this%fcs)

    if (allocated(this%edg))then
       call this%edg%free()
       deallocate(this%edg)
    end if
    allocate(manager_conn_obj_p4est_t::this%edg)

  end subroutine manager_conn_init_p4est

  !> Initialise connectivity information
  !! @param[in]    tdim    topological dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertex mapping
  !! @param[inout] fmap    element face mapping
  !! @param[inout] falgn   element face alignment
  !! @param[inout] emap    element edge mapping
  !! @param[inout] ealgn   element edge alignment
  !! @param[inout] hngel   element hanging flag
  !! @param[inout] hngfc   element face hanging flag
  !! @param[inout] hnged   element edge hanging flag
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_init_data_p4est(this, tdim, nel, vmap, fmap, falgn, &
       emap, ealgn, hngel, hngfc, hnged, ifsave)
    class(manager_conn_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    logical, allocatable, dimension(:), intent(inout) :: hngel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap, fmap, &
         falgn, emap, ealgn, hngfc, hnged
    logical, optional, intent(in) :: ifsave
    integer(i4) :: nvert, nface, nedge

    nvert = 2**tdim
    nface = 2 * tdim
    nedge = 12 * (tdim - 2)

    ! sanity check
    if ((nvert .ne. size(vmap, 1)) .or. (nel .ne. size(vmap, 2)) .or. &
         (nface .ne. size(fmap, 1)) .or. (nel .ne. size(fmap, 2)) .or. &
         (nface .ne. size(falgn, 1)) .or. (nel .ne. size(falgn, 2)) .or. &
         (nedge .ne. size(emap, 1)) .or. (nel .ne. size(emap, 2)) .or. &
         (nedge .ne. size(ealgn, 1)) .or. (nel .ne. size(ealgn, 2)) .or. &
         (nel .ne. size(hngel)) .or. &
         (nface .ne. size(hngfc, 1)) .or. (nel .ne. size(hngfc, 2)) .or. &
         (nedge .ne. size(hnged, 1)) .or. (nel .ne. size(hnged, 2))) &
         call neko_error('Inconsistent array sizes; p4est%conn')

    if (present(ifsave)) then
       call this%free_data(ifsave)
       call this%init_data_base(tdim, nel, vmap, fmap, emap, ifsave)
    else
       call this%free_data()
       call this%init_data_base(tdim, nel, vmap, fmap, emap)
    end if

    if (allocated(falgn)) call move_alloc(falgn, this%falgn)
    if (allocated(ealgn)) call move_alloc(ealgn, this%ealgn)
    if (allocated(hngel)) call move_alloc(hngel, this%hngel)
    if (allocated(hngfc)) call move_alloc(hngfc, this%hngfc)
    if (allocated(hnged)) call move_alloc(hnged, this%hnged)

  end subroutine manager_conn_init_data_p4est

  !> Initialise connectivity type based on another connectivity type
  !! @param[inout] conn   connectivity data
  subroutine manager_conn_init_type_p4est(this, conn)
    class(manager_conn_p4est_t), intent(inout) :: this
    class(manager_conn_t), intent(inout) :: conn

    call this%free_data()

    call this%init_type_base(conn)

    select type (conn)
    type is(manager_conn_p4est_t)
       if (allocated(conn%falgn)) call move_alloc(conn%falgn, this%falgn)
       if (allocated(conn%ealgn)) call move_alloc(conn%ealgn, this%ealgn)
       if (allocated(conn%hngel)) call move_alloc(conn%hngel, this%hngel)
       if (allocated(conn%hngfc)) call move_alloc(conn%hngfc, this%hngfc)
       if (allocated(conn%hnged)) call move_alloc(conn%hnged, this%hnged)
    end select

  end subroutine manager_conn_init_type_p4est

  !> Free connectivity information
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_free_data_p4est(this, ifsave)
    class(manager_conn_p4est_t), intent(inout) :: this
    logical, optional, intent(in) :: ifsave

    if (present(ifsave)) then
       call this%free_data_base(ifsave)
    else
       call this%free_data_base()
    end if

    if (allocated(this%falgn)) deallocate(this%falgn)
    if (allocated(this%ealgn)) deallocate(this%ealgn)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)

  end subroutine manager_conn_free_data_p4est

  !> Free connectivity information
  subroutine manager_conn_free_p4est(this)
    class(manager_conn_p4est_t), intent(inout) :: this

    call this%free_base()

    if (allocated(this%falgn)) deallocate(this%falgn)
    if (allocated(this%ealgn)) deallocate(this%ealgn)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)

  end subroutine manager_conn_free_p4est

end module manager_conn_p4est
