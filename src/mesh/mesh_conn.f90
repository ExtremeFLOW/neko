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
!> Implementation of the mesh connectivity type
module mesh_conn
  use num_types, only : i4, i8
  use utils, only : neko_error
  use comm, only : pe_size

  implicit none
  private

  !> Type for connectivity information regarding vertices, faces and edges.
  !! @details It contains global numbering of objects and element mapping
  !! information
  type, public :: mesh_conn_obj_t
     !> Number of local objects
     integer(i4) :: lnum
     !> Global number of objects
     integer(i8) :: gnum
     !> Global indexing of unique objects of given type
     integer(i8), allocatable, dimension(:) :: gidx
     !> Logical flag indicating sharing with other MPI ranks
     ! size of lnum
     logical, allocatable, dimension(:) :: lshare
     !> Number of shared objects
     integer(i4) :: nshare
     !> List of shared objects (local number)
     integer(i4), allocatable, dimension(:) :: sharelist
     !> Global ownership (MPI rank) of shared objects
     integer(i4), allocatable, dimension(:) :: shareown
     !> Number of MPI ranks sharing objects
     integer(i4) :: nrank
     !> List of ranks sharing objects
     integer(i4), allocatable, dimension(:) :: rank
     !> List of shared objects with respect to MPI rank
     integer(i4), allocatable, dimension(:) :: rankshare
     !> Mapping of the shared objects local number to the remote local number
     integer(i4), allocatable, dimension(:) :: sharemap
     !> Offset in the rankshare and sharemap lists
     ! size of nrank + 1
     integer(i4), allocatable, dimension(:) :: rankoff
     !> Local number of elements
     integer(i4) :: nel
     !> Number of objects per element
     integer(i4) :: nobj
     !> Element to object mapping
     integer(i4), allocatable, dimension(:,:) :: map
     !> Local mapping of object to element component
     ! inverse to map; related to the local numbering
     ! first dimension: 1- local element number, 2 - position in element
     integer(i4), allocatable, dimension(:,:) :: lmap
     !> Offset in the local mapping
     ! size of lnum + 1
     integer(i4), allocatable, dimension(:) :: lmapoff
     !> Global mapping of object to element component
     ! like lmap, but for shared objects; related to sharelist
     ! first dimension: 1 - MPI rank, 2- local element number,
     ! 3 - position in element
     integer(i4), allocatable, dimension(:,:) :: gmap
     !> Offset in the global mapping
     ! size of nshare + 1
     integer(i4), allocatable, dimension(:) :: gmapoff
     !> Is object alignment directly specified
     logical :: ifalgn
     !> Object alignment
     integer(i4), allocatable, dimension(:,:) :: algn
     !> Are hanging object directly specified
     logical :: ifhang_set
     !> Is there any hanging object
     logical :: ifhang
     !> Hanging object; position of hanging object, otherwise -1
     integer(i4), allocatable, dimension(:,:) :: hang
   contains
     procedure, pass(this) :: init => mesh_conn_obj_init
     procedure, pass(this) :: free => mesh_conn_obj_free
  end type mesh_conn_obj_t

  !> Type for element connectivity information
  type, public :: mesh_conn_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Local number of elements
     integer(i4) :: nel
     !> Connectivity information for vertices, edges and faces
     type(mesh_conn_obj_t) :: vrt, fcs, edg
     !> Are elements with hanging object directly specified
     logical :: ifhang_set
     !> Is there any element with hanging objects
     logical :: ifhang
     !> Elements with hanging object
     logical, allocatable, dimension(:) :: hang
   contains
     procedure, pass(this) :: init => mesh_conn_init
     procedure, pass(this) :: free => mesh_conn_free
  end type mesh_conn_t

contains

  !> Initialise object connectivity information
  !! @param[in]    gnum       global number of objects
  !! @param[in]    gidx       object global id
  !! @param[in]    map        element to objects mapping
  !! @param[in]    lmap       Local mapping of object to element component
  !! @param[in]    lmapoff    Offset in the local mapping
  !! @param[in]    sharelist  List of shared objects (local number)
  !! @param[in]    rank       Number of MPI ranks sharing objects
  !! @param[in]    rankshare  List of shared objects with respect to MPI rank
  !! @param[in]    sharemap   Mapping of shared objects lid to the remote lid
  !! @param[in]    rankoff    Offset in the rankshare and sharemap lists
  !! @param[in]    gmap       Global mapping of object to element component
  !! @param[in]    gmapoff    Offset in the global mapping
  !! @param[in]    algn       object alignment information
  !! @param[in]    hang       object hanging information
  subroutine mesh_conn_obj_init(this, gnum, gidx, map, lmap, lmapoff, &
       sharelist, rank, rankshare, sharemap, rankoff, gmap, gmapoff, algn, hang)
    class(mesh_conn_obj_t), intent(inout) :: this
    integer(i8), intent(in) :: gnum
    integer(i8), dimension(:), intent(in) :: gidx
    integer(i4), dimension(:, :), intent(in) :: map
    integer(i4), dimension(:), optional, intent(in) :: sharelist, rank, &
         rankshare, sharemap, rankoff, lmapoff, gmapoff
    integer(i4), dimension(:, :), optional, intent(in) :: lmap, gmap, algn, hang
    integer :: il, jl, itmp
    integer(i4), dimension(:), allocatable :: lown

    call this%free()

    this%gnum = gnum

    this%lnum = size(gidx)
    this%nel = size(map, 2)
    this%nobj = size(map, 1)
    allocate(this%gidx, source = gidx)
    allocate(this%map, source = map)
    ! basic sharing info
    allocate(this%lshare(this%lnum))
    this%lshare(:) = .false.

    if (present(lmap) .and. present(lmapoff)) then
       ! sanity check
       if ((this%lnum + 1) .ne. size(lmapoff) .or. 2 .ne. size(lmap, 1) &
            .or. (lmapoff(this%lnum + 1) - 1) .ne. size(lmap, 2)) &
            call neko_error('Inconsistent array sizes; conn_obj%lmap')
       allocate(this%lmap, source = lmap)
       allocate(this%lmapoff, source = lmapoff)
    end if

    if (present(sharelist)) then
       this%nshare = size(sharelist)
       allocate(this%sharelist, source = sharelist)
       ! sanity check
       il = minval(this%sharelist)
       jl = maxval(this%sharelist)
       if (il .le. 0 .or. jl .gt. this%lnum) &
            call neko_error('Shared node id out of range; conn_obj%sharelist')
       do il = 1, this%nshare
          this%lshare(this%sharelist(il)) = .true.
       end do

       if (present(gmap) .and. present(gmapoff)) then
          ! sanity check
          if ((this%nshare + 1) .ne. size(gmapoff) .or. 3 .ne. size(gmap, 1) &
               .or. (gmapoff(this%nshare + 1) - 1) .ne. size(gmap, 2)) &
               call neko_error('Inconsistent array sizes; conn_obj%gmap')
          allocate(this%gmap, source = gmap)
          allocate(this%gmapoff, source = gmapoff)
       end if
    else
       this%nshare = 0
    end if

    if (present(rank) .and. present(rankshare) .and. present(sharemap) &
         .and. present(rankoff)) then
       this%nrank = size(rank)
       ! sanity check
       if ((this%nrank + 1) .ne. size(rankoff) .or. &
            (rankoff(this%nrank + 1) - 1) .ne. size(rankshare) .or. &
            (rankoff(this%nrank + 1) - 1) .ne. size(sharemap)) &
            call neko_error('Inconsistent array sizes; conn_obj%rank')
       allocate(this%rank, source = rank)
       allocate(this%rankshare, source = rankshare)
       allocate(this%sharemap, source = sharemap)
       allocate(this%rankoff, source = rankoff)
    else
       this%nrank = 0
    end if

    if (present(algn)) then
       ! sanity check
       if (this%nel .ne. size(algn, 2) .or. this%nobj .ne. size(algn, 1)) &
            call neko_error('Inconsistent array sizes; conn_obj%algn')
       this%ifalgn = .true.
       allocate(this%algn, source = algn)
    end if

    if (present(hang)) then
       ! sanity check
       if (this%nel .ne. size(hang, 2) .or. this%nobj .ne. size(hang, 1)) &
            call neko_error('Inconsistent array sizes; conn_obj%hang')
       this%ifhang_set = .true.
       allocate(this%hang, source = hang)
       this%ifhang = maxval(hang) .ne. -1
    end if

    ! Get global object ownership
    if (allocated(this%sharelist).and. allocated(this%rankshare) .and. &
         allocated(this%rankoff)) then
       allocate(this%shareown(this%nshare))
       this%shareown(:) = -1
       ! find the unique global ownership (MPI rank) of the object
       ! I use the simplest method not requiring communication assigning
       ! ownership to the node with the smallest MPI rank. This is efficient,
       ! but does not guarantee the uniform object distribution among nodes.
       allocate(lown(this%lnum))
       lown(:) = pe_size
       do il = 1, this%nrank
          itmp = this%rank(il)
          do jl = this%rankoff(il), this%rankoff(il + 1) - 1
             lown(this%rankshare(jl)) = min(lown(this%rankshare(jl)), itmp)
          end do
       end do
       do il = 1, this%nshare
          if (lown(this%sharelist(il)) .ne. pe_size) then
             this%shareown(il) = lown(this%sharelist(il))
          else
             call neko_error('Unknown object ownership; conn_obj_t')
          end if
       end do
       deallocate(lown)
    end if

  end subroutine mesh_conn_obj_init

  !> Free object connectivity information
  subroutine mesh_conn_obj_free(this)
    class(mesh_conn_obj_t), intent(inout) :: this

    this%lnum = 0
    this%gnum = 0
    this%nshare = 0
    this%nrank = 0
    this%nel = 0
    this%nobj = 0
    this%ifalgn = .false.
    this%ifhang_set = .false.
    this%ifhang = .false.

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%lshare)) deallocate(this%lshare)
    if (allocated(this%sharelist)) deallocate(this%sharelist)
    if (allocated(this%shareown)) deallocate(this%shareown)
    if (allocated(this%rank)) deallocate(this%rank)
    if (allocated(this%rankshare)) deallocate(this%rankshare)
    if (allocated(this%sharemap)) deallocate(this%sharemap)
    if (allocated(this%rankoff)) deallocate(this%rankoff)
    if (allocated(this%map)) deallocate(this%map)
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%lmapoff)) deallocate(this%lmapoff)
    if (allocated(this%gmap)) deallocate(this%gmap)
    if (allocated(this%gmapoff)) deallocate(this%gmapoff)
    if (allocated(this%algn)) deallocate(this%algn)
    if (allocated(this%hang)) deallocate(this%hang)

  end subroutine mesh_conn_obj_free

  !> Initialise mesh connectivity information
  !! @param[in]    tdim     topological dimension
  !! @param[in]    nel      local number of elements
  !! @param[in]    hang     element hanging information
  subroutine mesh_conn_init(this, tdim, nel, hang)
    class(mesh_conn_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    logical, dimension(:), optional, intent(in) :: hang

    call this%free()

    this%tdim = tdim
    this%nel = nel

    if (present(hang)) then
       ! sanity check
       if (nel .ne. size(hang)) &
            call neko_error('Inconsistent array sizes; conn%hang')
       this%ifhang_set = .true.
       allocate(this%hang(nel))
       this%hang(:) = hang(:)
       this%ifhang = any(hang)
    end if

  end subroutine mesh_conn_init

  !> Free mesh connectivity information
  subroutine mesh_conn_free(this)
    class(mesh_conn_t), intent(inout) :: this

    call this%vrt%free()
    call this%fcs%free()
    call this%edg%free()

    this%tdim = 0
    this%nel = 0
    this%ifhang_set = .false.
    this%ifhang = .false.

    if (allocated(this%hang)) deallocate(this%hang)

  end subroutine mesh_conn_free

end module mesh_conn
