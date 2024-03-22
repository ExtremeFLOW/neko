! Copyright (c) 2018-2024, The Neko Authors
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
!> Abstract type for defining object ownership, mapping and sharing
!! @details These abstract types are bases for different type object lists
!! e.g, points, elements or topology objects.
module obj_sharing
  use num_types, only : i4
  use utils, only : neko_error
  use comm, only : NEKO_COMM, MPI_INTEGER, MPI_SUM, pe_rank, pe_size
  use stack, only : stack_i4_t

  implicit none
  private

  public :: obj_ownership_t, obj_mapping_t, obj_sharing_t

  !> Base type for the lists with defined object ownership
  type, abstract :: obj_ownership_t
     !> Number of local objects
     integer(i4), private :: lnum_ = -1
     !> Global number of objects
     integer(i4), private :: gnum_ = -1
     !> Global object offset based on local owned objects
     integer(i4), private :: goff_ = -1
     !> MPI id of the owner process
     integer(i4), dimension(:), allocatable :: owner_pid
     !> Number of local owned objects
     integer(i4), private :: lown_ = -1
     !> List of local numbers of local owned objects
     integer(i4), dimension(:), allocatable :: owned_list
     !> MPI ranks sharing owned object
     type(stack_i4_t), dimension(:), allocatable :: owned_ranks
     !> Global to local mapping
   contains
     !> Set space for ownership data
     procedure, pass(this) :: init_own => obj_ownership_init_own
     !> Free ownership data
     procedure, pass(this) :: free_own => obj_ownership_free_own
     !> Local number of objects
     procedure, pass(this) :: lnum => obj_ownership_lnum
     !> Global number of objects
     procedure, pass(this) :: gnum => obj_ownership_gnum
     !> Global object offset
     procedure, pass(this) :: goff => obj_ownership_goff
     !> Local number of owned objects
     procedure, pass(this) :: lown => obj_ownership_lown
     !> Global to local id mapping
     procedure, pass(this) :: glb_loc_id => obj_ownership_glb_loc_id
     !> Local to global id mapping
     procedure(obj_ownership_loc_glb_id), pass(this), deferred :: loc_glb_id
  end type obj_ownership_t

  !> Base type for the lists with defined object mapping
  type, extends(obj_ownership_t), abstract :: obj_mapping_t
   contains
     !> Free mapping data
     procedure, pass(this) :: free_map => obj_mapping_free_map
  end type obj_mapping_t

  !> Base type for the lists with defined object sharing
  type, extends(obj_mapping_t), abstract :: obj_sharing_t
   contains
     !> Free mapping data
     procedure, pass(this) :: free_share => obj_sharing_free_share
  end type obj_sharing_t

  !> Returns global object id based on it's local position
  !! @parameter[in]   lid   object local id position
  !! @return gid
  abstract interface
     function obj_ownership_loc_glb_id(this, lid) result(gid)
       import i4
       import obj_ownership_t
       class(obj_ownership_t), intent(in) :: this
       integer(i4), intent(in) :: lid
       integer(i4) :: gid
     end function obj_ownership_loc_glb_id
  end interface

contains

  !> Initialise ownership data
  !> @details The ownership data consists of list of global object ids and MPI
  !! rank id of the owning process. In addition global-to-local and
  !! local-to-global id mappings are added. There are possible simplifications
  !! of the general algorithm. If none of the objects is shared there is no need
  !! to identify owner MPI id. This situation is flagged by @a ifsimple
  !! argument. In some cases the ownership data can be given, so no additional
  !! operation is needed. In this situation data from @a owner_pid and
  !! (optionally) @a owned_ranks are used. In other case some addition
  !! constraint to the ownership data may be available (e.g, edges can be owned
  !! by ranks owning vertices only). In this case @a owner_cnst is used.
  !! @parameter[in]    gid_list    list of object global ids
  !! @parameter[in]    ifsimple    are objects shared among MPI ranks
  !! @parameter[inout] own_copy    available ownership to copy
  !! @parameter[inout] own_cnst    available ownership constraint
  subroutine obj_ownership_init_own(this, gid_list, ifsimple, owner_pid, &
       & owned_ranks, owner_cnst)
    class(obj_ownership_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: gid_list
    logical, intent(in) :: ifsimple
    integer(i4), dimension(size(gid_list)), optional, intent(inout) :: owner_pid
    type(stack_i4_t), dimension(:), allocatable, optional, intent(inout) :: &
         & owned_ranks
    integer(i4), dimension(:), optional, intent(inout) :: owner_cnst
    integer(i4) :: il, itmp, ierr

    call this%free_own()

    ! Simple distribution; no shared objects
    if (ifsimple .or. (pe_size == 1)) then
       this%lnum_ = size(gid_list)
       this%lown_ = this%lnum_
       allocate(this%owner_pid(this%lnum_), this%owned_list(this%lown_), &
            & this%owned_ranks(this%lown_))
       this%owner_pid(:) = pe_rank
       ! as little as possible
       itmp = 1
       do il = 1, this%lown_
          this%owned_list(il) = il
          call this%owned_ranks(il)%init(itmp)
       end do
    else if (present(owner_pid)) then
       this%lnum_ = size(gid_list)
       allocate(this%owner_pid(this%lnum_))
       this%owner_pid(:) = owner_pid(:)
       ! count owned objects
       itmp = 0
       do il = 1, this%lnum_
          if (owner_pid(il) == pe_rank) itmp = itmp + 1
       end do
       this%lown_ = itmp
       ! get owned object list
       if (this%lown_ > 0) then
          allocate(this%owned_list(this%lown_))
          itmp = 0
          do il = 1, this%lnum_
             if (owner_pid(il) == pe_rank) then
                itmp = itmp + 1
                this%owned_list(itmp) = il
             end if
          end do
          ! get list of MPI ranks sharing a given object
          if (present(owned_ranks)) then
             ! sanity check
             if (this%lown_ /= size(owned_ranks)) &
                  &call neko_error('Inconsistent owned_ranks size')
             call move_alloc(owned_ranks, this%owned_ranks)
          else
             allocate(this%owned_list(this%lown_))



             
          end if
       end if
    else if (present(owner_cnst)) then

    else

    end if

    ! global to local mapping


    ! get global object number (based on owned objects) and offset
    call MPI_Allreduce(this%lown_, this%gnum_, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    this%goff_ = 0
    call MPI_Exscan(this%lown_, this%goff_, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

  end subroutine obj_ownership_init_own

  !> Free ownership data
  subroutine obj_ownership_free_own(this)
    class(obj_ownership_t), intent(inout) :: this
    integer(i4) :: il
    this%lnum_ = -1
    this%gnum_ = -1
    this%goff_ = -1
    this%lown_ = -1
    if (allocated(this%owner_pid)) deallocate(this%owner_pid)
    if (allocated(this%owned_list)) deallocate(this%owned_list)
    if (allocated(this%owned_ranks)) then
       do il = 1, size(this%owned_ranks)
          call this%owned_ranks(il)%free()
       end do
       deallocate(this%owned_ranks)
    end if
  end subroutine obj_ownership_free_own

  !> Local number of objects
  !! @return itmp
  function obj_ownership_lnum(this) result(itmp)
    class(obj_ownership_t), intent(in) :: this
    integer(i4) :: itmp
    itmp = this%lnum_
  end function obj_ownership_lnum

  !> Global number of objects
  !! @return itmp
  function obj_ownership_gnum(this) result(itmp)
    class(obj_ownership_t), intent(in) :: this
    integer(i4) :: itmp
    itmp = this%gnum_
  end function obj_ownership_gnum

  !> Global object offset
  !! @return itmp
  function obj_ownership_goff(this) result(itmp)
    class(obj_ownership_t), intent(in) :: this
    integer(i4) :: itmp
    itmp = this%goff_
  end function obj_ownership_goff

  !> Local number of owned objects
  !! @return itmp
  function obj_ownership_lown(this) result(itmp)
    class(obj_ownership_t), intent(in) :: this
    integer(i4) :: itmp
    itmp = this%lown_
  end function obj_ownership_lown

  !> Returns local object position based on it's global id
  !! @parameter[in]   gid   object local id position
  !! @return lid
  function obj_ownership_glb_loc_id(this, gid) result(lid)
    class(obj_ownership_t), intent(in) :: this
    integer(i4), intent(in) :: gid
    integer(i4) :: lid
    lid = 0
  end function obj_ownership_glb_loc_id

  !> Free mapping data
  subroutine obj_mapping_free_map(this)
    class(obj_mapping_t), intent(inout) :: this

    call this%free_own()

  end subroutine obj_mapping_free_map

  !> Free sharing data
  subroutine obj_sharing_free_share(this)
    class(obj_sharing_t), intent(inout) :: this

    call this%free_map()

  end subroutine obj_sharing_free_share

end module obj_sharing
