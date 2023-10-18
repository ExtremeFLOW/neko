! Copyright (c) 2019-2021, The Neko Authors
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
! Implements a zone as a subset of GLL points in the mesh
module point_zone
  use stack, only: stack_i4_t
  use num_types, only: rp
  use utils, only: neko_error, nonlinear_index
  use dofmap, only: dofmap_t
  use json_module, only: json_file
  implicit none
  private

  type, public, abstract :: point_zone_t
     integer, allocatable :: mask(:)
     type(stack_i4_t), private :: scratch
     integer :: size
     logical, private :: finalized = .false.
     character(len=200) :: name
   contains
     procedure, pass(this) :: init_base => point_zone_init_base
     procedure, pass(this) :: free_base => point_zone_free_base
     procedure, pass(this) :: finalize => point_zone_finalize
     procedure, pass(this) :: add => point_zone_add
     procedure, pass(this) :: map => point_zone_map
     !> The common constructor using a JSON object.
     procedure(point_zone_init), pass(this), deferred :: init
     !> Destructor.
     procedure(point_zone_free), pass(this), deferred :: free
     procedure(point_zone_criterion), pass(this), deferred :: criterion
     
  end type point_zone_t

  abstract interface
     pure function point_zone_criterion(this, x, y, z, ix, iy, iz, ie) result(is_inside)
       import :: point_zone_t
       import :: rp
       class(point_zone_t), intent(in) :: this
       real(kind=rp), intent(in) :: x
       real(kind=rp), intent(in) :: y
       real(kind=rp), intent(in) :: z
       integer, intent(in) :: ix
       integer, intent(in) :: iy
       integer, intent(in) :: iz
       integer, intent(in) :: ie
       logical :: is_inside
     end function point_zone_criterion
  end interface

  abstract interface
     subroutine point_zone_init(this, json, dof)
       import :: point_zone_t
       import :: json_file
       import :: dofmap_t
       class(point_zone_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(dofmap_t), intent(inout) :: dof
     end subroutine point_zone_init
  end interface

  abstract interface
     subroutine point_zone_free(this)
       import :: point_zone_t
       class(point_zone_t), intent(inout) :: this
     end subroutine point_zone_free
  end interface

contains

  !> Initialize a point zone
  subroutine point_zone_init_base(this, size, name)
    class(point_zone_t), intent(inout) :: this
    integer, intent(in), optional :: size
    character(len=*), intent(in) :: name

    call point_zone_free_base(this)

    if (present(size)) then
       call this%scratch%init(size)
    else
       call this%scratch%init()
    end if

    this%name = trim(name)

  end subroutine point_zone_init_base

  !> Deallocate a point zone
  subroutine point_zone_free_base(this)
    class(point_zone_t), intent(inout) :: this
    if (allocated(this%mask)) then
       deallocate(this%mask)
    end if

    this%finalized = .false.
    this%size = 0

    call this%scratch%free()
    
  end subroutine point_zone_free_base

  !> Finalize a zone list
  !! @details Create a static list of integers
  subroutine point_zone_finalize(this)
    class(point_zone_t), intent(inout) :: this
    integer, pointer :: tp(:)
    integer :: i
    
    if (.not. this%finalized) then

       allocate(this%mask(this%scratch%size()))
       
       tp => this%scratch%array()
       do i = 1, this%scratch%size()
          this%mask(i) = tp(i)
       end do

       this%size = this%scratch%size()
       
       call this%scratch%clear()

       this%finalized = .true.
       
    end if
    
  end subroutine point_zone_finalize

  !> Add a point linear index to an unfinalized zone
  subroutine point_zone_add(this, idx)
    class(point_zone_t), intent(inout) :: this
    integer, intent(inout) :: idx
    
    if (this%finalized) then
       call neko_error('Point zone already finalized')
    end if
    
    call this%scratch%push(idx)
    
  end subroutine point_zone_add

  !> Map the GLL points which are in the point_zone, according to
  !! the given criterion (zone specific)
  !! @param dof Dofmap of points to go through
  subroutine point_zone_map(this, dof)
    class(point_zone_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dof

    integer :: i, ix, iy, iz, ie, nlindex(4), lx, idx
    real(kind=rp) :: x, y, z

    lx = dof%Xh%lx

    do i = 1, dof%size()
       nlindex = nonlinear_index(i, lx, lx, lx)
       x = dof%x(nlindex(1), nlindex(2), nlindex(3), nlindex(4))
       y = dof%y(nlindex(1), nlindex(2), nlindex(3), nlindex(4))
       z = dof%z(nlindex(1), nlindex(2), nlindex(3), nlindex(4))
       ix = nlindex(1)
       iy = nlindex(2)
       iz = nlindex(3)
       ie = nlindex(4)

       if (this%criterion(x, y, z, ix, iy, iz, ie)) then
          idx = i
          call this%add(idx)
       end if
    end do

  end subroutine point_zone_map

end module point_zone
