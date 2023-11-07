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
  use neko_config, only: NEKO_BCKND_DEVICE
  use device
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  private

  !> Base abstract type for point zones.
  type, public, abstract :: point_zone_t
     !> List of linear indices of the GLL points in the zone. 
     integer, allocatable :: mask(:)
     !> List of linear indices of the GLL points in the zone on the device.
     type(c_ptr) :: mask_d
     !> Scratch stack of integers to build the list mask.
     type(stack_i4_t), private :: scratch
     !> Size of the point zone mask.
     integer :: size = 0
     !> Flag to indicate if point_zone_finalize has been called and the mask
     !! has been built.
     logical, private :: finalized = .false.
     !> Name of the point zone (used for retrieval in the point_zone_registry).
     character(len=80) :: name
   contains
     !> Constructor for the point_zone_t base type.
     procedure, pass(this) :: init_base => point_zone_init_base
     !> Destructor for the point_zone_t base type.
     procedure, pass(this) :: free_base => point_zone_free_base
     !> Builds the mask from the scratch stack.
     procedure, pass(this) :: finalize => point_zone_finalize
     !> Adds a point's linear index to the scratch stack.
     procedure, pass(this) :: add => point_zone_add
     !> Maps the GLL points that verify a point_zone's `criterion` by adding
     !! them to the stack.
     procedure, pass(this) :: map => point_zone_map
     !> The common constructor using a JSON object.
     procedure(point_zone_init), pass(this), deferred :: init
     !> Destructor.
     procedure(point_zone_free), pass(this), deferred :: free
     !> Defines the criterion of selection of a GLL point to the point_zone.
     procedure(point_zone_criterion), pass(this), deferred :: criterion
  end type point_zone_t

  !> A helper type to build a list of polymorphic point_zones.
  type, public :: point_zone_wrapper_t
     class(point_zone_t), allocatable :: pz
  end type point_zone_wrapper_t

  abstract interface
     !> Defines the criterion of selection of a GLL point to the point_zone.
     !! @param x x-coordinate of the GLL point.
     !! @param y y-coordinate of the GLL point.
     !! @param z z-coordinate of the GLL point.
     !! @param j 1st nonlinear index of the GLL point.
     !! @param k 2nd nonlinear index of the GLL point.
     !! @param l 3rd nonlinear index of the GLL point.
     !! @param e element index of the GLL point.
     pure function point_zone_criterion(this, x, y, z, j, k, l, e) result(is_inside)
       import :: point_zone_t
       import :: rp
       class(point_zone_t), intent(in) :: this
       real(kind=rp), intent(in) :: x
       real(kind=rp), intent(in) :: y
       real(kind=rp), intent(in) :: z
       integer, intent(in) :: j
       integer, intent(in) :: k
       integer, intent(in) :: l
       integer, intent(in) :: e
       logical :: is_inside
     end function point_zone_criterion
  end interface

  abstract interface
     !> The common constructor using a JSON object.
     !! @param json Json object for the point zone.
     !! @param size Size with which to initialize the stack
     subroutine point_zone_init(this, json, size)
       import :: point_zone_t
       import :: json_file
       import :: dofmap_t
       class(point_zone_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       integer, intent(in) :: size
     end subroutine point_zone_init
  end interface

  abstract interface
     !> Destructor.
     subroutine point_zone_free(this)
       import :: point_zone_t
       class(point_zone_t), intent(inout) :: this
     end subroutine point_zone_free
  end interface

contains

  !> Constructor for the point_zone_t base type.
  !! @param size Size of the scratch stack.
  !! @param name Name of the point zone.
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

  !> Destructor for the point_zone_t base type.
  subroutine point_zone_free_base(this)
    class(point_zone_t), intent(inout) :: this
    if (allocated(this%mask)) then
       deallocate(this%mask)
    end if

    this%finalized = .false.
    this%size = 0

    call this%scratch%free()

    if (c_associated(this%mask_d)) then
       call device_free(this%mask_d)
    end if
    
  end subroutine point_zone_free_base

  !> Builds the mask from the scratch stack.
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

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%mask, this%mask_d, this%size)
          call device_memcpy(this%mask, this%mask_d, this%size, HOST_TO_DEVICE)
       end if

       this%finalized = .true.
       
    end if
    
  end subroutine point_zone_finalize

  !> Adds a point's linear index to the scratch stack.
  !! @param idx Linear index of the point to add.
  !! @note The linear index of a point `(j,k,l,e)` can be retrieved using the 
  !! subroutine `linear_index(j,k,l,e,lx)` in the `utils` module.
  subroutine point_zone_add(this, idx)
    class(point_zone_t), intent(inout) :: this
    integer, intent(inout) :: idx
    
    if (this%finalized) then
       call neko_error('Point zone already finalized')
    end if
    
    call this%scratch%push(idx)
    
  end subroutine point_zone_add

  !> Maps the GLL points that verify a point_zone's `criterion` by adding
  !! them to the stack.
  !! @param dof Dofmap of points to go through.
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
