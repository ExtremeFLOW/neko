! Copyright (c) 2023-2024, The Neko Authors
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
!
!> Defines a factory subroutine for point zones.
submodule (point_zone) point_zone_fctry
  use box_point_zone, only: box_point_zone_t
  use sphere_point_zone, only: sphere_point_zone_t
  use cylinder_point_zone, only: cylinder_point_zone_t
  use json_utils, only: json_get
  use utils, only : neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: POINTZ_KNOWN_TYPES(3) = [character(len=20) :: &
       "box", &
       "sphere", &
       "cylinder"]

contains

  !> Point zone factory. Constructs, initializes, and maps the
  !! point zone object.
  !! @param object The object allocated by the factory.
  !! @param json JSON object initializing the point zone.
  !! @param dof Dofmap from which to map the point zone.
  module subroutine point_zone_factory(object, json, dof)
    class(point_zone_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(inout), optional :: dof
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    call json_get(json, "geometry", type_name)

    call point_zone_allocator(object, type_name)

    if (present(dof)) then
       call object%init(json, dof%size())
       call object%map(dof)
       call object%finalize()
    else
       ! This is for initializing zones inside a combine point zone,
       ! as they don't need to be mapped
       call object%init(json, 1)
       call object%finalize()
    end if
  end subroutine point_zone_factory

  module subroutine point_zone_allocator(object, type_name)
    class(point_zone_t), allocatable, intent(inout) :: object
    character(len=:), allocatable, intent(in) :: type_name
    integer :: i

    if (allocated(object)) deallocate(object)

    select case (trim(type_name))
    case ('box')
       allocate(box_point_zone_t::object)
    case ('sphere')
       allocate(sphere_point_zone_t::object)
    case ('cylinder')
       allocate(cylinder_point_zone_t::object)
    case default
       do i = 1, point_zone_registry_size
          if (trim(type_name) .eq. trim(point_zone_registry(i)%type_name)) then
             call point_zone_registry(i)%allocator(object)
             return
          end if
       end do
       call neko_error("Unknown point zone type: " // trim(type_name))
    end select
  end subroutine point_zone_allocator

  !> Register a custom point zone allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_point_zone(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(point_zone_allocate), pointer, intent(in) :: allocator
    type(allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(POINTZ_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(POINTZ_KNOWN_TYPES(i))) then
          call neko_type_registration_error("point zone", type_name, .true.)
       end if
    end do

    do i = 1, point_zone_registry_size
       if (trim(type_name) .eq. trim(point_zone_registry(i)%type_name)) then
          call neko_type_registration_error("point zone", type_name, .false.)
       end if
    end do

    ! Expand registry
    if (point_zone_registry_size .eq. 0) then
       allocate(point_zone_registry(1))
    else
       allocate(temp(point_zone_registry_size + 1))
       temp(1:point_zone_registry_size) = point_zone_registry
       call move_alloc(temp, point_zone_registry)
    end if

    point_zone_registry_size = point_zone_registry_size + 1
    point_zone_registry(point_zone_registry_size)%type_name = type_name
    point_zone_registry(point_zone_registry_size)%allocator => allocator

  end subroutine register_point_zone

end submodule point_zone_fctry
