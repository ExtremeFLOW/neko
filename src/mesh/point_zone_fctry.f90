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
  use utils, only : concat_string_array
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

    if (trim(type_name) .eq. "box") then
       allocate(box_point_zone_t::object)
    else if (trim(type_name) .eq. "sphere") then
       allocate(sphere_point_zone_t::object)
    else if (trim(type_name) .eq. "cylinder") then
       allocate(cylinder_point_zone_t::object)
    else
       type_string =  concat_string_array(POINTZ_KNOWN_TYPES, &
            new_line('A') // "-  ", .true.)
       call neko_error("Unknown point zone type: " &
                       // trim(type_name) // ".  Known types are: " &
                       // type_string)
    end if

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

end submodule point_zone_fctry
