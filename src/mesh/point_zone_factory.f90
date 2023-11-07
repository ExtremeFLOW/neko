! Copyright (c) 2023, The Neko Authors
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
module point_zone_fctry
  use point_zone, only: point_zone_t
  use box_point_zone, only: box_point_zone_t
  use sphere_point_zone, only: sphere_point_zone_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use dofmap, only: dofmap_t
  use utils, only: neko_error
  implicit none
  private

  public :: point_zone_factory

  contains
    
    !> Point zone factory. Constructs, initializes, and maps the
    !! point zone object.
    !! @param json JSON object initializing the point zone.
    !! @param dof Dofmap from which to map the point zone.
    subroutine point_zone_factory(point_zone, json, dof)
      class(point_zone_t), allocatable, intent(inout) :: point_zone
      type(json_file), intent(inout) :: json
      type(dofmap_t), intent(inout) :: dof
      character(len=:), allocatable :: zone_type

      call json_get(json, "geometry", zone_type)

      if (trim(zone_type) .eq. "box") then
         allocate(box_point_zone_t::point_zone)
      else if (trim(zone_type) .eq. "sphere") then
         allocate(sphere_point_zone_t::point_zone)
      else
         call neko_error("Unknown source term "//trim(zone_type)//"! Valid &
              &source terms are 'box', 'sphere'.")
      end if

      call point_zone%init(json, dof%size())

      call point_zone%map(dof)
      call point_zone%finalize()

    end subroutine point_zone_factory

end module point_zone_fctry
