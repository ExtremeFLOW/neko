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
! Implements a sphere geometry subset

module sphere_point_zone
  use point_zone, only: point_zone_t
  use num_types, only: rp
  use dofmap, only: dofmap_t
  use json_utils, only: json_get
  use json_module, only: json_file
  use math, only: abscmp
  implicit none
  private

  type, public, extends(point_zone_t) :: sphere_point_zone_t
     real(kind=rp) :: x0
     real(kind=rp) :: y0
     real(kind=rp) :: z0
     real(kind=rp) :: radius
   contains
     procedure, pass(this) :: init => sphere_point_zone_init_from_json
     procedure, pass(this) :: free => sphere_point_zone_free
     procedure, pass(this) :: criterion => sphere_point_zone_criterion
  end type sphere_point_zone_t

contains

  subroutine sphere_point_zone_init_from_json(this, json, dof)
    class(sphere_point_zone_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(inout) :: dof

    character(len=:), allocatable :: str_read
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: value
    real(kind=rp) :: x0, y0, z0, radius

    call json_get(json, "center", values)
    x0 = values(1)
    y0 = values(2)
    z0 = values(3)
    call json_get(json, "radius", value)
    radius = value
    call json_get(json, "name", str_read)

    call sphere_point_zone_init_common(this, dof%size(), trim(str_read), x0, &
         y0, z0, radius)

    call this%map(dof)

    call this%finalize()

  end subroutine sphere_point_zone_init_from_json

  subroutine sphere_point_zone_init_common(this, size, name, x0, y0, z0, radius)
    class(sphere_point_zone_t), intent(inout) :: this
    integer, intent(in), optional :: size
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: x0
    real(kind=rp), intent(in) :: y0
    real(kind=rp), intent(in) :: z0
    real(kind=rp), intent(in) :: radius

    call this%init_base(size, name)

    this%x0 = x0
    this%y0 = y0
    this%z0 = z0
    this%radius = radius

  end subroutine sphere_point_zone_init_common

  subroutine sphere_point_zone_free(this)
    class(sphere_point_zone_t), intent(inout) :: this

    call this%free_base()

    this%x0 = 0.0_rp
    this%y0 = 0.0_rp
    this%z0 = 0.0_rp
    this%radius = 0.0_rp

  end subroutine sphere_point_zone_free

  pure function sphere_point_zone_criterion(this, x, y, z, ix, iy, iz, ie) result(is_inside)
    class(sphere_point_zone_t), intent(in) :: this
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    logical :: is_inside

    real(kind=rp) :: dist_from_center

    dist_from_center = sqrt( (x - this%x0)**2 + (y - this%y0)**2 + (z - this%z0)**2 )

    ! Inside if distance from center <= radius
    is_inside = (dist_from_center .lt. this%radius .or. &
                 abscmp(dist_from_center, this%radius))

  end function sphere_point_zone_criterion

end module sphere_point_zone
