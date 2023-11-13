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
! Implements a box geometry subset.
module box_point_zone
  use point_zone, only: point_zone_t
  use num_types, only: rp
  use json_utils, only: json_get
  use json_module, only: json_file
  use math, only: abscmp
  implicit none
  private

  !> A box-shaped point zone.
  !! @details As defined here, a box is described by its `x,y,z` bounds,
  !! specified in the json file as e.g. `"x_bounds": [<xmin>, <xmax>]"`, 
  !! etc for `y` and `z` coordinates.
  type, public, extends(point_zone_t) :: box_point_zone_t
     real(kind=rp) :: xmin
     real(kind=rp) :: xmax
     real(kind=rp) :: ymin
     real(kind=rp) :: ymax
     real(kind=rp) :: zmin
     real(kind=rp) :: zmax
   contains
     !> Constructor from json object file.
     procedure, pass(this) :: init => box_point_zone_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => box_point_zone_free
     !> Defines the criterion of selection of a GLL point in the box point zone.
     procedure, pass(this) :: criterion => box_point_zone_criterion
  end type box_point_zone_t

contains

  !> Constructor from json object file.
  !! @param json Json object file.
  !! @param size Size with which to initialize the stack
  subroutine box_point_zone_init_from_json(this, json, size)
    class(box_point_zone_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(in) :: size

    character(len=:), allocatable :: str_read
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: xmin, xmax, ymin, ymax, zmin, zmax

    call json_get(json, "x_bounds", values)
    xmin = values(1)
    xmax = values(2)
    call json_get(json, "y_bounds", values)
    ymin = values(1)
    ymax = values(2)
    call json_get(json, "z_bounds", values)
    zmin = values(1)
    zmax = values(2)
    call json_get(json, "name", str_read)

    call box_point_zone_init_common(this, size, trim(str_read), xmin, &
         xmax, ymin, ymax, zmin, zmax)

  end subroutine box_point_zone_init_from_json

  !> Initializes a box point zone from its coordinates.
  !! @param size Size of the scratch stack.
  !! @param name Name of the box point zone.
  !! @param xmin Lower x-bound of the box coordinates.
  !! @param xmax Upper x-bound of the box coordinates.
  !! @param ymin Lower y-bound of the box coordinates.
  !! @param ymax Upper y-bound of the box coordinates.
  !! @param zmin Lower z-bound of the box coordinates.
  !! @param zmax Upper z-bound of the box coordinates.
  subroutine box_point_zone_init_common(this, size, name, xmin, xmax, ymin, ymax, &
       zmin, zmax)
    class(box_point_zone_t), intent(inout) :: this
    integer, intent(in), optional :: size
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: xmin
    real(kind=rp), intent(in) :: xmax
    real(kind=rp), intent(in) :: ymin
    real(kind=rp), intent(in) :: ymax
    real(kind=rp), intent(in) :: zmin
    real(kind=rp), intent(in) :: zmax

    call this%init_base(size, name)

    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax
    this%zmin = zmin
    this%zmax = zmax

  end subroutine box_point_zone_init_common

  !> Destructor.
  subroutine box_point_zone_free(this)
    class(box_point_zone_t), intent(inout) :: this

    this%xmin = 0.0_rp
    this%xmax = 0.0_rp
    this%ymin = 0.0_rp
    this%ymax = 0.0_rp
    this%zmin = 0.0_rp
    this%zmax = 0.0_rp

    call this%free_base()

  end subroutine box_point_zone_free

  !> Defines the criterion of selection of a GLL point in the box point zone.
  !! In the case of a box point zone, an `x,y,z` GLL point is considered as 
  !! being inside the zone if:
  !! \f{eqnarray*}{
  !!    x_{min} \le x \le x_{max} \\
  !!    y_{min} \le y \le y_{max} \\
  !!    z_{min} \le z \le z_{max} \\
  !! \f}
  !! @param x x-coordinate of the GLL point.
  !! @param y y-coordinate of the GLL point.
  !! @param z z-coordinate of the GLL point.
  !! @param j 1st nonlinear index of the GLL point.
  !! @param k 2nd nonlinear index of the GLL point.
  !! @param l 3rd nonlinear index of the GLL point.
  !! @param e element index of the GLL point.
  pure function box_point_zone_criterion(this, x, y, z, j, k, l, e) result(is_inside)
    class(box_point_zone_t), intent(in) :: this
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    logical :: is_inside
    logical :: in_x, in_y, in_z

    ! inside x if xmin <= x <= xmax
    in_x = ( (x .gt. this%xmin .and. x .lt. this%xmax) .or. &
             (abscmp(x, this%xmin) .or. abscmp(x, this%xmax)))

    ! inside y if ymin <= y <= ymax
    in_y = ( (y .gt. this%ymin .and. y .lt. this%ymax) .or. &
             (abscmp(y, this%ymin) .or. abscmp(y, this%ymax)))

    ! inside z if zmin <= z <= zmax
    in_z = ( (z .gt. this%zmin .and. z .lt. this%zmax) .or. &
             (abscmp(z, this%zmin) .or. abscmp(z, this%zmax)))

    is_inside = in_x .and. in_y .and. in_z
  end function box_point_zone_criterion

end module box_point_zone
