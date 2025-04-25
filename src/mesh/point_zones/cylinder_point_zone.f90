! Copyright (c) 2024, The Neko Authors
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
!> Implements a cylinder geometry subset.
module cylinder_point_zone
  use point_zone, only: point_zone_t
  use num_types, only: rp
  use json_utils, only: json_get, json_get_or_default
  use json_module, only: json_file
  use utils, only: neko_error
  implicit none
  private

  !> A cylindrical point zone.
  !! @details As defined here, a cylinder is described by its two end points and
  !! its radius, specified in the json file
  !! as e.g. `"start": [<x0>, <y0>, <z0>]",
  !! "start": [<x1>, <y1>, <z1>]", "radius": <r>`.
  type, public, extends(point_zone_t) :: cylinder_point_zone_t
     real(kind=rp), dimension(3) :: p0
     real(kind=rp), dimension(3) :: p1
     real(kind=rp) :: radius
   contains
     !> Constructor from json object file.
     procedure, pass(this) :: init => cylinder_point_zone_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => cylinder_point_zone_free
     !> Defines the criterion of selection of a GLL point in the sphere point
     !! zone.
     procedure, pass(this) :: criterion => cylinder_point_zone_criterion
  end type cylinder_point_zone_t

contains

  !> Constructor from json object file.
  !! @param json Json object file.
  !! @param size Size with which to initialize the stack
  subroutine cylinder_point_zone_init_from_json(this, json, size)
    class(cylinder_point_zone_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(in) :: size

    character(len=:), allocatable :: name
    real(kind=rp), dimension(:), allocatable :: p0, p1
    real(kind=rp) :: radius
    logical :: invert

    call json_get(json, "name", name)
    call json_get(json, "start", p0)
    call json_get(json, "end", p1)
    call json_get(json, "radius", radius)

    ! Needed to use `shape` because of input name.
    if (all(shape(p0) .ne. (/3/))) then
       call neko_error("Cylinder point zone: invalid start point")
    end if

    if (all(shape(p1) .ne. (/3/))) then
       call neko_error("Cylinder point zone: invalid end point")
    end if

    if (radius .lt. 0.0_rp) then
       call neko_error("Cylinder point zone: invalid radius")
    end if

    call json_get_or_default(json, "invert", invert, .false.)

    call cylinder_point_zone_init_common(this, size, trim(name), invert, &
         p0, p1, radius)

  end subroutine cylinder_point_zone_init_from_json

  !> Initializes a cylinder point zone from its endpoint coordinates and radius.
  !! @param size Size of the scratch stack.
  !! @param name Name of the cylinder point zone.
  !! @param p0 Coordinates of the first endpoint.
  !! @param p1 Coordinates of the second endpoint.
  !! @param radius Sphere radius.
  subroutine cylinder_point_zone_init_common(this, size, name, invert, &
       p0, p1, radius)
    class(cylinder_point_zone_t), intent(inout) :: this
    integer, intent(in), optional :: size
    character(len=*), intent(in) :: name
    logical, intent(in) :: invert
    real(kind=rp), intent(in), dimension(3) :: p0
    real(kind=rp), intent(in), dimension(3) :: p1
    real(kind=rp), intent(in) :: radius

    call this%init_base(size, name, invert)

    this%p0 = p0
    this%p1 = p1
    this%radius = radius

  end subroutine cylinder_point_zone_init_common

  !> Destructor.
  subroutine cylinder_point_zone_free(this)
    class(cylinder_point_zone_t), intent(inout) :: this

    call this%free_base()

    this%p0 = 0.0_rp
    this%p1 = 0.0_rp
    this%radius = 0.0_rp

  end subroutine cylinder_point_zone_free

  !> @brief Defines the criterion of selection of a GLL point in the cylinder
  !! point zone.
  !! @details A GLL point of coordinates \f$ \vec{X} = (x, y, z) \f$ is
  !! considered as being inside the cylinder defined by endpoints
  !! \f$ \vec{p_0} \f$ and \f$ \vec{p_1} \f$ and radius \f$ r \f$ if it
  !! satisfies the following conditions:
  !! \f{eqnarray*}{
  !!   ||\vec{X} - \vec{X_0}|| &\le& r\\
  !!   0 &\le& t \le 1
  !! \f}
  !! where,
  !! \f{eqnarray*}{
  !!   t &=& (\vec{X} - \vec{p_0}) \cdot (\vec{p_1} - \vec{p_0})
  !!         / ||\vec{p_1} - \vec{p_0}|| \\
  !! \vec{X_0} &=& \vec{p_0} + t \cdot (\vec{p_1} - \vec{p_0})
  !! \f}
  !! @param x x-coordinate of the GLL point.
  !! @param y y-coordinate of the GLL point.
  !! @param z z-coordinate of the GLL point.
  !! @param j 1st nonlinear index of the GLL point.
  !! @param k 2nd nonlinear index of the GLL point.
  !! @param l 3rd nonlinear index of the GLL point.
  !! @param e element index of the GLL point.
  pure function cylinder_point_zone_criterion(this, x, y, z, j, k, l, e) &
       result(is_inside)
    class(cylinder_point_zone_t), intent(in) :: this
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    logical :: is_inside

    real(kind=rp), dimension(3) :: p
    real(kind=rp), dimension(3) :: centerline
    real(kind=rp), dimension(3) :: vec_p
    real(kind=rp) :: t
    real(kind=rp), dimension(3) :: projection
    real(kind=rp) :: distance

    p = [x, y, z]

    centerline = this%p1 - this%p0
    vec_p = p - this%p0
    t = dot_product(vec_p, centerline) / dot_product(centerline, centerline)

    projection = this%p0 + t * centerline
    distance = norm2(projection - p)

    is_inside = t >= 0.0_rp .and. t <= 1.0_rp .and. distance <= this%radius

  end function cylinder_point_zone_criterion

end module cylinder_point_zone
