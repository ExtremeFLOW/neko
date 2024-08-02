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

!> \submodule distance_primitives
!! This submodule contains the operators that are related to distance
!! calculations for the primitive geometric objects.
!!
!! @note Please note that this module only contains the implementations for the
!! Geometric Operator module.
submodule (geometric_operators) distance_primitives
  use num_types, only: dp
  implicit none

contains

  !> Distance to a point
  !! @param p Query point
  !! @param point Point
  !! @return Unsigned distance value
  module function distance_point_real(p, point)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: point
    real(kind=dp) :: distance_point_real

    distance_point_real = norm2(p - point)
  end function distance_point_real

  !> Distance to the infinite line defined by a point and a direction.
  !! @param p Query point
  !! @param point Point on the line
  !! @param direction Direction of the line
  !! @return Distance to the line
  module function distance_line(p, point, direction)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: point
    real(kind=dp), dimension(3), intent(in) :: direction
    real(kind=dp) :: distance_line

    real(kind=dp) :: t, projection(3)

    t = dot_product(p - point, direction) / norm2(direction)
    projection = point + t * direction

    distance_line = norm2(projection - p)
  end function distance_line

  !> Distance to the infinite ray starting at a point and going in a
  !! direction.
  !! @param p Query point
  !! @param point Point on the ray
  !! @param direction Direction of the ray
  !! @return Distance to the ray
  module function distance_line_ray(p, point, direction)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: point
    real(kind=dp), dimension(3), intent(in) :: direction
    real(kind=dp) :: distance_line_ray

    real(kind=dp) :: t, projection(3)

    ! Project the point onto the ray
    t = dot_product(p - point, direction) / norm2(direction)
    t = max(t, 0.0_dp)
    projection = point + t * direction

    distance_line_ray = norm2(projection - p)
  end function distance_line_ray

  !> Distance to a line segment defined by two points.
  !! @param p Query point
  !! @param point_0 First point of the line segment
  !! @param point_1 Second point of the line segment
  !! @return Distance to the line segment
  module function distance_line_segment(p, point_0, point_1)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: point_0
    real(kind=dp), dimension(3), intent(in) :: point_1
    real(kind=dp) :: distance_line_segment

    real(kind=dp) :: t, direction(3), projection(3)

    ! Compute direction of the segment
    direction = point_1 - point_0

    ! Project the point onto the line segment
    t = dot_product(p - point_0, direction) / norm2(direction)
    t = min(1.0_dp, max(0.0_dp, t))
    projection = point_0 + t * direction

    distance_line_segment = norm2(projection - p)
  end function distance_line_segment

  !> Distance to a plane defined by a point and a normal.
  !! @param p Query point
  !! @param point Point on the plane
  !! @param normal Normal of the plane
  !! @return Signed distance to the plane
  module function distance_plane(p, point, normal)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: point
    real(kind=dp), dimension(3), intent(in) :: normal
    real(kind=dp) :: distance_plane

    real(kind=dp) :: face_distance, projection(3)

    face_distance = dot_product(p - point, normal) / norm2(normal)
    projection = p - normal * face_distance / norm2(normal)

    distance_plane = sign(norm2(projection - p), face_distance)
  end function distance_plane

  !> Distance to a sphere defined by a center and a radius.
  !! @param p Query point
  !! @param sphere_center Center of the sphere
  !! @param sphere_radius Radius of the sphere
  !! @return Signed distance to the sphere
  module function distance_sphere(p, sphere_center, sphere_radius)
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), dimension(3), intent(in) :: sphere_center
    real(kind=dp), intent(in) :: sphere_radius
    real(kind=dp) :: distance_sphere

    distance_sphere = norm2(p - sphere_center) - sphere_radius
  end function distance_sphere

end submodule distance_primitives
