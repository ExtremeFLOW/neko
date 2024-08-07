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

!> \submodule distance_elements
!! This submodule contains the operators that are related to distance
!! calculations for mesh elements.
!!
!! @note Please note that this module only contains the implementations for the
!! Geometric Operator module.
submodule (geometric_operators) distance_elements
  use math, only: cross
  use utils, only: neko_error
  implicit none

contains

  !> @brief Distance function for a point
  !! @details This routine computes the distance between the query point and the
  !! given point.
  !!
  !! @param p Query point
  !! @param point Point
  !! @return Unsigned distance value
  real(kind=dp) module function distance_point_t(p, point)
    real(kind=dp), dimension(3), intent(in) :: p
    type(point_t), intent(in) :: point

    distance_point_t = norm2(p - point%x)

  end function distance_point_t

  !> @brief Distance function for a triangle.
  !! @details This routine computes the distance to a given triangle.
  !! We compute the barycentric coordinate of the point projected onto the
  !! triangle. If the projection is inside the triangle, the distance is the
  !! perpendicular distance to the plane. If the projection is outside the
  !! triangle, the distance is the distance to the nearest edge or vertex.
  !!
  !! @param p Query point
  !! @param triangle Triangle
  !! @return Unsigned distance value
  real(kind=dp) module function distance_triangle(p, triangle)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tri_t), intent(in) :: triangle

    real(kind=dp), dimension(3) :: v1, v2, v3
    real(kind=dp), dimension(3) :: normal
    real(kind=dp) :: normal_length
    real(kind=dp) :: b1, b2, b3

    real(kind=dp), dimension(3) :: projection
    real(kind=dp) :: tol = 1.0e-10_dp

    ! Get vertices and the normal vector
    v1 = triangle%pts(1)%p%x
    v2 = triangle%pts(2)%p%x
    v3 = triangle%pts(3)%p%x

    normal = cross(v2 - v1, v3 - v1)
    normal_length = norm2(normal)

    if (normal_length .lt. tol) then
       call neko_error('Triangle is degenerate')
    end if

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! triangular prism, off along an edge or by a vertex.
    projection = p - normal * dot_product(p - v1, normal) / normal_length**2
    b = barycentric_coordinate(projection, triangle)

    if (b(1) .le. tol) then
       distance_triangle = distance_line_segment(p, v3, v2)
    else if (b(2) .le. tol) then
       distance_triangle = distance_line_segment(p, v1, v3)
    else if (b(3) .le. tol) then
       distance_triangle = distance_line_segment(p, v2, v1)
    else
       distance_triangle = abs(distance_plane(p, v1, normal))
    end if

  end function distance_triangle

  !> @brief Distance function for a tetrahedron.
  !! @details This routine computes the distance to a given tetrahedron.
  !! We compute the barycentric coordinate of the point in relation to the
  !! tetrahedron. If any of the barycentric coordinates are less than or equal
  !! to zero, the point is off along a face, edge or by a vertex. Otherwise, the
  !! distance is the perpendicular distance to the plane.
  !!
  !! @param p Query point
  !! @param tetrahedron Tetrahedron
  !! @return Unsigned distance value
  real(kind=dp) module function distance_tetrahedron(p, tetrahedron)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tet_t), intent(in) :: tetrahedron

    type(point_t) :: v1, v2, v3, v4

    ! Variables for barycentric interpolation
    real(kind=dp), dimension(4) :: barycoord

    type(tri_t) :: face
    integer :: dummy_id

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! tetrahedron, off along a face, edge or by a vertex.
    barycoord = barycentric_coordinate(p, tetrahedron)

    ! Get vertices and the normal vector
    v1 = tetrahedron%pts(1)%p
    v2 = tetrahedron%pts(2)%p
    v3 = tetrahedron%pts(3)%p
    v4 = tetrahedron%pts(4)%p

    ! Check if the point is inside the tetrahedron
    select case (minloc(barycoord, 1))
      case (1)
       call face%init(dummy_id, v2, v3, v4)
      case (2)
       call face%init(dummy_id, v1, v3, v4)
      case (3)
       call face%init(dummy_id, v1, v2, v4)
      case (4)
       call face%init(dummy_id, v1, v2, v3)
    end select

    distance_tetrahedron = distance_triangle(p, face)
    if (all(barycoord .ge. 0.0_dp)) distance_tetrahedron = -distance_tetrahedron

  end function distance_tetrahedron

end submodule distance_elements
