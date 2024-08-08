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

!> \submodule linear coordinates
!! This submodule contains the operators that are related to linear
!! coordinates for mesh elements.
!!
!! @note Please note that this module only contains the implementations for the
!! Geometric Operator module.
submodule (geometric_operators) linear_coordinates
  use utils, only: neko_error
  use math, only: cross
  implicit none

contains

  !> @brief Barycentric coordinates for a point in a triangle
  !! @details This routine computes the barycentric coordinates of a point with
  !! respect to a triangle.
  !!
  !! @note The point is outside the triangle if the barycentric coordinates are
  !! negative, on an edge if one of the barycentric coordinates is zero, and on
  !! a vertex if two barycentric coordinates are zero.
  !!
  !! @param p Query point
  !! @param triangle Triangle
  !! @return Barycentric coordinates
  module function barycentric_coordinate_triangle(p, triangle) &
       result(barycoord)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tri_t), intent(in) :: triangle
    real(kind=dp), dimension(3) :: barycoord

    real(kind=dp), dimension(3) :: v1, v2, v3

    ! Variables for barycentric interpolation
    real(kind=dp), dimension(3, 3) :: barymatrix
    real(kind=dp), dimension(3) :: normal
    real(kind=dp) :: normal_length
    real(kind=dp), parameter :: tol = 1.0e-12_dp

    ! Get vertices and the normal vector
    v1 = triangle%pts(1)%p%x
    v2 = triangle%pts(2)%p%x
    v3 = triangle%pts(3)%p%x

    normal = cross(v2 - v1, v3 - v1)
    normal_length = norm2(normal)

    if (normal_length .lt. tol) then
       call neko_error('Triangle is degenerate')
    end if

    if (distance_plane(p, v1, normal) .gt. tol) then
       call neko_error('Point is not in the plane of the triangle')
    end if

    barycoord = [ &
         & dot_product(normal, cross(v3 - v2, p - v2)) / normal_length**2, &
         & dot_product(normal, cross(v1 - v3, p - v3)) / normal_length**2, &
         & dot_product(normal, cross(v2 - v1, p - v1)) / normal_length**2 &
         & ]

  end function barycentric_coordinate_triangle


  !> @brief Barycentric coordinates for a point in a tetrahedron
  !! @details This routine computes the barycentric coordinates of a point with
  !! respect to a tetrahedron.
  !!
  !! @note The point is outside the tetrahedron if the barycentric coordinates
  !! are negative, on a face if one of the barycentric coordinates is zero, on
  !! an edge if two of the barycentric coordinates are zero, and on a vertex if
  !! all three barycentric coordinates are zero.
  !!
  !! @param p Query point
  !! @param tetrahedron Tetrahedron
  !! @return Barycentric coordinates
  module function barycentric_coordinate_tetrahedron(p, tetrahedron) &
       result(barycoord)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tet_t), intent(in) :: tetrahedron
    real(kind=dp), dimension(4) :: barycoord

    real(kind=dp), dimension(3) :: v1, v2, v3, v4

    ! Variables for barycentric interpolation
    real(kind=dp), dimension(4, 4) :: barymatrix

    ! Lapack variables
    integer :: info
    integer, dimension(4) :: pivot_table

    ! Get vertices and the normal vector
    v1 = tetrahedron%pts(1)%p%x
    v2 = tetrahedron%pts(2)%p%x
    v3 = tetrahedron%pts(3)%p%x
    v4 = tetrahedron%pts(4)%p%x

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! tetrahedron, off along a face, edge or by a vertex.

    barycoord = reshape([1.0_dp, p], shape(barycoord))
    barymatrix = transpose(reshape( &
         & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         & v1(1), v2(1), v3(1), v4(1), &
         & v1(2), v2(2), v3(2), v4(2), &
         & v1(3), v2(3), v3(3), v4(3)], &
         & shape(barymatrix)))

    ! Solve the system of linear equations
    call dgesv(4, 1, barymatrix, 4, pivot_table, barycoord, 4, info )

    ! Check for the exact singularity.
    if (info .gt. 0) then
       call neko_error('Tetrahedron is degenerate')
    end if

  end function barycentric_coordinate_tetrahedron

end submodule linear_coordinates
