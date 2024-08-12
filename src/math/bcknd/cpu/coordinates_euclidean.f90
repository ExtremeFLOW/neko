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

!> \submodule coordinates euclidean
!! This submodule contains the operators that are related to computing euclidean
!! coordinates for mesh objects.
!!
!! @note Please note that this module only contains the implementations for the
!! Geometric Operator module.
submodule (geometric_operators) coordinates_euclidean
  use num_types, only: dp
  use utils, only: neko_error
  use math, only: cross
  use tri, only: tri_t
  use tet, only: tet_t
  use quad, only: quad_t
  use hex, only: hex_t
  implicit none

contains

  !> @brief Euclidean coordinates for a point in a triangle
  !! @details This routine computes the euclidean coordinates of a local
  !! barycentric point with respect to a triangle.
  !!
  !! @param bary Barycentric coordinates of the point
  !! @param triangle Triangle
  !! @return Euclidean coordinate
  module function euclidean_coordinate_triangle(bary, triangle) &
       result(point)
    real(kind=dp), dimension(3), intent(in) :: bary
    type(tri_t), intent(in) :: triangle
    real(kind=dp), dimension(3) :: point

    associate (p1 => triangle%pts(1)%p, &
               p2 => triangle%pts(2)%p, &
               p3 => triangle%pts(3)%p)
      point = bary(1) * p1%x + &
           bary(2) * p2%x + &
           bary(3) * p3%x
    end associate

  end function euclidean_coordinate_triangle

  !> @brief Euclidean coordinates for a point in a tetrahedron
  !! @details This routine computes the euclidean coordinates of a local
  !! barycentric point with respect to a tetrahedron.
  !!
  !! @param bary Barycentric coordinates of the point
  !! @param tetrahedron Tetrahedron
  !! @return Euclidean coordinates
  module function euclidean_coordinate_tetrahedron(bary, tetrahedron) &
       result(point)
    real(kind=dp), dimension(4), intent(in) :: bary
    type(tet_t), intent(in) :: tetrahedron
    real(kind=dp), dimension(3) :: point

    associate (p1 => tetrahedron%pts(1)%p, &
               p2 => tetrahedron%pts(2)%p, &
               p3 => tetrahedron%pts(3)%p, &
               p4 => tetrahedron%pts(4)%p)
      point = bary(1) * p1%x + &
           bary(2) * p2%x + &
           bary(3) * p3%x + &
           bary(4) * p4%x
    end associate

  end function euclidean_coordinate_tetrahedron

  !> @brief Euclidean coordinates for a point in a quadrilateral
  !! @details This routine computes the euclidean coordinates of a local
  !! bilinear point with respect to a quadrilateral.
  !!
  !! @param lin Bilinear coordinates of the point
  !! @param quadrilateral Quadrilateral
  !! @return Euclidean coordinates
  module function euclidean_coordinate_quadrilateral(lin, quadrilateral) &
       result(point)
    real(kind=dp), dimension(2), intent(in) :: lin
    type(quad_t), intent(in) :: quadrilateral
    real(kind=dp), dimension(3) :: point

    associate (p1 => quadrilateral%pts(1)%p, &
               p2 => quadrilateral%pts(2)%p, &
               p3 => quadrilateral%pts(3)%p, &
               p4 => quadrilateral%pts(4)%p)
      point = &
           (1.0_dp - lin(1)) * (1.0_dp - lin(2)) * p1%x + &
           (0.0_dp + lin(1)) * (1.0_dp - lin(2)) * p2%x + &
           (1.0_dp - lin(1)) * (0.0_dp + lin(2)) * p3%x + &
           (0.0_dp + lin(1)) * (0.0_dp + lin(2)) * p4%x
    end associate

  end function euclidean_coordinate_quadrilateral

  !> @brief Euclidean coordinates for a point in a hexahedron
  !! @details This routine computes the euclidean coordinates of a local
  !! trilinear point with respect to a hexahedron.
  !!
  !! @param lin Trilinear coordinates of the point
  !! @param hexahedron Hexahedron
  !! @return Euclidean coordinates
  module function euclidean_coordinate_hexahedron(lin, hexahedron) &
       result(point)
    real(kind=dp), dimension(3), intent(in) :: lin
    type(hex_t), intent(in) :: hexahedron
    real(kind=dp), dimension(3) :: point

    associate (p1 => hexahedron%pts(1)%p, &
               p2 => hexahedron%pts(2)%p, &
               p3 => hexahedron%pts(3)%p, &
               p4 => hexahedron%pts(4)%p, &
               p5 => hexahedron%pts(5)%p, &
               p6 => hexahedron%pts(6)%p, &
               p7 => hexahedron%pts(7)%p, &
               p8 => hexahedron%pts(8)%p)
      point = &
           (1.0_dp - lin(1)) * (1.0_dp - lin(2)) * (1.0_dp - lin(3)) * p1%x + &
           (0.0_dp + lin(1)) * (1.0_dp - lin(2)) * (1.0_dp - lin(3)) * p2%x + &
           (1.0_dp - lin(1)) * (0.0_dp + lin(2)) * (1.0_dp - lin(3)) * p3%x + &
           (0.0_dp + lin(1)) * (0.0_dp + lin(2)) * (1.0_dp - lin(3)) * p4%x + &
           (1.0_dp - lin(1)) * (1.0_dp - lin(2)) * (0.0_dp + lin(3)) * p5%x + &
           (0.0_dp + lin(1)) * (1.0_dp - lin(2)) * (0.0_dp + lin(3)) * p6%x + &
           (1.0_dp - lin(1)) * (0.0_dp + lin(2)) * (0.0_dp + lin(3)) * p7%x + &
           (0.0_dp + lin(1)) * (0.0_dp + lin(2)) * (0.0_dp + lin(3)) * p8%x
    end associate

  end function euclidean_coordinate_hexahedron

end submodule coordinates_euclidean
