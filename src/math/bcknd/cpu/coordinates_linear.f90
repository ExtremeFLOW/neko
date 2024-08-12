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

!> \submodule coordinates linear
!! This submodule contains the operators that are related to linear
!! coordinates for mesh elements.
!!
!! @note Please note that this module only contains the implementations for the
!! Geometric Operator module.
submodule (geometric_operators) coordinates_linear
  use num_types, only: dp
  use utils, only: neko_error
  use math, only: cross
  use tri, only: tri_t
  use tet, only: tet_t
  use quad, only: quad_t
  use hex, only: hex_t
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
       result(barycentric)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tri_t), intent(in) :: triangle
    real(kind=dp), dimension(3) :: barycentric

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

    barycentric = [ &
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
       result(barycentric)
    real(kind=dp), dimension(3), intent(in) :: p
    type(tet_t), intent(in) :: tetrahedron
    real(kind=dp), dimension(4) :: barycentric

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

    barycentric = reshape([1.0_dp, p], shape(barycentric))
    barymatrix = transpose(reshape( &
         & [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         & v1(1), v2(1), v3(1), v4(1), &
         & v1(2), v2(2), v3(2), v4(2), &
         & v1(3), v2(3), v3(3), v4(3)], &
         & shape(barymatrix)))

    ! Solve the system of linear equations
    call dgesv(4, 1, barymatrix, 4, pivot_table, barycentric, 4, info )

    ! Check for the exact singularity.
    if (info .gt. 0) then
       call neko_error('Tetrahedron is degenerate')
    end if

  end function barycentric_coordinate_tetrahedron

  !> @brief Linear coordinates for a point in a quadrilateral
  !! @details This routine computes the bi-linear coordinates of a point with
  !! respect to a quadrilateral.
  !!
  !! @param p Query point
  !! @param quadrilateral Quadrilateral
  !! @return Bi-linear coordinates
  module function bilinear_coordinate_quadrilateral(p, quadrilateral) &
       result(bilinear)
    real(kind=dp), dimension(3), intent(in) :: p
    type(quad_t), intent(in) :: quadrilateral
    real(kind=dp), dimension(2) :: bilinear

    real(kind=dp), dimension(3) :: v1, v2, v3, v4

    ! Variables for barycentric interpolation
    real(kind=dp), dimension(4, 2) :: barymatrix

    ! Lapack variables
    integer :: info
    integer, dimension(4) :: pivot_table

    ! Get vertices and the normal vector
    v1 = quadrilateral%pts(1)%p%x
    v2 = quadrilateral%pts(2)%p%x
    v3 = quadrilateral%pts(3)%p%x
    v4 = quadrilateral%pts(4)%p%x

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! quadrilateral, off along a face, edge or by a vertex.

    bilinear = 0.0_dp

  end function bilinear_coordinate_quadrilateral

  !> @brief Linear coordinates for a point in a hexahedron
  !! @details This routine computes the tri-linear coordinates of a point with
  !! respect to a hexahedron.
  !!
  !! @param p Query point
  !! @param hexahedron Hexahedron
  !! @return Tri-linear coordinates
  module function trilinear_coordinate_hexahedron(p, hexahedron) &
       result(trilinear)
    real(kind=dp), dimension(3), intent(in) :: p
    type(hex_t), intent(in) :: hexahedron
    real(kind=dp), dimension(3) :: trilinear

    real(kind=dp), dimension(3) :: p_est, residual

    ! Variables for barycentric interpolation
    real(kind=dp), dimension(3, 3) :: jacobian

    ! Parameters for the solver
    integer, parameter :: max_iter = 100
    real(kind=dp), parameter :: tol = 1.0e-12_dp
    integer :: iter

    ! Lapack variables
    integer :: info
    integer, dimension(3) :: pivot_table

    ! Initialize the estimate
    trilinear = 0.5_dp
    p_est = euclidean_coordinate(trilinear, hexahedron)
    residual = p_est - p

    do iter = 1, max_iter

       ! Compute the residual and check for convergence
       p_est = euclidean_coordinate(trilinear, hexahedron)
       residual = euclidean_coordinate(trilinear, hexahedron) - p

       if (norm2(residual) .lt. tol) exit

       ! Compute the Jacobian
       jacobian = jacobian_coordinate(trilinear, hexahedron)

       write (*,*) "trilinear: ", trilinear
       write (*,*) "p_est:     ", p_est
       write (*,*) "residual:  ", residual
       write (*,*) "Jacobian:  "
       write (*,'(F6.2, F6.2, F6.2)') jacobian(1, :)
       write (*,'(F6.2, F6.2, F6.2)') jacobian(2, :)
       write (*,'(F6.2, F6.2, F6.2)') jacobian(3, :)

       ! Solve the system of linear equations
       call dgesv(3, 1, jacobian, 3, pivot_table, residual, 3, info)

       ! Check for the exact singularity.
       !  if (info .gt. 0) call neko_error('Hexahedron is degenerate')

       ! Update the trilinear coordinates
       trilinear = trilinear - residual
    end do

    ! Check for convergence
    ! if (norm2(residual) .gt. tol) then
    !    call neko_error('Hexahedron is degenerate')
    ! end if

  end function trilinear_coordinate_hexahedron

end submodule coordinates_linear

! /** Hexahedral element.
! *
! * These hexahedral are static hexahedral. They know only their corners
! * and use those to compute local coordinates based on Linear Shape
! * Functions.
! *
! * These are an implementation of the Hex8 elements in Exodus.
! * https://sandialabs.github.io/seacas-docs/html/md_include_exodus_element_types.html#hex
! */
! class Hexahedron : public Polyhedron {

!    /** Compute the local coordinates of a point inside the cell.
!     *
!     * The local coordinates are computed using the trilinear map.
!     *
!     * @param p: The point to be tested.
!     * @return std::vector<double>: The local coordinates of the point.
!     */
!    std::vector<double> compute_local(CGLA::Vec3d p) const override {
!        const CGLA::Vec3d q = trilinear_inverse(p);
!        return S(q);
!    }

!  private:
!    // ------------------------------------------------------------------ //
!    // Element specific functions

!    /** Computes the Standard hex coordinates.
!     *
!     * https://stackoverflow.com/a/18332009
!     *
!     * Computes the inverse of the trilinear map from [-1,1]^3 to the box
!     * defined by corners c0,...,c7, where the corners are ordered
!     * consistent with the exodus format. Uses Gauss-Newton method.
!     *
!     * We compute the standardized position q given a point p in the
!     * generalized hex. "q" can then be used for trilinear interpolation.
!     *
!     * @todo We have only implemented a simple newton stepping method.
!     * Probably should update it to run as a Gauss-Newton as described
!     * above.
!     */
!    CGLA::Vec3d trilinear_inverse(CGLA::Vec3d p) const {

!        const double tol  = 1e-9 * this->diagonal().length();
!        const int    iter = 1000;

!        // Convert to eigen3 notation
!        const Eigen::Vector3d p_e(p[0], p[1], p[2]);

!        // Setup new computation
!        Eigen::Vector3d q_e(0.0, 0.0, 0.0);
!        Eigen::Vector3d residual = estimate(q_e, vertices) - p_e;
!        for (int k = 0; k < iter; k++) {

!            if (residual.norm() < tol) break;

!            const Eigen::Matrix3d J_e = jacobian(q_e, vertices);
!            q_e -= J_e.fullPivLu().solve(residual);
!            residual = estimate(q_e, vertices) - p_e;
!        }

!        // Return invalid point if we are too far away
!        if (tol < residual.norm()) return CGLA::Vec3d(-1e16);

!        const CGLA::Vec3d q(q_e[0], q_e[1], q_e[2]);
!        return q;
!    }

!    /** Evaluate the linear shape functions.
!     *
!     * @param q: The local coordinates to be evaluated.
!     * @return std::vector<double>: The shape functions evaluated at q.
!     */
!    template <typename T> std::vector<double> S(const T& q) const {
!        return {
!            0.125 * (1.0 - q[0]) * (1.0 - q[1]) * (1.0 - q[2]), //
!            0.125 * (1.0 + q[0]) * (1.0 - q[1]) * (1.0 - q[2]), //
!            0.125 * (1.0 + q[0]) * (1.0 + q[1]) * (1.0 - q[2]), //
!            0.125 * (1.0 - q[0]) * (1.0 + q[1]) * (1.0 - q[2]), //
!            0.125 * (1.0 - q[0]) * (1.0 - q[1]) * (1.0 + q[2]), //
!            0.125 * (1.0 + q[0]) * (1.0 - q[1]) * (1.0 + q[2]), //
!            0.125 * (1.0 + q[0]) * (1.0 + q[1]) * (1.0 + q[2]), //
!            0.125 * (1.0 - q[0]) * (1.0 + q[1]) * (1.0 + q[2])  //
!        };
!    }

!    /** Evaluate the derivative of the linear shape functions.
!     *
!     * @param q: The local coordinates to be evaluated.
!     * @return std::vector<std::vector<double>>: The derivatives of the
!     * shape functions evaluated at q.
!     */
!    template <typename T>
!    std::vector<std::vector<double>> dS(const T& q) const {

!        const std::vector<double> dSdu = {
!            -0.125 * (1.0 - q[1]) * (1.0 - q[2]), //
!            +0.125 * (1.0 - q[1]) * (1.0 - q[2]), //
!            +0.125 * (1.0 + q[1]) * (1.0 - q[2]), //
!            -0.125 * (1.0 + q[1]) * (1.0 - q[2]), //
!            -0.125 * (1.0 - q[1]) * (1.0 + q[2]), //
!            +0.125 * (1.0 - q[1]) * (1.0 + q[2]), //
!            +0.125 * (1.0 + q[1]) * (1.0 + q[2]), //
!            -0.125 * (1.0 + q[1]) * (1.0 + q[2])  //
!        };
!        const std::vector<double> dSdv = {
!            -0.125 * (1.0 - q[0]) * (1.0 - q[2]), //
!            -0.125 * (1.0 + q[0]) * (1.0 - q[2]), //
!            +0.125 * (1.0 + q[0]) * (1.0 - q[2]), //
!            +0.125 * (1.0 - q[0]) * (1.0 - q[2]), //
!            -0.125 * (1.0 - q[0]) * (1.0 + q[2]), //
!            -0.125 * (1.0 + q[0]) * (1.0 + q[2]), //
!            +0.125 * (1.0 + q[0]) * (1.0 + q[2]), //
!            +0.125 * (1.0 - q[0]) * (1.0 + q[2])  //
!        };
!        const std::vector<double> dSdw = {
!            -0.125 * (1.0 - q[0]) * (1.0 - q[1]), //
!            -0.125 * (1.0 + q[0]) * (1.0 - q[1]), //
!            -0.125 * (1.0 + q[0]) * (1.0 + q[1]), //
!            -0.125 * (1.0 - q[0]) * (1.0 + q[1]), //
!            +0.125 * (1.0 - q[0]) * (1.0 - q[1]), //
!            +0.125 * (1.0 + q[0]) * (1.0 - q[1]), //
!            +0.125 * (1.0 + q[0]) * (1.0 + q[1]), //
!            +0.125 * (1.0 - q[0]) * (1.0 + q[1])  //
!        };

!        return {dSdu, dSdv, dSdw};
!    }

!    /** Estimate the position of a point given local coordinates.
!     *
!     * @param q: The local coordinates to be evaluated.
!     * @param C: The vertex positions of the current polyhedron.
!     * @return CGLA::Vec3d: The position of the point.
!     */
!    template <typename T>
!    T estimate(const T& q, const std::vector<CGLA::Vec3d>& C) const {

!        const std::vector<double> N   = S(q);
!        const CGLA::Vec3d         xyz = std::inner_product(
!            begin(N), end(N), begin(C), CGLA::Vec3d{0.0});

!        return T(xyz[0], xyz[1], xyz[2]);
!    }

!    /** Estimate the Jacobian of a point given local coordinates.
!     *
!     * @param q: The local coordinates to be evaluated.
!     * @param C: The vertex positions of the current polyhedron.
!     * @return Eigen::Matrix3d: The Jacobian of the point.
!     */
!    template <typename T>
!    Eigen::Matrix3d
!        jacobian(const T& q, const std::vector<CGLA::Vec3d>& C) const {

!        const std::vector<std::vector<double>> dN = dS(q);

!        const CGLA::Vec3d Ju = std::inner_product(
!            begin(dN[0]), end(dN[0]), begin(C), CGLA::Vec3d{0.0});
!        const CGLA::Vec3d Jv = std::inner_product(
!            begin(dN[1]), end(dN[1]), begin(C), CGLA::Vec3d{0.0});
!        const CGLA::Vec3d Jw = std::inner_product(
!            begin(dN[2]), end(dN[2]), begin(C), CGLA::Vec3d{0.0});

!        Eigen::Matrix3d J;
!        J << Ju[0], Jv[0], Jw[0], //
!            Ju[1], Jv[1], Jw[1],  //
!            Ju[2], Jv[2], Jw[2];  //

!        return J;
!    }
! };
