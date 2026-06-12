! Copyright (c) 2026, The Neko Authors
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
!> Helpers for LPT wall-collision handling.
module lpt_wall_collision
  use num_types, only : rp
  use utils, only : neko_error
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use global_interpolation, only : global_interpolation_t
  use matrix, only : matrix_t
  use vector, only : vector_t
  use point, only : point_t
  use point_interpolator, only : point_interpolator_t
  implicit none
  private

  public :: lpt_handle_elastic_wall_collisions
  public :: lpt_init_wall_facet_mask

contains

  subroutine lpt_init_wall_facet_mask(wall_facet_mask, msh, wall_zone_indices)
    logical, allocatable, intent(inout) :: wall_facet_mask(:, :)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: wall_zone_indices(:)
    integer :: i
    integer :: j
    integer :: facet
    integer :: el
    integer :: zone_id

    if (allocated(wall_facet_mask)) deallocate(wall_facet_mask)
    allocate(wall_facet_mask(2 * msh%gdim, msh%nelv))
    wall_facet_mask = .false.

    do i = 1, size(wall_zone_indices)
       zone_id = wall_zone_indices(i)
       if (zone_id .lt. 1 .or. zone_id .gt. size(msh%labeled_zones)) then
          call neko_error("lpt wall_zone_indices contains an invalid zone id")
       end if

       do j = 1, msh%labeled_zones(zone_id)%size
          facet = msh%labeled_zones(zone_id)%facet_el(j)%x(1)
          el = msh%labeled_zones(zone_id)%facet_el(j)%x(2)
          wall_facet_mask(facet, el) = .true.
       end do
    end do
  end subroutine lpt_init_wall_facet_mask

  !> Restore particles that crossed a configured wall zone and reflect the
  !! remaining trajectory together with the current particle velocity and
  !! velocity history for a purely elastic collision.
  subroutine lpt_handle_elastic_wall_collisions(global_interp, msh, dm_Xh, &
       coef, wall_facet_mask, xyz, vel, vel_lag, vel_rhs, lag_len)
    type(global_interpolation_t), intent(inout) :: global_interp
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(inout) :: dm_Xh
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: wall_facet_mask(:, :)
    real(kind=rp), intent(inout) :: xyz(:, :)
    real(kind=rp), intent(inout) :: vel(:, :)
    real(kind=rp), intent(inout) :: vel_lag(:, :, :)
    real(kind=rp), intent(inout) :: vel_rhs(:, :)
    integer, intent(in) :: lag_len
    type(matrix_t) :: rst_new
    type(vector_t) :: x_t
    type(vector_t) :: y_t
    type(vector_t) :: z_t
    type(vector_t) :: resx
    type(vector_t) :: resy
    type(vector_t) :: resz
    integer, allocatable :: el_list(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: facet
    integer :: el
    integer :: el_mesh
    real(kind=rp) :: normal(3)

    n = size(xyz, 2)
    if (n .eq. 0) return

    allocate(el_list(n))
    call rst_new%init(3, n)
    call x_t%init(n)
    call y_t%init(n)
    call z_t%init(n)
    call resx%init(n)
    call resy%init(n)
    call resz%init(n)

    do i = 1, n
       el_list(i) = global_interp%el_owner0_local(i)
       x_t%x(i) = xyz(1, i)
       y_t%x(i) = xyz(2, i)
       z_t%x(i) = xyz(3, i)
    end do

    call global_interp%rst_finder%find(rst_new, x_t, y_t, z_t, el_list, n, &
         resx, resy, resz)

    do i = 1, n
       el = el_list(i)
       if (el .lt. 0) cycle
       el_mesh = el + 1
       if (el_mesh .gt. msh%nelv) cycle

       facet = identify_wall_facet(wall_facet_mask, global_interp%rst_local(:, i), &
            rst_new%x(:, i), el_mesh, msh%gdim)
       if (facet .eq. 0) cycle

       call wall_facet_normal(coef, el_mesh, facet, normal)
       if (norm2(normal) .le. epsilon(1.0_rp)) cycle

       call reflect_position(dm_Xh, xyz(:, i), global_interp%rst_local(:, i), &
            rst_new%x(:, i), el_mesh, facet)
       call reflect_vector(vel_rhs(:, i), normal)
       do j = 1, lag_len
          call reflect_vector(vel_lag(:, j, i), normal)
       end do
       vel(:, i) = vel_rhs(:, i)
    end do

    if (allocated(el_list)) deallocate(el_list)
    call rst_new%free()
    call x_t%free()
    call y_t%free()
    call z_t%free()
    call resx%free()
    call resy%free()
    call resz%free()
  end subroutine lpt_handle_elastic_wall_collisions

  !> Return the wall facet crossed by a particle, or 0 if none matched.
  integer function identify_wall_facet(wall_facet_mask, rst_old, rst_new, el, &
       gdim) result(facet)
    logical, intent(in) :: wall_facet_mask(:, :)
    real(kind=rp), intent(in) :: rst_old(3)
    real(kind=rp), intent(in) :: rst_new(3)
    integer, intent(in) :: el
    integer, intent(in) :: gdim
    real(kind=rp), parameter :: tol = 1.0e-8_rp
    real(kind=rp) :: severity
    real(kind=rp) :: best_severity

    facet = 0
    best_severity = 0.0_rp

    severity = max((-1.0_rp - rst_new(1)), 0.0_rp)
    if (wall_facet_mask(1, el) .and. severity .gt. tol .and. &
         rst_new(1) .le. rst_old(1) + tol .and. severity .gt. best_severity) then
       facet = 1
       best_severity = severity
    end if

    severity = max((rst_new(1) - 1.0_rp), 0.0_rp)
    if (wall_facet_mask(2, el) .and. severity .gt. tol .and. &
         rst_new(1) .ge. rst_old(1) - tol .and. severity .gt. best_severity) then
       facet = 2
       best_severity = severity
    end if

    severity = max((-1.0_rp - rst_new(2)), 0.0_rp)
    if (wall_facet_mask(3, el) .and. severity .gt. tol .and. &
         rst_new(2) .le. rst_old(2) + tol .and. severity .gt. best_severity) then
       facet = 3
       best_severity = severity
    end if

    severity = max((rst_new(2) - 1.0_rp), 0.0_rp)
    if (wall_facet_mask(4, el) .and. severity .gt. tol .and. &
         rst_new(2) .ge. rst_old(2) - tol .and. severity .gt. best_severity) then
       facet = 4
       best_severity = severity
    end if

    if (gdim .eq. 3) then
       severity = max((-1.0_rp - rst_new(3)), 0.0_rp)
       if (wall_facet_mask(5, el) .and. severity .gt. tol .and. &
            rst_new(3) .le. rst_old(3) + tol .and. &
            severity .gt. best_severity) then
          facet = 5
          best_severity = severity
       end if

       severity = max((rst_new(3) - 1.0_rp), 0.0_rp)
       if (wall_facet_mask(6, el) .and. severity .gt. tol .and. &
            rst_new(3) .ge. rst_old(3) - tol .and. &
            severity .gt. best_severity) then
          facet = 6
       end if
    end if
  end function identify_wall_facet

  !> Use the face-center SEM normal as the reflection normal.
  subroutine wall_facet_normal(coef, el, facet, normal)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: normal(3)
    integer :: ic
    integer :: jc
    integer :: kc

    ic = max(1, (coef%Xh%lx + 1) / 2)
    jc = max(1, (coef%Xh%ly + 1) / 2)
    kc = max(1, (coef%Xh%lz + 1) / 2)

    select case (facet)
    case (1, 2)
       normal = coef%get_normal(1, jc, kc, el, facet)
    case (3, 4)
       normal = coef%get_normal(ic, 1, kc, el, facet)
    case (5, 6)
       normal = coef%get_normal(ic, jc, 1, el, facet)
    case default
       normal = 0.0_rp
    end select
  end subroutine wall_facet_normal

  !> Reflect the post-collision position by mirroring the overshoot in rst and
  !! interpolating the reflected point back to xyz on the same element.
  subroutine reflect_position(dm_Xh, xyz, rst_old, rst_new, el, facet)
    type(dofmap_t), intent(inout) :: dm_Xh
    real(kind=rp), intent(inout) :: xyz(3)
    real(kind=rp), intent(in) :: rst_old(3)
    real(kind=rp), intent(in) :: rst_new(3)
    integer, intent(in) :: el
    integer, intent(in) :: facet
    type(point_interpolator_t) :: interp
    type(point_t) :: rst_point(1)
    type(point_t), allocatable :: xyz_ref(:)
    real(kind=rp) :: collision_rst(3)
    real(kind=rp) :: reflected_rst(3)
    real(kind=rp) :: drst(3)
    real(kind=rp) :: alpha
    real(kind=rp) :: boundary_value
    integer :: dim

    reflected_rst = rst_new
    call facet_rst_info(facet, dim, boundary_value)
    drst = rst_new - rst_old

    if (abs(drst(dim)) .le. epsilon(1.0_rp)) return

    alpha = (boundary_value - rst_old(dim)) / drst(dim)
    alpha = max(0.0_rp, min(1.0_rp, alpha))
    collision_rst = rst_old + alpha * drst
    reflected_rst = collision_rst + (rst_new - collision_rst)
    reflected_rst(dim) = 2.0_rp * boundary_value - rst_new(dim)

    call rst_point(1)%init(real(reflected_rst, kind=kind(rst_point(1)%x)))
    call interp%init(dm_Xh%Xh)
    xyz_ref = interp%interpolate(rst_point, dm_Xh%x(:, :, :, el), &
         dm_Xh%y(:, :, :, el), dm_Xh%z(:, :, :, el))
    xyz = real(xyz_ref(1)%x, kind=rp)
    call interp%free()
    if (allocated(xyz_ref)) deallocate(xyz_ref)
  end subroutine reflect_position

  pure subroutine facet_rst_info(facet, dim, boundary_value)
    integer, intent(in) :: facet
    integer, intent(out) :: dim
    real(kind=rp), intent(out) :: boundary_value

    select case (facet)
    case (1)
       dim = 1
       boundary_value = -1.0_rp
    case (2)
       dim = 1
       boundary_value = 1.0_rp
    case (3)
       dim = 2
       boundary_value = -1.0_rp
    case (4)
       dim = 2
       boundary_value = 1.0_rp
    case (5)
       dim = 3
       boundary_value = -1.0_rp
    case (6)
       dim = 3
       boundary_value = 1.0_rp
    case default
       dim = 1
       boundary_value = 0.0_rp
    end select
  end subroutine facet_rst_info

  pure subroutine reflect_vector(vec, normal)
    real(kind=rp), intent(inout) :: vec(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp) :: nhat(3)
    real(kind=rp) :: nmag
    real(kind=rp) :: vn

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) return

    nhat = normal / nmag
    vn = dot_product(vec, nhat)
    vec = vec - 2.0_rp * vn * nhat
  end subroutine reflect_vector

end module lpt_wall_collision
