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
       coef, wall_facet_mask, x_old, y_old, z_old, x, y, z, d, u, v, w, &
       u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
       acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag, u_old, &
       v_old, w_old, acc_x, acc_y, acc_z, lag_len)
    type(global_interpolation_t), intent(inout) :: global_interp
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: wall_facet_mask(:, :)
    real(kind=rp), intent(in) :: x_old(:), y_old(:), z_old(:)
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    real(kind=rp), intent(in) :: d(:)
    real(kind=rp), intent(inout) :: u(:), v(:), w(:)
    real(kind=rp), intent(inout) :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), intent(inout) :: u_laglag(:), v_laglag(:), w_laglag(:)
    real(kind=rp), intent(inout) :: acc_xlag(:), acc_ylag(:), acc_zlag(:)
    real(kind=rp), intent(inout) :: acc_xlaglag(:), acc_ylaglag(:)
    real(kind=rp), intent(inout) :: acc_zlaglag(:)
    real(kind=rp), intent(inout) :: u_old(:), v_old(:), w_old(:)
    real(kind=rp), intent(inout) :: acc_x(:), acc_y(:), acc_z(:)
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
    integer :: n
    integer :: facet
    integer :: candidate
    integer :: el
    integer :: el_mesh
    integer :: hit_count
    integer :: hit_idx
    integer :: hit_facets(6)
    real(kind=rp) :: normal(3)
    real(kind=rp) :: wall_point(3)
    real(kind=rp) :: xyz_old(3)
    real(kind=rp) :: xyz_new(3)
    real(kind=rp) :: radius

    n = size(x)
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
       x_t%x(i) = x(i)
       y_t%x(i) = y(i)
       z_t%x(i) = z(i)
    end do

    call global_interp%rst_finder%find(rst_new, x_t, y_t, z_t, el_list, n, &
         resx, resy, resz)

    do i = 1, n
       el = el_list(i)
       if (el .lt. 0) cycle
       el_mesh = el + 1
       if (el_mesh .gt. msh%nelv) cycle

       radius = 0.5_rp * d(i)
       xyz_old = [x_old(i), y_old(i), z_old(i)]
       xyz_new = [x(i), y(i), z(i)]

       hit_count = 0
       do candidate = 1, 2 * msh%gdim
          if (wall_facet_is_hit(wall_facet_mask, dm_Xh, coef, xyz_old, &
               xyz_new, radius, el_mesh, candidate, msh%gdim)) then
             hit_count = hit_count + 1
             hit_facets(hit_count) = candidate
          end if
       end do
       if (hit_count .eq. 0) cycle

       do hit_idx = 1, hit_count
          facet = hit_facets(hit_idx)
          call wall_facet_normal(coef, el_mesh, facet, normal)
          if (norm2(normal) .le. epsilon(1.0_rp)) cycle

          call wall_facet_center(dm_Xh, el_mesh, facet, wall_point)
          call reflect_position(xyz_new, wall_point, normal, radius)
          call reflect_vector_components(u(i), v(i), w(i), normal)
          call reflect_vector_components(u_old(i), v_old(i), w_old(i), &
               normal)
          call reflect_vector_components(acc_x(i), acc_y(i), acc_z(i), &
               normal)
          if (lag_len .ge. 1) then
             call reflect_vector_components(u_lag(i), v_lag(i), w_lag(i), &
                  normal)
             call reflect_vector_components(acc_xlag(i), acc_ylag(i), &
                  acc_zlag(i), normal)
          end if
          if (lag_len .ge. 2) then
             call reflect_vector_components(u_laglag(i), v_laglag(i), &
                  w_laglag(i), normal)
             call reflect_vector_components(acc_xlaglag(i), acc_ylaglag(i), &
                  acc_zlaglag(i), normal)
          end if
       end do

       x(i) = xyz_new(1)
       y(i) = xyz_new(2)
       z(i) = xyz_new(3)
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

  !> Return true if the particle surface reaches a given wall facet.
  logical function wall_facet_is_hit(wall_facet_mask, dm_Xh, coef, xyz_old, &
       xyz_new, radius, el, facet, gdim) result(is_hit)
    logical, intent(in) :: wall_facet_mask(:, :)
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: xyz_old(3)
    real(kind=rp), intent(in) :: xyz_new(3)
    real(kind=rp), intent(in) :: radius
    integer, intent(in) :: el
    integer, intent(in) :: facet
    integer, intent(in) :: gdim
    real(kind=rp), parameter :: tol = 1.0e-8_rp
    real(kind=rp) :: normal(3)
    real(kind=rp) :: wall_point(3)
    real(kind=rp) :: dist_old
    real(kind=rp) :: dist_new
    real(kind=rp) :: penetration

    is_hit = .false.

    if (facet .lt. 1 .or. facet .gt. 2 * gdim) return
    if (.not. wall_facet_mask(facet, el)) return

    call wall_facet_normal(coef, el, facet, normal)
    if (norm2(normal) .le. epsilon(1.0_rp)) return
    call wall_facet_center(dm_Xh, el, facet, wall_point)
    dist_old = signed_plane_distance(xyz_old, wall_point, normal)
    dist_new = signed_plane_distance(xyz_new, wall_point, normal)
    penetration = dist_new + radius

    if (penetration .le. tol) return
    if (dist_new .le. dist_old + tol) return

    is_hit = .true.
  end function wall_facet_is_hit

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

  !> Return the facet-center position that matches the reflection normal.
  subroutine wall_facet_center(dm_Xh, el, facet, wall_point)
    type(dofmap_t), intent(in) :: dm_Xh
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: wall_point(3)
    integer :: ic
    integer :: jc
    integer :: kc

    ic = max(1, (dm_Xh%Xh%lx + 1) / 2)
    jc = max(1, (dm_Xh%Xh%ly + 1) / 2)
    kc = max(1, (dm_Xh%Xh%lz + 1) / 2)

    select case (facet)
    case (1)
       wall_point = [dm_Xh%x(1, jc, kc, el), dm_Xh%y(1, jc, kc, el), &
            dm_Xh%z(1, jc, kc, el)]
    case (2)
       wall_point = [dm_Xh%x(dm_Xh%Xh%lx, jc, kc, el), &
            dm_Xh%y(dm_Xh%Xh%lx, jc, kc, el), &
            dm_Xh%z(dm_Xh%Xh%lx, jc, kc, el)]
    case (3)
       wall_point = [dm_Xh%x(ic, 1, kc, el), dm_Xh%y(ic, 1, kc, el), &
            dm_Xh%z(ic, 1, kc, el)]
    case (4)
       wall_point = [dm_Xh%x(ic, dm_Xh%Xh%ly, kc, el), &
            dm_Xh%y(ic, dm_Xh%Xh%ly, kc, el), &
            dm_Xh%z(ic, dm_Xh%Xh%ly, kc, el)]
    case (5)
       wall_point = [dm_Xh%x(ic, jc, 1, el), dm_Xh%y(ic, jc, 1, el), &
            dm_Xh%z(ic, jc, 1, el)]
    case (6)
       wall_point = [dm_Xh%x(ic, jc, dm_Xh%Xh%lz, el), &
            dm_Xh%y(ic, jc, dm_Xh%Xh%lz, el), &
            dm_Xh%z(ic, jc, dm_Xh%Xh%lz, el)]
    case default
       wall_point = 0.0_rp
    end select
  end subroutine wall_facet_center

  !> Reflect the particle center across the contact plane located one radius
  !! inward from the wall plane.
  subroutine reflect_position(point, wall_point, normal, radius)
    real(kind=rp), intent(inout) :: point(3)
    real(kind=rp), intent(in) :: wall_point(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp), intent(in) :: radius
    real(kind=rp) :: nhat(3)
    real(kind=rp) :: nmag
    real(kind=rp) :: signed_contact_distance

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) return

    nhat = normal / nmag
    signed_contact_distance = dot_product(point - wall_point, nhat) + radius
    if (signed_contact_distance .le. 0.0_rp) return

    point = point - 2.0_rp * signed_contact_distance * nhat
  end subroutine reflect_position

  pure real(kind=rp) function signed_plane_distance(point, wall_point, normal)
    real(kind=rp), intent(in) :: point(3)
    real(kind=rp), intent(in) :: wall_point(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp) :: nmag

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) then
       signed_plane_distance = -huge(1.0_rp)
       return
    end if

    signed_plane_distance = dot_product(point - wall_point, normal / nmag)
  end function signed_plane_distance

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

  pure subroutine reflect_vector_components(x, y, z, normal)
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(inout) :: y
    real(kind=rp), intent(inout) :: z
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp) :: vec(3)

    vec = [x, y, z]
    call reflect_vector(vec, normal)
    x = vec(1)
    y = vec(2)
    z = vec(3)
  end subroutine reflect_vector_components

end module lpt_wall_collision
