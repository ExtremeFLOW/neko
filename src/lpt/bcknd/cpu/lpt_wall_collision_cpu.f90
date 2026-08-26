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
!> CPU implementation for LPT wall-collision handling.
module lpt_wall_collision_cpu
  use num_types, only : rp
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use vector, only : vector_t
  implicit none
  private

  public :: lpt_handle_elastic_wall_collisions_cpu

contains

  !> Reflect inertial particles that hit configured wall facets on the CPU.
  !! @param msh Mesh containing wall facets.
  !! @param dm_Xh Coordinate dofmap used for facet locations.
  !! @param coef Coefficients used for facet normals.
  !! @param wall_facet_mask Mask of facets that behave as elastic walls.
  !! @param el_list Local owner element for each particle.
  !! @param n Number of local particles.
  subroutine lpt_handle_elastic_wall_collisions_cpu(msh, dm_Xh, coef, &
       wall_facet_mask, el_list, x_old, y_old, z_old, x, y, z, d, u, v, w, &
       u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
       acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag, u_old, &
       v_old, w_old, acc_x, acc_y, acc_z, lag_len, n)
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), target, intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: wall_facet_mask(:, :)
    integer, intent(in) :: el_list(:)
    type(vector_t), intent(in) :: x_old, y_old, z_old
    type(vector_t), intent(inout) :: x, y, z
    type(vector_t), intent(in) :: d
    type(vector_t), intent(inout) :: u, v, w
    type(vector_t), intent(inout) :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout) :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout) :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout) :: acc_xlaglag, acc_ylaglag
    type(vector_t), intent(inout) :: acc_zlaglag
    type(vector_t), intent(inout) :: u_old, v_old, w_old
    type(vector_t), intent(inout) :: acc_x, acc_y, acc_z
    integer, intent(in) :: lag_len
    integer, intent(in) :: n
    integer :: i
    integer :: facet
    integer :: candidate
    integer :: el
    integer :: el_mesh
    integer :: hit_count
    integer :: hit_idx
    integer :: hit_facets(6)
    real(kind=rp) :: normal(3)
    real(kind=rp) :: wall_point(3)
    real(kind=rp) :: radius
    real(kind=rp), pointer, dimension(:,:,:,:) :: dm_x, dm_y, dm_z
    integer :: lx
    integer :: ly
    integer :: lz

    ! Once per call rather than at setup: wall_facet_normal() below is only
    ! reached from here, and this routine has no init to hang the check on
    call coef%require_facets('lpt wall collisions')

    dm_x => dm_Xh%x
    dm_y => dm_Xh%y
    dm_z => dm_Xh%z
    lx = dm_Xh%Xh%lx
    ly = dm_Xh%Xh%ly
    lz = dm_Xh%Xh%lz

    do i = 1, n
       el = el_list(i)
       if (el .lt. 0) cycle
       el_mesh = el + 1
       if (el_mesh .gt. msh%nelv) cycle

       radius = 0.5_rp * d%x(i)

       hit_count = 0
       do candidate = 1, 2 * msh%gdim
          if (wall_facet_is_hit(wall_facet_mask, dm_x, dm_y, dm_z, &
               lx, ly, lz, coef, &
               x_old%x(i), y_old%x(i), z_old%x(i), x%x(i), y%x(i), &
               z%x(i), radius, el_mesh, candidate, msh%gdim)) then
             hit_count = hit_count + 1
             hit_facets(hit_count) = candidate
          end if
       end do
       if (hit_count .eq. 0) cycle

       do hit_idx = 1, hit_count
          facet = hit_facets(hit_idx)
          call wall_facet_normal(coef, lx, ly, lz, el_mesh, facet, normal)
          if (norm2(normal) .le. epsilon(1.0_rp)) cycle

          call wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el_mesh, &
               facet, wall_point)
          call reflect_position(x%x(i), y%x(i), z%x(i), wall_point, &
               normal, radius)
          call reflect_vector_components(u%x(i), v%x(i), w%x(i), normal)
          call reflect_vector_components(u_old%x(i), v_old%x(i), &
               w_old%x(i), normal)
          call reflect_vector_components(acc_x%x(i), acc_y%x(i), &
               acc_z%x(i), normal)
          if (lag_len .ge. 1) then
             call reflect_vector_components(u_lag%x(i), v_lag%x(i), &
                  w_lag%x(i), normal)
             call reflect_vector_components(acc_xlag%x(i), acc_ylag%x(i), &
                  acc_zlag%x(i), normal)
          end if
          if (lag_len .ge. 2) then
             call reflect_vector_components(u_laglag%x(i), v_laglag%x(i), &
                  w_laglag%x(i), normal)
             call reflect_vector_components(acc_xlaglag%x(i), &
                  acc_ylaglag%x(i), acc_zlaglag%x(i), normal)
          end if
       end do
    end do
  end subroutine lpt_handle_elastic_wall_collisions_cpu

  !> Return true if the particle surface reaches a given wall facet.
  logical function wall_facet_is_hit(wall_facet_mask, dm_x, dm_y, dm_z, &
       lx, ly, lz, coef, x_old, y_old, z_old, x, y, z, radius, el, facet, &
       gdim) result(is_hit)
    logical, intent(in) :: wall_facet_mask(:, :)
    real(kind=rp), pointer, dimension(:,:,:,:), intent(in) :: dm_x
    real(kind=rp), pointer, dimension(:,:,:,:), intent(in) :: dm_y
    real(kind=rp), pointer, dimension(:,:,:,:), intent(in) :: dm_z
    integer, intent(in) :: lx
    integer, intent(in) :: ly
    integer, intent(in) :: lz
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: x_old
    real(kind=rp), intent(in) :: y_old
    real(kind=rp), intent(in) :: z_old
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
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

    call wall_facet_normal(coef, lx, ly, lz, el, facet, normal)
    if (norm2(normal) .le. epsilon(1.0_rp)) return
    call wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el, facet, wall_point)
    dist_old = signed_plane_distance(x_old, y_old, z_old, wall_point, normal)
    dist_new = signed_plane_distance(x, y, z, wall_point, normal)
    penetration = dist_new + radius

    if (penetration .le. tol) return
    if (dist_new .le. dist_old + tol) return

    is_hit = .true.
  end function wall_facet_is_hit

  !> Use the face-center SEM normal as the reflection normal.
  !! @param normal Output normal vector.
  subroutine wall_facet_normal(coef, lx, ly, lz, el, facet, normal)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: lx
    integer, intent(in) :: ly
    integer, intent(in) :: lz
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: normal(3)
    integer :: ic
    integer :: jc
    integer :: kc

    ic = max(1, (lx + 1) / 2)
    jc = max(1, (ly + 1) / 2)
    kc = max(1, (lz + 1) / 2)

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
  !! @param wall_point Output point on the wall plane.
  subroutine wall_facet_center(dm_x, dm_y, dm_z, lx, ly, lz, el, facet, &
       wall_point)
    real(kind=rp), pointer, dimension(:,:,:,:), intent(in) :: dm_x, dm_y, dm_z
    integer, intent(in) :: lx, ly, lz
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: wall_point(3)
    integer :: ic
    integer :: jc
    integer :: kc

    ic = max(1, (lx + 1) / 2)
    jc = max(1, (ly + 1) / 2)
    kc = max(1, (lz + 1) / 2)

    select case (facet)
    case (1)
       wall_point = [dm_x(1, jc, kc, el), dm_y(1, jc, kc, el), &
            dm_z(1, jc, kc, el)]
    case (2)
       wall_point = [dm_x(lx, jc, kc, el), &
            dm_y(lx, jc, kc, el), &
            dm_z(lx, jc, kc, el)]
    case (3)
       wall_point = [dm_x(ic, 1, kc, el), dm_y(ic, 1, kc, el), &
            dm_z(ic, 1, kc, el)]
    case (4)
       wall_point = [dm_x(ic, ly, kc, el), &
            dm_y(ic, ly, kc, el), &
            dm_z(ic, ly, kc, el)]
    case (5)
       wall_point = [dm_x(ic, jc, 1, el), dm_y(ic, jc, 1, el), &
            dm_z(ic, jc, 1, el)]
    case (6)
       wall_point = [dm_x(ic, jc, lz, el), &
            dm_y(ic, jc, lz, el), &
            dm_z(ic, jc, lz, el)]
    case default
       wall_point = 0.0_rp
    end select
  end subroutine wall_facet_center

  !> Reflect the particle center across the contact plane located one radius
  !! inward from the wall plane.
  !! @param x Particle x coordinate updated in place.
  !! @param y Particle y coordinate updated in place.
  !! @param z Particle z coordinate updated in place.
  !! @param wall_point Point on the wall plane.
  !! @param normal Wall normal.
  !! @param radius Particle radius.
  subroutine reflect_position(x, y, z, wall_point, normal, radius)
    real(kind=rp), intent(inout) :: x, y, z
    real(kind=rp), intent(in) :: wall_point(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp), intent(in) :: radius
    real(kind=rp) :: nhat(3)
    real(kind=rp) :: nmag
    real(kind=rp) :: signed_contact_distance

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) return

    nhat = normal / nmag
    signed_contact_distance = dot_product([x, y, z] - wall_point, nhat) + &
         radius
    if (signed_contact_distance .le. 0.0_rp) return

    x = x - 2.0_rp * signed_contact_distance * nhat(1)
    y = y - 2.0_rp * signed_contact_distance * nhat(2)
    z = z - 2.0_rp * signed_contact_distance * nhat(3)
  end subroutine reflect_position

  !> Return signed distance from a point to a wall plane.
  !! @param wall_point Point on the wall plane.
  !! @param normal Wall normal.
  pure real(kind=rp) function signed_plane_distance(x, y, z, wall_point, &
       normal)
    real(kind=rp), intent(in) :: x, y, z
    real(kind=rp), intent(in) :: wall_point(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp) :: nmag

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) then
       signed_plane_distance = -huge(1.0_rp)
       return
    end if

    signed_plane_distance = dot_product([x, y, z] - wall_point, normal / nmag)
  end function signed_plane_distance

  !> Reflect a vector across a plane with the supplied normal.
  !! @param vec Vector updated in place.
  !! @param normal Plane normal.
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

  !> Reflect vector components across a plane with the supplied normal.
  !! @param x X component updated in place.
  !! @param y Y component updated in place.
  !! @param z Z component updated in place.
  !! @param normal Plane normal.
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

end module lpt_wall_collision_cpu
