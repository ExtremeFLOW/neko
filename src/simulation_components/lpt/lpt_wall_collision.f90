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
  use utils, only : neko_error
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use global_interpolation, only : global_interpolation_t
  use lpt_wall_collision_cpu, only : lpt_handle_elastic_wall_collisions_cpu
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
    type(matrix_t) :: rst_new
      !  type(vector_t) :: x_t
      !  type(vector_t) :: y_t
      !  type(vector_t) :: z_t
    type(vector_t) :: resx
    type(vector_t) :: resy
    type(vector_t) :: resz
    integer, allocatable :: el_list(:)
    integer :: i
    integer :: n

    n = x%size()
    if (n .eq. 0) return

    allocate(el_list(n))
    call rst_new%init(3, n)
    call resx%init(n)
    call resy%init(n)
    call resz%init(n)

    do i = 1, n
       el_list(i) = global_interp%el_owner0_local(i)
    end do

    call global_interp%rst_finder%find(rst_new, x, y, z, el_list, n, &
         resx, resy, resz)

    call lpt_handle_elastic_wall_collisions_cpu(msh, dm_Xh, coef, &
         wall_facet_mask, el_list, x_old, y_old, z_old, x, y, z, d, u, v, &
         w, u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
         acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag, u_old, &
         v_old, w_old, acc_x, acc_y, acc_z, lag_len, n)

    if (allocated(el_list)) deallocate(el_list)
    call rst_new%free()
    call resx%free()
    call resy%free()
    call resz%free()
  end subroutine lpt_handle_elastic_wall_collisions

end module lpt_wall_collision
