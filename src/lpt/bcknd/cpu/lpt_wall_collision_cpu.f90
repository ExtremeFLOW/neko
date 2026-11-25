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
  use vector, only : vector_t
  use lpt_wall_geometry, only : lpt_wall_geometry_t, LPT_WALL_MAX_BOUNCES, &
       LPT_WALL_RESTITUTION
  implicit none
  private

  public :: lpt_handle_elastic_wall_collisions_cpu

contains

  !> Detect and resolve collisions of inertial particles with the curved
  !! high-order wall facets on the CPU.
  !! The position is corrected by lpt_wall_geometry_t::collide, which
  !! intersects the particle path with the offset wall surface; every
  !! contact normal is then applied to the current velocity, the velocity
  !! and acceleration used by the time integrator in this step and their
  !! histories, so that the multistep scheme continues on the reflected
  !! trajectory.
  !! @param wall_geom High-order wall geometry and candidate facets.
  !! @param el_list Local owner element (0-based) of each particle at the
  !! start of the step.
  !! @param rst_old Reference coordinates of each particle in that element
  !! at the start of the step, used as initial guesses.
  !! @param lag_len Number of history levels in use.
  !! @param n Number of local particles.
  subroutine lpt_handle_elastic_wall_collisions_cpu(wall_geom, el_list, &
       rst_old, x_old, y_old, z_old, x, y, z, d, u, v, w, u_lag, v_lag, w_lag, &
       u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag, u_old, v_old, w_old, acc_x, &
       acc_y, acc_z, lag_len, n)
    type(lpt_wall_geometry_t), intent(inout) :: wall_geom
    integer, intent(in) :: el_list(:)
    real(kind=rp), intent(in) :: rst_old(:,:)
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
    real(kind=rp) :: normals(3, LPT_WALL_MAX_BOUNCES)
    real(kind=rp) :: p1(3), nrm(3)
    integer :: i, k, n_hits

    do i = 1, n
       if (el_list(i) .lt. 0) cycle

       p1 = [x%x(i), y%x(i), z%x(i)]
       call wall_geom%collide(el_list(i) + 1, &
            [x_old%x(i), y_old%x(i), z_old%x(i)], p1, 0.5_rp * d%x(i), &
            LPT_WALL_RESTITUTION, n_hits, normals, rst_old(:, i))
       x%x(i) = p1(1)
       y%x(i) = p1(2)
       z%x(i) = p1(3)

       do k = 1, n_hits
          nrm = normals(:, k)
          call reflect(u%x(i), v%x(i), w%x(i), nrm)
          call reflect(u_old%x(i), v_old%x(i), w_old%x(i), nrm)
          call reflect(acc_x%x(i), acc_y%x(i), acc_z%x(i), nrm)
          if (lag_len .ge. 1) then
             call reflect(u_lag%x(i), v_lag%x(i), w_lag%x(i), nrm)
             call reflect(acc_xlag%x(i), acc_ylag%x(i), acc_zlag%x(i), nrm)
          end if
          if (lag_len .ge. 2) then
             call reflect(u_laglag%x(i), v_laglag%x(i), w_laglag%x(i), nrm)
             call reflect(acc_xlaglag%x(i), acc_ylaglag%x(i), &
                  acc_zlaglag%x(i), nrm)
          end if
       end do
    end do

  end subroutine lpt_handle_elastic_wall_collisions_cpu

  !> Apply \f$ v^+ = v - (1 + e_n)(v \cdot n) n \f$ to vector components.
  !! @param nrm Unit contact normal.
  pure subroutine reflect(vx, vy, vz, nrm)
    real(kind=rp), intent(inout) :: vx, vy, vz
    real(kind=rp), intent(in) :: nrm(3)
    real(kind=rp) :: vn

    vn = (1.0_rp + LPT_WALL_RESTITUTION) * &
         (vx * nrm(1) + vy * nrm(2) + vz * nrm(3))
    vx = vx - vn * nrm(1)
    vy = vy - vn * nrm(2)
    vz = vz - vn * nrm(3)

  end subroutine reflect

end module lpt_wall_collision_cpu
