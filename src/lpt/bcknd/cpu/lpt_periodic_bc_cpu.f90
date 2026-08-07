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
!> CPU kernels for LPT periodic boundary-condition wrapping.
module lpt_periodic_bc_cpu
  use num_types, only : rp
  use vector, only : vector_t
  implicit none
  private

  real(kind=rp), parameter :: LPT_PERIODIC_TOL = 1.0e-8_rp

  public :: lpt_periodic_bc_wrap_rotational_cpu
  public :: lpt_periodic_bc_wrap_translational_cpu

contains

  !> Wrap particles through a rotational periodic sector on the CPU.
  !! @param x Particle x coordinates.
  !! @param y Particle y coordinates.
  !! @param z Particle z coordinates.
  !! @param n Number of local particles.
  !! @param theta_min Minimum sector angle.
  !! @param theta_max Maximum sector angle.
  !! @param theta_len Angular sector length.
  subroutine lpt_periodic_bc_wrap_rotational_cpu(x, y, z, n, theta_min, &
       theta_max, theta_len, u, v, w, u_lag, v_lag, w_lag, u_laglag, &
       v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, &
       acc_ylaglag, acc_zlaglag)
    type(vector_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: theta_min
    real(kind=rp), intent(in) :: theta_max
    real(kind=rp), intent(in) :: theta_len
    type(vector_t), intent(inout), optional :: u, v, w
    type(vector_t), intent(inout), optional :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout), optional :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout), optional :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout), optional :: acc_xlaglag
    type(vector_t), intent(inout), optional :: acc_ylaglag
    type(vector_t), intent(inout), optional :: acc_zlaglag
    integer :: i
    real(kind=rp) :: radius
    real(kind=rp) :: theta_old
    real(kind=rp) :: theta
    real(kind=rp) :: dtheta
    real(kind=rp) :: pi

    pi = acos(-1.0_rp)
    do i = 1, n
       radius = sqrt(x%x(i) * x%x(i) + y%x(i) * y%x(i))
       theta_old = modulo(atan2(y%x(i), x%x(i)) + 2.0_rp * pi, &
            2.0_rp * pi)
       theta = theta_old

       do while (theta .lt. theta_min - LPT_PERIODIC_TOL)
          theta = theta + theta_len
       end do

       do while (theta .gt. theta_max + LPT_PERIODIC_TOL)
          theta = theta - theta_len
       end do

       dtheta = theta - theta_old
       x%x(i) = radius * cos(theta)
       y%x(i) = radius * sin(theta)
       if (abs(dtheta) .le. LPT_PERIODIC_TOL) cycle

       if (present(u) .and. present(v)) &
            call lpt_rotate_xy(u%x(i), v%x(i), dtheta)
       if (present(u_lag) .and. present(v_lag)) &
            call lpt_rotate_xy(u_lag%x(i), v_lag%x(i), dtheta)
       if (present(u_laglag) .and. present(v_laglag)) &
            call lpt_rotate_xy(u_laglag%x(i), v_laglag%x(i), dtheta)
       if (present(acc_xlag) .and. present(acc_ylag)) &
            call lpt_rotate_xy(acc_xlag%x(i), acc_ylag%x(i), dtheta)
       if (present(acc_xlaglag) .and. present(acc_ylaglag)) &
            call lpt_rotate_xy(acc_xlaglag%x(i), acc_ylaglag%x(i), dtheta)
    end do
  end subroutine lpt_periodic_bc_wrap_rotational_cpu

  !> Wrap particle coordinates through translational periodic directions.
  !! @param x Particle x coordinates.
  !! @param y Particle y coordinates.
  !! @param z Particle z coordinates.
  !! @param n Number of local particles.
  subroutine lpt_periodic_bc_wrap_translational_cpu(x, y, z, n, &
       periodic_enabled, n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, &
       periodic_dir_z1, periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, &
       periodic_dir_x3, periodic_dir_y3, periodic_dir_z3, periodic_min1, &
       periodic_min2, periodic_min3, periodic_max1, periodic_max2, &
       periodic_max3, periodic_shift_x1, periodic_shift_y1, &
       periodic_shift_z1, periodic_shift_x2, periodic_shift_y2, &
       periodic_shift_z2, periodic_shift_x3, periodic_shift_y3, &
       periodic_shift_z3, periodic_len1, periodic_len2, periodic_len3)
    type(vector_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    logical, intent(in) :: periodic_enabled
    integer, intent(in) :: n_periodic_dirs
    real(kind=rp), intent(in) :: periodic_dir_x1, periodic_dir_y1
    real(kind=rp), intent(in) :: periodic_dir_z1, periodic_dir_x2
    real(kind=rp), intent(in) :: periodic_dir_y2, periodic_dir_z2
    real(kind=rp), intent(in) :: periodic_dir_x3, periodic_dir_y3
    real(kind=rp), intent(in) :: periodic_dir_z3
    real(kind=rp), intent(in) :: periodic_min1, periodic_min2, periodic_min3
    real(kind=rp), intent(in) :: periodic_max1, periodic_max2, periodic_max3
    real(kind=rp), intent(in) :: periodic_shift_x1, periodic_shift_y1
    real(kind=rp), intent(in) :: periodic_shift_z1, periodic_shift_x2
    real(kind=rp), intent(in) :: periodic_shift_y2, periodic_shift_z2
    real(kind=rp), intent(in) :: periodic_shift_x3, periodic_shift_y3
    real(kind=rp), intent(in) :: periodic_shift_z3
    real(kind=rp), intent(in) :: periodic_len1, periodic_len2, periodic_len3
    integer :: i
    integer :: j
    real(kind=rp) :: point(3)
    real(kind=rp) :: dir(3)
    real(kind=rp) :: shift(3)
    real(kind=rp) :: coord
    real(kind=rp) :: periodic_min
    real(kind=rp) :: periodic_max
    real(kind=rp) :: periodic_len

    if (.not. periodic_enabled) return

    do i = 1, n
       point = [x%x(i), y%x(i), z%x(i)]
       do j = 1, n_periodic_dirs
          call lpt_periodic_bc_get_translational_slot(j, &
               periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
               periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, &
               periodic_dir_x3, periodic_dir_y3, periodic_dir_z3, &
               periodic_min1, periodic_min2, periodic_min3, periodic_max1, &
               periodic_max2, periodic_max3, periodic_shift_x1, &
               periodic_shift_y1, periodic_shift_z1, periodic_shift_x2, &
               periodic_shift_y2, periodic_shift_z2, periodic_shift_x3, &
               periodic_shift_y3, periodic_shift_z3, periodic_len1, &
               periodic_len2, periodic_len3, dir, periodic_min, &
               periodic_max, shift, periodic_len)

          coord = dot_product(point, dir)

          do while (coord .lt. periodic_min - LPT_PERIODIC_TOL)
             point = point + shift
             coord = coord + periodic_len
          end do

          do while (coord .gt. periodic_max + LPT_PERIODIC_TOL)
             point = point - shift
             coord = coord - periodic_len
          end do
       end do
       x%x(i) = point(1)
       y%x(i) = point(2)
       z%x(i) = point(3)
    end do
  end subroutine lpt_periodic_bc_wrap_translational_cpu

  !> Rotate a two-component vector in the x-y plane.
  !! @param x X component updated in place.
  !! @param y Y component updated in place.
  !! @param theta Rotation angle.
  subroutine lpt_rotate_xy(x, y, theta)
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(inout) :: y
    real(kind=rp), intent(in) :: theta
    real(kind=rp) :: x_old
    real(kind=rp) :: y_old
    real(kind=rp) :: cos_theta
    real(kind=rp) :: sin_theta

    x_old = x
    y_old = y
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    x = cos_theta * x_old - sin_theta * y_old
    y = sin_theta * x_old + cos_theta * y_old
  end subroutine lpt_rotate_xy

  !> Extract one translational periodic slot from scalar storage.
  !! @param idx Direction slot to extract.
  !! @param dir Periodic direction.
  !! @param periodic_min Minimum projected coordinate.
  !! @param periodic_max Maximum projected coordinate.
  !! @param shift Periodic translation vector.
  !! @param periodic_len Periodic length along `dir`.
  pure subroutine lpt_periodic_bc_get_translational_slot(idx, &
       periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, periodic_dir_x2, &
       periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, periodic_dir_y3, &
       periodic_dir_z3, periodic_min1, periodic_min2, periodic_min3, &
       periodic_max1, periodic_max2, periodic_max3, periodic_shift_x1, &
       periodic_shift_y1, periodic_shift_z1, periodic_shift_x2, &
       periodic_shift_y2, periodic_shift_z2, periodic_shift_x3, &
       periodic_shift_y3, periodic_shift_z3, periodic_len1, periodic_len2, &
       periodic_len3, dir, periodic_min, periodic_max, shift, periodic_len)
    integer, intent(in) :: idx
    real(kind=rp), intent(in) :: periodic_dir_x1, periodic_dir_y1
    real(kind=rp), intent(in) :: periodic_dir_z1, periodic_dir_x2
    real(kind=rp), intent(in) :: periodic_dir_y2, periodic_dir_z2
    real(kind=rp), intent(in) :: periodic_dir_x3, periodic_dir_y3
    real(kind=rp), intent(in) :: periodic_dir_z3
    real(kind=rp), intent(in) :: periodic_min1, periodic_min2, periodic_min3
    real(kind=rp), intent(in) :: periodic_max1, periodic_max2, periodic_max3
    real(kind=rp), intent(in) :: periodic_shift_x1, periodic_shift_y1
    real(kind=rp), intent(in) :: periodic_shift_z1, periodic_shift_x2
    real(kind=rp), intent(in) :: periodic_shift_y2, periodic_shift_z2
    real(kind=rp), intent(in) :: periodic_shift_x3, periodic_shift_y3
    real(kind=rp), intent(in) :: periodic_shift_z3
    real(kind=rp), intent(in) :: periodic_len1, periodic_len2, periodic_len3
    real(kind=rp), intent(out) :: dir(3)
    real(kind=rp), intent(out) :: periodic_min
    real(kind=rp), intent(out) :: periodic_max
    real(kind=rp), intent(out) :: shift(3)
    real(kind=rp), intent(out) :: periodic_len

    select case (idx)
    case (1)
       dir = [periodic_dir_x1, periodic_dir_y1, periodic_dir_z1]
       periodic_min = periodic_min1
       periodic_max = periodic_max1
       shift = [periodic_shift_x1, periodic_shift_y1, periodic_shift_z1]
       periodic_len = periodic_len1
    case (2)
       dir = [periodic_dir_x2, periodic_dir_y2, periodic_dir_z2]
       periodic_min = periodic_min2
       periodic_max = periodic_max2
       shift = [periodic_shift_x2, periodic_shift_y2, periodic_shift_z2]
       periodic_len = periodic_len2
    case (3)
       dir = [periodic_dir_x3, periodic_dir_y3, periodic_dir_z3]
       periodic_min = periodic_min3
       periodic_max = periodic_max3
       shift = [periodic_shift_x3, periodic_shift_y3, periodic_shift_z3]
       periodic_len = periodic_len3
    case default
       dir = 0.0_rp
       periodic_min = 0.0_rp
       periodic_max = 0.0_rp
       shift = 0.0_rp
       periodic_len = 0.0_rp
    end select
  end subroutine lpt_periodic_bc_get_translational_slot

end module lpt_periodic_bc_cpu
