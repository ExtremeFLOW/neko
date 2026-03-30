! Copyright (c) 2025-2026 The Neko Authors
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

!> Defines data structures and algorithms for configuring, calculating,
!> and time-integrating the rigid-body motion (translation and rotation)
!> of objects in an ALE simulation.
!> CPU-only module.
module ale_rigid_kinematics
  use num_types, only : rp
  use math, only : pi, math_dstepf
  use time_state, only : time_state_t
  use ab_time_scheme, only : ab_time_scheme_t
  use math, only : rzero
  use logger, only : neko_log

  implicit none
  private

  public :: compute_body_kinematics_built_in
  public :: update_pivot_location
  public :: init_pivot_state
  public :: advance_point_tracker
  public :: ab_integrate_point_pos

  !> Stiff Geometry
  type, public :: stiffness_geometry_t
     !> 'sphere', 'cylinder', 'cheap_dist', 'box (not implemented)'
     character(len=32) :: type = 'cheap_dist'
     !> 'gaussian' or 'tanh'
     character(len=32) :: decay_profile = 'gaussian'
     !> Stiffness Gain
     real(kind=rp) :: gain
     !> Center for geometric shapes
     real(kind=rp) :: center(3)
     !> Box dimensions (ToDo)
     real(kind=rp) :: dims(3) = 1.0_rp
     !> Box rotation angle CW (ToDo)
     real(kind=rp) :: rot_angle_deg_cw = 0.0_rp
     !> Acts as stiffness_decay for shapes
     real(kind=rp) :: radius
     !> Acts as stiffness_decay for cheap_dist
     real(kind=rp) :: stiff_dist
     !> Coefficient defining the cutoff sharpness (Default: 9.0 for Gaussian,
     !> 3.5 for Tanh)
     real(kind=rp) :: cutoff_coef
  end type stiffness_geometry_t

  !> Configuration for a single moving body
  type, public :: ale_body_t
     integer :: id
     character(len=32) :: name = 'body'
     type(stiffness_geometry_t) :: stiff_geom
     !> Oscillation (x, y, z)
     real(kind=rp) :: osc_amp(3) = 0.0_rp
     real(kind=rp) :: osc_freq(3) = 0.0_rp
     !> Rotation Control
     character(len=32) :: rotation_center_type = 'relative'
     character(len=32) :: rotation_type
     real(kind=rp) :: rot_amp_degree(3) = 0.0_rp
     real(kind=rp) :: rot_freq(3) = 0.0_rp
     real(kind=rp) :: rot_center(3) = 0.0_rp
     !> Smooth Step Rotation Parameters
     !> Single vector of 4 times (applied only to rotation_axis)
     !> 1=t0, 2=t1, 3=t2, 4=t3
     real(kind=rp) :: step_control_times(4) = 0.0_rp
     real(kind=rp) :: target_rot_angle_deg = 0.0_rp
     integer :: rotation_axis = 3
     !> Ramp Parameters
     real(kind=rp) :: ramp_omega0(3) = 0.0_rp
     real(kind=rp) :: ramp_t0(3) = 1.0_rp
     !> Boundary Zone Indices associated with this body
     integer, allocatable :: zone_indices(:)
  end type ale_body_t

  !> Global ALE Configuration
  type, public :: ale_config_t
     !> Stiffness Control.
     !> A placeholder for future options.
     character(len=32) :: stiffness_type = 'built-in'
     logical :: if_output_phi = .true.
     logical :: if_output_stiffness = .false.
     !> Array of Moving Bodies
     type(ale_body_t), allocatable :: bodies(:)
     integer :: nbodies = 0
  end type ale_config_t

  !> Calculated Kinematics for a body at current time
  type, public :: body_kinematics_t
     real(kind=rp) :: vel_trans(3) = 0.0_rp
     real(kind=rp) :: vel_ang(3) = 0.0_rp
     real(kind=rp) :: center(3) = 0.0_rp
  end type body_kinematics_t

  !> State history for time-integration of pivots
  type, public :: pivot_state_t
     real(kind=rp) :: pos(3) = 0.0_rp
     real(kind=rp) :: vel_lag(3, 3) = 0.0_rp
     real(kind=rp) :: vel(3) = 0.0_rp

  end type pivot_state_t

  !> Type for a tracked point linked to a body
  type, public :: point_tracker_t
     real(kind=rp) :: pos(3) = 0.0_rp
     real(kind=rp) :: vel_lag(3, 3) = 0.0_rp
     integer :: body_id = 0
  end type point_tracker_t

contains

  !> Initialize pivot state
  subroutine init_pivot_state(pivot, body_conf)
    type(pivot_state_t), intent(out) :: pivot
    type(ale_body_t), intent(in) :: body_conf
    pivot%pos = body_conf%rot_center
    pivot%vel_lag = 0.0_rp
    pivot%vel = 0.0_rp
  end subroutine init_pivot_state

  !> Advance a single point position (x,y,z) from the point's velocity
  !> using AB time-integration
  subroutine ab_integrate_point_pos(pos, vel_lag, current_vel, time, nadv)
    real(kind=rp), intent(inout) :: pos(3)
    real(kind=rp), intent(inout) :: vel_lag(3, 3)
    real(kind=rp), intent(in) :: current_vel(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv

    type(ab_time_scheme_t) :: ab_scheme_obj
    real(kind=rp) :: ab_coeffs(4), dt_history(10)
    integer :: j

    call rzero(ab_coeffs, 4)
    dt_history(1) = time%dt
    dt_history(2) = time%dtlag(1)
    dt_history(3) = time%dtlag(2)
    call ab_scheme_obj%compute_coeffs(ab_coeffs, dt_history, nadv)

    pos(1) = pos(1) + time%dt * ab_coeffs(1) * current_vel(1)
    pos(2) = pos(2) + time%dt * ab_coeffs(1) * current_vel(2)
    pos(3) = pos(3) + time%dt * ab_coeffs(1) * current_vel(3)

    do j = 2, nadv
       pos(1) = pos(1) + time%dt * ab_coeffs(j) * vel_lag(1, j - 1)
       pos(2) = pos(2) + time%dt * ab_coeffs(j) * vel_lag(2, j - 1)
       pos(3) = pos(3) + time%dt * ab_coeffs(j) * vel_lag(3, j - 1)
    end do

    ! Shift history
    do j = 3, 2, -1
       vel_lag(1, j) = vel_lag(1, j - 1)
       vel_lag(2, j) = vel_lag(2, j - 1)
       vel_lag(3, j) = vel_lag(3, j - 1)
    end do

    ! Store current velocity
    vel_lag(1, 1) = current_vel(1)
    vel_lag(2, 1) = current_vel(2)
    vel_lag(3, 1) = current_vel(3)

  end subroutine ab_integrate_point_pos

  !> Updates the point tracker's position and velocity history using
  !> AB time integration based on the current velocity.
  subroutine advance_point_tracker(tracker, current_vel, time, nadv)
    type(point_tracker_t), intent(inout) :: tracker
    real(kind=rp), intent(in) :: current_vel(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv

    call ab_integrate_point_pos(tracker%pos, tracker%vel_lag, &
          current_vel, time, nadv)

  end subroutine advance_point_tracker

  !> Compute built-in kinematics for a body. Uses inputs from JSON.
  !> CPU-only.
  subroutine compute_body_kinematics_built_in(kinematics, body_conf, time)
    type(body_kinematics_t), intent(out) :: kinematics
    type(ale_body_t), intent(in) :: body_conf
    type(time_state_t), intent(in) :: time

    real(kind=rp) :: w_scalar, amp_scalar, ex_term
    real(kind=rp) :: t_curr, t0, t1, t2, t3, target_rad, tau, current_omega
    integer :: i

    ! Oscillation
    kinematics%vel_trans = 0.0_rp
    do i = 1, 3
       if (abs(body_conf%osc_amp(i)) > 0.0_rp) then
          w_scalar = body_conf%osc_freq(i) * 2.0_rp * pi
          kinematics%vel_trans(i) = body_conf%osc_amp(i) * w_scalar * &
                cos(w_scalar * time%t)
       end if
    end do

    ! Rotation
    kinematics%vel_ang = 0.0_rp
    select case (trim(body_conf%rotation_type))
    case ('smooth_step')
       if (abs(body_conf%target_rot_angle_deg) > 0.0_rp) then
          t_curr = time%t
          t0 = body_conf%step_control_times(1)
          t1 = body_conf%step_control_times(2)
          t2 = body_conf%step_control_times(3)
          t3 = body_conf%step_control_times(4)
          target_rad = body_conf%target_rot_angle_deg * (pi / 180.0_rp)
          current_omega = 0.0_rp

          if (t_curr < t0) then
             current_omega = 0.0_rp
          elseif (t_curr < t1) then ! Rise
             if (t1 > t0) then
                tau = (t_curr - t0) / (t1 - t0)
                current_omega = target_rad * math_dstepf(tau) * &
                      (1.0_rp / (t1 - t0))
             end if
          elseif (t_curr < t2) then ! Hold
             current_omega = 0.0_rp
          elseif (t_curr < t3) then ! Fall
             if (t3 > t2) then
                tau = (t_curr - t2) / (t3 - t2)
                current_omega = -1.0_rp * target_rad * math_dstepf(tau) * &
                      (1.0_rp / (t3 - t2))
             end if
          end if
          kinematics%vel_ang(body_conf%rotation_axis) = current_omega
       end if


    case ('ramp')
       do i = 1, 3
          if (body_conf%ramp_t0(i) > 0.0_rp .and. &
                abs(body_conf%ramp_omega0(i)) > 0.0_rp) then
             ex_term = exp(-4.6_rp * time%t / body_conf%ramp_t0(i))
             kinematics%vel_ang(i) = body_conf%ramp_omega0(i) * &
                   (1.0_rp - ex_term)
          end if
       end do

    case ('harmonic')
       do i = 1, 3
          if (abs(body_conf%rot_amp_degree(i)) > 0.0_rp) then
             amp_scalar = body_conf%rot_amp_degree(i) * (pi / 180.0_rp)
             w_scalar = body_conf%rot_freq(i) * 2.0_rp * pi
             kinematics%vel_ang(i) = w_scalar * amp_scalar * &
                   cos(w_scalar * time%t)
          end if
       end do
    end select

  end subroutine compute_body_kinematics_built_in

  !> Updates pivot location
  subroutine update_pivot_location(pivot, pivot_loc, pivot_vel, time, &
       nadv, body_conf)
    type(pivot_state_t), intent(inout) :: pivot
    real(kind=rp), intent(out) :: pivot_loc(3)
    real(kind=rp), intent(in) :: pivot_vel(3)
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    type(ale_body_t), intent(in) :: body_conf
    real(kind=rp) :: omega_osc(3), trans_disp(3)

    select case (trim(body_conf%rotation_center_type))

       ! This mode uses vel_trans to move the pivot by integrating the velocity
    case ('relative')

       if (time%tstep > 0) then
          call ab_integrate_point_pos(pivot%pos, pivot%vel_lag, &
               pivot_vel, time, nadv)
       end if
       pivot_loc = pivot%pos

       ! Mostly for validation. It's too restrictive.
       ! maybe I remove it totally in future.
    case ('relative_sin')
       omega_osc = body_conf%osc_freq * 2.0_rp * pi
       trans_disp(1) = body_conf%osc_amp(1) * sin(omega_osc(1) * time%t)
       trans_disp(2) = body_conf%osc_amp(2) * sin(omega_osc(2) * time%t)
       trans_disp(3) = body_conf%osc_amp(3) * sin(omega_osc(3) * time%t)
       pivot_loc = body_conf%rot_center + trans_disp
       pivot%pos = pivot_loc
    end select
  end subroutine update_pivot_location

end module ale_rigid_kinematics
