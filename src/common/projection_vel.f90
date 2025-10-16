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
!> Project x onto X , the space of old solutions and back again
!> Couple projections for velocity
module projection_vel
  use num_types, only : rp, c_rp
  use math, only : add2, copy
  use coefs, only : coef_t
  use ax_product, only : ax_t
  use bc_list, only : bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_get_ptr
  use device_math, only : device_copy, device_add2
  use profiler, only : profiler_start_region, profiler_end_region
  use, intrinsic :: iso_c_binding
  use time_step_controller, only : time_step_controller_t
  use projection, only : projection_t, proj_ortho

  implicit none
  private

  type, public :: projection_vel_t
     type(projection_t) :: proj_u, proj_v, proj_w
     integer :: activ_step
     integer :: L
   contains
     procedure, pass(this) :: init => projection_init_vel
     procedure, pass(this) :: free => projection_free_vel
     procedure, pass(this) :: pre_solving => projection_pre_solving_vel
     procedure, pass(this) :: post_solving => projection_post_solving_vel
     procedure, pass(this) :: project_back => bcknd_project_back_vel
  end type projection_vel_t

contains

  subroutine projection_init_vel(this, n, L, activ_step)
    class(projection_vel_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, optional, intent(in) :: L, activ_step
    integer :: i
    integer(c_size_t) :: ptr_size
    type(c_ptr) :: ptr
    real(c_rp) :: dummy

    call this%free()

    call this%proj_u%init(n, L, activ_step)
    call this%proj_v%init(n, L, activ_step)
    call this%proj_w%init(n, L, activ_step)

    this%L = L
    this%activ_step = activ_step


  end subroutine projection_init_vel

  subroutine projection_free_vel(this)
    class(projection_vel_t), intent(inout) :: this
    integer :: i

    call this%proj_u%free()
    call this%proj_v%free()
    call this%proj_w%free()

  end subroutine projection_free_vel

  subroutine projection_pre_solving_vel(this, b_u, b_v, b_w, tstep, coef, n, &
       dt_controller, &
       stringx, stringy, stringz)
    class(projection_vel_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: b_u, b_v, b_w
    integer, intent(in) :: tstep
    class(coef_t), intent(inout) :: coef
    type(time_step_controller_t), intent(in) :: dt_controller
    character(len=*), optional :: stringx, stringy, stringz

    call this%proj_u%pre_solving(b_u, tstep, coef, n, dt_controller, stringx)
    call this%proj_v%pre_solving(b_v, tstep, coef, n, dt_controller, stringy)
    call this%proj_w%pre_solving(b_w, tstep, coef, n, dt_controller, stringz)

  end subroutine projection_pre_solving_vel

  subroutine projection_post_solving_vel(this, x_u, x_v, x_w, Ax, coef, &
       bclst_u, bclst_v, bclst_w, &
       gs_h, n, tstep, dt_controller)
    class(projection_vel_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(ax_t), intent(inout) :: Ax
    class(coef_t), intent(inout) :: coef
    class(bc_list_t), intent(inout) :: bclst_u, bclst_v, bclst_w
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x_u, x_v, x_w
    integer, intent(in) :: tstep
    type(time_step_controller_t), intent(in) :: dt_controller

    ! Here we assume the projection space sizes and activate steps
    ! for all three velocity equations are the same
    if (tstep .gt. this%activ_step .and. this%L .gt. 0) then
       if ((.not. dt_controller%is_variable_dt) .or. &
            (dt_controller%dt_last_change .gt. this%activ_step - 1)) then
          call this%project_back(x_u, x_v, x_w, Ax, coef, bclst_u, bclst_v, &
               bclst_w, gs_h, n)
       end if
    end if

  end subroutine projection_post_solving_vel

  subroutine bcknd_project_back_vel(this, x_u, x_v, x_w, &
       Ax, coef, bclst_u, bclst_v, bclst_w, gs_h, n)
    class(projection_vel_t) :: this
    integer, intent(inout) :: n
    class(ax_t), intent(inout) :: Ax
    class(coef_t), intent(inout) :: coef
    class(bc_list_t), intent(inout) :: bclst_u, bclst_v, bclst_w
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x_u, x_v, x_w
    type(c_ptr) :: x_u_d, x_v_d, x_w_d

    call profiler_start_region('Project back', 17)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_u_d = device_get_ptr(x_u)
       x_v_d = device_get_ptr(x_v)
       x_w_d = device_get_ptr(x_w)
       if (this%proj_u%m .gt. 0) then ! Restore desired solution
          call device_add2(x_u_d, this%proj_u%xbar_d, n)
       end if
       if (this%proj_v%m .gt. 0) then ! Restore desired solution
          call device_add2(x_v_d, this%proj_v%xbar_d, n)
       end if
       if (this%proj_w%m .gt. 0) then ! Restore desired solution
          call device_add2(x_w_d, this%proj_w%xbar_d, n)
       end if

       if (this%proj_u%m .eq. this%proj_u%L) then
          this%proj_u%m = 1
       else
          this%proj_u%m = min(this%proj_u%m + 1, this%proj_u%L)
       end if
       if (this%proj_v%m .eq. this%proj_v%L) then
          this%proj_v%m = 1
       else
          this%proj_v%m = min(this%proj_v%m + 1, this%proj_v%L)
       end if
       if (this%proj_w%m .eq. this%proj_w%L) then
          this%proj_w%m = 1
       else
          this%proj_w%m = min(this%proj_w%m + 1, this%proj_w%L)
       end if

       call device_copy(this%proj_u%xx_d(this%proj_u%m), &
            x_u_d,n) ! Update (X,B)
       call device_copy(this%proj_v%xx_d(this%proj_u%m), &
            x_v_d,n) ! Update (X,B)
       call device_copy(this%proj_w%xx_d(this%proj_u%m), &
            x_w_d,n) ! Update (X,B)

    else
       if (this%proj_u%m .gt. 0) then
          call add2(x_u, this%proj_u%xbar, n) ! Restore desired solution
       end if
       if (this%proj_v%m .gt. 0) then
          call add2(x_v, this%proj_v%xbar, n) ! Restore desired solution
       end if
       if (this%proj_w%m .gt. 0) then
          call add2(x_w, this%proj_w%xbar, n) ! Restore desired solution
       end if

       if (this%proj_u%m .eq. this%proj_u%L) then
          this%proj_u%m = 1
       else
          this%proj_u%m = min(this%proj_u%m + 1, this%proj_u%L)
       end if
       if (this%proj_v%m .eq. this%proj_v%L) then
          this%proj_v%m = 1
       else
          this%proj_v%m = min(this%proj_v%m + 1, this%proj_v%L)
       end if
       if (this%proj_w%m .eq. this%proj_w%L) then
          this%proj_w%m = 1
       else
          this%proj_w%m = min(this%proj_w%m + 1, this%proj_w%L)
       end if

       call copy(this%proj_u%xx(1, this%proj_u%m), x_u, n) ! Update (X,B)
       call copy(this%proj_v%xx(1, this%proj_v%m), x_v, n) ! Update (X,B)
       call copy(this%proj_w%xx(1, this%proj_w%m), x_w, n) ! Update (X,B)
    end if

    call Ax%compute_vector(this%proj_u%bb(1, this%proj_u%m), &
         this%proj_v%bb(1, this%proj_v%m), &
         this%proj_w%bb(1, this%proj_w%m), x_u, x_v, x_w, &
         coef, coef%msh, coef%Xh)

    call gs_h%gs_op_vector(this%proj_u%bb(1, this%proj_u%m), n, GS_OP_ADD)
    call gs_h%gs_op_vector(this%proj_v%bb(1, this%proj_v%m), n, GS_OP_ADD)
    call gs_h%gs_op_vector(this%proj_w%bb(1, this%proj_w%m), n, GS_OP_ADD)

    call bclst_u%apply_scalar(this%proj_u%bb(1, this%proj_u%m), n)
    call bclst_v%apply_scalar(this%proj_v%bb(1, this%proj_v%m), n)
    call bclst_w%apply_scalar(this%proj_w%bb(1, this%proj_w%m), n)

    call proj_ortho(this%proj_u, coef, n)
    call proj_ortho(this%proj_v, coef, n)
    call proj_ortho(this%proj_w, coef, n)
    call profiler_end_region('Project back', 17)
  end subroutine bcknd_project_back_vel
end module projection_vel
