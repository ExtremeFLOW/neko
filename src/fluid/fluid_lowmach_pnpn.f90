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
!> Low-Mach Pn-Pn fluid scheme.
!!
!! Extends fluid_pnpn_t and overrides step() with a body copied verbatim
!! from fluid_pnpn_step. This override currently contains no low-Mach
!! physics and is expected to produce bit-identical output to the standard
!! Pn-Pn solver. It is the starting point for incremental addition of
!! variable-density momentum, Q_T source in the pressure Poisson equation,
!! and full-stress variable-viscosity viscous terms.
module fluid_lowmach_pnpn
  use num_types, only : rp
  use fluid_pnpn, only : fluid_pnpn_t
  use krylov, only : ksp_monitor_t
  use bc, only : bc_t
  use non_normal, only : non_normal_t
  use file, only : file_t
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use profiler, only : profiler_start_region, profiler_end_region
  use field_math, only : field_add2
  use operators, only : ortho, rotate_cyc
  use opr_device, only : device_ortho
  use device, only : device_event_sync, glb_cmd_event
  use device_mathops, only : device_opadd2cm
  use mathops, only : opadd2cm
  use fluid_aux, only : fluid_step_info
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  implicit none
  private

  type, public, extends(fluid_pnpn_t) :: fluid_lowmach_pnpn_t
   contains
     procedure, pass(this) :: step => fluid_lowmach_pnpn_step
  end type fluid_lowmach_pnpn_t

contains

  !> Advance the solution by one time step.
  !! Body copied verbatim from fluid_pnpn_step. Low-Mach physics will be
  !! introduced here in subsequent commits.
  subroutine fluid_lowmach_pnpn_step(this, time, dt_controller)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    ! number of degrees of freedom
    integer :: n
    ! Solver results monitors (pressure + 3 velocity)
    type(ksp_monitor_t) :: ksp_results(4)

    type(file_t) :: dump_file
    class(bc_t), pointer :: bc_i
    type(non_normal_t), pointer :: bc_j

    if (this%freeze) return

    n = this%dm_Xh%size()

    call profiler_start_region('Fluid', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_res => this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax_vel => this%Ax_vel, Ax_prs => this%Ax_prs, &
         Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, prs_res => this%prs_res, &
         source_term => this%source_term, vel_res => this%vel_res, &
         sumab => this%sumab, makeoifs => this%makeoifs, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         vel_projection_dim => this%vel_projection_dim, &
         pr_projection_dim => this%pr_projection_dim, &
         oifs => this%oifs, &
         rho => this%rho, mu_tot => this%mu_tot, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
         t => time%t, tstep => time%tstep, dt => time%dt, &
         ext_bdf => this%ext_bdf, event => glb_cmd_event)

      ! Extrapolate the velocity if it's not done in nut_field estimation
      call sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
           ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)

      ! Compute the source terms
      call this%source_term%compute(time)

      ! Add Neumann bc contributions to the RHS
      call this%bcs_vel%apply_vector(f_x%x, f_y%x, f_z%x, &
           this%dm_Xh%size(), time, strong = .false.)

      if (oifs) then
         ! Add the advection operators to the right-hand-side.
         call this%adv%compute(u, v, w, &
              this%advx, this%advy, this%advz, &
              Xh, this%c_Xh, dm_Xh%size(), dt)

         ! At this point the RHS contains the sum of the advection operator and
         ! additional source terms, evaluated using the velocity field from the
         ! previous time-step. Now, this value is used in the explicit time
         ! scheme to advance both terms in time.

         call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
              this%abx2, this%aby2, this%abz2, &
              f_x%x, f_y%x, f_z%x, &
              rho%x(1,1,1,1), ext_bdf%advection_coeffs%x, n)

         ! Now, the source terms from the previous time step are added to the RHS.
         call makeoifs%compute_fluid(this%advx%x, this%advy%x, this%advz%x, &
              f_x%x, f_y%x, f_z%x, &
              rho%x(1,1,1,1), dt, n)
      else
         ! Add the advection operators to the right-hand-side.
         call this%adv%compute(u, v, w, &
              f_x, f_y, f_z, &
              Xh, this%c_Xh, dm_Xh%size())

         ! At this point the RHS contains the sum of the advection operator and
         ! additional source terms, evaluated using the velocity field from the
         ! previous time-step. Now, this value is used in the explicit time
         ! scheme to advance both terms in time.
         call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
              this%abx2, this%aby2, this%abz2, &
              f_x%x, f_y%x, f_z%x, &
              rho%x(1,1,1,1), ext_bdf%advection_coeffs%x, n)

         ! Add the RHS contributions coming from the BDF scheme.
         call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
              u, v, w, c_Xh%B, rho%x(1,1,1,1), dt, &
              ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
      end if

      call ulag%update()
      call vlag%update()
      call wlag%update()

      call this%bc_apply_vel(time, strong = .true.)
      call this%bc_apply_prs(time)

      ! Update material properties if necessary
      call this%update_material_properties(time)

      ! Compute pressure residual.
      call profiler_start_region('Pressure_residual', 18)

      call prs_res%compute(p, p_res,&
           u, v, w, &
           u_e, v_e, w_e, &
           f_x, f_y, f_z, &
           c_Xh, gs_Xh, &
           this%bc_prs_surface, this%bc_sym_surface,&
           Ax_prs, ext_bdf%diffusion_coeffs%x(1), dt, &
           mu_tot, rho, event)

      ! De-mean the pressure residual when no strong pressure boundaries present
      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p_res%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p_res%x, this%glb_n_points, n)
      end if

      call gs_Xh%op(p_res, GS_OP_ADD, event)
      call device_event_sync(event)

      ! Set the residual to zero at strong pressure boundaries.
      call this%bclst_dp%apply_scalar(p_res%x, p%dof%size(), time)


      call profiler_end_region('Pressure_residual', 18)

      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
           Ax=Ax_prs, gs_h=gs_Xh, bclst=this%bclst_dp, string='Pressure')

      call this%pc_prs%update()

      call profiler_start_region('Pressure_solve', 3)

      ! Solve for the pressure increment.
      ksp_results(1) = &
           this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, &
           this%bclst_dp, gs_Xh)
      ksp_results(1)%name = 'Pressure'


      call profiler_end_region('Pressure_solve', 3)

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
           this%bclst_dp, gs_Xh, n, tstep, dt_controller)

      ! Update the pressure with the increment. Demean if necessary.
      call field_add2(p, dp, n)
      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p%x, this%glb_n_points, n)
      end if

      ! Compute velocity residual.
      call profiler_start_region('Velocity_residual', 19)
      call vel_res%compute(Ax_vel, u, v, w, &
           u_res, v_res, w_res, &
           p, &
           f_x, f_y, f_z, &
           c_Xh, msh, Xh, &
           mu_tot, rho, ext_bdf%diffusion_coeffs%x(1), &
           dt, dm_Xh%size())

      call rotate_cyc(u_res%x, v_res%x, w_res%x, 1, c_Xh)
      call gs_Xh%op(u_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(v_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(w_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 0, c_Xh)

      ! Set residual to zero at strong velocity boundaries.
      call this%bclst_vel_res%apply(u_res, v_res, w_res, time)


      call profiler_end_region('Velocity_residual', 19)

      call this%proj_vel%pre_solving(u_res%x, v_res%x, w_res%x, &
           tstep, c_Xh, n, dt_controller, 'Velocity')

      call this%pc_vel%update()

      call profiler_start_region("Velocity_solve", 4)
      ksp_results(2:4) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh, &
           this%ksp_vel%max_iter)
      call profiler_end_region("Velocity_solve", 4)
      if (this%full_stress_formulation) then
         ksp_results(2)%name = 'Momentum'
      else
         ksp_results(2)%name = 'X-Velocity'
         ksp_results(3)%name = 'Y-Velocity'
         ksp_results(4)%name = 'Z-Velocity'
      end if

      call this%proj_vel%post_solving(du%x, dv%x, dw%x, Ax_vel, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh, n, tstep, &
           dt_controller)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if

      if (this%forced_flow_rate) then
         ! Horrible mu hack?!
         call this%vol_flow%adjust( u, v, w, p, u_res, v_res, w_res, p_res, &
              c_Xh, gs_Xh, ext_bdf, rho%x(1,1,1,1), mu_tot, &
              dt, time, this%bclst_dp, this%bclst_du, this%bclst_dv, &
              this%bclst_dw, this%bclst_vel_res, Ax_vel, Ax_prs, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(time, ksp_results, &
           this%full_stress_formulation, this%strict_convergence, &
           this%allow_stabilization)

    end associate
    call profiler_end_region('Fluid', 1)
  end subroutine fluid_lowmach_pnpn_step

end module fluid_lowmach_pnpn
