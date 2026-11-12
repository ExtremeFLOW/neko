! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
! Copyright (c) 2025-2026, The Neko Authors
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
module fluid_volflow
  use operators, only : opgrad, cdtp, rotate_cyc
  use num_types, only : rp
  use mathops, only : opchsign
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use coefs, only : coef_t
  use time_state, only : time_state_t
  use time_scheme_controller, only : time_scheme_controller_t
  use math, only : copy, glsc2, glsc3, glmin, glmax, add2, addcol3, cmult, &
       abscmp, NEKO_EPS
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math, only : device_cfill, device_rzero, device_copy, &
       device_add2, device_add2s2, device_glsc2, device_glsc3, &
       device_addcol3, device_cmult
  use device, only : HOST_TO_DEVICE
  use utils, only : neko_error
  use device_mathops, only : device_opchsign
  use gather_scatter, only : gs_t, GS_OP_ADD
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_get_or_lookup
  use scratch_registry, only : neko_scratch_registry
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t
  use ax_product, only : ax_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Flow direction value selecting the azimuthal flow about the cyclic axis
  integer, public, parameter :: VOLFLOW_DIR_CYCLIC = 4

  !> Defines volume flow
  !! @details For flow_dir = 1, 2, 3 the flow is driven along x, y or z by a
  !! uniform force, and the flow rate is \f$ \int u_i dV / L_i \f$ with
  !! \f$ L_i \f$ the extent of the domain. For flow_dir = VOLFLOW_DIR_CYCLIC
  !! the domain is a sector of angle \f$ \Theta \f$ about the cyclic axis,
  !! e.g. a bent pipe. The force is \f$ g = e_\theta / r \f$, which is
  !! the gradient of the angle and therefore equivalent to a uniform
  !! pressure drop per radian, and the flow rate through a cross-section is
  !! \f$ \int u \cdot g dV / \Theta \f$.
  type, public :: fluid_volflow_t
     integer :: flow_dir !< these two should be moved to params
     logical :: avflow
     logical :: log = .true.
     real(kind=rp) :: flow_rate
     real(kind=rp) :: dtlag = 0d0
     real(kind=rp) :: bdlag = 0d0 !< Pointless since we do not vary the timestep
     type(field_t) :: u_vol, v_vol, w_vol, p_vol
     !> Forcing direction \f$ e_\theta / r \f$ for the cyclic case
     type(field_t) :: g_x, g_y, g_z
     !> Extent of the domain in the flow direction, the sector angle for the
     !! cyclic case
     real(kind=rp) :: domain_length
     !> Cross-section area, used with avflow
     real(kind=rp) :: xsec
     !> coef%metrics_version the cyclic geometry was built for, -1 if never
     integer :: geom_version = -1
     real(kind=rp) :: base_flow
   contains
     procedure, pass(this) :: init => fluid_vol_flow_init
     procedure, pass(this) :: free => fluid_vol_flow_free
     procedure, pass(this) :: adjust => fluid_vol_flow
     procedure, private, pass(this) :: compute => fluid_vol_flow_compute
  end type fluid_volflow_t

contains

  subroutine fluid_vol_flow_init(this, dm_Xh, params)
    class(fluid_volflow_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dm_Xh
    type(json_file), intent(inout) :: params
    logical average, log_output, cyclic
    integer :: direction
    real(kind=rp) :: rate

    call this%free()

    !Initialize vol_flow (if there is a forced volume flow)
    call json_get_or_lookup(params, 'case.fluid.flow_rate_force.direction', &
         direction)
    call json_get_or_lookup(params, 'case.fluid.flow_rate_force.value', rate)
    call json_get(params, 'case.fluid.flow_rate_force.use_averaged_flow',&
         average)
    call json_get_or_default(params, 'case.fluid.flow_rate_force.log', &
         log_output, .true.)

    this%flow_dir = direction
    this%avflow = average
    this%log = log_output
    this%flow_rate = rate

    if (this%flow_dir .lt. 0 .or. this%flow_dir .gt. VOLFLOW_DIR_CYCLIC) then
       call neko_error("flow_rate_force.direction must be 0, 1, 2, 3 or 4")
    end if

    if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC) then
       call json_get_or_default(params, 'case.fluid.cyclic', cyclic, .false.)
       if (.not. cyclic) then
          call neko_error("flow_rate_force.direction = 4 (cyclic) " // &
               "requires case.fluid.cyclic = true")
       end if
       call this%g_x%init(dm_Xh, 'volflow_g_x')
       call this%g_y%init(dm_Xh, 'volflow_g_y')
       call this%g_z%init(dm_Xh, 'volflow_g_z')
    end if

    if (this%flow_dir .ne. 0) then
       call this%u_vol%init(dm_Xh, 'u_vol')
       call this%v_vol%init(dm_Xh, 'v_vol')
       call this%w_vol%init(dm_Xh, 'w_vol')
       call this%p_vol%init(dm_Xh, 'p_vol')
    end if

  end subroutine fluid_vol_flow_init

  subroutine fluid_vol_flow_free(this)
    class(fluid_volflow_t), intent(inout) :: this

    call this%u_vol%free()
    call this%v_vol%free()
    call this%w_vol%free()
    call this%p_vol%free()
    call this%g_x%free()
    call this%g_y%free()
    call this%g_z%free()
    this%geom_version = -1

  end subroutine fluid_vol_flow_free

  !> Build the cyclic forcing direction \f$ g = e_\theta / r =
  !! a \times (x - x_0) / |r_\perp|^2 \f$ and set the sector angle as the
  !! domain length and \f$ \int r^{-1} dV / \Theta \f$ as the cross-section.
  subroutine fluid_vol_flow_cyclic_geometry(this, c_Xh)
    class(fluid_volflow_t), intent(inout) :: this
    type(coef_t), intent(in) :: c_Xh
    real(kind=rp) :: a(3), d(3), t(3), r2, sum_inv_r
    integer :: i, n, ierr

    n = c_Xh%dof%size()
    a = c_Xh%cyc_axis

    sum_inv_r = 0.0_rp
    do i = 1, n
       d = [c_Xh%dof%x%x(i,1,1,1), c_Xh%dof%y%x(i,1,1,1), &
            c_Xh%dof%z%x(i,1,1,1)] - c_Xh%cyc_origin
       t = [a(2) * d(3) - a(3) * d(2), &
            a(3) * d(1) - a(1) * d(3), &
            a(1) * d(2) - a(2) * d(1)]
       ! |a x d| = |r_perp|, since a is a unit vector
       r2 = dot_product(t, t)
       if (r2 .le. NEKO_EPS) then
          call neko_error("Cyclic flow rate forcing is singular on the " // &
               "rotation axis, which must lie outside the domain.")
       end if
       this%g_x%x(i,1,1,1) = t(1) / r2
       this%g_y%x(i,1,1,1) = t(2) / r2
       this%g_z%x(i,1,1,1) = t(3) / r2
       sum_inv_r = sum_inv_r + c_Xh%B(i,1,1,1) / sqrt(r2)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, sum_inv_r, 1, MPI_REAL_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%g_x%copy_from(HOST_TO_DEVICE, .false.)
       call this%g_y%copy_from(HOST_TO_DEVICE, .false.)
       call this%g_z%copy_from(HOST_TO_DEVICE, .true.)
    end if

    this%domain_length = c_Xh%cyc_angle
    this%xsec = sum_inv_r / c_Xh%cyc_angle

  end subroutine fluid_vol_flow_cyclic_geometry

  !> Flow rate of (u, v, w) in the forced direction
  function fluid_vol_flow_rate(this, u, v, w, c_Xh) result(rate)
    class(fluid_volflow_t), intent(in) :: this
    type(field_t), intent(in) :: u, v, w
    type(coef_t), intent(in) :: c_Xh
    real(kind=rp) :: rate
    integer :: n

    n = c_Xh%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       select case (this%flow_dir)
       case (1)
          rate = device_glsc2(u%x_d, c_Xh%B_d, n)
       case (2)
          rate = device_glsc2(v%x_d, c_Xh%B_d, n)
       case (3)
          rate = device_glsc2(w%x_d, c_Xh%B_d, n)
       case (VOLFLOW_DIR_CYCLIC)
          rate = device_glsc3(u%x_d, c_Xh%B_d, this%g_x%x_d, n) &
               + device_glsc3(v%x_d, c_Xh%B_d, this%g_y%x_d, n) &
               + device_glsc3(w%x_d, c_Xh%B_d, this%g_z%x_d, n)
       end select
    else
       select case (this%flow_dir)
       case (1)
          rate = glsc2(u%x, c_Xh%B, n)
       case (2)
          rate = glsc2(v%x, c_Xh%B, n)
       case (3)
          rate = glsc2(w%x, c_Xh%B, n)
       case (VOLFLOW_DIR_CYCLIC)
          rate = glsc3(u%x, c_Xh%B, this%g_x%x, n) &
               + glsc3(v%x, c_Xh%B, this%g_y%x, n) &
               + glsc3(w%x, c_Xh%B, this%g_z%x, n)
       end select
    end if

    rate = rate / this%domain_length

  end function fluid_vol_flow_rate

  !> Compute flow adjustment
  !! @brief Compute pressure and velocity using fractional step method.
  !! (Tombo splitting scheme).
  subroutine fluid_vol_flow_compute(this, u_res, v_res, w_res, p_res, &
       ext_bdf, gs_Xh, c_Xh, rho, mu, bd, dt, &
       dp_projector, vel_projector, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, &
       vel_max_iter)
    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(scalar_bc_projector_t), intent(inout) :: dp_projector
    class(vector_bc_projector_t), intent(inout) :: vel_projector
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: rho, dt
    type(field_t) :: mu
    integer, intent(in) :: vel_max_iter, prs_max_iter
    integer :: n, i
    real(kind=rp) :: xlmin, xlmax
    real(kind=rp) :: ylmin, ylmax
    real(kind=rp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_results(4)
    type(field_t), pointer :: ta1, ta2, ta3, ta4
    integer :: temp_indices(4)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(ta4, temp_indices(4), .false.)


    associate(msh => c_Xh%msh, p_vol => this%p_vol, &
         u_vol => this%u_vol, v_vol => this%v_vol, w_vol => this%w_vol)

      n = c_Xh%dof%size()
      xlmin = glmin(c_Xh%dof%x%x, n)
      xlmax = glmax(c_Xh%dof%x%x, n)
      ylmin = glmin(c_Xh%dof%y%x, n) !  for Y!
      ylmax = glmax(c_Xh%dof%y%x, n)
      zlmin = glmin(c_Xh%dof%z%x, n) !  for Z!
      zlmax = glmax(c_Xh%dof%z%x, n)
      if (this%flow_dir .eq. 1) then
         this%domain_length = xlmax - xlmin
      end if
      if (this%flow_dir .eq. 2) then
         this%domain_length = ylmax - ylmin
      end if
      if (this%flow_dir .eq. 3) then
         this%domain_length = zlmax - zlmin
      end if
      ! The geometry only changes when the mesh moves
      if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC .and. &
           this%geom_version .ne. c_Xh%metrics_version) then
         call fluid_vol_flow_cyclic_geometry(this, c_Xh)
         this%geom_version = c_Xh%metrics_version
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, 1.0_rp/rho, n)
         call device_rzero(c_Xh%h2_d, n)
      else
         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp parallel do
         do i = 1, n
            c_Xh%h1(i,1,1,1) = 1.0_rp / rho
            c_Xh%h2(i,1,1,1) = 0.0_rp
         end do
         !$omp end parallel do
      end if
      c_Xh%ifh2 = .false.

      !   Compute pressure

      if (this%flow_dir .eq. 1) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
      end if

      if (this%flow_dir .eq. 2) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
      end if

      if (this%flow_dir .eq. 3) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
      end if

      if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC) then
         ! Weak divergence of g / rho, component by component
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_copy(ta1%x_d, this%g_x%x_d, n)
            call device_copy(ta2%x_d, this%g_y%x_d, n)
            call device_copy(ta3%x_d, this%g_z%x_d, n)
            call device_cmult(ta1%x_d, 1.0_rp / rho, n)
            call device_cmult(ta2%x_d, 1.0_rp / rho, n)
            call device_cmult(ta3%x_d, 1.0_rp / rho, n)
         else
            call copy(ta1%x, this%g_x%x, n)
            call copy(ta2%x, this%g_y%x, n)
            call copy(ta3%x, this%g_z%x, n)
            call cmult(ta1%x, 1.0_rp / rho, n)
            call cmult(ta2%x, 1.0_rp / rho, n)
            call cmult(ta3%x, 1.0_rp / rho, n)
         end if
         call cdtp(p_res%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
         call cdtp(ta4%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_add2(p_res%x_d, ta4%x_d, n)
         else
            call add2(p_res%x, ta4%x, n)
         end if
         call cdtp(ta4%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_add2(p_res%x_d, ta4%x_d, n)
         else
            call add2(p_res%x, ta4%x, n)
         end if
      end if

      call gs_Xh%op(p_res, GS_OP_ADD)
      call dp_projector%apply(p_res%x, n)
      call pc_prs%update()
      ksp_results(1) = ksp_prs%solve(Ax_prs, p_vol, p_res%x, n, &
           c_Xh, dp_projector, gs_Xh, prs_max_iter)

      !   Compute velocity

      call opgrad(u_res%x, v_res%x, w_res%x, p_vol%x, c_Xh)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opchsign(u_res%x_d, v_res%x_d, w_res%x_d, msh%gdim, n)
         call device_copy(ta1%x_d, c_Xh%B_d, n)
         call device_copy(ta2%x_d, c_Xh%B_d, n)
         call device_copy(ta3%x_d, c_Xh%B_d, n)
      else
         call opchsign(u_res%x, v_res%x, w_res%x, msh%gdim, n)
         call copy(ta1%x, c_Xh%B, n)
         call copy(ta2%x, c_Xh%B, n)
         call copy(ta3%x, c_Xh%B, n)
      end if
      call vel_projector%apply(ta1%x, ta2%x, ta3%x, n)

      ! add forcing

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            call device_add2(u_res%x_d, ta1%x_d, n)
         else if (this%flow_dir .eq. 2) then
            call device_add2(v_res%x_d, ta2%x_d, n)
         else if (this%flow_dir .eq. 3) then
            call device_add2(w_res%x_d, ta3%x_d, n)
         else if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC) then
            call device_addcol3(u_res%x_d, ta1%x_d, this%g_x%x_d, n)
            call device_addcol3(v_res%x_d, ta2%x_d, this%g_y%x_d, n)
            call device_addcol3(w_res%x_d, ta3%x_d, this%g_z%x_d, n)
         end if
      else
         if (this%flow_dir .eq. 1) then
            call add2(u_res%x, ta1%x, n)
         else if (this%flow_dir .eq. 2) then
            call add2(v_res%x, ta2%x, n)
         else if (this%flow_dir .eq. 3) then
            call add2(w_res%x, ta3%x, n)
         else if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC) then
            call addcol3(u_res%x, ta1%x, this%g_x%x, n)
            call addcol3(v_res%x, ta2%x, this%g_y%x, n)
            call addcol3(w_res%x, ta3%x, this%g_z%x, n)
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_copy(c_Xh%h1_d, mu%x_d, n)
         call device_cfill(c_Xh%h2_d, rho * (bd / dt), n)
      else
         call copy(c_Xh%h1, mu%x, n)
         c_Xh%h2 = rho * (bd / dt)
      end if
      c_Xh%ifh2 = .true.

      call rotate_cyc(u_res, v_res, w_res, 1, c_Xh)
      call gs_Xh%op(u_res%x, v_res%x, w_res%x, n, GS_OP_ADD)
      call rotate_cyc(u_res, v_res, w_res, 0, c_Xh)

      call vel_projector%apply(u_res%x, v_res%x, w_res%x, n)
      call pc_vel%update()

      ksp_results(2:4) = ksp_vel%solve_coupled(Ax_vel, &
           u_vol, v_vol, w_vol, &
           u_res%x, v_res%x, w_res%x, &
           n, c_Xh, &
           vel_projector, &
           gs_Xh, vel_max_iter)

      this%base_flow = fluid_vol_flow_rate(this, u_vol, v_vol, w_vol, c_Xh)
    end associate

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine fluid_vol_flow_compute

  !> Adjust flow volume
  !! @brief  Adjust flow volume at end of time step to keep flow rate fixed by
  !! adding an appropriate multiple of the linear solution to the Stokes
  !! problem arising from a unit forcing in the X-direction.  This assumes
  !! that the flow rate in the X-direction is to be fixed (as opposed to Y-
  !! or Z-) *and* that the periodic boundary conditions in the X-direction
  !! occur at the extreme left and right ends of the mesh.
  !!
  !! pff 6/28/98
  subroutine fluid_vol_flow(this, u, v, w, p, u_res, v_res, w_res, p_res, &
       c_Xh, gs_Xh, ext_bdf, rho, mu, dt, time, &
       dp_projector, vel_projector, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, &
       vel_max_iter)

    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w, p
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(in) :: rho, dt
    type(field_t) :: mu
    type(scalar_bc_projector_t), intent(inout) :: dp_projector
    class(vector_bc_projector_t), intent(inout) :: vel_projector
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    integer, intent(in) :: prs_max_iter, vel_max_iter
    real(kind=rp) :: ifcomp, flow_rate
    real(kind=rp) :: current_flow, delta_flow, scale
    integer :: n, ierr, i
    character(len=5) :: flow_dir_label
    character(len=12) :: step_str

    character(len=200) :: log_buf

    associate(u_vol => this%u_vol, v_vol => this%v_vol, &
         w_vol => this%w_vol, p_vol => this%p_vol)

      n = c_Xh%dof%size()

      ! If either dt or the backwards difference coefficient change,
      ! then recompute base flow solution corresponding to unit forcing:

      ifcomp = 0.0_rp

      if ((.not. abscmp(dt, this%dtlag)) .or. &
           (.not. abscmp(ext_bdf%diffusion_coeffs%x(1), this%bdlag))) then
         ifcomp = 1.0_rp
      end if

      this%dtlag = dt
      this%bdlag = ext_bdf%diffusion_coeffs%x(1)

      call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      if (ifcomp .gt. 0d0) then
         call this%compute(u_res, v_res, w_res, p_res, &
              ext_bdf, gs_Xh, c_Xh, rho, mu, ext_bdf%diffusion_coeffs%x(1), dt, &
              dp_projector, vel_projector, &
              Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, &
              vel_max_iter)
      end if

      current_flow = fluid_vol_flow_rate(this, u, v, w, c_Xh)

      if (this%avflow) then
         if (this%flow_dir .ne. VOLFLOW_DIR_CYCLIC) then
            this%xsec = c_Xh%volume / this%domain_length
         end if
         flow_rate = this%flow_rate * this%xsec
      else
         flow_rate = this%flow_rate
      end if

      delta_flow = flow_rate - current_flow
      scale = delta_flow / this%base_flow

      if (this%log .and. pe_rank .eq. 0) then
         if (this%flow_dir .eq. 1) then
            flow_dir_label = '    x'
         else if (this%flow_dir .eq. 2) then
            flow_dir_label = '    y'
         else if (this%flow_dir .eq. 3) then
            flow_dir_label = '    z'
         else if (this%flow_dir .eq. VOLFLOW_DIR_CYCLIC) then
            flow_dir_label = '  cyc'
         end if
         write(step_str, '(I12)') time%tstep
         step_str = adjustl(step_str)
         write(log_buf, '(A,A3,A5,1X,5A18)') &
              'Flow rate   ', ' | ', 'Dir.:', 'Time:', 'Scale:', 'Rate:', &
              'Current:', 'Base:'
         call neko_log%message(log_buf)
         write(log_buf, '(A12,A3,A5,1X,5E18.9)') step_str, ' | ', &
              flow_dir_label, time%t, scale, flow_rate, current_flow, &
              this%base_flow
         call neko_log%message(log_buf)
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(u%x_d, u_vol%x_d, scale, n)
         call device_add2s2(v%x_d, v_vol%x_d, scale, n)
         call device_add2s2(w%x_d, w_vol%x_d, scale, n)
         call device_add2s2(p%x_d, p_vol%x_d, scale, n)
      else
         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp parallel do
         do i = 1, n
            u%x(i,1,1,1) = u%x(i,1,1,1) + scale * u_vol%x(i,1,1,1)
            v%x(i,1,1,1) = v%x(i,1,1,1) + scale * v_vol%x(i,1,1,1)
            w%x(i,1,1,1) = w%x(i,1,1,1) + scale * w_vol%x(i,1,1,1)
            p%x(i,1,1,1) = p%x(i,1,1,1) + scale * p_vol%x(i,1,1,1)
         end do
         !$omp end parallel do
      end if
    end associate

  end subroutine fluid_vol_flow

end module fluid_volflow
