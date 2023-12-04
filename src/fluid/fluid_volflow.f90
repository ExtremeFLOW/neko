! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
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
  use operators
  use num_types
  use mathops    
  use krylov, only : ksp_t, ksp_monitor_t
  use precon
  use dofmap
  use field
  use coefs
  use time_scheme_controller
  use math
  use comm
  use neko_config
  use device_math
  use device_mathops
  use gather_scatter, only : gs_t, GS_OP_ADD
  use json_module, only : json_file
  use json_utils, only: json_get
  use scratch_registry, only : scratch_registry_t
  use bc, only : bc_list_t, bc_list_apply, bc_list_apply_vector, &
                 bc_list_apply_scalar
  use ax_product, only : ax_t
  implicit none
  private
  
  !> Defines volume flow
  type, public :: fluid_volflow_t
     integer :: flow_dir !< these two should be moved to params
     logical :: avflow 
     real(kind=rp) :: flow_rate 
     real(kind=rp) :: dtlag = 0d0
     real(kind=rp) :: bdlag = 0d0 !< Really quite pointless since we do not vary the timestep
     type(field_t) :: u_vol, v_vol, w_vol, p_vol
     real(kind=rp) :: domain_length, base_flow
     !> Manager for temporary fields
     type(scratch_registry_t) :: scratch
   contains
     procedure, pass(this) :: init => fluid_vol_flow_init
     procedure, pass(this) :: free => fluid_vol_flow_free
     procedure, pass(this) :: adjust => fluid_vol_flow
     procedure, private, pass(this) :: compute => fluid_vol_flow_compute
  end type fluid_volflow_t

contains

  subroutine fluid_vol_flow_init(this, dm_Xh, params)
    class(fluid_volflow_t), intent(inout) :: this
    type(dofmap_t), intent(inout) :: dm_Xh
    type(json_file), intent(inout) :: params
    logical average, found
    integer :: direction
    real(kind=rp) :: rate

    call this%free()

    !Initialize vol_flow (if there is a forced volume flow)
    call json_get(params, 'case.fluid.flow_rate_force.direction', direction)
    call json_get(params, 'case.fluid.flow_rate_force.value', rate)
    call json_get(params, 'case.fluid.flow_rate_force.use_averaged_flow',&
                  average)

    this%flow_dir = direction
    this%avflow = average
    this%flow_rate = rate

    if (this%flow_dir .ne. 0) then
       call this%u_vol%init(dm_Xh, 'u_vol')
       call this%v_vol%init(dm_Xh, 'v_vol')
       call this%w_vol%init(dm_Xh, 'w_vol')
       call this%p_vol%init(dm_Xh, 'p_vol')
    end if

    this%scratch = scratch_registry_t(dm_Xh, 3, 1)
    
  end subroutine fluid_vol_flow_init

  subroutine fluid_vol_flow_free(this)
    class(fluid_volflow_t), intent(inout) :: this

    call this%u_vol%free()
    call this%v_vol%free()
    call this%w_vol%free()
    call this%p_vol%free()

    call this%scratch%free()
        
  end subroutine fluid_vol_flow_free

  !> Compute flow adjustment
  !! @brief Compute pressure and velocity using fractional step method.
  !! (Tombo splitting scheme).
  subroutine fluid_vol_flow_compute(this, u_res, v_res, w_res, p_res, &
       ext_bdf, gs_Xh, c_Xh, rho, mu, bd, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)
    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(inout) :: Ax
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    real(kind=rp), intent(inout) :: bd
    real(kind=rp), intent(in) :: rho, mu, dt
    integer, intent(in) :: vel_max_iter, prs_max_iter
    integer :: n, i
    real(kind=rp) :: xlmin, xlmax
    real(kind=rp) :: ylmin, ylmax
    real(kind=rp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_result
    type(field_t), pointer :: ta1, ta2, ta3
    integer :: temp_indices(3)
    
    call this%scratch%request_field(ta1, temp_indices(1))
    call this%scratch%request_field(ta2, temp_indices(2))
    call this%scratch%request_field(ta3, temp_indices(3))
    

    associate(msh => c_Xh%msh, p_vol => this%p_vol, &
         u_vol => this%u_vol, v_vol => this%v_vol, w_vol => this%w_vol)

      n = c_Xh%dof%size()
      xlmin = glmin(c_Xh%dof%x, n)
      xlmax = glmax(c_Xh%dof%x, n)
      ylmin = glmin(c_Xh%dof%y, n)          !  for Y!
      ylmax = glmax(c_Xh%dof%y, n)
      zlmin = glmin(c_Xh%dof%z, n)          !  for Z!
      zlmax = glmax(c_Xh%dof%z, n)
      if (this%flow_dir.eq.1) this%domain_length = xlmax - xlmin
      if (this%flow_dir.eq.2) this%domain_length = ylmax - ylmin
      if (this%flow_dir.eq.3) this%domain_length = zlmax - zlmin

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, 1.0_rp/rho, n)
         call device_rzero(c_Xh%h2_d, n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = 1.0_rp / rho
            c_Xh%h2(i,1,1,1) = 0.0_rp
         end do
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

      call gs_Xh%op(p_res, GS_OP_ADD) 
      call bc_list_apply_scalar(bclst_dp, p_res%x, n)
      call pc_prs%update()
      ksp_result = ksp_prs%solve(Ax, p_vol, p_res%x, n, &
           c_Xh, bclst_dp, gs_Xh, prs_max_iter)    

      !   Compute velocity
      
      call opgrad(u_res%x, v_res%x, w_res%x, p_vol%x, c_Xh)
         
      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
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
      call bc_list_apply_vector(bclst_vel_res,&
           ta1%x, ta2%x, ta3%x, n)
       
      ! add forcing

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            call device_add2(u_res%x_d, ta1%x_d, n) 
         else if (this%flow_dir .eq. 2) then
            call device_add2(v_res%x_d, ta2%x_d, n)
         else if (this%flow_dir .eq. 3) then
            call device_add2(w_res%x_d, ta3%x_d, n)
         end if
      else
         if (this%flow_dir .eq. 1) then
            call add2(u_res%x, ta1%x, n) 
         else if (this%flow_dir .eq. 2) then
            call add2(v_res%x, ta2%x, n)
         else if (this%flow_dir .eq. 3) then
            call add2(w_res%x, ta3%x, n)
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, mu, n)
         call device_cfill(c_Xh%h2_d, rho * (bd / dt), n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = mu
            c_Xh%h2(i,1,1,1) = rho * (bd / dt)
         end do
      end if
      c_Xh%ifh2 = .true.

       call gs_Xh%op(u_res, GS_OP_ADD) 
       call gs_Xh%op(v_res, GS_OP_ADD) 
       call gs_Xh%op(w_res, GS_OP_ADD) 

       call bc_list_apply_vector(bclst_vel_res,&
            u_res%x, v_res%x, w_res%x, n)
       call pc_vel%update()

       ksp_result = ksp_vel%solve(Ax, u_vol, u_res%x, n, &
            c_Xh, bclst_du, gs_Xh, vel_max_iter)
       ksp_result = ksp_vel%solve(Ax, v_vol, v_res%x, n, &
            c_Xh, bclst_dv, gs_Xh, vel_max_iter)
       ksp_result = ksp_vel%solve(Ax, w_vol, w_res%x, n, &
            c_Xh, bclst_dw, gs_Xh, vel_max_iter)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            this%base_flow = &
                 device_glsc2(u_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if
         
         if (this%flow_dir .eq. 2) then
            this%base_flow = &
                 device_glsc2(v_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if
         
         if (this%flow_dir .eq. 3) then
            this%base_flow = &
                 device_glsc2(w_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if
      else
         if (this%flow_dir .eq. 1) then
            this%base_flow = glsc2(u_vol%x, c_Xh%B, n) / this%domain_length
         end if
         
         if (this%flow_dir .eq. 2) then
            this%base_flow = glsc2(v_vol%x, c_Xh%B, n) / this%domain_length
         end if
         
         if (this%flow_dir .eq. 3) then
            this%base_flow = glsc2(w_vol%x, c_Xh%B, n) / this%domain_length
         end if
      end if
     end associate

    call this%scratch%relinquish_field(temp_indices)
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
       c_Xh, gs_Xh, ext_bdf, rho, mu, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)

    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w, p
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    real(kind=rp), intent(in) :: rho, mu, dt
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(inout) :: Ax
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    integer, intent(in) :: prs_max_iter, vel_max_iter
    real(kind=rp) :: ifcomp, flow_rate, xsec
    real(kind=rp) :: current_flow, delta_flow, base_flow, scale
    integer :: n, ierr
    type(field_t), pointer :: ta1, ta2, ta3
    integer :: temp_indices(3)
    
    associate(u_vol => this%u_vol, v_vol => this%v_vol, &
         w_vol => this%w_vol, p_vol => this%p_vol)
      
      n = c_Xh%dof%size()

      ! If either dt or the backwards difference coefficient change,
      ! then recompute base flow solution corresponding to unit forcing:
      
      ifcomp = 0.0_rp

      if (dt .ne. this%dtlag .or. ext_bdf%diffusion_coeffs(1) .ne. this%bdlag) then
         ifcomp = 1.0_rp
      end if
      
      this%dtlag = dt
      this%bdlag = ext_bdf%diffusion_coeffs(1)

      call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    
      if (ifcomp .gt. 0d0) then
         call this%compute(u_res, v_res, w_res, p_res, &
              ext_bdf, gs_Xh, c_Xh, rho, mu, ext_bdf%diffusion_coeffs(1), dt, &
              bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
              Ax, ksp_vel, ksp_prs, pc_prs, pc_vel, prs_max_iter, vel_max_iter)
      end if
      
      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            current_flow = &
                 device_glsc2(u%x_d, c_Xh%B_d, n) / this%domain_length  ! for X
         else if (this%flow_dir .eq. 2) then
            current_flow = &
                 device_glsc2(v%x_d, c_Xh%B_d, n) / this%domain_length  ! for Y
         else if (this%flow_dir .eq. 3) then
            current_flow = &
                 device_glsc2(w%x_d, c_Xh%B_d, n) / this%domain_length  ! for Z
         end if
      else
         if (this%flow_dir .eq. 1) then
            current_flow = glsc2(u%x, c_Xh%B, n) / this%domain_length  ! for X
         else if (this%flow_dir .eq. 2) then
            current_flow = glsc2(v%x, c_Xh%B, n) / this%domain_length  ! for Y
         else if (this%flow_dir .eq. 3) then
            current_flow = glsc2(w%x, c_Xh%B, n) / this%domain_length  ! for Z
         end if
      end if

      if (this%avflow) then
         xsec = c_Xh%volume / this%domain_length
         flow_rate = this%flow_rate*xsec
      endif
      
      delta_flow = flow_rate - current_flow            
      scale = delta_flow / this%base_flow
      
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(u%x_d, u_vol%x_d, scale, n)
         call device_add2s2(v%x_d, v_vol%x_d, scale, n)
         call device_add2s2(w%x_d, w_vol%x_d, scale, n)
         call device_add2s2(p%x_d, p_vol%x_d, scale, n)
      else
         call add2s2(u%x, u_vol%x, scale, n)
         call add2s2(v%x, v_vol%x, scale, n)
         call add2s2(w%x, w_vol%x, scale, n)
         call add2s2(p%x, p_vol%x, scale, n)
      end if
    end associate
    
  end subroutine fluid_vol_flow


end module fluid_volflow
