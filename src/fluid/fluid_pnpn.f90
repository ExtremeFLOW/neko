! Copyright (c) 2022, The Neko Authors
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
!> Modular version of the Classic Nek5000 Pn/Pn formulation for fluids
module fluid_pnpn
  use pnpn_res_fctry
  use ax_helm_fctry
  use rhs_maker_fctry
  use fluid_volflow
  use fluid_scheme
  use field_series
  use facet_normal
  use device_math
  use device_mathops
  use fluid_aux
  use time_scheme_controller
  use projection
  use device
  use logger
  use advection
  use profiler
  use json_utils, only : json_get
  use json_module, only : json_file
  use material_properties, only : material_properties_t
  implicit none
  private


  type, public, extends(fluid_scheme_t) :: fluid_pnpn_t
     type(field_t) :: p_res, u_res, v_res, w_res

     type(field_t) :: dp, du, dv, dw

     class(ax_t), allocatable :: Ax

     type(projection_t) :: proj_prs
     type(projection_t) :: proj_u
     type(projection_t) :: proj_v
     type(projection_t) :: proj_w

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(facet_normal_t) :: bc_sym_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_res   !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_dp   !< Dirichlet condition vel. res.
     type(non_normal_t) :: bc_vel_res_non_normal   !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_res
     type(bc_list_t) :: bclst_du
     type(bc_list_t) :: bclst_dv
     type(bc_list_t) :: bclst_dw
     type(bc_list_t) :: bclst_dp

     class(advection_t), allocatable :: adv

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2

     !> Pressure residual
     class(pnpn_prs_res_t), allocatable :: prs_res

     !> Velocity residual
     class(pnpn_vel_res_t), allocatable :: vel_res

     !> Summation of AB/BDF contributions
     class(rhs_maker_sumab_t), allocatable :: sumab

     !> Contributions to kth order extrapolation scheme
     class(rhs_maker_ext_t), allocatable :: makeabf

     !> Contributions to F from lagged BD terms
     class(rhs_maker_bdf_t), allocatable :: makebdf

     !> Adjust flow volume
     type(fluid_volflow_t) :: vol_flow

   contains
     procedure, pass(this) :: init => fluid_pnpn_init
     procedure, pass(this) :: free => fluid_pnpn_free
     procedure, pass(this) :: step => fluid_pnpn_step
     procedure, pass(this) :: restart => fluid_pnpn_restart
  end type fluid_pnpn_t

contains

  subroutine fluid_pnpn_init(this, msh, lx, params, user, material_properties)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), intent(in) :: user
    type(material_properties_t), intent(inout) :: material_properties
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    logical :: found, logical_val
    integer :: integer_val
    real(kind=rp) :: real_val

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, .true., .true., scheme, user, &
                          material_properties)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax)

    ! Setup backend dependent prs residual routines
    call pnpn_prs_res_factory(this%prs_res)

    ! Setup backend dependent vel residual routines
    call pnpn_vel_res_factory(this%vel_res)

    ! Setup backend dependent summation of AB/BDF
    call rhs_maker_sumab_fctry(this%sumab)

    ! Setup backend dependent summation of extrapolation scheme
    call rhs_maker_ext_fctry(this%makeabf)

    ! Setup backend depenent contributions to F from lagged BD terms
    call rhs_maker_bdf_fctry(this%makebdf)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%p_res%init(dm_Xh, "p_res")
      call this%u_res%init(dm_Xh, "u_res")
      call this%v_res%init(dm_Xh, "v_res")
      call this%w_res%init(dm_Xh, "w_res")
      call this%abx1%init(dm_Xh, "abx1")
      call this%aby1%init(dm_Xh, "aby1")
      call this%abz1%init(dm_Xh, "abz1")
      call this%abx2%init(dm_Xh, "abx2")
      call this%aby2%init(dm_Xh, "aby2")
      call this%abz2%init(dm_Xh, "abz2")
      this%abx1 = 0.0_rp
      this%aby1 = 0.0_rp
      this%abz1 = 0.0_rp
      this%abx2 = 0.0_rp
      this%aby2 = 0.0_rp
      this%abz2 = 0.0_rp

      call this%du%init(dm_Xh, 'du')
      call this%dv%init(dm_Xh, 'dv')
      call this%dw%init(dm_Xh, 'dw')
      call this%dp%init(dm_Xh, 'dp')

    end associate

    ! Initialize velocity surface terms in pressure rhs
    call this%bc_prs_surface%init(this%dm_Xh)
    call this%bc_prs_surface%mark_zone(msh%inlet)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
                                                 'v', this%bc_labels)
    call this%bc_prs_surface%finalize()
    call this%bc_prs_surface%set_coef(this%c_Xh)
    ! Initialize symmetry surface terms in pressure rhs
    call this%bc_sym_surface%init(this%dm_Xh)
    call this%bc_sym_surface%mark_zone(msh%sympln)
    call this%bc_sym_surface%mark_zones_from_list(msh%labeled_zones,&
                                                 'sym', this%bc_labels)
    call this%bc_sym_surface%finalize()
    call this%bc_sym_surface%set_coef(this%c_Xh)
    ! Initialize dirichlet bcs for velocity residual
    call this%bc_vel_res_non_normal%init(this%dm_Xh)
    call this%bc_vel_res_non_normal%mark_zone(msh%outlet_normal)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
                                                         'on', this%bc_labels)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
                                                         'on+dong', &
                                                         this%bc_labels)
    call this%bc_vel_res_non_normal%finalize()
    call this%bc_vel_res_non_normal%init_msk(this%c_Xh)

    call this%bc_dp%init(this%dm_Xh)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones, 'on+dong', &
                                         this%bc_labels)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones, &
                                         'o+dong', this%bc_labels)
    call this%bc_dp%finalize()
    call this%bc_dp%set_g(0.0_rp)
    call bc_list_init(this%bclst_dp)
    call bc_list_add(this%bclst_dp, this%bc_dp)
    !Add 0 prs bcs
    call bc_list_add(this%bclst_dp, this%bc_prs)

    call this%bc_vel_res%init(this%dm_Xh)
    call this%bc_vel_res%mark_zone(msh%inlet)
    call this%bc_vel_res%mark_zone(msh%wall)
    call this%bc_vel_res%mark_zones_from_list(msh%labeled_zones, &
                                              'v', this%bc_labels)
    call this%bc_vel_res%mark_zones_from_list(msh%labeled_zones, &
                                              'w', this%bc_labels)
    call this%bc_vel_res%finalize()
    call this%bc_vel_res%set_g(0.0_rp)
    call bc_list_init(this%bclst_vel_res)
    call bc_list_add(this%bclst_vel_res, this%bc_vel_res)
    call bc_list_add(this%bclst_vel_res, this%bc_vel_res_non_normal)
    call bc_list_add(this%bclst_vel_res, this%bc_sym)

    !Initialize bcs for u, v, w velocity components
    call bc_list_init(this%bclst_du)
    call bc_list_add(this%bclst_du,this%bc_sym%bc_x)
    call bc_list_add(this%bclst_du,this%bc_vel_res_non_normal%bc_x)
    call bc_list_add(this%bclst_du, this%bc_vel_res)

    call bc_list_init(this%bclst_dv)
    call bc_list_add(this%bclst_dv,this%bc_sym%bc_y)
    call bc_list_add(this%bclst_dv,this%bc_vel_res_non_normal%bc_y)
    call bc_list_add(this%bclst_dv, this%bc_vel_res)

    call bc_list_init(this%bclst_dw)
    call bc_list_add(this%bclst_dw,this%bc_sym%bc_z)
    call bc_list_add(this%bclst_dw,this%bc_vel_res_non_normal%bc_z)
    call bc_list_add(this%bclst_dw, this%bc_vel_res)

    !Intialize projection space thingy
    if (this%pr_projection_dim .gt. 0) then
       call this%proj_prs%init(this%dm_Xh%size(), this%pr_projection_dim, &
                              this%pr_projection_activ_step)
    end if

    if (this%vel_projection_dim .gt. 0) then
       call this%proj_u%init(this%dm_Xh%size(), this%vel_projection_dim, &
                              this%vel_projection_activ_step)
       call this%proj_v%init(this%dm_Xh%size(), this%vel_projection_dim, &
                              this%vel_projection_activ_step)
       call this%proj_w%init(this%dm_Xh%size(), this%vel_projection_dim, &
                              this%vel_projection_activ_step)
    end if

    ! Add lagged term to checkpoint
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)

    call json_get(params, 'case.numerics.dealias', logical_val)
    call params%get('case.numerics.dealiased_polynomial_order', integer_val, &
                    found)
    if (.not. found) then
       call json_get(params, 'case.numerics.polynomial_order', integer_val)
       integer_val =  3.0_rp / 2.0_rp * (integer_val + 1) - 1
    end if
    ! an extra +1 below to go from poly order to space size
    call advection_factory(this%adv, this%c_Xh, logical_val, integer_val + 1)

    if (params%valid_path('case.fluid.flow_rate_force')) then
       call this%vol_flow%init(this%dm_Xh, params)
    end if

  end subroutine fluid_pnpn_init

  subroutine fluid_pnpn_restart(this,dtlag, tlag)
    class(fluid_pnpn_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
    type(field_t) :: u_temp, v_temp, w_temp
    integer :: i, n

    n = this%u%dof%size()
    ! Make sure that continuity is maintained (important for interpolation) 
    ! Do not do this for lagged rhs (derivatives are not necessairly coninous across elements)
    call col2(this%u%x,this%c_Xh%mult,this%u%dof%size())
    call col2(this%v%x,this%c_Xh%mult,this%u%dof%size())
    call col2(this%w%x,this%c_Xh%mult,this%u%dof%size())
    call col2(this%p%x,this%c_Xh%mult,this%u%dof%size())
    do i = 1, this%ulag%size()
       call col2(this%ulag%lf(i)%x,this%c_Xh%mult,this%u%dof%size())
       call col2(this%vlag%lf(i)%x,this%c_Xh%mult,this%u%dof%size())
       call col2(this%wlag%lf(i)%x,this%c_Xh%mult,this%u%dof%size())
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u=>this%u, v=>this%v, w=>this%w, &
            ulag=>this%ulag, vlag=>this%vlag, wlag=>this%wlag,&
            p=>this%p)
         call device_memcpy(u%x, u%x_d, u%dof%size(), &
                            HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(v%x, v%x_d, v%dof%size(), &
                            HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(w%x, w%x_d, w%dof%size(), &
                            HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(p%x, p%x_d, p%dof%size(), &
                            HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
                            u%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
                            u%dof%size(), HOST_TO_DEVICE, sync=.false.)

         call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
                            v%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
                            v%dof%size(), HOST_TO_DEVICE, sync=.false.)

         call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%abx1%x, this%abx1%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%abx2%x, this%abx2%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%aby1%x, this%aby1%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%aby2%x, this%aby2%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%abz1%x, this%abz1%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
         call device_memcpy(this%abz2%x, this%abz2%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync=.false.)
       end associate
    end if


    call this%gs_Xh%op(this%u,GS_OP_ADD)
    call this%gs_Xh%op(this%v,GS_OP_ADD)
    call this%gs_Xh%op(this%w,GS_OP_ADD)
    call this%gs_Xh%op(this%p,GS_OP_ADD)

    do i = 1, this%ulag%size()
       call this%gs_Xh%op(this%ulag%lf(i),GS_OP_ADD)
       call this%gs_Xh%op(this%vlag%lf(i),GS_OP_ADD)
       call this%gs_Xh%op(this%wlag%lf(i),GS_OP_ADD)
    end do

    !! If we would decide to only restart from lagged fields instead of asving abx1, aby1 etc.
    !! Observe that one also needs to recompute the focing at the old time steps
    !u_temp = this%ulag%lf(2)
    !v_temp = this%vlag%lf(2)
    !w_temp = this%wlag%lf(2)
    !! Compute the source terms
    !call this%source_term%compute(tlag(2), -1)
    !
    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !! Add the advection operators to the right-hand-side.
    !call this%adv%compute(u_temp, v_temp, w_temp, &
    !                      this%f_x%x, this%f_y%x, this%f_z%x, &
    !                      this%Xh, this%c_Xh, this%dm_Xh%size())
    !this%abx2 = this%f_x
    !this%aby2 = this%f_y
    !this%abz2 = this%f_z
    !
    !u_temp = this%ulag%lf(1)
    !v_temp = this%vlag%lf(1)
    !w_temp = this%wlag%lf(1)
    !call this%source_term%compute(tlag(1), 0)

    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !call this%adv%compute(u_temp, v_temp, w_temp, &
    !                      this%f_x%x, this%f_y%x, this%f_z%x, &
    !                      this%Xh, this%c_Xh, this%dm_Xh%size())
    !this%abx1 = this%f_x
    !this%aby1 = this%f_y
    !this%abz1 = this%f_z

  end subroutine fluid_pnpn_restart

  subroutine fluid_pnpn_free(this)
    class(fluid_pnpn_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call this%bc_prs_surface%free()
    call this%bc_sym_surface%free()
    call bc_list_free(this%bclst_vel_res)
    call bc_list_free(this%bclst_dp)
    call this%proj_prs%free()
    call this%proj_u%free()
    call this%proj_v%free()
    call this%proj_w%free()

    call this%p_res%free()
    call this%u_res%free()
    call this%v_res%free()
    call this%w_res%free()

    call this%du%free()
    call this%dv%free()
    call this%dw%free()
    call this%dp%free()

    call this%abx1%free()
    call this%aby1%free()
    call this%abz1%free()

    call this%abx2%free()
    call this%aby2%free()
    call this%abz2%free()

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    if (allocated(this%prs_res)) then
       deallocate(this%prs_res)
    end if

    if (allocated(this%vel_res)) then
       deallocate(this%vel_res)
    end if

    if (allocated(this%sumab)) then
       deallocate(this%sumab)
    end if

    if (allocated(this%makeabf)) then
       deallocate(this%makeabf)
    end if

    if (allocated(this%makebdf)) then
       deallocate(this%makebdf)
    end if

    call this%vol_flow%free()

  end subroutine fluid_pnpn_free

  !> Advance fluid simulation in time.
  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  subroutine fluid_pnpn_step(this, t, tstep, dt, ext_bdf, if_variable_dt, dt_last_change)
    class(fluid_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    ! number of degrees of freedom
    integer :: n
    ! Solver results monitors (pressure + 3 velocity)
    type(ksp_monitor_t) :: ksp_results(4)
    ! Extrapolated velocity for the pressure residual
    type(field_t), pointer :: u_e, v_e, w_e
    ! Indices for tracking temporary fields
    integer :: temp_indices(3)
    ! Counter
    integer :: i
    ! time step controller
    logical, intent(in) :: if_variable_dt
    integer, intent(in) :: dt_last_change ! 0 at the step where dt changes

    if (this%freeze) return

    n = this%dm_Xh%size()

    call profiler_start_region('Fluid', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, prs_res => this%prs_res, &
         source_term => this%source_term, &
         vel_res => this%vel_res, sumab => this%sumab, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         vel_projection_dim => this%vel_projection_dim, &
         pr_projection_dim => this%pr_projection_dim, &
         rho => this%rho, mu => this%mu, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z)

      ! Get temporary arrays
      call this%scratch%request_field(u_e, temp_indices(1))
      call this%scratch%request_field(v_e, temp_indices(2))
      call this%scratch%request_field(w_e, temp_indices(3))
      call sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
           ulag, vlag, wlag, ext_bdf%advection_coeffs, ext_bdf%nadv)

      ! Compute the source terms
      call this%source_term%compute(t, tstep)

      ! Pre-multiply the source terms with the mass matrix.
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opcolv(f_x%x_d, f_y%x_d, f_z%x_d, c_Xh%B_d, msh%gdim, n)
      else
         call opcolv(f_x%x, f_y%x, f_z%x, c_Xh%B, msh%gdim, n)
      end if

      ! Add the advection operators to the right-hand-side.
      call this%adv%compute(u, v, w, &
                            f_x%x, f_y%x, f_z%x, &
                            Xh, this%c_Xh, dm_Xh%size())

      ! At this point the RHS contains the sum of the advection operator and
      ! additional source terms, evaluated using the velocity field from the
      ! previous time-step. Now, this value is used in the explicit time
      ! scheme to advance both terms in time.
      call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
                           this%abx2, this%aby2, this%abz2, &
                           f_x%x, f_y%x, f_z%x, &
                           rho, ext_bdf%advection_coeffs, n)

      ! Add the RHS contributions coming from the BDF scheme.
      call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
                           u, v, w, c_Xh%B, rho, dt, &
                           ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)

      call ulag%update()
      call vlag%update()
      call wlag%update()

      !> We assume that no change of boundary conditions
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply dirichlet
      call this%bc_apply_vel(t, tstep)
      call this%bc_apply_prs(t, tstep)

      ! Compute pressure.
      call profiler_start_region('Pressure residual', 18)
      call prs_res%compute(p, p_res, u, v, w, u_e, v_e, w_e, &
                           f_x, f_y, f_z, c_Xh, gs_Xh, this%bc_prs_surface, &
                           this%bc_sym_surface, Ax, ext_bdf%diffusion_coeffs(1), &
                           dt, mu, rho)

      call gs_Xh%op(p_res, GS_OP_ADD)
      call bc_list_apply_scalar(this%bclst_dp, p_res%x, p%dof%size(), t, tstep)
      call profiler_end_region

      if( tstep .gt. this%proj_prs%activ_step .and. pr_projection_dim .gt. 0) then
         if (if_variable_dt) then
            if (dt_last_change .eq. 0) then ! the time step at which dt is changed
               call this%proj_prs%clear(n) 
            else if (dt_last_change .gt. this%proj_prs%activ_step - 1) then
               ! activate projection some steps after dt is changed
               ! note that dt_last_change start from 0
               call this%proj_prs%project_on(p_res%x, c_Xh, n)
               call this%proj_prs%log_info('Pressure')
            end if
         else
            call this%proj_prs%project_on(p_res%x, c_Xh, n)
            call this%proj_prs%log_info('Pressure')
         end if
      end if

      call this%pc_prs%update()
      call profiler_start_region('Pressure solve', 3)
      ksp_results(1) = &
         this%ksp_prs%solve(Ax, dp, p_res%x, n, c_Xh,  this%bclst_dp, gs_Xh)

      call profiler_end_region

      if( tstep .gt. this%proj_prs%activ_step .and. pr_projection_dim .gt. 0) then
         if (.not.(if_variable_dt) .or. &
            (dt_last_change .gt. this%proj_prs%activ_step - 1)) then

            call this%proj_prs%project_back(dp%x, Ax, c_Xh, &
                                         this%bclst_dp, gs_Xh, n)
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2(p%x_d, dp%x_d,n)
      else
         call add2(p%x, dp%x,n)
      end if


      ! Compute velocity.
      call profiler_start_region('Velocity residual', 19)
      call vel_res%compute(Ax, u, v, w, &
                           u_res, v_res, w_res, &
                           p, &
                           f_x, f_y, f_z, &
                           c_Xh, msh, Xh, &
                           mu, rho, ext_bdf%diffusion_coeffs(1), &
                           dt, dm_Xh%size())

      call gs_Xh%op(u_res, GS_OP_ADD)
      call gs_Xh%op(v_res, GS_OP_ADD)
      call gs_Xh%op(w_res, GS_OP_ADD)

      call bc_list_apply_vector(this%bclst_vel_res,&
                                u_res%x, v_res%x, w_res%x, dm_Xh%size(),&
                                t, tstep)

      call profiler_end_region

      if (tstep .gt. this%proj_u%activ_step .and. vel_projection_dim .gt. 0) then         
         if (if_variable_dt) then
            if (dt_last_change .eq. 0) then ! the time step at which dt is changed
               call this%proj_u%clear(n)
               call this%proj_v%clear(n)
               call this%proj_w%clear(n) 
            else if (dt_last_change .gt. this%proj_u%activ_step - 1) then
               ! activate projection some steps after dt is changed
               ! note that dt_last_change start from 0
               call this%proj_u%project_on(u_res%x, c_Xh, n)
               call this%proj_v%project_on(v_res%x, c_Xh, n)
               call this%proj_w%project_on(w_res%x, c_Xh, n)
            end if
         else
            call this%proj_u%project_on(u_res%x, c_Xh, n)
            call this%proj_v%project_on(v_res%x, c_Xh, n)
            call this%proj_w%project_on(w_res%x, c_Xh, n)
         end if
      end if

      call this%pc_vel%update()

      call profiler_start_region("Velocity solve", 4)
      ksp_results(2) = this%ksp_vel%solve(Ax, du, u_res%x, n, &
           c_Xh, this%bclst_du, gs_Xh)
      ksp_results(3) = this%ksp_vel%solve(Ax, dv, v_res%x, n, &
           c_Xh, this%bclst_dv, gs_Xh)
      ksp_results(4) = this%ksp_vel%solve(Ax, dw, w_res%x, n, &
           c_Xh, this%bclst_dw, gs_Xh)
      call profiler_end_region

      if (tstep .gt. this%proj_u%activ_step .and. vel_projection_dim .gt. 0) then
         if (.not.(if_variable_dt) .or. &
            (dt_last_change .gt. this%proj_u%activ_step - 1)) then
            call this%proj_u%project_back(du%x, Ax, c_Xh, &
                                    this%bclst_du, gs_Xh, n)
            call this%proj_v%project_back(dv%x, Ax, c_Xh, &
                                    this%bclst_dv, gs_Xh, n)
            call this%proj_w%project_back(dw%x, Ax, c_Xh, &
                                    this%bclst_dw, gs_Xh, n)
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if

      if (this%forced_flow_rate) then
         call this%vol_flow%adjust( u, v, w, p, u_res, v_res, w_res, p_res, &
              c_Xh, gs_Xh, ext_bdf, rho, mu,&
              dt, this%bclst_dp, this%bclst_du, this%bclst_dv, &
              this%bclst_dw, this%bclst_vel_res, Ax, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(tstep, t, dt, ksp_results)

      call this%scratch%relinquish_field(temp_indices)

    end associate
    call profiler_end_region
  end subroutine fluid_pnpn_step


end module fluid_pnpn
