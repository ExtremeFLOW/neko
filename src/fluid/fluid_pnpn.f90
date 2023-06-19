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
  use logger
  use advection
  use profiler
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file, json_value
  implicit none
  private

  
  type, public, extends(fluid_scheme_t) :: fluid_pnpn_t
     type(field_t) :: p_res, u_res, v_res, w_res

     type(field_series_t) :: ulag, vlag, wlag

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
  end type fluid_pnpn_t

contains
  
  subroutine fluid_pnpn_init(this, msh, lx, params)    
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    character(len=20), dimension(:), allocatable :: bc_labels
    ! Variables for retrieving json parameters
    logical :: found, logical_val
    real(kind=rp) :: real_val
    integer :: integer_val
    character(len=:), allocatable :: string_val1, string_val2

    call this%free()
    
    ! Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, params, .true., .true., scheme)

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

      call field_init(this%p_res, dm_Xh, "p_res")
      call field_init(this%u_res, dm_Xh, "u_res")
      call field_init(this%v_res, dm_Xh, "v_res")
      call field_init(this%w_res, dm_Xh, "w_res")            
      call field_init(this%abx1, dm_Xh, "abx1")
      call field_init(this%aby1, dm_Xh, "aby1")
      call field_init(this%abz1, dm_Xh, "abz1")

      call field_init(this%abx2, dm_Xh, "abx2")
      call field_init(this%aby2, dm_Xh, "aby2")
      call field_init(this%abz2, dm_Xh, "abz2")
                  
      call field_init(this%du, dm_Xh, 'du')
      call field_init(this%dv, dm_Xh, 'dv')
      call field_init(this%dw, dm_Xh, 'dw')
      call field_init(this%dp, dm_Xh, 'dp')

      call this%ulag%init(this%u, 2)
      call this%vlag%init(this%v, 2)
      call this%wlag%init(this%w, 2)
      
    end associate
    
    call params%get('case.fluid.boundary_types', bc_labels, found) 
    if (.not. found) then
       if (allocated(bc_labels)) then
          deallocate(bc_labels)
       end if
       allocate(bc_labels(NEKO_MSH_MAX_ZLBLS))
       bc_labels = "not"
    end if

    ! Initialize velocity surface terms in pressure rhs
    call this%bc_prs_surface%init(this%dm_Xh)
    call this%bc_prs_surface%mark_zone(msh%inlet)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
                                                 'v', bc_labels)
    call this%bc_prs_surface%finalize()
    call this%bc_prs_surface%set_coef(this%c_Xh)
    ! Initialize symmetry surface terms in pressure rhs
    call this%bc_sym_surface%init(this%dm_Xh)
    call this%bc_sym_surface%mark_zone(msh%sympln)
    call this%bc_sym_surface%mark_zones_from_list(msh%labeled_zones,&
                                                 'sym', bc_labels)
    call this%bc_sym_surface%finalize()
    call this%bc_sym_surface%set_coef(this%c_Xh)
    ! Initialize dirichlet bcs for velocity residual
    call this%bc_vel_res_non_normal%init(this%dm_Xh)
    call this%bc_vel_res_non_normal%mark_zone(msh%outlet_normal)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
                                                         'on', bc_labels)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
                                                        'on+dong', bc_labels)
    call this%bc_vel_res_non_normal%finalize()
    call this%bc_vel_res_non_normal%init_msk(this%c_Xh)    

    call this%bc_dp%init(this%dm_Xh)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones, 'on+dong', &
                                         bc_labels)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones, 'o+dong', bc_labels)
    call this%bc_dp%finalize()
    call this%bc_dp%set_g(0.0_rp)
    call bc_list_init(this%bclst_dp)
    call bc_list_add(this%bclst_dp, this%bc_dp)
    !Add 0 prs bcs
    call bc_list_add(this%bclst_dp, this%bc_prs)

    call this%bc_vel_res%init(this%dm_Xh)
    call this%bc_vel_res%mark_zone(msh%inlet)
    call this%bc_vel_res%mark_zone(msh%wall)
    call this%bc_vel_res%mark_zones_from_list(msh%labeled_zones, 'v', bc_labels)
    call this%bc_vel_res%mark_zones_from_list(msh%labeled_zones, 'w', bc_labels)
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
    call json_get_or_default(params,&
                            'case.fluid.pressure_solver.projection_space_size',&
                             integer_val, 0)
    if (integer_val .gt. 0) then
       call this%proj_prs%init(this%dm_Xh%size(), integer_val)
    end if
    
    call json_get_or_default(params,&
                            'case.fluid.velocity_solver.projection_space_size',&
                             integer_val, 0)
    if (integer_val .gt. 0) then
       call this%proj_u%init(this%dm_Xh%size(), integer_val)
       call this%proj_v%init(this%dm_Xh%size(), integer_val)
       call this%proj_w%init(this%dm_Xh%size(), integer_val)
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
   
    call field_free(this%p_res)        
    call field_free(this%u_res)
    call field_free(this%v_res)
    call field_free(this%w_res)
    
    call field_free(this%du)
    call field_free(this%dv)
    call field_free(this%dw)
    call field_free(this%dp)
    
    call field_free(this%abx1)
    call field_free(this%aby1)
    call field_free(this%abz1)

    call field_free(this%abx2)
    call field_free(this%aby2)
    call field_free(this%abz2)
    
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
    
    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()
    
  end subroutine fluid_pnpn_free

  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  subroutine fluid_pnpn_step(this, t, tstep, dt, ext_bdf)
    class(fluid_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    integer :: n
    type(ksp_monitor_t) :: ksp_results(4)
    ! Variables for retrieving json parameters
    type(json_value), pointer :: json_val
    logical :: found, logical_val
    real(kind=rp) :: real_val, rho, mu, Re
    integer :: integer_val, ksp_vel_maxiter, ksp_pr_maxiter
    type(field_t), pointer :: u_e, v_e, w_e
    integer :: temp_indices(3)

    n = this%dm_Xh%size()

    call profiler_start_region('Fluid')
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         params => this%params, msh => this%msh, prs_res => this%prs_res, &
         vel_res => this%vel_res, sumab => this%sumab, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         rho => this%rho, Re => this%Re, mu => this%mu)


      call json_get(params, 'case.fluid.velocity_solver.max_iterations', &
                    ksp_vel_maxiter)
      call json_get(params, 'case.fluid.pressure_solver.max_iterations', &
                    ksp_pr_maxiter)
         
      ! Get temporary arrays
      call this%scratch%request_field(u_e, temp_indices(1))
      call this%scratch%request_field(v_e, temp_indices(2))
      call this%scratch%request_field(w_e, temp_indices(3))

      call sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
           ulag, vlag, wlag, ext_bdf%advection_coeffs, ext_bdf%nadv)
     
      call f_Xh%eval(t)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opcolv(f_Xh%u_d, f_Xh%v_d, f_Xh%w_d, c_Xh%B_d, msh%gdim, n)
      else
         call opcolv(f_Xh%u, f_Xh%v, f_Xh%w, c_Xh%B, msh%gdim, n)
      end if

      call this%adv%apply(this%u, this%v, this%w, &
                          f_Xh%u, f_Xh%v, f_Xh%w, &
                          Xh, this%c_Xh, dm_Xh%size())

      call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
                           this%abx2, this%aby2, this%abz2, &
                           f_Xh%u, f_Xh%v, f_Xh%w,&
                           rho, ext_bdf%advection_coeffs, n)

      call makebdf%compute_fluid(ulag, vlag, wlag, f_Xh%u, f_Xh%v, f_Xh%w, &
                           u, v, w, c_Xh%B, rho, dt, &
                           ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)

      call ulag%update()
      call vlag%update()
      call wlag%update()
      !> We assume that no change of boundary conditions 
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply dirichlet
      call this%bc_apply_vel()
      call this%bc_apply_prs()

      ! compute pressure
      call profiler_start_region('Pressure residual')
      call prs_res%compute(p, p_res, u, v, w, u_e, v_e, w_e, &
                           f_Xh, c_Xh, gs_Xh, this%bc_prs_surface, &
                           this%bc_sym_surface, Ax, ext_bdf%diffusion_coeffs(1), &
                           dt, Re, rho)
      
      call gs_op(gs_Xh, p_res, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_dp, p_res%x, p%dof%size())
      call profiler_end_region

      call params%get('case.fluid.pressure_solver.projection_space_size', &
                      integer_val, found)
     
      if( tstep .gt. 5 .and. found .and. integer_val .gt. 0) then
         call this%proj_prs%project_on(p_res%x, c_Xh, n)
         call this%proj_prs%log_info('Pressure')
      end if
      
      call this%pc_prs%update()
      call profiler_start_region('Pressure solve')
      ksp_results(1) = this%ksp_prs%solve(Ax, dp, p_res%x, n, c_Xh, &
                                          this%bclst_dp, gs_Xh, ksp_pr_maxiter)
      call profiler_end_region

      if( tstep .gt. 5 .and. found .and. integer_val .gt. 0) then
         call this%proj_prs%project_back(dp%x, Ax, c_Xh, &
                                         this%bclst_dp, gs_Xh, n)
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2(p%x_d, dp%x_d,n)
      else
         call add2(p%x, dp%x,n)
      end if
      

      ! compute velocity
      call profiler_start_region('Velocity residual')
      call vel_res%compute(Ax, u, v, w, &
                           u_res, v_res, w_res, &
                           p, &
                           f_Xh, c_Xh, msh, Xh, &
                           Re, rho, ext_bdf%diffusion_coeffs(1), &
                           dt, dm_Xh%size())
      
      call gs_op(gs_Xh, u_res, GS_OP_ADD) 
      call gs_op(gs_Xh, v_res, GS_OP_ADD) 
      call gs_op(gs_Xh, w_res, GS_OP_ADD) 

      call bc_list_apply_vector(this%bclst_vel_res,&
                                u_res%x, v_res%x, w_res%x, dm_Xh%size())
      call profiler_end_region
      
      call params%get('case.fluid.velocity_solver.projection_space_size', &
                      integer_val, found)
      if (tstep .gt. 5 .and. found .and. integer_val .gt. 0) then 
         call this%proj_u%project_on(u_res%x, c_Xh, n)
         call this%proj_v%project_on(v_res%x, c_Xh, n)
         call this%proj_w%project_on(w_res%x, c_Xh, n)
      end if

      call this%pc_vel%update()

      call profiler_start_region("Velocity solve")
      ksp_results(2) = this%ksp_vel%solve(Ax, du, u_res%x, n, &
           c_Xh, this%bclst_du, gs_Xh, ksp_vel_maxiter)
      ksp_results(3) = this%ksp_vel%solve(Ax, dv, v_res%x, n, &
           c_Xh, this%bclst_dv, gs_Xh, ksp_vel_maxiter)
      ksp_results(4) = this%ksp_vel%solve(Ax, dw, w_res%x, n, &
           c_Xh, this%bclst_dw, gs_Xh, ksp_vel_maxiter)
      call profiler_end_region

      if (tstep .gt. 5 .and. found .and. integer_val .gt. 0) then 
         call this%proj_u%project_back(du%x, Ax, c_Xh, &
                                  this%bclst_du, gs_Xh, n)
         call this%proj_v%project_back(dv%x, Ax, c_Xh, &
                                  this%bclst_dv, gs_Xh, n)
         call this%proj_w%project_back(dw%x, Ax, c_Xh, &
                                  this%bclst_dw, gs_Xh, n)
      end if
      
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if

      call params%get('case.fluid.flow_rate_force.rate', real_val, found)
      if (found) then
         call this%vol_flow%adjust( u, v, w, p, u_res, v_res, w_res, p_res, &
              c_Xh, gs_Xh, ext_bdf, rho, Re,&
              dt, this%bclst_dp, this%bclst_du, this%bclst_dv, &
              this%bclst_dw, this%bclst_vel_res, Ax, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, ksp_pr_maxiter, &
              ksp_vel_maxiter)
      end if
      
      call fluid_step_info(tstep, t, dt, ksp_results)
      
      call this%scratch%relinquish_field(temp_indices)
      
    end associate
    call profiler_end_region
  end subroutine fluid_pnpn_step

  
end module fluid_pnpn
