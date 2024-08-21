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
!> Linearized Navier Stokes solver based on Pn/Pn formulation.
module adjoint_pnpn
  use num_types, only : rp, dp
  use krylov, only : ksp_monitor_t
  use pnpn_res_fctry, only : pnpn_prs_res_factory, pnpn_vel_res_factory
  use pnpn_res_stress_fctry, only : pnpn_prs_res_stress_factory, &
       pnpn_vel_res_stress_factory
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t
  use ax_helm_fctry, only : ax_helm_factory
  use rhs_maker_fctry, only : rhs_maker_sumab_fctry, rhs_maker_bdf_fctry, &
       rhs_maker_ext_fctry
  use rhs_maker, only : rhs_maker_sumab_t, rhs_maker_bdf_t, rhs_maker_ext_t
  use fluid_volflow, only : fluid_volflow_t
  use adjoint_scheme, only : adjoint_scheme_t
  use field_series, only : field_series_t
  use device_math, only : device_add2, device_col2, device_cmult
  use device_mathops, only : device_opcolv, device_opadd2cm
  use fluid_aux, only : fluid_step_info
  use time_scheme_controller, only : time_scheme_controller_t
  use projection, only : projection_t
  use device, only : device_memcpy, HOST_TO_DEVICE
  use logger, only : neko_log, NEKO_LOG_DEBUG
  use advection, only : advection_lin_t
  use profiler, only : profiler_start_region, profiler_end_region
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use material_properties, only : material_properties_t
  use advection_fctry, only : advection_lin_factory
  use ax_product, only : ax_t
  use field, only : field_t
  use dirichlet, only : dirichlet_t
  use facet_normal, only : facet_normal_t
  use non_normal, only : non_normal_t
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use coefs, only : coef_t
  use time_step_controller, only : time_step_controller_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : col2, add2, copy, glsc2, vdot3, cmult
  use mathops, only : opadd2cm, opcolv
  use bc, only: bc_list_t, bc_list_init, bc_list_add, bc_list_free, &
       bc_list_apply_scalar, bc_list_apply_vector
  use utils, only : neko_error
  use field_math, only : field_add2
  use flow_ic, only : set_flow_ic
  use file, only : file_t, file_free, fld_file_data_t
  use field_registry, only : neko_field_registry
  use operators, only : curl
  use fld_file_output, only : fld_file_output_t
  use comm, only : pe_rank
  use vector, only : vector_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private


  type, public, extends(adjoint_scheme_t) :: adjoint_pnpn_t
     type(field_t) :: p_res, u_res, v_res, w_res

     type(field_t) :: dp, du, dv, dw

     ! Coupled Helmholz operator for velocity
     class(ax_t), allocatable :: Ax_vel
     ! Helmholz operator for pressure
     class(ax_t), allocatable :: Ax_prs

     type(field_t), pointer :: u_b => null() !< x-component of baseflow velocity
     type(field_t), pointer :: v_b => null() !< y-component of baseflow Velocity
     type(field_t), pointer :: w_b => null() !< z-component of baseflow Velocity
     type(field_t), pointer :: p_b => null() !< Baseflow pressure

     type(projection_t) :: proj_prs
     type(projection_t) :: proj_u
     type(projection_t) :: proj_v
     type(projection_t) :: proj_w

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(facet_normal_t) :: bc_sym_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_res !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_field_dirichlet_p !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_field_dirichlet_u !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_field_dirichlet_v !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_field_dirichlet_w !< Dirichlet condition vel. res.
     type(non_normal_t) :: bc_vel_res_non_normal !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_res
     type(bc_list_t) :: bclst_du
     type(bc_list_t) :: bclst_dv
     type(bc_list_t) :: bclst_dw
     type(bc_list_t) :: bclst_dp

     class(advection_lin_t), allocatable :: adv

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

     ! ======================================================================= !
     ! Addressable attributes

     real(kind=rp) :: norm_scaling !< Constant for the norm of the velocity field.
     real(kind=rp) :: norm_target !< Target norm for the velocity field.
     real(kind=rp) :: norm_tolerance !< Tolerance for when to rescale the flow.

     ! ======================================================================= !
     ! Definition of shorthands and local variables

     !> Norm of the base field.
     real(kind=rp) :: norm_l2_base

     !> Upper limit for the norm
     real(kind=rp) :: norm_l2_upper
     !> Lower limit for the norm
     real(kind=rp) :: norm_l2_lower

     !> Output file
     type(file_t) :: file_output

   contains
     procedure, pass(this) :: init => adjoint_pnpn_init
     procedure, pass(this) :: free => adjoint_pnpn_free
     procedure, pass(this) :: step => adjoint_pnpn_step
     procedure, pass(this) :: restart => adjoint_pnpn_restart

     !> Compute the power_iterations field.
     procedure, public, pass(this) :: PW_compute_ => power_iterations_compute
  end type adjoint_pnpn_t

contains

  subroutine adjoint_pnpn_init(this, msh, lx, params, user, material_properties)
    class(adjoint_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), intent(in) :: user
    type(material_properties_t), target, intent(inout) :: material_properties
    character(len=20), parameter :: scheme = 'Perturbation (Pn/Pn)'
    type(file_t) :: field_file, out_file
    type(fld_file_data_t) :: field_data
    integer :: n

    ! Temporary field pointers
    real(kind=rp) :: norm_l2_base
    character(len=:), allocatable :: file_name
    character(len=256) :: header_line

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, .true., .true., scheme, user, &
         material_properties)

    if (this%variable_material_properties .eqv. .true.) then
       ! Setup backend dependent Ax routines
       call ax_helm_factory(this%Ax_vel, full_formulation = .true.)

       ! Setup backend dependent prs residual routines
       call pnpn_prs_res_stress_factory(this%prs_res)

       ! Setup backend dependent vel residual routines
       call pnpn_vel_res_stress_factory(this%vel_res)
    else
       ! Setup backend dependent Ax routines
       call ax_helm_factory(this%Ax_vel, full_formulation = .false.)

       ! Setup backend dependent prs residual routines
       call pnpn_prs_res_factory(this%prs_res)

       ! Setup backend dependent vel residual routines
       call pnpn_vel_res_factory(this%vel_res)
    end if

    ! Setup Ax for the pressure
    call ax_helm_factory(this%Ax_prs, full_formulation = .false.)


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
    call this%bc_prs_surface%init_base(this%c_Xh)
    call this%bc_prs_surface%mark_zone(msh%inlet)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
         'v', this%bc_labels)
    !This impacts the rhs of the pressure, need to check what is correct to add here
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_u', this%bc_labels)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_v', this%bc_labels)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_w', this%bc_labels)
    call this%bc_prs_surface%finalize()
    ! Initialize symmetry surface terms in pressure rhs
    call this%bc_sym_surface%init_base(this%c_Xh)
    call this%bc_sym_surface%mark_zone(msh%sympln)
    call this%bc_sym_surface%mark_zones_from_list(msh%labeled_zones,&
         'sym', this%bc_labels)
    ! Same here, should du, dv, dw be marked here?
    call this%bc_sym_surface%finalize()
    ! Initialize dirichlet bcs for velocity residual
    call this%bc_vel_res_non_normal%init_base(this%c_Xh)
    call this%bc_vel_res_non_normal%mark_zone(msh%outlet_normal)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
         'on', this%bc_labels)
    call this%bc_vel_res_non_normal%mark_zones_from_list(msh%labeled_zones,&
         'on+dong', &
         this%bc_labels)
    call this%bc_vel_res_non_normal%finalize()
    call this%bc_vel_res_non_normal%init(this%c_Xh)

    call this%bc_field_dirichlet_p%init_base(this%c_Xh)
    call this%bc_field_dirichlet_p%mark_zones_from_list(msh%labeled_zones, &
         'on+dong', this%bc_labels)
    call this%bc_field_dirichlet_p%mark_zones_from_list(msh%labeled_zones, &
         'o+dong', this%bc_labels)
    call this%bc_field_dirichlet_p%mark_zones_from_list(msh%labeled_zones, &
         'd_pres', this%bc_labels)
    call this%bc_field_dirichlet_p%finalize()
    call this%bc_field_dirichlet_p%set_g(0.0_rp)
    call bc_list_init(this%bclst_dp)
    call bc_list_add(this%bclst_dp, this%bc_field_dirichlet_p)
    !Add 0 prs bcs
    call bc_list_add(this%bclst_dp, this%bc_prs)

    call this%bc_field_dirichlet_u%init_base(this%c_Xh)
    call this%bc_field_dirichlet_u%mark_zones_from_list( &
         msh%labeled_zones, 'd_vel_u', this%bc_labels)
    call this%bc_field_dirichlet_u%finalize()
    call this%bc_field_dirichlet_u%set_g(0.0_rp)

    call this%bc_field_dirichlet_v%init_base(this%c_Xh)
    call this%bc_field_dirichlet_v%mark_zones_from_list(msh%labeled_zones, &
         'd_vel_v', &
         this%bc_labels)
    call this%bc_field_dirichlet_v%finalize()
    call this%bc_field_dirichlet_v%set_g(0.0_rp)

    call this%bc_field_dirichlet_w%init_base(this%c_Xh)
    call this%bc_field_dirichlet_w%mark_zones_from_list(msh%labeled_zones, &
         'd_vel_w', &
         this%bc_labels)
    call this%bc_field_dirichlet_w%finalize()
    call this%bc_field_dirichlet_w%set_g(0.0_rp)

    call this%bc_vel_res%init_base(this%c_Xh)
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
    call bc_list_add(this%bclst_du, this%bc_sym%bc_x)
    call bc_list_add(this%bclst_du, this%bc_vel_res_non_normal%bc_x)
    call bc_list_add(this%bclst_du, this%bc_vel_res)
    call bc_list_add(this%bclst_du, this%bc_field_dirichlet_u)

    call bc_list_init(this%bclst_dv)
    call bc_list_add(this%bclst_dv, this%bc_sym%bc_y)
    call bc_list_add(this%bclst_dv, this%bc_vel_res_non_normal%bc_y)
    call bc_list_add(this%bclst_dv, this%bc_vel_res)
    call bc_list_add(this%bclst_dv, this%bc_field_dirichlet_v)

    call bc_list_init(this%bclst_dw)
    call bc_list_add(this%bclst_dw, this%bc_sym%bc_z)
    call bc_list_add(this%bclst_dw, this%bc_vel_res_non_normal%bc_z)
    call bc_list_add(this%bclst_dw, this%bc_vel_res)
    call bc_list_add(this%bclst_dw, this%bc_field_dirichlet_w)

    !Intialize projection space thingy

    if (this%variable_material_properties .and. &
         this%vel_projection_dim .gt. 0) then
       call neko_error("Velocity projection not available for full stress &
            &formulation")
    end if


    call this%proj_prs%init(this%dm_Xh%size(), this%pr_projection_dim, &
         this%pr_projection_activ_step)

    call this%proj_u%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)
    call this%proj_v%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)
    call this%proj_w%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)


    ! Add lagged term to checkpoint
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)

    call advection_lin_factory(this%adv, params, this%c_Xh)

    if (params%valid_path('case.fluid.flow_rate_force')) then
       call this%vol_flow%init(this%dm_Xh, params)
    end if

    ! ------------------------------------------------------------------------ !
    ! Handling the rescaling and baseflow

    ! Read the norm scaling from the json file
    call json_get_or_default(params, 'norm_scaling', &
         this%norm_scaling, 0.5_rp)

    ! Read the baseflow from file or user input.
    call neko_field_registry%add_field(this%dm_Xh, 'u_b')
    call neko_field_registry%add_field(this%dm_Xh, 'v_b')
    call neko_field_registry%add_field(this%dm_Xh, 'w_b')
    call neko_field_registry%add_field(this%dm_Xh, 'p_b')
    this%u_b => neko_field_registry%get_field('u_b')
    this%v_b => neko_field_registry%get_field('v_b')
    this%w_b => neko_field_registry%get_field('w_b')
    this%p_b => neko_field_registry%get_field('p_b')

    ! Read the json file
    call json_get_or_default(params, 'norm_target', &
         this%norm_target, -1.0_rp)
    call json_get_or_default(params, 'norm_tolerance', &
         this%norm_tolerance, 10.0_rp)

    ! Build the header
    call json_get_or_default(params, 'output_file', &
         file_name, 'power_iterations.csv')
    this%file_output = file_t(trim(file_name))
    write(header_line, '(A)') 'Time, Norm, Scaling'
    call this%file_output%set_header(header_line)

  end subroutine adjoint_pnpn_init

  subroutine adjoint_pnpn_restart(this, dtlag, tlag)
    class(adjoint_pnpn_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: i, n

    n = this%u_adj%dof%size()
    ! Make sure that continuity is maintained (important for interpolation)
    ! Do not do this for lagged rhs
    ! (derivatives are not necessairly coninous across elements)
    call col2(this%u_adj%x, this%c_Xh%mult, this%u_adj%dof%size())
    call col2(this%v_adj%x, this%c_Xh%mult, this%u_adj%dof%size())
    call col2(this%w_adj%x, this%c_Xh%mult, this%u_adj%dof%size())
    call col2(this%p_adj%x, this%c_Xh%mult, this%u_adj%dof%size())
    do i = 1, this%ulag%size()
       call col2(this%ulag%lf(i)%x, this%c_Xh%mult, this%u_adj%dof%size())
       call col2(this%vlag%lf(i)%x, this%c_Xh%mult, this%u_adj%dof%size())
       call col2(this%wlag%lf(i)%x, this%c_Xh%mult, this%u_adj%dof%size())
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u => this%u_adj, v => this%v_adj, w => this%w_adj, &
            ulag => this%ulag, vlag => this%vlag, wlag => this%wlag,&
            p => this%p_adj)
         call device_memcpy(u%x, u%x_d, u%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(v%x, v%x_d, v%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(w%x, w%x_d, w%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(p%x, p%x_d, p%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
              u%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
              u%dof%size(), HOST_TO_DEVICE, sync = .false.)

         call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
              v%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
              v%dof%size(), HOST_TO_DEVICE, sync = .false.)

         call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abx1%x, this%abx1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abx2%x, this%abx2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%aby1%x, this%aby1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%aby2%x, this%aby2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abz1%x, this%abz1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abz2%x, this%abz2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
       end associate
    end if


    call this%gs_Xh%op(this%u_adj, GS_OP_ADD)
    call this%gs_Xh%op(this%v_adj, GS_OP_ADD)
    call this%gs_Xh%op(this%w_adj, GS_OP_ADD)
    call this%gs_Xh%op(this%p_adj, GS_OP_ADD)

    do i = 1, this%ulag%size()
       call this%gs_Xh%op(this%ulag%lf(i), GS_OP_ADD)
       call this%gs_Xh%op(this%vlag%lf(i), GS_OP_ADD)
       call this%gs_Xh%op(this%wlag%lf(i), GS_OP_ADD)
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
    !   call device_opcolv(this%f_adj_x%x_d, this%f_adj_y%x_d, this%f_adj_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_adj_x%x, this%f_adj_y%x, this%f_adj_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !! Add the advection operators to the right-hand-side.
    !call this%adv%compute(u_temp, v_temp, w_temp, &
    !                      this%f_adj_x%x, this%f_adj_y%x, this%f_adj_z%x, &
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
    !   call device_opcolv(this%f_adj_x%x_d, this%f_adj_y%x_d, this%f_adj_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_adj_x%x, this%f_adj_y%x, this%f_adj_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_adj_x%x_d, this%f_adj_y%x_d, this%f_adj_z%x_d, this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_adj_x%x, this%f_adj_y%x, this%f_adj_z%x, this%c_Xh%B, this%msh%gdim, n)
    !end if

    !call this%adv%compute(u_temp, v_temp, w_temp, &
    !                      this%f_adj_x%x, this%f_adj_y%x, this%f_adj_z%x, &
    !                      this%Xh, this%c_Xh, this%dm_Xh%size())
    !this%abx1 = this%f_x
    !this%aby1 = this%f_y
    !this%abz1 = this%f_z

  end subroutine adjoint_pnpn_restart

  subroutine adjoint_pnpn_free(this)
    class(adjoint_pnpn_t), intent(inout) :: this

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

    if (allocated(this%Ax_vel)) then
       deallocate(this%Ax_vel)
    end if

    if (allocated(this%Ax_prs)) then
       deallocate(this%Ax_prs)
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

  end subroutine adjoint_pnpn_free

  !> Advance fluid simulation in time.
  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  !! @param dt_controller timestep controller
  subroutine adjoint_pnpn_step(this, t, tstep, dt, ext_bdf, dt_controller)
    class(adjoint_pnpn_t), target, intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    ! number of degrees of freedom
    integer :: n
    ! Solver results monitors (pressure + 3 velocity)
    type(ksp_monitor_t) :: ksp_results(4)
    ! Extrapolated velocity for the pressure residual
    type(field_t), pointer :: u_e, v_e, w_e
    ! Indices for tracking temporary fields
    integer :: temp_indices(3)

    if (this%freeze) return

    n = this%dm_Xh%size()

    call profiler_start_region('Fluid', 1)
    associate(u => this%u_adj, v => this%v_adj, w => this%w_adj, p => this%p_adj, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_b => this%u_b, v_b => this%v_b, w_b => this%w_b, &
         u_res => this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax_vel => this%Ax_vel, Ax_prs => this%Ax_prs, &
         Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, prs_res => this%prs_res, &
         source_term => this%source_term, &
         vel_res => this%vel_res, sumab => this%sumab, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         vel_projection_dim => this%vel_projection_dim, &
         pr_projection_dim => this%pr_projection_dim, &
         rho => this%rho, mu => this%mu, &
         rho_field => this%rho_field, mu_field => this%mu_field, &
         f_x => this%f_adj_x, f_y => this%f_adj_y, f_z => this%f_adj_z, &
         if_variable_dt => dt_controller%if_variable_dt, &
         dt_last_change => dt_controller%dt_last_change)

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

      call this%adv%compute_adjoint(u, v, w, u_b, v_b, w_b, &
           f_x, f_y, f_z, &
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
      !> Apply the user dirichlet boundary condition
      call this%user_field_bc_vel%update(this%user_field_bc_vel%field_list, &
           this%user_field_bc_vel%bc_list, this%c_Xh, t, tstep, "fluid")

      call this%bc_apply_vel(t, tstep)
      call this%bc_apply_prs(t, tstep)

      ! Update material properties if necessary
      call this%update_material_properties()

      ! Compute pressure.
      call profiler_start_region('Pressure residual', 18)
      call prs_res%compute(p, p_res,&
           u, v, w, &
           u_e, v_e, w_e, &
           f_x, f_y, f_z, &
           c_Xh, gs_Xh, &
           this%bc_prs_surface, this%bc_sym_surface,&
           Ax_prs, ext_bdf%diffusion_coeffs(1), dt, &
           mu_field, rho_field)

      call gs_Xh%op(p_res, GS_OP_ADD)
      call bc_list_apply_scalar(this%bclst_dp, p_res%x, p%dof%size(), t, tstep)
      call profiler_end_region

      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
           'Pressure')

      call this%pc_prs%update()
      call profiler_start_region('Pressure solve', 3)
      ksp_results(1) = &
           this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, this%bclst_dp, gs_Xh)

      call profiler_end_region

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
           this%bclst_dp, gs_Xh, n, tstep, dt_controller)

      call field_add2(p, dp, n)

      ! Compute velocity.
      call profiler_start_region('Velocity residual', 19)
      call vel_res%compute(Ax_vel, u, v, w, &
           u_res, v_res, w_res, &
           p, &
           f_x, f_y, f_z, &
           c_Xh, msh, Xh, &
           mu_field, rho_field, ext_bdf%diffusion_coeffs(1), &
           dt, dm_Xh%size())

      call gs_Xh%op(u_res, GS_OP_ADD)
      call gs_Xh%op(v_res, GS_OP_ADD)
      call gs_Xh%op(w_res, GS_OP_ADD)

      call bc_list_apply_vector(this%bclst_vel_res,&
           u_res%x, v_res%x, w_res%x, dm_Xh%size(),&
           t, tstep)

      !We should implement a bc that takes three field_bcs and implements vector_apply
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call this%bc_field_dirichlet_u%apply_scalar_dev(u_res%x_d, t, tstep)
         call this%bc_field_dirichlet_v%apply_scalar_dev(v_res%x_d, t, tstep)
         call this%bc_field_dirichlet_w%apply_scalar_dev(w_res%x_d, t, tstep)
      else
         call this%bc_field_dirichlet_u%apply_scalar(u_res%x, n, t, tstep)
         call this%bc_field_dirichlet_v%apply_scalar(v_res%x, n, t, tstep)
         call this%bc_field_dirichlet_w%apply_scalar(w_res%x, n, t, tstep)
      end if

      call profiler_end_region

      call this%proj_u%pre_solving(u_res%x, tstep, c_Xh, n, dt_controller)
      call this%proj_v%pre_solving(v_res%x, tstep, c_Xh, n, dt_controller)
      call this%proj_w%pre_solving(w_res%x, tstep, c_Xh, n, dt_controller)

      call this%pc_vel%update()

      call profiler_start_region("Velocity solve", 4)
      ksp_results(2:4) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh)
      call profiler_end_region

      call this%proj_u%post_solving(du%x, Ax_vel, c_Xh, &
           this%bclst_du, gs_Xh, n, tstep, dt_controller)
      call this%proj_v%post_solving(dv%x, Ax_vel, c_Xh, &
           this%bclst_dv, gs_Xh, n, tstep, dt_controller)
      call this%proj_w%post_solving(dw%x, Ax_vel, c_Xh, &
           this%bclst_dw, gs_Xh, n, tstep, dt_controller)

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
              this%bclst_dw, this%bclst_vel_res, Ax_vel, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(tstep, t, dt, ksp_results)

      call this%scratch%relinquish_field(temp_indices)

    end associate
    call profiler_end_region

    ! Compute the norm of the field and determine if we should do a rescale.
    call this%PW_compute_(t, tstep)

  end subroutine adjoint_pnpn_step

  subroutine rescale_fluid(fluid_data, scale)
    use neko_config, only : NEKO_BCKND_DEVICE
    implicit none

    !> Fluid data
    class(adjoint_pnpn_t), intent(inout) :: fluid_data
    !> Scaling factor
    real(kind=rp), intent(in) :: scale

    ! Local variables
    integer :: i

    ! Scale the velocity fields
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%u_adj%x_d, scale, fluid_data%u_adj%size())
       call device_cmult(fluid_data%v_adj%x_d, scale, fluid_data%v_adj%size())
       call device_cmult(fluid_data%w_adj%x_d, scale, fluid_data%w_adj%size())
    else
       call cmult(fluid_data%u_adj%x, scale, fluid_data%u_adj%size())
       call cmult(fluid_data%v_adj%x, scale, fluid_data%v_adj%size())
       call cmult(fluid_data%w_adj%x, scale, fluid_data%w_adj%size())
    end if

    ! Scale the right hand sides
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%f_adj_x%x_d, scale, fluid_data%f_adj_x%size())
       call device_cmult(fluid_data%f_adj_y%x_d, scale, fluid_data%f_adj_y%size())
       call device_cmult(fluid_data%f_adj_z%x_d, scale, fluid_data%f_adj_z%size())
       ! HARRY
       ! maybe the abx's too
       call device_cmult(fluid_data%abx1%x_d, scale, fluid_data%abx1%size())
       call device_cmult(fluid_data%aby1%x_d, scale, fluid_data%aby1%size())
       call device_cmult(fluid_data%abz1%x_d, scale, fluid_data%abz1%size())
       call device_cmult(fluid_data%abx2%x_d, scale, fluid_data%abx2%size())
       call device_cmult(fluid_data%aby2%x_d, scale, fluid_data%aby2%size())
       call device_cmult(fluid_data%abz2%x_d, scale, fluid_data%abz2%size())

    else
       call cmult(fluid_data%f_adj_x%x, scale, fluid_data%f_adj_x%size())
       call cmult(fluid_data%f_adj_y%x, scale, fluid_data%f_adj_y%size())
       call cmult(fluid_data%f_adj_z%x, scale, fluid_data%f_adj_z%size())

       call cmult(fluid_data%abx1%x, scale, fluid_data%abx1%size())
       call cmult(fluid_data%aby1%x, scale, fluid_data%aby1%size())
       call cmult(fluid_data%abz1%x, scale, fluid_data%abz1%size())

       call cmult(fluid_data%abx2%x, scale, fluid_data%abx2%size())
       call cmult(fluid_data%aby2%x, scale, fluid_data%aby2%size())
       call cmult(fluid_data%abz2%x, scale, fluid_data%abz2%size())
    end if

    ! Scale the lag terms
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, fluid_data%ulag%size()
          call device_cmult(fluid_data%ulag%lf(i)%x_d, &
               scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call device_cmult(fluid_data%vlag%lf(i)%x_d, &
               scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call device_cmult(fluid_data%wlag%lf(i)%x_d, &
               scale, fluid_data%wlag%lf(i)%size())
       end do
    else
       do i = 1, fluid_data%ulag%size()
          call cmult(fluid_data%ulag%lf(i)%x, &
               scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call cmult(fluid_data%vlag%lf(i)%x, &
               scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call cmult(fluid_data%wlag%lf(i)%x, &
               scale, fluid_data%wlag%lf(i)%size())
       end do
    end if

  end subroutine rescale_fluid

  function norm(x, y, z, B, volume, n)
    use mpi_f08, only : MPI_SUM, MPI_COMM_WORLD, &
         MPI_IN_PLACE, mpi_allreduce
    use comm, only : MPI_REAL_PRECISION
    use math, only : vlsc3

    real(kind=rp), dimension(n), intent(in) :: x, y, z
    real(kind=rp), dimension(n), intent(in) :: B
    real(kind=rp), intent(in) :: volume
    integer, intent(in) :: n

    real(kind=rp) :: norm

    norm = vlsc3(x, x, B, n) + vlsc3(y, y, B, n) + vlsc3(z, z, B, n)

    call mpi_allreduce(MPI_IN_PLACE, norm, 1, &
         MPI_REAL_PRECISION, MPI_SUM, MPI_COMM_WORLD)

    norm = sqrt(norm / volume)
  end function norm

  function device_norm(x_d, y_d, z_d, B_d, volume, n)
    use neko_config, only : NEKO_BCKND_DEVICE
    use device_math, only : device_vlsc3
    use comm, only : MPI_REAL_PRECISION
    use mpi_f08, only : MPI_SUM, MPI_COMM_WORLD, &
         MPI_IN_PLACE, mpi_allreduce

    implicit none

    type(c_ptr), intent(in) :: x_d, y_d, z_d
    type(c_ptr), intent(in) :: B_d
    real(kind=rp), intent(in) :: volume
    integer, intent(in) :: n

    real(kind=rp) :: device_norm

    device_norm = device_vlsc3(x_d, x_d, B_d, n) + &
         device_vlsc3(y_d, y_d, B_d, n) + &
         device_vlsc3(z_d, z_d, B_d, n)

    call mpi_allreduce(MPI_IN_PLACE, device_norm, 1, &
         MPI_REAL_PRECISION, MPI_SUM, MPI_COMM_WORLD)

    device_norm = sqrt(device_norm / volume)

  end function device_norm

  !> Compute the power_iterations field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine power_iterations_compute(this, t, tstep)
    class(adjoint_pnpn_t), target, intent(inout) :: this

    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Local variables
    real(kind=rp) :: scaling_factor
    real(kind=rp) :: norm_l2, norm_l2_base
    real(kind=rp) :: lambda
    real(kind=rp) :: dt
    character(len=256) :: log_message
    type(vector_t) :: data_line
    integer :: n

    n = this%c_Xh%dof%size()
    if (tstep .eq. 1) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          norm_l2_base = device_norm(this%u_adj%x_d, this%v_adj%x_d, this%w_adj%x_d, &
               this%c_Xh%B_d, this%c_Xh%volume, n)
       else
          norm_l2_base = this%norm_scaling * norm(this%u_adj%x, this%v_adj%x, this%w_adj%x, &
               this%c_Xh%B, this%c_Xh%volume, n)
       end if
       if (this%norm_target .lt. 0.0_rp) then
          this%norm_target = norm_l2_base
       end if

       this%norm_l2_upper = this%norm_tolerance * this%norm_target
       this%norm_l2_lower = this%norm_target / this%norm_tolerance

    end if

    ! Compute the norm of the velocity field and eigenvalue estimate
    if (NEKO_BCKND_DEVICE .eq. 1) then
       norm_l2 = device_norm(this%u_adj%x_d, this%v_adj%x_d, this%w_adj%x_d, &
            this%c_Xh%B_d, this%c_Xh%volume, n)
    else
       norm_l2 = norm(this%u_adj%x, this%v_adj%x, this%w_adj%x, &
            this%c_Xh%B, this%c_Xh%volume, n)
    end if
    norm_l2 = sqrt(this%norm_scaling) * norm_l2
    scaling_factor = 1.0_rp

    ! Rescale the flow if necessary
    if (norm_l2 .gt. this%norm_l2_upper &
         .or. norm_l2 .lt. this%norm_l2_lower) then
       scaling_factor = this%norm_target / norm_l2
       call rescale_fluid(this, scaling_factor)
       norm_l2 = this%norm_target

       if (tstep .eq. 1) then
          scaling_factor = 1.0_rp
       end if
    end if

    ! Log the results
    call neko_log%section('Power Iterations', lvl=NEKO_LOG_DEBUG)

    write (log_message, '(A7,E20.14)') 'Norm: ', norm_l2
    call neko_log%message(log_message, lvl=NEKO_LOG_DEBUG)
    write (log_message, '(A7,E20.14)') 'Scaling: ', scaling_factor
    call neko_log%message(log_message, lvl=NEKO_LOG_DEBUG)

    ! Save to file
    call data_line%init(2)
    data_line%x = [norm_l2, scaling_factor]
    call this%file_output%write(data_line, t)

    call neko_log%end_section('Power Iterations', lvl=NEKO_LOG_DEBUG)
  end subroutine power_iterations_compute


end module adjoint_pnpn
