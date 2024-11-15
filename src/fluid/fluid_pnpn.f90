! Copyright (c) 2022-2024, The Neko Authors
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
  use comm
  use coefs, only : coef_t
  use symmetry, only : symmetry_t
  use field_registry, only : neko_field_registry
  use logger, only: neko_log, LOG_SIZE
  use num_types, only : rp
  use krylov, only : ksp_monitor_t
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t, &
       pnpn_prs_res_factory, pnpn_vel_res_factory, &
       pnpn_prs_res_stress_factory, pnpn_vel_res_stress_factory
  use rhs_maker, only : rhs_maker_sumab_t, rhs_maker_bdf_t, rhs_maker_ext_t, &
       rhs_maker_oifs_t, rhs_maker_sumab_fctry, rhs_maker_bdf_fctry, &
       rhs_maker_ext_fctry, rhs_maker_oifs_fctry
  use fluid_volflow, only : fluid_volflow_t
  use fluid_scheme, only : fluid_scheme_t
  use device_mathops, only : device_opcolv, device_opadd2cm
  use fluid_aux, only : fluid_step_info
  use time_scheme_controller, only : time_scheme_controller_t
  use projection, only : projection_t
  use device, only : device_memcpy, HOST_TO_DEVICE
  use advection, only : advection_t, advection_factory
  use profiler, only : profiler_start_region, profiler_end_region
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file
  use ax_product, only : ax_t, ax_helm_factory
  use field, only : field_t
  use dirichlet, only : dirichlet_t
  use facet_normal, only : facet_normal_t
  use non_normal, only : non_normal_t
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use time_step_controller, only : time_step_controller_t
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : col2, glsum
  use mathops, only : opadd2cm, opcolv
  use bc_list, only: bc_list_t
  use zero_dirichlet, only : zero_dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use utils, only : neko_error
  use field_math, only : field_add2, field_copy
  use bc, only : bc_t
  implicit none
  private


  type, public, extends(fluid_scheme_t) :: fluid_pnpn_t

     ! The right-hand sides in the linear solves.
     type(field_t) :: p_res, u_res, v_res, w_res

     !> The unknowns in the linear solves.
     type(field_t) :: dp, du, dv, dw

     !
     ! Implicit operators, i.e. the left-hand-side of the Helmholz problem.
     !

     ! Coupled Helmholz operator for velocity
     class(ax_t), allocatable :: Ax_vel
     ! Helmholz operator for pressure
     class(ax_t), allocatable :: Ax_prs

     !
     ! Projections for solver speed-up
     !

     !> Pressure projection
     type(projection_t) :: proj_prs
     !> X velocity projection
     type(projection_t) :: proj_u
     !> Y velocity projection
     type(projection_t) :: proj_v
     !> Z velocity projection
     type(projection_t) :: proj_w

     !
     ! Special Karniadakis scheme boundary conditions in the pressure equation
     !

     !> Surface term in pressure rhs. Masks all strong velocity bcs.
     type(facet_normal_t) :: bc_prs_surface

     !> Surface term in pressure rhs. Masks symmetry bcs.
     type(facet_normal_t) :: bc_sym_surface

     type(dong_outflow_t) :: bc_dong           !< Dong outflow condition


     type(zero_dirichlet_t) :: bc_vel_res   !< Dirichlet condition vel. res.
     type(zero_dirichlet_t) :: bc_field_dirichlet_p   !< Dirichlet condition vel. res.
     type(zero_dirichlet_t) :: bc_field_dirichlet_u   !< Dirichlet condition vel. res.
     type(zero_dirichlet_t) :: bc_field_dirichlet_v   !< Dirichlet condition vel. res.
     type(zero_dirichlet_t) :: bc_field_dirichlet_w   !< Dirichlet condition vel. res.
     type(non_normal_t) :: bc_vel_res_non_normal   !< Dirichlet condition vel. res.


     !
     ! Boundary condition lists
     !
     type(bc_list_t) :: bclst_vel_res
     type(bc_list_t) :: bclst_du
     type(bc_list_t) :: bclst_dv
     type(bc_list_t) :: bclst_dw
     type(bc_list_t) :: bclst_dp

     ! TEMP
     type(bc_list_t) :: bcs_prs

     class(advection_t), allocatable :: adv

     ! Time interpolation scheme
     logical :: oifs

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2
     ! Advection terms for the oifs method
     type(field_t) :: advx, advy, advz

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

     !> Contributions to the RHS from the OIFS method
     class(rhs_maker_oifs_t), allocatable :: makeoifs

     !> Adjust flow volume
     type(fluid_volflow_t) :: vol_flow

   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_pnpn_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_pnpn_free
     !> Perform a single time-step of the scheme.
     procedure, pass(this) :: step => fluid_pnpn_step
     !> Restart from a previous solution.
     procedure, pass(this) :: restart => fluid_pnpn_restart
     !> Set up boundary conditions.
     procedure, pass(this) :: pnpn_setup_bcs => fluid_pnpn_setup_bcs
  end type fluid_pnpn_t

contains

  subroutine fluid_pnpn_init(this, msh, lx, params, user, time_scheme)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    integer :: i
    real(kind=rp) :: abs_tol
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr, integer_val, solver_maxiter
    character(len=:), allocatable :: solver_type, precon_type
    logical :: monitor

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, .true., .true., scheme, user)

    ! Add pressure field to the registery. For this scheme it is in the same
    ! Xh as the velocity
    call neko_field_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_field_registry%get_field('p')

    !
    ! Select governing equations via associated residual and Ax types
    !

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

    ! Setup backend dependent summations of the OIFS method
    call rhs_maker_oifs_fctry(this%makeoifs)

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
      call this%advx%init(dm_Xh, "advx")
      call this%advy%init(dm_Xh, "advy")
      call this%advz%init(dm_Xh, "advz")
      this%abx1 = 0.0_rp
      this%aby1 = 0.0_rp
      this%abz1 = 0.0_rp
      this%abx2 = 0.0_rp
      this%aby2 = 0.0_rp
      this%abz2 = 0.0_rp
      this%advx = 0.0_rp
      this%advy = 0.0_rp
      this%advz = 0.0_rp
    end associate

    call this%du%init(this%dm_Xh, 'du')
    call this%dv%init(this%dm_Xh, 'dv')
    call this%dw%init(this%dm_Xh, 'dw')
    call this%dp%init(this%dm_Xh, 'dp')


    !
    ! Boundary conditions
    !

    ! Initialize velocity surface terms in pressure rhs. Masks all strong
    ! velocity bcs.
    call this%bc_prs_surface%init(this%c_Xh, params)
    do i = 1, this%bcs%size()
       if (this%bcs%strong(i) .eqv. .true.) then
          call this%bc_prs_surface%mark_facets( &
               this%bcs%items(i)%obj%marked_facet)
       end if
    end do
    call this%bc_prs_surface%finalize()

    ! Initialize symmetry surface terms in pressure rhs. Masks symmetry bcs.
    call this%bc_sym_surface%init(this%c_Xh, params)
    do i = 1, this%bcs%size()
       select type (vel_bc => this%bcs%items(i)%obj)
       type is (symmetry_t)
          write(*,*) "MARKING PRESSURE SYMMETRY"
          call this%bc_sym_surface%mark_facets( &
               this%bcs%items(i)%obj%marked_facet)
       end select
    end do
    call this%bc_sym_surface%finalize()

    !call this%bc_prs_surface%mark_zones_from_list('v', this%bc_labels)

    ! This impacts the rhs of the pressure, need to check what is correct to add
    ! here

    !call this%bc_prs_surface%mark_zones_from_list('d_vel_u', this%bc_labels)
    !call this%bc_prs_surface%mark_zones_from_list('d_vel_v', this%bc_labels)
    !call this%bc_prs_surface%mark_zones_from_list('d_vel_w', this%bc_labels)
    !call this%bc_prs_surface%finalize()

    ! Generate pressure boundary condition list
    call this%pnpn_setup_bcs(user)

    ! Add all strong bcs to the list

!    call this%bc_prs%mark_zones_from_list('o', this%bc_labels)
!    call this%bc_prs%mark_zones_from_list('on', this%bc_labels)

    ! Field dirichlet pressure bc
    call this%user_field_bc_prs%init_base(this%c_Xh)
!    call this%user_field_bc_prs%mark_zones_from_list('d_pres', this%bc_labels)
!    call this%user_field_bc_prs%finalize()
!    call MPI_Allreduce(this%user_field_bc_prs%msk(0), integer_val, 1, &
!         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

!    call this%bc_field_dirichlet_p%init(this%c_Xh, params)
!    call this%bc_field_dirichlet_p%mark_zones_from_list('on+dong', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_p%mark_zones_from_list('o+dong', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_p%mark_zones_from_list('d_pres', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_p%finalize()

!    if (integer_val .gt. 0)  call this%user_field_bc_prs%init_field('d_pres')
!    call this%bclst_prs%append(this%user_field_bc_prs)
!    call this%user_field_bc_vel%bc_list%append(this%user_field_bc_prs)

!    call this%bclst_prs%append(this%bc_prs)


!    call this%bc_dong%init(this%c_Xh, params)
!    call this%bc_dong%mark_zones_from_list('o+dong', this%bc_labels)
!    call this%bc_dong%mark_zones_from_list('on+dong', this%bc_labels)
!    call this%bc_dong%finalize()



    !call this%bc_sym_surface%mark_zones_from_list('sym', this%bc_labels)
    ! Same here, should du, dv, dw be marked here?
    !call this%bc_sym_surface%finalize()


    ! Initialize dirichlet bcs for velocity residual
    call this%bc_vel_res_non_normal%init(this%c_Xh, params)
!    call this%bc_vel_res_non_normal%mark_zones_from_list('on', this%bc_labels)
!    call this%bc_vel_res_non_normal%mark_zones_from_list('on+dong', &
!                                                         this%bc_labels)
    call this%bc_vel_res_non_normal%finalize()



    ! Mark Dirichlet bcs for pressure
    call this%bclst_dp%init()

    call this%bclst_prs%init()
    do i = 1, this%bcs_prs%size()
       if (this%bcs_prs%strong(i) .eqv. .true.) then
          call this%bclst_dp%append(this%bcs_prs%items(i)%obj)
       end if
    end do

!    call this%bclst_dp%append(this%bc_field_dirichlet_p)
    !Add 0 prs bcs
!    call this%bclst_dp%append(this%bc_prs)

!    call this%bc_field_dirichlet_u%init(this%c_Xh, params)
!    call this%bc_field_dirichlet_u%mark_zones_from_list('d_vel_u', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_u%finalize()

!    call this%bc_field_dirichlet_v%init(this%c_Xh, params)
!    call this%bc_field_dirichlet_v%mark_zones_from_list('d_vel_v', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_v%finalize()

!    call this%bc_field_dirichlet_w%init(this%c_Xh, params)
!    call this%bc_field_dirichlet_w%mark_zones_from_list('d_vel_w', &
!                                         this%bc_labels)
!    call this%bc_field_dirichlet_w%finalize()

!    call this%bc_vel_res%init(this%c_Xh, params)
!    call this%bc_vel_res%mark_zones_from_list('v', this%bc_labels)
!    call this%bc_vel_res%mark_zones_from_list('w', this%bc_labels)
!    call this%bc_vel_res%finalize()

    ! Populate lists for the velocity residual and solution increments
    call this%bclst_vel_res%init()
    call this%bclst_du%init()
    call this%bclst_dv%init()
    call this%bclst_dw%init()

    ! Add all strong velocity bcs.
    do i = 1, this%bcs%size()
       if (this%bcs%strong(i) .eqv. .true.) then

         ! We need to treat bcs with nested bcs in a delicate way
         select type (vel_bc => this%bcs%items(i)%obj)
         type is (symmetry_t)
            write(*,*) "MARKING SYMMETRY IN VELOCITY BLISTS"
            call this%bclst_vel_res%append(vel_bc)
            call this%bclst_du%append(vel_bc%bc_x)
            call this%bclst_dv%append(vel_bc%bc_y)
            call this%bclst_dw%append(vel_bc%bc_z)
         class default
            call this%bclst_vel_res%append(this%bcs%items(i)%obj)
            call this%bclst_du%append(this%bcs%items(i)%obj)
            call this%bclst_dv%append(this%bcs%items(i)%obj)
            call this%bclst_dw%append(this%bcs%items(i)%obj)
         end select
       end if
    end do

    !call this%bclst_vel_res%append(this%bc_vel_res)
    !call this%bclst_vel_res%append(this%bc_vel_res_non_normal)

!    call this%bclst_vel_res%append(this%bc_sym)
!    call this%bclst_vel_res%append(this%bc_sh%symmetry)
!    call this%bclst_vel_res%append(this%bc_wallmodel%symmetry)

    !Initialize bcs for u, v, w velocity components
!    call this%bclst_du%append(this%bc_sym%bc_x)
!    call this%bclst_du%append(this%bc_sh%symmetry%bc_x)
!    call this%bclst_du%append(this%bc_wallmodel%symmetry%bc_x)
!    call this%bclst_du%append(this%bc_vel_res_non_normal%bc_x)
!    call this%bclst_du%append(this%bc_vel_res)
!    call this%bclst_du%append(this%bc_field_dirichlet_u)

!    call this%bclst_dv%append(this%bc_sym%bc_y)
!    call this%bclst_dv%append(this%bc_sh%symmetry%bc_y)
!    call this%bclst_dv%append(this%bc_wallmodel%symmetry%bc_y)
!    call this%bclst_dv%append(this%bc_vel_res_non_normal%bc_y)
!    call this%bclst_dv%append(this%bc_vel_res)
!    call this%bclst_dv%append(this%bc_field_dirichlet_v)

!    call this%bclst_dw%append(this%bc_sym%bc_z)
!    call this%bclst_dw%append(this%bc_sh%symmetry%bc_z)
!    call this%bclst_dw%append(this%bc_wallmodel%symmetry%bc_z)
!    call this%bclst_dw%append(this%bc_vel_res_non_normal%bc_z)
!    call this%bclst_dw%append(this%bc_vel_res)
!    call this%bclst_dw%append(this%bc_field_dirichlet_w)


    ! Intialize projection space

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

    ! Determine the time-interpolation scheme
    call json_get_or_default(params, 'case.numerics.oifs', this%oifs, .false.)

    ! Initialize the advection factory
    call advection_factory(this%adv, params, this%c_Xh, &
                           this%ulag, this%vlag, this%wlag, &
                           this%chkp%dtlag, this%chkp%tlag, time_scheme)

    if (params%valid_path('case.fluid.flow_rate_force')) then
       call this%vol_flow%init(this%dm_Xh, params)
    end if

    ! Setup pressure solver
    call neko_log%section("Pressure solver")

    call json_get_or_default(params, &
                            'case.fluid.pressure_solver.max_iterations', &
                            solver_maxiter, 800)
    call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
    call json_get(params, 'case.fluid.pressure_solver.preconditioner', &
         precon_type)
    call json_get(params, 'case.fluid.pressure_solver.absolute_tolerance', &
         abs_tol)
     call json_get_or_default(params, 'case.fluid.velocity_solver.monitor', &
          monitor, .false.)
    call neko_log%message('Type       : ('// trim(solver_type) // &
          ', ' // trim(precon_type) // ')')
    write(log_buf, '(A,ES13.6)') 'Abs tol    :',  abs_tol
    call neko_log%message(log_buf)

    call this%solver_factory(this%ksp_prs, this%dm_Xh%size(), &
         solver_type, solver_maxiter, abs_tol, monitor)
    call this%precon_factory(this%pc_prs, this%ksp_prs, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs, precon_type)

    call neko_log%end_section()

  end subroutine fluid_pnpn_init

  subroutine fluid_pnpn_restart(this, dtlag, tlag)
    class(fluid_pnpn_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
    type(field_t) :: u_temp, v_temp, w_temp
    integer :: i, n

    n = this%u%dof%size()
    if (allocated(this%chkp%previous_mesh%elements) .or. &
        this%chkp%previous_Xh%lx .ne. this%Xh%lx) then
       call col2(this%u%x, this%c_Xh%mult, this%u%dof%size())
       call col2(this%v%x, this%c_Xh%mult, this%u%dof%size())
       call col2(this%w%x, this%c_Xh%mult, this%u%dof%size())
       call col2(this%p%x, this%c_Xh%mult, this%u%dof%size())
       do i = 1, this%ulag%size()
          call col2(this%ulag%lf(i)%x, this%c_Xh%mult, this%u%dof%size())
          call col2(this%vlag%lf(i)%x, this%c_Xh%mult, this%u%dof%size())
          call col2(this%wlag%lf(i)%x, this%c_Xh%mult, this%u%dof%size())
       end do
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u => this%u, v => this%v, w => this%w, &
            ulag => this%ulag, vlag => this%vlag, wlag => this%wlag,&
            p => this%p)
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
         call device_memcpy(this%advx%x, this%advx%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%advy%x, this%advy%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%advz%x, this%advz%x_d, &
                            w%dof%size(), HOST_TO_DEVICE, sync = .false.)
       end associate
    end if
    ! Make sure that continuity is maintained (important for interpolation)
    ! Do not do this for lagged rhs
    ! (derivatives are not necessairly coninous across elements)

    if (allocated(this%chkp%previous_mesh%elements) &
         .or. this%chkp%previous_Xh%lx .ne. this%Xh%lx) then
       call this%gs_Xh%op(this%u, GS_OP_ADD)
       call this%gs_Xh%op(this%v, GS_OP_ADD)
       call this%gs_Xh%op(this%w, GS_OP_ADD)
       call this%gs_Xh%op(this%p, GS_OP_ADD)

       do i = 1, this%ulag%size()
          call this%gs_Xh%op(this%ulag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%vlag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%wlag%lf(i), GS_OP_ADD)
       end do
    end if

    !! If we would decide to only restart from lagged fields instead of saving
    !! abx1, aby1 etc.
    !! Observe that one also needs to recompute the focing at the old time steps
    !u_temp = this%ulag%lf(2)
    !v_temp = this%vlag%lf(2)
    !w_temp = this%wlag%lf(2)
    !! Compute the source terms
    !call this%source_term%compute(tlag(2), -1)
    !
    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, &
    !                      this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, &
    !               this%c_Xh%B, this%msh%gdim, n)
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
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, &
    !                      this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, &
    !               this%c_Xh%B, this%msh%gdim, n)
    !end if

    !! Pre-multiply the source terms with the mass matrix.
    !if (NEKO_BCKND_DEVICE .eq. 1) then
    !   call device_opcolv(this%f_x%x_d, this%f_y%x_d, this%f_z%x_d, &
    !                      this%c_Xh%B_d, this%msh%gdim, n)
    !else
    !   call opcolv(this%f_x%x, this%f_y%x, this%f_z%x, &
    !               this%c_Xh%B, this%msh%gdim, n)
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
    call this%bclst_vel_res%free()
    call this%bclst_dp%free()
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

    call this%advx%free()
    call this%advy%free()
    call this%advz%free()

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

    if (allocated(this%makeoifs)) then
       deallocate(this%makeoifs)
    end if

    call this%vol_flow%free()

  end subroutine fluid_pnpn_free

  !> Advance fluid simulation in time.
  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  !! @param dt_controller timestep controller
  subroutine fluid_pnpn_step(this, t, tstep, dt, ext_bdf, dt_controller)
    class(fluid_pnpn_t), target, intent(inout) :: this
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
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
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
         rho => this%rho, mu => this%mu, oifs => this%oifs, &
         rho_field => this%rho_field, mu_field => this%mu_field, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
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

      ! Add Neumann bc contributions to the RHS
      !call this%bclst_vel_neumann%apply_vector(f_x%x, f_y%x, f_z%x, &
      !     this%dm_Xh%size(), t, tstep)

      ! Compute the grandient jump penalty term
      if (this%if_gradient_jump_penalty .eqv. .true.) then
         call this%gradient_jump_penalty_u%compute(u, v, w, u)
         call this%gradient_jump_penalty_v%compute(u, v, w, v)
         call this%gradient_jump_penalty_w%compute(u, v, w, w)
         call this%gradient_jump_penalty_u%perform(f_x)
         call this%gradient_jump_penalty_v%perform(f_y)
         call this%gradient_jump_penalty_w%perform(f_z)
      end if

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
                                    rho, ext_bdf%advection_coeffs, n)

         ! Now, the source terms from the previous time step are added to the RHS.
         call makeoifs%compute_fluid(this%advx%x, this%advy%x, this%advz%x, &
                                     f_x%x, f_y%x, f_z%x, &
                                     rho, dt, n)
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
                              rho, ext_bdf%advection_coeffs, n)

         ! Add the RHS contributions coming from the BDF scheme.
         call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
                              u, v, w, c_Xh%B, rho, dt, &
                              ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
      end if

      call ulag%update()
      call vlag%update()
      call wlag%update()

      !> We assume that no change of boundary conditions
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply the user dirichlet boundary condition
      !call this%user_field_bc_vel%update(this%user_field_bc_vel%field_list, &
      !        this%user_field_bc_vel%bc_list, this%c_Xh, t, tstep, "fluid")

      call this%bc_apply_vel(t, tstep)
      call this%bc_apply_prs(t, tstep)

      ! Update material properties if necessary
      call this%update_material_properties()

      ! Compute pressure.
      call profiler_start_region('Pressure_residual', 18)
      call prs_res%compute(p, p_res,&
                           u, v, w, &
                           u_e, v_e, w_e, &
                           f_x, f_y, f_z, &
                           c_Xh, gs_Xh, &
                           this%bc_prs_surface, this%bc_sym_surface,&
                           Ax_prs, ext_bdf%diffusion_coeffs(1), dt, &
                           mu_field, rho_field)

      call gs_Xh%op(p_res, GS_OP_ADD)
      call this%bclst_dp%apply_scalar(p_res%x, p%dof%size(), t, tstep)
      call profiler_end_region('Pressure_residual', 18)

      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
                                     'Pressure')

      call this%pc_prs%update()
      call profiler_start_region('Pressure_solve', 3)
      ksp_results(1) = &
         this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, this%bclst_dp, gs_Xh)

      call profiler_end_region('Pressure_solve', 3)

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
                                 this%bclst_dp, gs_Xh, n, tstep, dt_controller)

      call field_add2(p, dp, n)

      ! Compute velocity.
      call profiler_start_region('Velocity_residual', 19)
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

      call this%bclst_vel_res%apply_vector(u_res%x, v_res%x, w_res%x, &
           dm_Xh%size(), t, tstep)

      ! We should implement a bc that takes three field_bcs and implements
      ! vector_apply
!      if (NEKO_BCKND_DEVICE .eq. 1) then
!         call this%bc_field_dirichlet_u%apply_scalar_dev(u_res%x_d, t, tstep)
!         call this%bc_field_dirichlet_v%apply_scalar_dev(v_res%x_d, t, tstep)
!         call this%bc_field_dirichlet_w%apply_scalar_dev(w_res%x_d, t, tstep)
!      else
!         call this%bc_field_dirichlet_u%apply_scalar(u_res%x, n, t, tstep)
!         call this%bc_field_dirichlet_v%apply_scalar(v_res%x, n, t, tstep)
!         call this%bc_field_dirichlet_w%apply_scalar(w_res%x, n, t, tstep)
!      end if

      call profiler_end_region('Velocity_residual', 19)

      call this%proj_u%pre_solving(u_res%x, tstep, c_Xh, n, dt_controller)
      call this%proj_v%pre_solving(v_res%x, tstep, c_Xh, n, dt_controller)
      call this%proj_w%pre_solving(w_res%x, tstep, c_Xh, n, dt_controller)

      call this%pc_vel%update()

      call profiler_start_region("Velocity_solve", 4)
      ksp_results(2:4) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh, &
           this%ksp_vel%max_iter)
      call profiler_end_region("Velocity_solve", 4)

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
              c_Xh, gs_Xh, ext_bdf, rho, mu, dt, &
              this%bclst_dp, this%bclst_du, this%bclst_dv, &
              this%bclst_dw, this%bclst_vel_res, Ax_vel, Ax_prs, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(tstep, t, dt, ksp_results)

      call this%scratch%relinquish_field(temp_indices)

    end associate
    call profiler_end_region('Fluid', 1)
  end subroutine fluid_pnpn_step

  subroutine fluid_pnpn_setup_bcs(this, user)
    class(fluid_pnpn_t), intent(inout) :: this
    type(user_t), target, intent(in) :: user
    integer :: i, j, n_bcs
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found

    if (this%params%valid_path('case.fluid.boundary_conditions')) then
       call this%params%info('case.fluid.boundary_conditions', n_children=n_bcs)
       call this%params%get_core(core)
       call this%params%get('case.fluid.boundary_conditions', bc_object, found)

       call this%bcs_prs%init(n_bcs)

       write(*,*) "PRESSURE BCS"
       j = 1
       do i=1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)

          write(*,*) "i", i
          call pressure_bc_factory(this%bcs_prs%items(j)%obj, bc_subdict, &
               this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (allocated(this%bcs_prs%items(j)%obj)) then
             write(*,*) "Allocated", j
             if (this%bcs_prs%strong(j)) then
!                this%n_strong = this%n_strong + 1
             end if
             j = j + 1
             this%bcs_prs%size_ = this%bcs_prs%size_ + 1

          end if

          write(*,*) "Done", i
       end do
       write(*,*) "N PRESSURE BCS", j-1, this%bcs_prs%size()
    end if
  end subroutine

  subroutine pressure_bc_factory(object, json, coef, user)
    class(bc_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: zone_index

    call json_get(json, "type", type)

    if (trim(type) .eq. "outflow") then
       allocate(zero_dirichlet_t::object)
    else if (trim(type) .eq. "dong_outflow") then
       allocate(dong_outflow_t::object)
!    else if (trim(type) .eq. "no_slip") then
!       allocate(zero_dirichlet_t::object)
!    else if (trim(type) .eq. "normal_outlet") then
!       allocate(non_normal_t::object)
!    else if (trim(type) .eq. "blasius_profile") then
!       allocate(blasius_t::object)
!    else
!       call neko_error("Unknown boundary condition for the fluid.")
    else
      return
    end if

    call json_get(json, "zone_index", zone_index)
    call object%init(coef, json)
    call object%mark_zone(coef%msh%labeled_zones(zone_index))
    call object%finalize()
  end subroutine

end module fluid_pnpn
