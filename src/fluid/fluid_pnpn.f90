! Copyright (c) 2022-2025, The Neko Authors
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
  use, intrinsic :: iso_fortran_env, only : error_unit
  use coefs, only : coef_t
  use registry, only : neko_registry
  use logger, only : neko_log, LOG_SIZE
  use num_types, only : rp
  use krylov, only : ksp_monitor_t
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t, &
       pnpn_prs_res_factory, pnpn_vel_res_factory, &
       pnpn_prs_res_stress_factory, pnpn_vel_res_stress_factory
  use rhs_maker, only : rhs_maker_sumab_t, rhs_maker_bdf_t, rhs_maker_ext_t, &
       rhs_maker_oifs_t, rhs_maker_sumab_fctry, rhs_maker_bdf_fctry, &
       rhs_maker_ext_fctry, rhs_maker_oifs_fctry
  use fluid_volflow, only : fluid_volflow_t
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use scratch_registry, only : neko_scratch_registry
  use device_mathops, only : device_opcolv, device_opadd2cm
  use fluid_aux, only : fluid_step_info
  use projection, only : projection_t
  use projection_vel, only : projection_vel_t
  use device, only : device_memcpy, HOST_TO_DEVICE, device_event_sync, &
       glb_cmd_event
  use advection, only : advection_t, advection_factory
  use profiler, only : profiler_start_region, profiler_end_region
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use json_module, only : json_file
  use ax_product, only : ax_t, ax_helm_factory
  use field, only : field_t
  use dirichlet, only : dirichlet_t
  use shear_stress, only : shear_stress_t
  use wall_model_bc, only : wall_model_bc_t
  use facet_normal, only : facet_normal_t
  use non_normal_aligned, only : non_normal_aligned_t
  use non_normal, only : non_normal_t
  use symmetry_aligned, only : symmetry_aligned_t
  use symmetry, only : symmetry_t
  use checkpoint, only : chkp_t
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use time_step_controller, only : time_step_controller_t
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use mathops, only : opadd2cm, opcolv
  use zero_dirichlet, only : zero_dirichlet_t
  use utils, only : neko_error, neko_type_error
  use field_math, only : field_add2, field_copy
  use bc, only : bc_t, mixed_bc_t, BC_TYPES
  use bc_resolver, only : scalar_bc_resolver_t, vector_bc_resolver_t, &
       segregated_vector_bc_resolver_t, coupled_vector_bc_resolver_t
  use file, only : file_t
  use operators, only : ortho, rotate_cyc
  use opr_device, only : device_ortho
  use time_state, only : time_state_t
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MAX, MPI_LOR, &
       MPI_INTEGER, MPI_LOGICAL
  implicit none
  private



  type, public, extends(fluid_scheme_incompressible_t) :: fluid_pnpn_t

     !> The right-hand sides in the linear solves.
     type(field_t) :: p_res, u_res, v_res, w_res

     !> The unknowns in the linear solves, i.e. the solution increments with
     !! respect to the previous time-step.
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
     type(projection_vel_t) :: proj_vel

     !
     ! Special Karniadakis scheme boundary conditions in the pressure equation
     !

     !> Surface term in pressure rhs. Masks all strong velocity bcs.
     type(facet_normal_t) :: bc_prs_surface

     !> Surface term in pressure rhs. Masks symmetry bcs.
     type(facet_normal_t) :: bc_sym_surface

     !> Boundary conditions resolver for velocity constraints.
     class(vector_bc_resolver_t), allocatable :: bcs_vel_resolver
     !> Boundary conditions resolver for pressure constraints.
     type(scalar_bc_resolver_t) :: bcs_prs_resolver


     ! Checker for wether we have a strong pressure bc. If not, the pressure
     ! is demeaned at every time step.
     logical :: prs_dirichlet = .false.


     ! The advection operator.
     class(advection_t), allocatable :: adv

     ! Time OIFS interpolation scheme for advection.
     logical :: oifs

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2

     ! Advection terms for the oifs method
     type(field_t) :: advx, advy, advz

     !> Pressure residual equation for computing `p_res`.
     class(pnpn_prs_res_t), allocatable :: prs_res

     !> Velocity residual equation for computing `u_res`, `v_res`, `w_res`.
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

     !> Whether to use the full formulation of the viscous stress term
     logical :: full_stress_formulation = .false.

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
     procedure, pass(this) :: setup_bcs => fluid_pnpn_setup_bcs
     !> Write a field with boundary condition specifications.
     procedure, pass(this) :: write_boundary_conditions => &
          fluid_pnpn_write_boundary_conditions
  end type fluid_pnpn_t

  interface
     !> Boundary condition factory for pressure.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `scalar_pnpn` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine pressure_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_pnpn_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), target, intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine pressure_bc_factory
  end interface

  interface
     !> Boundary condition factory for velocity
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `scalar_pnpn` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine velocity_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_pnpn_t), intent(inout) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), target, intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine velocity_bc_factory
  end interface

contains

  subroutine fluid_pnpn_init(this, msh, lx, params, user, chkp)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    integer :: i
    class(bc_t), pointer :: bc_i, vel_bc
    real(kind=rp) :: abs_tol
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr, integer_val, solver_maxiter
    character(len=:), allocatable :: solver_type, precon_type
    logical :: monitor, found
    logical :: advection
    type(json_file) :: numerics_params, precon_params

    call this%free()

    ! Initialize base class
    call this%init_base(msh, lx, params, scheme, user, .true.)

    ! Add pressure field to the registry. For this scheme it is in the same
    ! Xh as the velocity
    call neko_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_registry%get_field('p')

    !
    ! Select governing equations via associated residual and Ax types
    !

    call json_get_or_lookup(params, 'case.numerics.time_order', integer_val)
    allocate(this%ext_bdf)
    call this%ext_bdf%init(integer_val)

    call json_get_or_default(params, "case.fluid.full_stress_formulation", &
         this%full_stress_formulation, .false.)


    call json_get_or_default(params, "case.fluid.cyclic", this%c_Xh%cyclic, &
         .false.)
    call this%c_Xh%generate_cyclic_bc()

    if (this%full_stress_formulation) then
       ! Setup backend dependent Ax routines
       call ax_helm_factory(this%Ax_vel, full_formulation = .true.)

       ! Setup backend dependent prs residual routines
       call pnpn_prs_res_stress_factory(this%prs_res)

       ! Setup backend dependent vel residual routines
       call pnpn_vel_res_stress_factory(this%vel_res)

       ! Allocate coupled resolver for velocity boundary conditions
       allocate(coupled_vector_bc_resolver_t :: this%bcs_vel_resolver)
    else
       ! Setup backend dependent Ax routines
       call ax_helm_factory(this%Ax_vel, full_formulation = .false.)

       ! Setup backend dependent prs residual routines
       call pnpn_prs_res_factory(this%prs_res)

       ! Setup backend dependent vel residual routines
       call pnpn_vel_res_factory(this%vel_res)

       ! Allocate segregated resolver for velocity boundary conditions
       allocate(segregated_vector_bc_resolver_t :: this%bcs_vel_resolver)
    end if

    ! Initialize the velocity bc resolver
    call this%bcs_vel_resolver%init(this%c_Xh)

    if (params%valid_path('case.fluid.nut_field')) then
       if (.not. this%full_stress_formulation) then
          call neko_error("You need to set full_stress_formulation to " // &
               "true for the fluid to have a spatially varying " // &
               "viscocity field.")
       end if
       call json_get(params, 'case.fluid.nut_field', this%nut_field_name)
    else
       this%nut_field_name = ""
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
    end associate

    call this%du%init(this%dm_Xh, 'du')
    call this%dv%init(this%dm_Xh, 'dv')
    call this%dw%init(this%dm_Xh, 'dw')
    call this%dp%init(this%dm_Xh, 'dp')

    call neko_log%section("Fluid boundary conditions")
    ! Set up boundary conditions
    call this%setup_bcs(user, params)

    ! Check if we need to output boundaries
    call json_get_or_default(params, 'case.output_boundary', found, .false.)
    if (found) call this%write_boundary_conditions()
    call neko_log%end_section()

    call this%proj_prs%init(this%dm_Xh%size(), this%pr_projection_dim, &
         this%pr_projection_activ_step, this%pr_projection_reorthogonalize_basis)

    call this%proj_vel%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)



    ! Determine the time-interpolation scheme
    call json_get_or_default(params, 'case.numerics.oifs', this%oifs, .false.)
    if (params%valid_path('case.fluid.flow_rate_force')) then
       call this%vol_flow%init(this%dm_Xh, params)
    end if

    ! Setup pressure solver
    call neko_log%section("Pressure solver")

    call json_get_or_lookup_or_default(params, &
         'case.fluid.pressure_solver.max_iterations', &
         solver_maxiter, 800)
    call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
    call json_get(params, 'case.fluid.pressure_solver.preconditioner.type', &
         precon_type)
    call json_get(params, &
         'case.fluid.pressure_solver.preconditioner', precon_params)
    call json_get_or_lookup(params, &
         'case.fluid.pressure_solver.absolute_tolerance', &
         abs_tol)
    call json_get_or_default(params, 'case.fluid.pressure_solver.monitor', &
         monitor, .false.)
    call neko_log%message('Type       : ('// trim(solver_type) // &
         ', ' // trim(precon_type) // ')')
    write(log_buf, '(A,ES13.6)') 'Abs tol    :', abs_tol
    call neko_log%message(log_buf)

    call this%solver_factory(this%ksp_prs, this%dm_Xh%size(), &
         solver_type, solver_maxiter, abs_tol, monitor)
    call this%precon_factory_(this%pc_prs, this%ksp_prs, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_prs, &
         precon_type, precon_params)
    call neko_log%end_section()

    ! Initialize the advection factory
    call json_get_or_default(params, 'case.fluid.advection', advection, .true.)
    call json_get(params, 'case.numerics', numerics_params)
    call advection_factory(this%adv, numerics_params, this%c_Xh, &
         this%ulag, this%vlag, this%wlag, &
         chkp%dtlag, chkp%tlag, this%ext_bdf, &
         .not. advection)
    ! Should be in init_base maybe?
    this%chkp => chkp
    ! This is probably scheme specific
    call this%chkp%add_fluid(this%u, this%v, this%w, this%p)

    this%chkp%abx1 => this%abx1
    this%chkp%abx2 => this%abx2
    this%chkp%aby1 => this%aby1
    this%chkp%aby2 => this%aby2
    this%chkp%abz1 => this%abz1
    this%chkp%abz2 => this%abz2
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)



    call neko_log%end_section()

  end subroutine fluid_pnpn_init

  subroutine fluid_pnpn_restart(this, chkp)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    real(kind=rp) :: dtlag(10), tlag(10)
    type(field_t) :: u_temp, v_temp, w_temp
    integer :: i, j, n

    dtlag = chkp%dtlag
    tlag = chkp%tlag

    n = this%u%dof%size()
    if (allocated(chkp%previous_mesh%elements) .or. &
         chkp%previous_Xh%lx .ne. this%Xh%lx) then
       associate(u => this%u, v => this%v, w => this%w, p => this%p, &
            c_Xh => this%c_Xh, ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag)
         do concurrent (j = 1:n)
            u%x(j,1,1,1) = u%x(j,1,1,1) * c_Xh%mult(j,1,1,1)
            v%x(j,1,1,1) = v%x(j,1,1,1) * c_Xh%mult(j,1,1,1)
            w%x(j,1,1,1) = w%x(j,1,1,1) * c_Xh%mult(j,1,1,1)
            p%x(j,1,1,1) = p%x(j,1,1,1) * c_Xh%mult(j,1,1,1)
         end do
         do i = 1, this%ulag%size()
            do concurrent (j = 1:n)
               ulag%lf(i)%x(j,1,1,1) = ulag%lf(i)%x(j,1,1,1) &
                    * c_Xh%mult(j,1,1,1)
               vlag%lf(i)%x(j,1,1,1) = vlag%lf(i)%x(j,1,1,1) &
                    * c_Xh%mult(j,1,1,1)
               wlag%lf(i)%x(j,1,1,1) = wlag%lf(i)%x(j,1,1,1) &
                    * c_Xh%mult(j,1,1,1)
            end do
         end do
       end associate
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

    if (allocated(chkp%previous_mesh%elements) &
         .or. chkp%previous_Xh%lx .ne. this%Xh%lx) then

       call rotate_cyc(this%u%x, this%v%x, this%w%x, 1, this%c_Xh)
       call this%gs_Xh%op(this%u, GS_OP_ADD)
       call this%gs_Xh%op(this%v, GS_OP_ADD)
       call this%gs_Xh%op(this%w, GS_OP_ADD)
       call this%gs_Xh%op(this%p, GS_OP_ADD)
       call rotate_cyc(this%u%x, this%v%x, this%w%x, 0, this%c_Xh)

       do i = 1, this%ulag%size()
          call rotate_cyc(this%ulag%lf(i)%x, this%vlag%lf(i)%x, &
               this%wlag%lf(i)%x, 1, this%c_Xh)
          call this%gs_Xh%op(this%ulag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%vlag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%wlag%lf(i), GS_OP_ADD)
          call rotate_cyc(this%ulag%lf(i)%x, this%vlag%lf(i)%x, &
               this%wlag%lf(i)%x, 0, this%c_Xh)
       end do
    end if

  end subroutine fluid_pnpn_restart

  subroutine fluid_pnpn_free(this)
    class(fluid_pnpn_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    if (allocated(this%ext_bdf)) then
       call this%ext_bdf%free()
       deallocate(this%ext_bdf)
    end if

    call this%bc_prs_surface%free()
    call this%bc_sym_surface%free()
    if (allocated(this%bcs_vel_resolver)) then
       call this%bcs_vel_resolver%free()
       deallocate(this%bcs_vel_resolver)
    end if
    call this%bcs_prs_resolver%free()
    call this%proj_prs%free()
    call this%proj_vel%free()

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

    if (allocated(this%adv)) then
       call this%adv%free()
       deallocate(this%adv)
    end if

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

    if (allocated(this%ext_bdf)) then
       deallocate(this%ext_bdf)
    end if

    call this%vol_flow%free()

  end subroutine fluid_pnpn_free

  !> Advance fluid simulation in time.
  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  !! @param dt_controller timestep controller
  subroutine fluid_pnpn_step(this, time, dt_controller)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    ! number of degrees of freedom
    integer :: n
    ! Solver results monitors (pressure + 3 velocity)
    type(ksp_monitor_t) :: ksp_results(4)

    type(file_t) :: dump_file
    class(bc_t), pointer :: bc_i

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

      ! Apply Dirichlet boundary conditions to the velocity.
      call this%bc_apply_vel(time, strong = .true.)
      ! Apply Dirichlet boundary conditions to the pressure.
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
      call this%bcs_prs_resolver%apply(p_res%x, p%dof%size())


      call profiler_end_region('Pressure_residual', 18)

      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
           Ax=Ax_prs, gs_h=gs_Xh, bclst=this%bcs_prs_resolver, &
           string='Pressure')

      call this%pc_prs%update()

      call profiler_start_region('Pressure_solve', 3)

      ! Solve for the pressure increment.
      ksp_results(1) = &
           this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, &
           this%bcs_prs_resolver, gs_Xh)
      ksp_results(1)%name = 'Pressure'


      call profiler_end_region('Pressure_solve', 3)

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
           this%bcs_prs_resolver, gs_Xh, n, tstep, dt_controller)

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
      call this%bcs_vel_resolver%apply(u_res%x, v_res%x, w_res%x, dm_Xh%size())


      call profiler_end_region('Velocity_residual', 19)

      call this%proj_vel%pre_solving(u_res%x, v_res%x, w_res%x, &
           tstep, c_Xh, n, dt_controller, 'Velocity')

      call this%pc_vel%update()

      call profiler_start_region("Velocity_solve", 4)
      ksp_results(2:4) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bcs_vel_resolver, gs_Xh, &
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
           this%bcs_vel_resolver, gs_Xh, n, tstep, &
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
              dt, time, this%bcs_prs_resolver, this%bcs_vel_resolver, &
              Ax_vel, Ax_prs, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(time, ksp_results, &
           this%full_stress_formulation, this%strict_convergence, &
           this%allow_stabilization)

    end associate
    call profiler_end_region('Fluid', 1)
  end subroutine fluid_pnpn_step

  !> Sets up the boundary condition for the scheme.
  !! @param user The user interface.
  subroutine fluid_pnpn_setup_bcs(this, user, params)
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: i, n_bcs, zone_index, j, zone_size, global_zone_size, ierr
    class(bc_t), pointer :: bc_i
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found
    ! Monitor which boundary zones have been marked
    logical, allocatable :: marked_zones(:)
    integer, allocatable :: zone_indices(:)

    ! Special PnPn boundary conditions for pressure
    call this%bc_prs_surface%init_from_components(this%c_Xh)
    call this%bc_sym_surface%init_from_components(this%c_Xh)

    ! Populate bcs_vel and bcs_prs based on the case file
    if (params%valid_path('case.fluid.boundary_conditions')) then
       call params%info('case.fluid.boundary_conditions', n_children = n_bcs)
       call params%get_core(core)
       call params%get('case.fluid.boundary_conditions', bc_object, found)

       !
       ! Velocity bcs
       !
       call this%bcs_vel%init(n_bcs)

       allocate(marked_zones(size(this%msh%labeled_zones)))
       marked_zones = .false.

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)

          call json_get_or_lookup(bc_subdict, "zone_indices", zone_indices)

          ! Check that we are not trying to assing a bc to zone, for which one
          ! has already been assigned and that the zone has more than 0 size
          ! in the mesh.
          do j = 1, size(zone_indices)
             zone_size = this%msh%labeled_zones(zone_indices(j))%size
             call MPI_Allreduce(zone_size, global_zone_size, 1, &
                  MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

             if (global_zone_size .eq. 0) then
                write(error_unit, '(A, A, I0, A, A, I0, A)') "*** ERROR ***: ",&
                     "Zone index ", zone_indices(j), &
                     " is invalid as this zone has 0 size, meaning it ", &
                     "is not in the mesh. Check fluid boundary condition ", &
                     i, "."
                error stop
             end if

             if (marked_zones(zone_indices(j))) then
                write(error_unit, '(A, A, I0, A, A, A, A)') "*** ERROR ***: ", &
                     "Zone with index ", zone_indices(j), &
                     " has already been assigned a boundary condition. ", &
                     "Please check your boundary_conditions entry for the ", &
                     "fluid and make sure that each zone index appears only ", &
                     "in a single boundary condition."
                error stop
             else
                marked_zones(zone_indices(j)) = .true.
             end if
          end do

          bc_i => null()
          call velocity_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for velocity in particular,
          ! so we check.
          if (associated(bc_i)) then

             select type (bc_i)
             type is (symmetry_aligned_t)
                ! In this case we need to tell the segregated resolver where
                ! we have the dirichlet dofs component-wise. This is stored
                ! in the nested bcs. Of course, we rely on axis-alignment of
                ! the geometry.
                call this%bcs_vel_resolver%mark(bc_i%bc_x, component='x')
                call this%bcs_vel_resolver%mark(bc_i%bc_y, component='y')
                call this%bcs_vel_resolver%mark(bc_i%bc_z, component='z')
                call this%bcs_vel%append(bc_i)
                call this%bc_sym_surface%mark_facets(bc_i%marked_facet)
             type is (symmetry_t)
                ! In this case we add the bc itself to the resolver, which
                ! should be coupled.
                if (.not. this%full_stress_formulation) then
                   call neko_error("The symmetry boundary condition " // &
                        "requires the full stress formulation to be enabled.")
                end if
                call this%bcs_vel_resolver%mark(bc_i)
                call this%bcs_vel%append(bc_i)
                call this%bc_sym_surface%mark_facets(bc_i%marked_facet)
             type is (non_normal_aligned_t)
                ! The situation is the same as symmetry.
                call this%bcs_vel_resolver%mark(bc_i%bc_x, component='x')
                call this%bcs_vel_resolver%mark(bc_i%bc_y, component='y')
                call this%bcs_vel_resolver%mark(bc_i%bc_z, component='z')
                call this%bcs_vel%append(bc_i)
             type is (non_normal_t)
                call this%bcs_vel_resolver%mark(bc_i)
                call this%bcs_vel%append(bc_i)
             type is (shear_stress_t)
                if (.not. this%full_stress_formulation) then
                   call neko_error("The shear_stress boundary condition " // &
                        "requires the full stress formulation to be enabled.")
                end if
                call this%bcs_vel_resolver%mark(bc_i)
                call this%bcs_vel%append(bc_i)
             type is (wall_model_bc_t)
                if (.not. this%full_stress_formulation) then
                   call neko_error("The wall_model boundary condition " // &
                        "requires the full stress formulation to be enabled.")
                end if
                call this%bcs_vel_resolver%mark(bc_i)
                call this%bcs_vel%append(bc_i)
             class default

                ! Additionally we mark the special PnPn pressure bc.
                if (bc_i%bc_type .eq. BC_TYPES%DIRICHLET) then
                   call this%bc_prs_surface%mark_labeled_zones( &
                        bc_i%zone_indices)
                   call this%bcs_vel_resolver%mark(bc_i, component='x')
                   call this%bcs_vel_resolver%mark(bc_i, component='y')
                   call this%bcs_vel_resolver%mark(bc_i, component='z')
                end if

                call this%bcs_vel%append(bc_i)
             end select
          end if
       end do

       ! Make sure all labeled zones with non-zero size have been marked
       do i = 1, size(this%msh%labeled_zones)
          if ((this%msh%labeled_zones(i)%size .gt. 0) .and. &
               (.not. marked_zones(i))) then
             write(error_unit, '(A, A, I0)') "*** ERROR ***: ", &
                  "No fluid boundary condition assigned to zone ", i
             error stop
          end if
       end do

       !
       ! Pressure bcs
       !
       call this%bcs_prs%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call pressure_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_prs%append(bc_i)

             ! Mark strong pressure bcs in the resolver to force zero change.
             if (bc_i%bc_type .eq. BC_TYPES%DIRICHLET) then
                call this%bcs_prs_resolver%mark(bc_i)
             end if

          end if

       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             call neko_error("No boundary_conditions entry in the case file!")
          end if
       end do

       ! For a pure periodic case, we still need to initilise the bc lists
       ! to a zero size to avoid issues with apply() in step()
       call this%bcs_vel%init()
       call this%bcs_prs%init()

    end if

    call this%bc_prs_surface%finalize()
    call this%bc_sym_surface%finalize()
    call this%bcs_vel_resolver%finalize()

!    do i = 1, this%bcs_vel%size()
!       bc_i => this%bcs_vel%get(i)
!       select type (bci => bc_i)
!         class is (mixed_bc_t)
!            call bci%debug_output()
!         class default
!       end select
!    end do

    ! If we have no strong pressure bcs, we will demean the pressure
    this%prs_dirichlet = this%bcs_prs_resolver%dof_mask%is_set()
    call MPI_Allreduce(MPI_IN_PLACE, this%prs_dirichlet, 1, &
         MPI_LOGICAL, MPI_LOR, NEKO_COMM)


    if (allocated(marked_zones)) then
       deallocate(marked_zones)
    end if

    if (allocated(zone_indices)) then
       deallocate(zone_indices)
    end if

  end subroutine fluid_pnpn_setup_bcs

  !> Write a field with boundary condition specifications
  subroutine fluid_pnpn_write_boundary_conditions(this)
    use inflow, only : inflow_t
    use field_dirichlet, only : field_dirichlet_t
    use blasius, only : blasius_t
    use field_dirichlet_vector, only : field_dirichlet_vector_t
    use dong_outflow, only : dong_outflow_t
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(dirichlet_t) :: bdry_mask
    type(field_t), pointer :: bdry_field
    type(file_t) :: bdry_file
    integer :: temp_index, i
    class(bc_t), pointer :: bci
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A)') 'Marking using integer keys in bdry0.f00000'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') 'Condition-value pairs: '
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  no_slip                         = 1'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  velocity_value                  = 2'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  outflow, normal_outflow (+dong) = 3'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  symmetry                        = 4'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  periodic                        = 6'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  user_velocity                   = 7'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  user_pressure                   = 8'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  shear_stress                    = 9'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  wall_modelling                  = 10'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  blasius_profile                 = 11'
    call neko_log%message(log_buf)

    call neko_scratch_registry%request_field(bdry_field, temp_index, .true.)



    call bdry_mask%init_from_components(this%c_Xh, 6.0_rp)
    call bdry_mask%mark_zone(this%msh%periodic)
    call bdry_mask%finalize()
    call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
    call bdry_mask%free()

    do i = 1, this%bcs_prs%size()
       bci => this%bcs_prs%get(i)
       select type (bc => bci)
       type is (zero_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 3.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (dong_outflow_t)
          call bdry_mask%init_from_components(this%c_Xh, 3.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (field_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 8.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       end select
    end do

    do i = 1, this%bcs_vel%size()
       bci => this%bcs_vel%get(i)
       select type (bc => bci)
       type is (zero_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 1.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (inflow_t)
          call bdry_mask%init_from_components(this%c_Xh, 2.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (symmetry_aligned_t)
          call bdry_mask%init_from_components(this%c_Xh, 4.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (symmetry_t)
          call bdry_mask%init_from_components(this%c_Xh, 4.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (field_dirichlet_vector_t)
          call bdry_mask%init_from_components(this%c_Xh, 7.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (shear_stress_t)
          call bdry_mask%init_from_components(this%c_Xh, 9.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (wall_model_bc_t)
          call bdry_mask%init_from_components(this%c_Xh, 10.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (blasius_t)
          call bdry_mask%init_from_components(this%c_Xh, 11.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       end select
    end do


    call bdry_file%init('bdry.fld')
    call bdry_file%write(bdry_field)

    call neko_scratch_registry%relinquish_field(temp_index)
  end subroutine fluid_pnpn_write_boundary_conditions

end module fluid_pnpn
