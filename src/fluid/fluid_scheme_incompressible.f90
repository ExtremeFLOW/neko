! Copyright (c) 2020-2025, The Neko Authors
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
!> Fluid formulations
module fluid_scheme_incompressible
  use fluid_scheme_base, only : fluid_scheme_base_t
  use gather_scatter, only : gs_t, GS_OP_MIN, GS_OP_MAX
  use neko_config, only : NEKO_BCKND_DEVICE
  use checkpoint, only : chkp_t
  use num_types, only : rp, i8
  use fluid_source_term, only: fluid_source_term_t
  use field, only : field_t
  use space, only : GLL
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t, krylov_solver_factory, KSP_MAX_ITER
  use coefs, only: coef_t
  use dirichlet, only : dirichlet_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use phmg, only : phmg_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use fluid_stats, only : fluid_stats_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use mesh, only : mesh_t
  use math, only : glsum
  use operators, only : cfl
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use scratch_registry, only : scratch_registry_t
  use user_intf, only : user_t, dummy_user_material_properties, &
       user_material_properties_intf
  use utils, only : neko_error
  use field_math, only : field_cfill, field_addcol3, field_copy
  use device, only : device_event_sync, glb_cmd_event, DEVICE_TO_HOST, &
       device_memcpy
  use time_state, only : time_state_t
  implicit none
  private

  !> Base type of all fluid formulations
  type, abstract, extends(fluid_scheme_base_t) :: fluid_scheme_incompressible_t
     !> The source term for the momentum equation.
     type(fluid_source_term_t) :: source_term
     class(ksp_t), allocatable :: ksp_vel !< Krylov solver for velocity
     class(ksp_t), allocatable :: ksp_prs !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs !< Velocity Preconditioner
     integer :: vel_projection_dim !< Size of the projection space for ksp_vel
     integer :: pr_projection_dim !< Size of the projection space for ksp_pr
     integer :: vel_projection_activ_step !< Steps to activate projection for ksp_vel
     integer :: pr_projection_activ_step !< Steps to activate projection for ksp_pr
     logical :: strict_convergence !< Strict convergence for the velocity solver
     !> Extrapolation velocity fields for LES
     type(field_t), pointer :: u_e => null() !< Extrapolated x-Velocity
     type(field_t), pointer :: v_e => null() !< Extrapolated y-Velocity
     type(field_t), pointer :: w_e => null() !< Extrapolated z-Velocity

     type(fluid_stats_t) :: stats !< Fluid statistics
     logical :: forced_flow_rate = .false. !< Is the flow rate forced?

     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name

     ! The total viscosity field
     type(field_t), pointer :: mu_tot => null()

     !> Global number of GLL points for the fluid (not unique)
     integer(kind=i8) :: glb_n_points
     !> Global number of GLL points for the fluid (unique)
     integer(kind=i8) :: glb_unique_points
     type(scratch_registry_t) :: scratch !< Manager for temporary fields
   contains
     !> Constructor for the base type
     procedure, pass(this) :: init_base => fluid_scheme_init_base
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => fluid_scheme_validate
     !> Apply pressure boundary conditions
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     !> Apply velocity boundary conditions
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => fluid_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: set_material_properties => &
          fluid_scheme_set_material_properties

     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          fluid_scheme_update_material_properties
     !> Linear solver factory, wraps a KSP constructor
     procedure, nopass :: solver_factory => fluid_scheme_solver_factory
     !> Preconditioner factory
     procedure, pass(this) :: precon_factory_ => fluid_scheme_precon_factory
  end type fluid_scheme_incompressible_t

  interface
     !> Initialise a fluid scheme
     module subroutine fluid_scheme_factory(object, type_name)
       class(fluid_scheme_base_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine fluid_scheme_factory
  end interface

  public :: fluid_scheme_incompressible_t, fluid_scheme_factory

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_base(this, msh, lx, params, scheme, user, &
       kspv_init)
    implicit none
    class(fluid_scheme_incompressible_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical, intent(in) :: kspv_init
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), allocatable :: real_vec(:)
    real(kind=rp) :: real_val, kappa, B, z0
    logical :: logical_val
    integer :: integer_val, ierr
    type(json_file) :: wm_json
    character(len=:), allocatable :: string_val1, string_val2
    real(kind=rp) :: GJP_param_a, GJP_param_b
    type(json_file) :: json_subdict

    !
    ! SEM simulation fundamentals
    !

    this%msh => msh

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Local scratch registry
    call this%scratch%init(this%dm_Xh, 10, 2)

    ! Assign a name
    call json_get_or_default(params, 'case.fluid.name', this%name, "fluid")

    !
    ! First section of fluid log
    !

    call neko_log%section('Fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, 'u')
    call neko_field_registry%add_field(this%dm_Xh, 'v')
    call neko_field_registry%add_field(this%dm_Xh, 'w')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    !
    ! Material properties
    !
    call this%set_material_properties(params, user)

    ! Projection spaces
    call json_get_or_default(params, &
         'case.fluid.velocity_solver.projection_space_size', &
         this%vel_projection_dim, 0)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_space_size', &
         this%pr_projection_dim, 0)
    call json_get_or_default(params, &
         'case.fluid.velocity_solver.projection_hold_steps', &
         this%vel_projection_activ_step, 5)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_hold_steps', &
         this%pr_projection_activ_step, 5)


    call json_get_or_default(params, 'case.fluid.freeze', this%freeze, .false.)

    if (params%valid_path("case.fluid.flow_rate_force")) then
       this%forced_flow_rate = .true.
    end if


    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'Poly order : ', lx-1
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'Poly order : ', lx-1
    else
       write(log_buf, '(A, I3)') 'Poly order : ', lx-1
    end if
    call neko_log%message(log_buf)
    this%glb_n_points = int(this%msh%glb_nelv, i8)*int(this%Xh%lxyz, i8)
    this%glb_unique_points = int(glsum(this%c_Xh%mult, this%dm_Xh%size()), i8)

    write(log_buf, '(A, I0)') 'GLL points : ', this%glb_n_points
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') 'Unique pts.: ', this%glb_unique_points
    call neko_log%message(log_buf)


    call json_get(params, 'case.numerics.dealias', logical_val)
    write(log_buf, '(A, L1)') 'Dealias    : ', logical_val
    call neko_log%message(log_buf)


    call json_get_or_default(params, 'case.output_boundary', logical_val, &
         .false.)
    write(log_buf, '(A, L1)') 'Save bdry  : ', logical_val
    call neko_log%message(log_buf)

    call json_get_or_default(params, "case.fluid.full_stress_formulation", &
         logical_val, .false.)
    write(log_buf, '(A, L1)') 'Full stress: ', logical_val
    call neko_log%message(log_buf)


    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_x)
    allocate(this%f_y)
    allocate(this%f_z)
    call this%f_x%init(this%dm_Xh, fld_name = "fluid_rhs_x")
    call this%f_y%init(this%dm_Xh, fld_name = "fluid_rhs_y")
    call this%f_z%init(this%dm_Xh, fld_name = "fluid_rhs_z")

    ! Initialize velocity solver
    if (kspv_init) then
       call neko_log%section("Velocity solver")
       call json_get_or_default(params, &
            'case.fluid.velocity_solver.max_iterations', &
            integer_val, KSP_MAX_ITER)
       call json_get(params, 'case.fluid.velocity_solver.type', string_val1)
       call json_get(params, 'case.fluid.velocity_solver.preconditioner.type', &
            string_val2)
       call json_get(params, &
            'case.fluid.velocity_solver.preconditioner', json_subdict)
       call json_get(params, 'case.fluid.velocity_solver.absolute_tolerance', &
            real_val)
       call json_get_or_default(params, &
            'case.fluid.velocity_solver.monitor', &
            logical_val, .false.)

       call neko_log%message('Type       : ('// trim(string_val1) // &
            ', ' // trim(string_val2) // ')')

       write(log_buf, '(A,ES13.6)') 'Abs tol    :', real_val
       call neko_log%message(log_buf)
       call this%solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            string_val1, integer_val, real_val, logical_val)
       call this%precon_factory_(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_vel, &
            string_val2, json_subdict)
       call neko_log%end_section()
    end if

    ! Strict convergence for the velocity solver
    call json_get_or_default(params, 'case.fluid.strict_convergence', &
         this%strict_convergence, .false.)


    !! Initialize time-lag fields
    call this%ulag%init(this%u, 2)
    call this%vlag%init(this%v, 2)
    call this%wlag%init(this%w, 2)

    call neko_field_registry%add_field(this%dm_Xh, 'u_e')
    call neko_field_registry%add_field(this%dm_Xh, 'v_e')
    call neko_field_registry%add_field(this%dm_Xh, 'w_e')
    this%u_e => neko_field_registry%get_field('u_e')
    this%v_e => neko_field_registry%get_field('v_e')
    this%w_e => neko_field_registry%get_field('w_e')

    ! Initialize the source term
    call neko_log%section('Fluid Source term')
    call this%source_term%init(this%f_x, this%f_y, this%f_z, this%c_Xh, user, &
         this%name)
    call this%source_term%add(params, 'case.fluid.source_terms')
    call neko_log%end_section()

  end subroutine fluid_scheme_init_base

  subroutine fluid_scheme_free(this)
    class(fluid_scheme_incompressible_t), intent(inout) :: this

    call this%Xh%free()

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if

    if (allocated(this%pc_vel)) then
       call precon_destroy(this%pc_vel)
       deallocate(this%pc_vel)
    end if

    if (allocated(this%pc_prs)) then
       call precon_destroy(this%pc_prs)
       deallocate(this%pc_prs)
    end if

    call this%source_term%free()

    call this%gs_Xh%free()

    call this%c_Xh%free()

    call this%scratch%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    nullify(this%u_e)
    nullify(this%v_e)
    nullify(this%w_e)

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()


    if (associated(this%f_x)) then
       call this%f_x%free()
       deallocate(this%f_x)
    end if

    if (associated(this%f_y)) then
       call this%f_y%free()
       deallocate(this%f_y)
    end if

    if (associated(this%f_z)) then
       call this%f_z%free()
       deallocate(this%f_z)
    end if

    nullify(this%f_x)
    nullify(this%f_y)
    nullify(this%f_z)
    nullify(this%rho)
    nullify(this%mu)
    nullify(this%mu_tot)

    call this%dm_Xh%free()
    call this%Xh%free()
    nullify(this%msh)

  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_incompressible_t), target, intent(inout) :: this
    ! Variables for retrieving json parameters
    logical :: logical_val

    if ( (.not. associated(this%u)) .or. &
         (.not. associated(this%v)) .or. &
         (.not. associated(this%w)) .or. &
         (.not. associated(this%p))) then
       call neko_error('Fields are not registered')
    end if

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%p%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp_vel)) then
       call neko_error('No Krylov solver for velocity defined')
    end if

    if (.not. allocated(this%ksp_prs)) then
       call neko_error('No Krylov solver for pressure defined')
    end if

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! Here we perform additional gs operations to take care of
  !! shared points between elements that have different BCs, as done in Nek5000.
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this, time, strong)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical, intent(in) :: strong
    integer :: i
    class(bc_t), pointer :: b
    b => null()

    call this%bcs_vel%apply_vector(&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size(), time, strong)
    call this%gs_Xh%op(this%u, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)


    call this%bcs_vel%apply_vector(&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size(), time, strong)
    call this%gs_Xh%op(this%u, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)

    do i = 1, this%bcs_vel%size()
       b => this%bcs_vel%get(i)
       b%updated = .false.
    end do
    nullify(b)

  end subroutine fluid_scheme_bc_apply_vel

  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this, time)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i
    class(bc_t), pointer :: b
    b => null()

    call this%bcs_prs%apply(this%p, time)
    call this%gs_Xh%op(this%p, GS_OP_MIN, glb_cmd_event)
    call device_event_sync(glb_cmd_event)

    call this%bcs_prs%apply(this%p, time)
    call this%gs_Xh%op(this%p, GS_OP_MAX, glb_cmd_event)
    call device_event_sync(glb_cmd_event)

    do i = 1, this%bcs_prs%size()
       b => this%bcs_prs%get(i)
       b%updated = .false.
    end do
    nullify(b)

  end subroutine fluid_scheme_bc_apply_prs

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver, &
       max_iter, abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*), intent(in) :: solver
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, abstol, &
         monitor = monitor)

  end subroutine fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine fluid_scheme_precon_factory(this, pc, ksp, coef, dof, gs, bclst, &
       pctype, pcparams)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype
    type(json_file), intent(inout) :: pcparams

    call precon_factory(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (hsmg_t)
       call pcp%init(coef, bclst, pcparams)
    type is (phmg_t)
       call pcp%init(coef, bclst, pcparams)
    end select

    call ksp%set_pc(pc)

  end subroutine fluid_scheme_precon_factory

  !> Compute CFL
  function fluid_compute_cfl(this, dt) result(c)
    class(fluid_scheme_incompressible_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u%x, this%v%x, this%w%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)

  end function fluid_compute_cfl


  !> Call user material properties routine and update the values of `mu`
  !! if necessary.
  !! @param t Time value.
  !! @param tstep Current time step.
  subroutine fluid_scheme_update_material_properties(this, time)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: nut

    call this%user_material_properties(this%name, this%material_properties, &
         time)

    if (len(trim(this%nut_field_name)) > 0) then
       nut => neko_field_registry%get_field(this%nut_field_name)
       ! Copy material property
       call field_copy(this%mu_tot, this%mu)
       ! Add turbulent contribution
       call field_addcol3(this%mu_tot, nut, this%rho)
    end if

    ! Since mu, rho is a field_t, and we use the %x(1,1,1,1)
    ! host array data to pass constant density and viscosity
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, this%rho%size(), &
            DEVICE_TO_HOST, sync=.false.)
    end if
  end subroutine fluid_scheme_update_material_properties

  !> Sets rho and mu
  !! @param params The case paramter file.
  !! @param user The user interface.
  subroutine fluid_scheme_set_material_properties(this, params, user)
    class(fluid_scheme_incompressible_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties_intf), pointer :: dummy_mp_ptr
    logical :: nondimensional
    real(kind=rp) :: dummy_lambda, dummy_cp
    real(kind=rp) :: const_mu, const_rho
    type(time_state_t) :: dummy_time_state


    dummy_mp_ptr => dummy_user_material_properties

    call neko_field_registry%add_field(this%dm_Xh, this%name // "_mu")
    call neko_field_registry%add_field(this%dm_Xh, this%name // "_mu_tot")
    call neko_field_registry%add_field(this%dm_Xh, this%name // "_rho")
    this%mu => neko_field_registry%get_field(this%name // "_mu")
    this%mu_tot => neko_field_registry%get_field(this%name // "_mu_tot")
    this%rho => neko_field_registry%get_field(this%name // "_rho")

    call this%material_properties%init(2)
    call this%material_properties%assign(1, this%rho)
    call this%material_properties%assign(2, this%mu)

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') 'Material properties must be set in the user' // &
            ' file!'
       call neko_log%message(log_buf)
       this%user_material_properties => user%material_properties

       call user%material_properties(this%name, this%material_properties, &
            dummy_time_state)

    else
       this%user_material_properties => dummy_user_material_properties
       ! Incorrect user input
       if (params%valid_path('case.fluid.Re') .and. &
            (params%valid_path('case.fluid.mu') .or. &
            params%valid_path('case.fluid.rho'))) then
          call neko_error("To set the material properties for the fluid, " // &
               "either provide Re OR mu and rho in the case file.")

       else if (params%valid_path('case.fluid.Re')) then
          ! Non-dimensional case
          write(log_buf, '(A)') 'Non-dimensional fluid material properties &
          & input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Density will be set to 1, dynamic viscosity to&
          & 1/Re.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Re into mu for further manipulation.
          call json_get(params, 'case.fluid.Re', const_mu)
          write(log_buf, '(A)') 'Read non-dimensional material properties'
          call neko_log%message(log_buf)
          write(log_buf, '(A,ES13.6)') 'Re         :', const_mu
          call neko_log%message(log_buf)

          ! Set rho to 1 since the setup is non-dimensional.
          const_rho = 1.0_rp
          ! Invert the Re to get viscosity.
          const_mu = 1.0_rp/const_mu
       else
          ! Dimensional case
          call json_get(params, 'case.fluid.mu', const_mu)
          call json_get(params, 'case.fluid.rho', const_rho)
       end if
    end if

    ! We need to fill the fields based on the parsed const values
    ! if the user routine is not used.
    if (associated(user%material_properties, dummy_mp_ptr)) then
       ! Fill mu and rho field with the physical value
       call field_cfill(this%mu, const_mu)
       call field_cfill(this%mu_tot, const_mu)
       call field_cfill(this%rho, const_rho)


       write(log_buf, '(A,ES13.6)') 'rho        :', const_rho
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'mu         :', const_mu
       call neko_log%message(log_buf)
    end if

    ! Copy over material property to the total one
    call field_copy(this%mu_tot, this%mu)

    ! Since mu, rho is a field_t, and we use the %x(1,1,1,1)
    ! host array data to pass constant density and viscosity
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, this%rho%size(), &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%mu%x, this%mu%x_d, this%mu%size(), &
            DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(this%mu_tot%x, this%mu_tot%x_d, this%mu%size(), &
            DEVICE_TO_HOST, sync = .false.)
    end if
  end subroutine fluid_scheme_set_material_properties

end module fluid_scheme_incompressible
