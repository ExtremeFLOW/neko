! Copyright (c) 2022-2026, The Neko Authors
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
!> Contains the scalar_scheme_t type.

module scalar_scheme
  use gather_scatter, only : gs_t
  use checkpoint, only : chkp_t
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use space, only : space_t
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t, krylov_solver_factory, KSP_MAX_ITER, ksp_monitor_t
  use coefs, only : coef_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use hsmg, only : hsmg_t
  use bc_list, only : bc_list_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use mesh, only : mesh_t
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use registry, only : neko_registry
  use json_utils, only : json_get, json_get_or_default, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use json_module, only : json_file
  use user_intf, only : user_t, dummy_user_material_properties, &
       user_material_properties_intf
  use utils, only : neko_error
  use scalar_source_term, only : scalar_source_term_t
  use field_series, only : field_series_t
  use field_math, only : field_cmult2, field_col3, field_cfill, field_add3, &
       field_copy, field_col2
  use neko_config, only : NEKO_BCKND_DEVICE
  use time_step_controller, only : time_step_controller_t
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use device, only : device_memcpy, DEVICE_TO_HOST
  implicit none

  !> Base type for a scalar advection-diffusion solver.
  type, abstract :: scalar_scheme_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name
     !> x-component of Velocity
     type(field_t), pointer :: u
     !> y-component of Velocity
     type(field_t), pointer :: v
     !> z-component of Velocity
     type(field_t), pointer :: w
     !> The scalar.
     type(field_t), pointer :: s
     !> Lag arrays, i.e. solutions at previous timesteps.
     type(field_series_t) :: slag
     !> Function space \f$ X_h \f$.
     type(space_t), pointer :: Xh
     !> Dofmap associated with \f$ X_h \f$.
     type(dofmap_t), pointer :: dm_Xh
     !> Gather-scatter associated with \f$ X_h \f$.
     type(gs_t), pointer :: gs_Xh
     !> Coefficients associated with \f$ X_h \f$.
     type(coef_t), pointer :: c_Xh
     !> Right-hand side.
     type(field_t), pointer :: f_Xh => null()
     !> The source term for equation.
     type(scalar_source_term_t) :: source_term
     !> Krylov solver.
     class(ksp_t), allocatable :: ksp
     !> Max iterations in the Krylov solver.
     integer :: ksp_maxiter
     !> Projection space size.
     integer :: projection_dim
     !< Steps to activate projection for ksp
     integer :: projection_activ_step
     !> Preconditioner.
     class(pc_t), allocatable :: pc
     !> List of boundary conditions, including the user one.
     type(bc_list_t) :: bcs
     !> Case paramters.
     type(json_file), pointer :: params
     !> Mesh.
     type(mesh_t), pointer :: msh => null()
     !> Checkpoint for restarts.
     type(chkp_t), pointer :: chkp => null()
     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name
     !> The turbulent diffusivity field name
     character(len=:), allocatable :: alphat_field_name
     !> Density.
     type(field_t), pointer :: rho => null()
     !> Thermal diffusivity.
     type(field_t), pointer :: lambda => null()
     !> Specific heat capacity.
     type(field_t), pointer :: cp => null()
     !> Total diffusivity.
     type(field_t), pointer :: lambda_tot => null()
     !> Turbulent Prandtl number.
     real(kind=rp) :: urb
     !> Field list with cp and lambda
     type(field_list_t) :: material_properties
     procedure(user_material_properties_intf), nopass, pointer :: &
          user_material_properties => null()
     !> Freeze the scheme, i.e. do nothing in step()
     logical :: freeze = .false.
   contains
     !> Constructor for the base type.
     procedure, pass(this) :: scheme_init => scalar_scheme_init
     !> Destructor for the base type.
     procedure, pass(this) :: scheme_free => scalar_scheme_free
     !> Validate successful initialization.
     procedure, pass(this) :: validate => scalar_scheme_validate
     !> Set lambda and cp
     procedure, pass(this) :: set_material_properties => &
          scalar_scheme_set_material_properties
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          scalar_scheme_update_material_properties
     !> Constructor.
     procedure(scalar_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor.
     procedure(scalar_scheme_free_intrf), pass(this), deferred :: free
     !> Solve for the current timestep.
     procedure(scalar_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint.
     procedure(scalar_scheme_restart_intrf), pass(this), deferred :: restart
  end type scalar_scheme_t

  !> Abstract interface to initialize a scalar formulation
  abstract interface
     subroutine scalar_scheme_init_intrf(this, msh, coef, gs, params, &
          numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
       import scalar_scheme_t
       import json_file
       import coef_t
       import gs_t
       import mesh_t
       import user_t
       import field_series_t, field_t
       import time_scheme_controller_t
       import rp
       import chkp_t
       class(scalar_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(in) :: msh
       type(coef_t), target, intent(in) :: coef
       type(gs_t), target, intent(inout) :: gs
       type(json_file), target, intent(inout) :: params
       type(json_file), target, intent(inout) :: numerics_params
       type(user_t), target, intent(in) :: user
       type(chkp_t), target, intent(inout) :: chkp
       type(field_series_t), target, intent(in) :: ulag, vlag, wlag
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
       type(field_t), target, intent(in) :: rho
     end subroutine scalar_scheme_init_intrf
  end interface

  !> Abstract interface to restart a scalar formulation
  abstract interface
     subroutine scalar_scheme_restart_intrf(this, chkp)
       import scalar_scheme_t
       import chkp_t
       import rp
       class(scalar_scheme_t), target, intent(inout) :: this
       type(chkp_t), intent(inout) :: chkp
     end subroutine scalar_scheme_restart_intrf
  end interface

  !> Abstract interface to dealocate a scalar formulation
  abstract interface
     subroutine scalar_scheme_free_intrf(this)
       import scalar_scheme_t
       class(scalar_scheme_t), intent(inout) :: this
     end subroutine scalar_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine scalar_scheme_step_intrf(this, time, ext_bdf, dt_controller, &
          ksp_results)
       import scalar_scheme_t
       import time_state_t
       import time_scheme_controller_t
       import time_step_controller_t
       import ksp_monitor_t
       class(scalar_scheme_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
       type(time_scheme_controller_t), intent(in) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
       type(ksp_monitor_t), intent(inout) :: ksp_results
     end subroutine scalar_scheme_step_intrf
  end interface

  ! ========================================================================== !
  ! Helper functions and types for scalar_scheme_t
  ! ========================================================================== !

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: scalar_scheme_wrapper_t
     class(scalar_scheme_t), allocatable :: scalar
   contains
     !> Constructor. Just allocates the object.
     procedure, pass(this) :: init => scalar_scheme_wrapper_init
     !> Destructor. Just deallocates the object.
     procedure, pass(this) :: free => scalar_scheme_wrapper_free
     !> Move operator for the wrapper, needed for storing schemes
     !! in lists and arrays.
     procedure, pass(this) :: move_from => &
          scalar_scheme_wrapper_move_from
     !> Return allocation status.
     procedure, pass(this) :: is_allocated => &
          scalar_scheme_wrapper_is_allocated
  end type scalar_scheme_wrapper_t


  interface
     !> Scalar scheme factory.
     !! Both constructs and initializes the object.
     !! @param object The object to be created and initialized.
     !! @param msh The mesh.
     !! @param coef The coefficients.
     !! @param gs The gather-scatter.
     !! @param params The parameter dictionary in json.
     !! @param numerics_params The numerical parameter dictionary in json.
     !! @param user Type with user-defined procedures.
     !! @param chkp Checkpoint for restarts.
     !! @param ulag, vlag, wlag The lagged velocity fields.
     !! @param time_scheme The time scheme controller.
     !! @param rho The density field.
     module subroutine scalar_scheme_factory(object, msh, coef, gs, params, &
          numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
       class(scalar_scheme_t), allocatable, intent(inout) :: object
       type(mesh_t), target, intent(in) :: msh
       type(coef_t), target, intent(in) :: coef
       type(gs_t), target, intent(inout) :: gs
       type(json_file), target, intent(inout) :: params
       type(json_file), target, intent(inout) :: numerics_params
       type(user_t), target, intent(in) :: user
       type(chkp_t), target, intent(inout) :: chkp
       type(field_series_t), target, intent(in) :: ulag, vlag, wlag
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
       type(field_t), target, intent(in) :: rho
     end subroutine scalar_scheme_factory
  end interface

  interface
     !> Scalar scheme allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the scalar scheme type.
     module subroutine scalar_scheme_allocator(object, type_name)
       class(scalar_scheme_t), allocatable, intent(inout) :: object
       character(len=*), intent(in):: type_name
     end subroutine scalar_scheme_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine scalar_scheme_allocate(obj)
       import scalar_scheme_t
       class(scalar_scheme_t), allocatable, intent(inout) :: obj
     end subroutine scalar_scheme_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_scalar_scheme(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(scalar_scheme_allocate), pointer, intent(in) :: allocator
     end subroutine register_scalar_scheme
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type scalar_scheme_allocator_entry
     character(len=20) :: type_name
     procedure(scalar_scheme_allocate), pointer, nopass :: allocator
  end type scalar_scheme_allocator_entry

  !> Registry of allocators for user-defined types
  type(scalar_scheme_allocator_entry), allocatable, private :: &
       scalar_scheme_registry(:)

  !> The size of the `scalar_scheme_registry`
  integer, private :: scalar_scheme_registry_size = 0

contains

  !> Initialize all related components of the current scheme
  !! @param msh The mesh.
  !! @param c_Xh The coefficients.
  !! @param gs_Xh The gather-scatter.
  !! @param params The parameter dictionary in json.
  !! @param scheme The name of the scalar scheme.
  !! @param user Type with user-defined procedures.
  !! @param rho The density of the fluid.
  subroutine scalar_scheme_init(this, msh, c_Xh, gs_Xh, params, scheme, user, &
       rho)
    class(scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: c_Xh
    type(gs_t), target, intent(inout) :: gs_Xh
    type(json_file), target, intent(inout) :: params
    character(len=*), intent(in) :: scheme
    type(user_t), target, intent(in) :: user
    type(field_t), target, intent(in) :: rho
    ! IO buffer for log output
    character(len=LOG_SIZE) :: log_buf
    ! Variables for retrieving json parameters
    logical :: logical_val
    real(kind=rp) :: real_val, solver_abstol
    integer :: integer_val, ierr
    character(len=:), allocatable :: solver_type, solver_precon
    type(json_file) :: precon_params
    type(json_file) :: json_subdict
    logical :: nut_dependency

    this%u => neko_registry%get_field('u')
    this%v => neko_registry%get_field('v')
    this%w => neko_registry%get_field('w')
    this%rho => rho

    ! Assign a name
    ! Note that the keyword is added by `scalars_t`, so there is always a
    ! default.
    call json_get(params, 'name', this%name)

    ! Set the freeze flag
    call json_get_or_default(params, 'freeze', this%freeze, .false.)

    call neko_log%section('Scalar')
    call json_get(params, 'solver.type', solver_type)
    call json_get(params, 'solver.preconditioner.type', &
         solver_precon)
    call json_get(params, 'solver.preconditioner', precon_params)
    call json_get_or_lookup(params, 'solver.absolute_tolerance', &
         solver_abstol)

    call json_get_or_lookup_or_default(params, &
         'solver.projection_space_size', &
         this%projection_dim, 0)
    call json_get_or_lookup_or_default(params, &
         'solver.projection_hold_steps', &
         this%projection_activ_step, 5)


    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)
    call neko_log%message('Ksp scalar : ('// trim(solver_type) // &
         ', ' // trim(solver_precon) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :', solver_abstol
    call neko_log%message(log_buf)

    this%Xh => this%u%Xh
    this%dm_Xh => this%u%dof
    this%params => params
    this%msh => msh

    call neko_registry%add_field(this%dm_Xh, this%name, &
         ignore_existing = .true.)

    this%s => neko_registry%get_field(this%name)

    call this%slag%init(this%s, 2)

    this%gs_Xh => gs_Xh
    this%c_Xh => c_Xh

    !
    ! Material properties
    !
    call this%set_material_properties(params, user)


    !
    ! Turbulence modelling
    !
    this%alphat_field_name = ""
    this%nut_field_name = ""
    if (params%valid_path('alphat')) then
       call json_get(this%params, 'alphat', json_subdict)
       call json_get(json_subdict, 'nut_dependency', nut_dependency)
       if (nut_dependency) then
          call json_get_or_lookup(json_subdict, 'Pr_t', this%pr_turb)
          call json_get(json_subdict, 'nut_field', this%nut_field_name)
       else
          call json_get(json_subdict, 'alphat_field', this%alphat_field_name)
       end if
    end if

    !
    ! Setup right-hand side field.
    !
    allocate(this%f_Xh)
    call this%f_Xh%init(this%dm_Xh, fld_name = "scalar_rhs")

    ! Initialize the source term
    call this%source_term%init(this%f_Xh, this%c_Xh, user, this%name)
    call this%source_term%add(params, 'source_terms')

    ! todo parameter file ksp tol should be added
    call json_get_or_lookup_or_default(params, &
         'solver.max_iterations', &
         integer_val, KSP_MAX_ITER)
    call json_get_or_default(params, &
         'solver.monitor', &
         logical_val, .false.)
    call scalar_scheme_solver_factory(this%ksp, this%dm_Xh%size(), &
         solver_type, integer_val, solver_abstol, logical_val)
    call scalar_scheme_precon_factory(this%pc, this%ksp, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs, &
         solver_precon, precon_params)

    call neko_log%end_section()

  end subroutine scalar_scheme_init


  !> Deallocate a scalar formulation
  subroutine scalar_scheme_free(this)
    class(scalar_scheme_t), intent(inout) :: this

    nullify(this%Xh)
    nullify(this%dm_Xh)
    nullify(this%gs_Xh)
    nullify(this%c_Xh)
    nullify(this%params)

    if (allocated(this%ksp)) then
       call this%ksp%free()
       deallocate(this%ksp)
    end if

    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if

    if (allocated(this%name)) then
       deallocate(this%name)
    end if

    call this%source_term%free()

    call this%bcs%free()
    call this%slag%free()

    nullify(this%cp)
    nullify(this%lambda)
    nullify(this%lambda_tot)

  end subroutine scalar_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine scalar_scheme_validate(this)
    class(scalar_scheme_t), target, intent(inout) :: this

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%s%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp)) then
       call neko_error('No Krylov solver for velocity defined')
    end if

    if (.not. associated(this%Xh)) then
       call neko_error('No function space defined')
    end if

    if (.not. associated(this%dm_Xh)) then
       call neko_error('No dofmap defined')
    end if

    if (.not. associated(this%c_Xh)) then
       call neko_error('No coefficients defined')
    end if

    if (.not. associated(this%f_Xh)) then
       call neko_error('No rhs allocated')
    end if

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

    if (.not. associated(this%rho)) then
       call neko_error('No density field defined')
    end if

  end subroutine scalar_scheme_validate

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine scalar_scheme_solver_factory(ksp, n, solver, max_iter, &
       abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    integer, intent(in) :: max_iter
    character(len=*), intent(in) :: solver
    real(kind=rp) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, &
         abstol, monitor = monitor)

  end subroutine scalar_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine scalar_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, &
       pctype, pcparams)
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
    end select

    call ksp%set_pc(pc)

  end subroutine scalar_scheme_precon_factory

  !> Call user material properties routine and update the values of `lambda`
  !! if necessary.
  !! @param t Time value.
  !! @param tstep Current time step.
  subroutine scalar_scheme_update_material_properties(this, time)
    class(scalar_scheme_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: nut, alphat
    integer :: index
    ! Factor to transform nu_t to lambda_t
    type(field_t), pointer :: lambda_factor

    call this%user_material_properties(this%name, this%material_properties, &
         time)

    ! factor = rho * cp / pr_turb
    if (len_trim(this%nut_field_name) .gt. 0 &
         .and. len_trim(this%alphat_field_name) .eq. 0 ) then
       nut => neko_registry%get_field(this%nut_field_name)

       ! lambda_tot = lambda + rho * cp * nut / pr_turb
       call neko_scratch_registry%request_field(lambda_factor, index, .false.)
       call field_cmult2(lambda_factor, nut, 1.0_rp / this%pr_turb)
       call field_col2(lambda_factor, this%cp)
       call field_col2(lambda_factor, this%rho)
       call field_add3(this%lambda_tot, this%lambda, lambda_factor)
       call neko_scratch_registry%relinquish_field(index)

    else if (len_trim(this%alphat_field_name) .gt. 0 &
         .and. len_trim(this%nut_field_name) .eq. 0 ) then
       alphat => neko_registry%get_field(this%alphat_field_name)

       ! lambda_tot = lambda + rho * cp * alphat
       call neko_scratch_registry%request_field(lambda_factor, index, .false.)
       call field_col3(lambda_factor, this%cp, alphat)
       call field_col2(lambda_factor, this%rho)
       call field_add3(this%lambda_tot, this%lambda, lambda_factor)
       call neko_scratch_registry%relinquish_field(index)

    else if (len_trim(this%alphat_field_name) .gt. 0 &
         .and. len_trim(this%nut_field_name) .gt. 0 ) then
       call neko_error("Conflicting definition of eddy diffusivity " // &
            "for the scalar equation")
    end if

    ! Since cp is a fields and we use the %x(1,1,1,1) of the
    ! host array data to pass constant  material properties
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%cp%x, this%cp%x_d, this%cp%size(), &
            DEVICE_TO_HOST, sync = .false.)
    end if

  end subroutine scalar_scheme_update_material_properties

  !> Set lamdba and cp.
  !! @param params The case file configuration dictionary.
  !! @param user The user interface.
  subroutine scalar_scheme_set_material_properties(this, params, user)
    class(scalar_scheme_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties_intf), pointer :: dummy_mp_ptr
    real(kind=rp) :: const_cp, const_lambda
    ! Dummy time state set to 0
    type(time_state_t) :: time

    dummy_mp_ptr => dummy_user_material_properties

    ! Fill lambda field with the physical value

    call neko_registry%add_field(this%dm_Xh, this%name // "_lambda")
    call neko_registry%add_field(this%dm_Xh, this%name // "_lambda_tot")
    call neko_registry%add_field(this%dm_Xh, this%name // "_cp")
    this%lambda => neko_registry%get_field(this%name // "_lambda")
    this%lambda_tot => neko_registry%get_field(this%name // "_lambda_tot")
    this%cp => neko_registry%get_field(this%name // "_cp")

    call this%material_properties%init(2)
    call this%material_properties%assign(1, this%cp)
    call this%material_properties%assign(2, this%lambda)

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user " // &
            "file!"
       call neko_log%message(log_buf)
       this%user_material_properties => user%material_properties

       call user%material_properties(this%name, this%material_properties, time)
    else
       this%user_material_properties => dummy_user_material_properties
       if (params%valid_path('Pe') .and. &
            (params%valid_path('lambda') .or. &
            params%valid_path('cp'))) then
          call neko_error("To set the material properties for the scalar, " // &
               "either provide Pe OR lambda and cp in the case file.")
          ! Non-dimensional case
       else if (params%valid_path('Pe')) then
          write(log_buf, '(A)') 'Non-dimensional scalar material properties' //&
               ' input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Specific heat capacity will be set to 1, '
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'conductivity to 1/Pe. Assumes density is 1.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Pe into lambda for further manipulation.
          call json_get_or_lookup(params, 'Pe', const_lambda)
          write(log_buf, '(A,ES13.6)') 'Pe         :', const_lambda
          call neko_log%message(log_buf)

          ! Set cp and rho to 1 since the setup is non-dimensional.
          const_cp = 1.0_rp
          ! Invert the Pe to get conductivity
          const_lambda = 1.0_rp/const_lambda
          ! Dimensional case
       else
          call json_get_or_lookup(params, 'lambda', const_lambda)
          call json_get_or_lookup(params, 'cp', const_cp)
       end if
    end if
    ! We need to fill the fields based on the parsed const values
    ! if the user routine is not used.
    if (associated(user%material_properties, dummy_mp_ptr)) then
       ! Fill mu and rho field with the physical value
       call field_cfill(this%lambda, const_lambda)
       call field_cfill(this%cp, const_cp)

       write(log_buf, '(A,ES13.6)') 'lambda     :', const_lambda
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'cp         :', const_cp
       call neko_log%message(log_buf)
    end if

    ! Copy over material property to the total one
    call field_copy(this%lambda_tot, this%lambda)

    ! Since cp is a field and we use the %x(1,1,1,1) of the
    ! host array data to pass constant material properties
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%cp%x, this%cp%x_d, this%cp%size(), &
            DEVICE_TO_HOST, sync = .false.)
    end if
  end subroutine scalar_scheme_set_material_properties

  ! ========================================================================== !
  ! Scalar scheme wrapper type methods

  !> Constructor. Initializes the object.
  subroutine scalar_scheme_wrapper_init(this, msh, coef, gs, params, &
       numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalar_scheme_wrapper_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_t), target, intent(in) :: rho

    call this%free()
    call scalar_scheme_factory(this%scalar, msh, coef, gs, params, &
         numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)

  end subroutine scalar_scheme_wrapper_init

  !> Destructor. Just deallocates the pointer.
  subroutine scalar_scheme_wrapper_free(this)
    class(scalar_scheme_wrapper_t), intent(inout) :: this

    if (allocated(this%scalar)) then
       call this%scalar%free()
       deallocate(this%scalar)
    end if

  end subroutine scalar_scheme_wrapper_free

  !> Move assignment operator for the wrapper, needed for storing schemes
  !! in lists and arrays.
  !! @param this The wrapper to move to.
  !! @param other The other wrapper to move from. Will be deallocated.
  subroutine scalar_scheme_wrapper_move_from(this, other)
    class(scalar_scheme_wrapper_t), intent(inout) :: this
    class(scalar_scheme_wrapper_t), intent(inout) :: other

    ! Move the pointer
    call move_alloc(other%scalar, this%scalar)

  end subroutine scalar_scheme_wrapper_move_from

  !> Return allocation status.
  !! @param this The wrapper to check.
  function scalar_scheme_wrapper_is_allocated(this) result(is_alloc)
    class(scalar_scheme_wrapper_t), intent(in) :: this
    logical :: is_alloc
    is_alloc = allocated(this%scalar)
  end function scalar_scheme_wrapper_is_allocated

end module scalar_scheme
