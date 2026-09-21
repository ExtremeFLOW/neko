!> @file adjoint_scalar_scheme.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Contains the adjoint_scalar_scheme_t type.

module adjoint_scalar_scheme
  use gather_scatter, only : gs_t
  use checkpoint, only : chkp_t
  use num_types, only: rp
  use field, only : field_t
  use field_list, only: field_list_t
  use space, only : space_t
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t, krylov_solver_factory, KSP_MAX_ITER, ksp_monitor_t
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use hsmg, only : hsmg_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use precon, only : pc_t, precon_allocator, precon_destroy
  use field_dirichlet, only: field_dirichlet_t, field_dirichlet_update
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use facet_zone, only : facet_zone_t
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use registry, only : neko_registry
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file
  use user_intf, only : user_t, dummy_user_material_properties, &
       user_material_properties_intf
  use utils, only : neko_error, neko_warning
  use comm, only: NEKO_COMM
  use scalar_source_term, only : scalar_source_term_t
  use field_series, only : field_series_t
  use math, only : cfill, add2s2
  use field_math, only : field_cmult, field_col3, field_cfill, field_add2, &
       field_col2
  use device_math, only : device_cfill, device_add2s2
  use neko_config, only : NEKO_BCKND_DEVICE
  use field_series, only : field_series_t
  use time_step_controller, only : time_step_controller_t
  use json_utils_ext, only: json_key_fallback
  use scalar_scheme, only: scalar_scheme_precon_factory, &
       scalar_scheme_solver_factory
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use device, only : device_memcpy, DEVICE_TO_HOST
  use field_math, only : field_col3, field_cmult2, field_add2, field_cfill
  use mpi_f08, only: MPI_INTEGER, MPI_SUM
  implicit none

  !> Base type for a scalar advection-diffusion solver.
  type, abstract :: adjoint_scalar_scheme_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name
     !> The name of the corresponding primal scalar
     character(len=:), allocatable :: primal_name
     !> x-component of Velocity
     type(field_t), pointer :: u
     !> y-component of Velocity
     type(field_t), pointer :: v
     !> z-component of Velocity
     type(field_t), pointer :: w
     !> The forward scalar.
     type(field_t), pointer :: s
     !> The adjoint scalar.
     type(field_t), pointer :: s_adj
     !> Lag arrays, i.e. solutions at previous timesteps.
     type(field_series_t) :: s_adj_lag
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
     !> Density.
     type(field_t), pointer :: rho => null()
     !> Thermal diffusivity.
     type(field_t) :: lambda
     !> Specific heat capacity.
     type(field_t) :: cp
     !> Turbulent Prandtl number.
     real(kind=rp) :: pr_turb
     !> Field list with cp and lambda
     type(field_list_t) :: material_properties
     !> Is lambda varying in time? Currently only due to LES models.
     logical :: variable_material_properties = .false.
     ! Lag arrays for the RHS.
     type(field_t) :: abx1, abx2
     procedure(user_material_properties_intf), nopass, pointer :: &
          user_material_properties => null()
     !> Freeze the scheme, i.e. do nothing in step()
     logical :: freeze = .false.
   contains
     !> Constructor for the base type.
     procedure, pass(this) :: scheme_init => adjoint_scalar_scheme_init
     !> Destructor for the base type.
     procedure, pass(this) :: scheme_free => adjoint_scalar_scheme_free
     !> Validate successful initialization.
     procedure, pass(this) :: validate => adjoint_scalar_scheme_validate
     !> Set lambda and cp
     procedure, pass(this) :: set_material_properties => &
          adjoint_scalar_scheme_set_material_properties
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          adjoint_scalar_scheme_update_material_properties
     !> Constructor.
     procedure(adjoint_scalar_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor.
     procedure(adjoint_scalar_scheme_free_intrf), pass(this), deferred :: free
     !> Solve for the current timestep.
     procedure(adjoint_scalar_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint.
     procedure(adjoint_scalar_scheme_restart_intrf), pass(this), deferred :: &
          restart
  end type adjoint_scalar_scheme_t

  !> Abstract interface to initialize a scalar formulation
  abstract interface
     subroutine adjoint_scalar_scheme_init_intrf(this, msh, coef, gs, &
          params_adjoint, params_primal, numerics_params, user, chkp, ulag, &
          vlag, wlag, time_scheme, rho)
       import adjoint_scalar_scheme_t
       import json_file
       import coef_t
       import gs_t
       import mesh_t
       import user_t
       import field_series_t, field_t
       import time_scheme_controller_t
       import rp
       import chkp_t
       class(adjoint_scalar_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(in) :: msh
       type(coef_t), target, intent(in) :: coef
       type(gs_t), target, intent(inout) :: gs
       type(json_file), target, intent(inout) :: params_adjoint
       type(json_file), target, intent(inout) :: params_primal
       type(json_file), target, intent(inout) :: numerics_params
       type(user_t), target, intent(in) :: user
       type(chkp_t), target, intent(inout) :: chkp
       type(field_series_t), target, intent(in) :: ulag, vlag, wlag
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
       type(field_t), target, intent(in) :: rho
     end subroutine adjoint_scalar_scheme_init_intrf
  end interface

  !> Abstract interface to restart a scalar formulation
  abstract interface
     subroutine adjoint_scalar_scheme_restart_intrf(this, chkp)
       import adjoint_scalar_scheme_t
       import chkp_t
       import rp
       class(adjoint_scalar_scheme_t), target, intent(inout) :: this
       type(chkp_t), intent(inout) :: chkp
     end subroutine adjoint_scalar_scheme_restart_intrf
  end interface

  !> Abstract interface to dealocate a scalar formulation
  abstract interface
     subroutine adjoint_scalar_scheme_free_intrf(this)
       import adjoint_scalar_scheme_t
       class(adjoint_scalar_scheme_t), intent(inout) :: this
     end subroutine adjoint_scalar_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine adjoint_scalar_scheme_step_intrf(this, time, ext_bdf, &
          dt_controller, ksp_results)
       import adjoint_scalar_scheme_t
       import time_state_t
       import time_scheme_controller_t
       import time_step_controller_t
       import ksp_monitor_t
       class(adjoint_scalar_scheme_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
       type(time_scheme_controller_t), intent(in) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
       type(ksp_monitor_t), intent(inout) :: ksp_results
     end subroutine adjoint_scalar_scheme_step_intrf
  end interface

contains

  !> Initialize all related components of the current scheme
  !! @param[inout] this The object.
  !! @param msh The mesh.
  !! @param c_Xh The coefficients.
  !! @param gs_Xh The gather-scatter.
  !! @param params_adjoint The parameter dictionary in json for the adjoint.
  !! @param params_primal The parameter dictionary in json for the primal.
  !! @param scheme The name of the scalar scheme.
  !! @param user Type with user-defined procedures.
  !! @param rho The density of the fluid.
  subroutine adjoint_scalar_scheme_init(this, msh, c_Xh, gs_Xh, &
       params_adjoint, params_primal, scheme, user, rho)
    class(adjoint_scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: c_Xh
    type(gs_t), target, intent(inout) :: gs_Xh
    type(json_file), target, intent(inout) :: params_primal, params_adjoint
    character(len=*), intent(in) :: scheme
    type(user_t), target, intent(in) :: user
    type(field_t), target, intent(in) :: rho
    type(json_file), pointer :: params_selected
    ! IO buffer for log output
    character(len=LOG_SIZE) :: log_buf
    ! Variables for retrieving json parameters
    logical :: logical_val
    real(kind=rp) :: solver_abstol
    integer :: integer_val
    character(len=:), allocatable :: solver_type, solver_precon
    type(json_file) :: precon_params

    this%u => neko_registry%get_field('u')
    this%v => neko_registry%get_field('v')
    this%w => neko_registry%get_field('w')
    this%rho => rho

    ! get the primal adjoint's name
    call json_get_or_default(params_adjoint, 'primal_name', this%primal_name, &
         's')
    ! Assign a name
    call json_get_or_default(params_adjoint, 'name', this%name, &
         this%primal_name // '_adj')

    call neko_log%section('Adjoint scalar')
    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.type')
    call json_get(params_selected, 'solver.type', solver_type)

    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.preconditioner.type')
    call json_get(params_selected, 'solver.preconditioner.type', solver_precon)

    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.preconditioner')
    call json_get(params_selected, 'solver.preconditioner', &
         precon_params)

    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.absolute_tolerance')
    call json_get(params_selected, 'solver.absolute_tolerance', &
         solver_abstol)

    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.projection_space_size')
    call json_get_or_default(params_selected, &
         'solver.projection_space_size', &
         this%projection_dim, 0)

    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.projection_hold_steps')
    call json_get_or_default(params_selected, &
         'solver.projection_hold_steps', &
         this%projection_activ_step, 5)


    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)
    call neko_log%message('Ksp adjoint scalar : ('// trim(solver_type) // &
         ', ' // trim(solver_precon) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :', solver_abstol
    call neko_log%message(log_buf)

    this%Xh => this%u%Xh
    this%dm_Xh => this%u%dof
    this%params => params_adjoint
    this%msh => msh

    if (.not. neko_registry%field_exists(this%name)) then
       call neko_registry%add_field(this%dm_Xh, this%name)
    end if

    this%s_adj => neko_registry%get_field(this%name)

    call this%s_adj_lag%init(this%s_adj, 2)

    this%s => neko_registry%get_field(this%primal_name)

    this%gs_Xh => gs_Xh
    this%c_Xh => c_Xh

    !
    ! Material properties
    !
    call this%set_material_properties(params_primal, user)


    !
    ! Turbulence modelling and variable material properties
    !
    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'variable_material_properties')
    if (params_selected%valid_path('variable_material_properties')) then
       call neko_error('variable material properties no supported for adjoint')
    end if

    write(log_buf, '(A,L1)') 'LES        : ', this%variable_material_properties
    call neko_log%message(log_buf)

    !
    ! Setup right-hand side field.
    !
    allocate(this%f_Xh)
    call this%f_Xh%init(this%dm_Xh, fld_name = "adjoint_scalar_rhs")

    ! Initialize the source term
    call this%source_term%init(this%f_Xh, this%c_Xh, user, this%name)
    ! We should ONLY read the adjoint
    call this%source_term%add(params_primal, 'source_terms')

    ! todo parameter file ksp tol should be added
    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.max_iterations')
    call json_get_or_default(params_selected, &
         'solver.max_iterations', &
         integer_val, KSP_MAX_ITER)
    params_selected => json_key_fallback(params_adjoint, params_primal, &
         'solver.monitor')
    call json_get_or_default(params_selected, &
         'solver.monitor', &
         logical_val, .false.)
    call adjoint_scalar_scheme_solver_factory(this%ksp, this%dm_Xh%size(), &
         solver_type, integer_val, solver_abstol, logical_val)
    call scalar_scheme_precon_factory(this%pc, this%ksp, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs, &
         solver_precon, precon_params)

    call neko_log%end_section()

  end subroutine adjoint_scalar_scheme_init


  !> Deallocate a scalar formulation
  subroutine adjoint_scalar_scheme_free(this)
    class(adjoint_scalar_scheme_t), intent(inout) :: this
    class(bc_t), pointer :: bc
    integer :: i

    bc => null()

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

    call this%source_term%free()

    if (associated(this%f_Xh)) then
       call this%f_Xh%free()
       deallocate(this%f_Xh)
       nullify(this%f_Xh)
    end if

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)
       if (associated(bc)) then
          call bc%free()
          deallocate(bc)
       end if
    end do

    call this%bcs%free()

    call this%cp%free()
    call this%lambda%free()
    call this%s_adj_lag%free()

    nullify(bc)

  end subroutine adjoint_scalar_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine adjoint_scalar_scheme_validate(this)
    class(adjoint_scalar_scheme_t), target, intent(inout) :: this

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%s%x)) .or. &
         (.not. allocated(this%s_adj%x))) then
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

  end subroutine adjoint_scalar_scheme_validate

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine adjoint_scalar_scheme_solver_factory(ksp, n, solver, max_iter, &
       abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    integer, intent(in) :: max_iter
    character(len=*), intent(in) :: solver
    real(kind=rp) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, &
         abstol, monitor = monitor)

  end subroutine adjoint_scalar_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine adjoint_scalar_scheme_precon_factory(pc, ksp, coef, dof, gs, &
       bclst, pctype, pcparams)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype
    type(json_file), intent(inout) :: pcparams

    call precon_allocator(pc, pctype)

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

  end subroutine adjoint_scalar_scheme_precon_factory

  !> Call user material properties routine and update the values of `lambda`
  !! if necessary.
  !! @param[inout] this The object.
  !! @param time The time state.
  subroutine adjoint_scalar_scheme_update_material_properties(this, time)
    class(adjoint_scalar_scheme_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: nut
    integer :: index
    ! Factor to transform nu_t to lambda_t
    type(field_t), pointer :: lambda_factor

    call this%user_material_properties(this%name, this%material_properties, &
         time)

    ! factor = rho * cp / pr_turb
    if (this%variable_material_properties .and. &
         len(trim(this%nut_field_name)) > 0) then
       nut => neko_registry%get_field(this%nut_field_name)

       ! lambda = lambda + rho * cp * nut / pr_turb
       call neko_scratch_registry%request_field(lambda_factor, index, .false.)

       call field_col3(lambda_factor, this%cp, this%rho)
       call field_col2(lambda_factor, nut)
       call field_cmult(lambda_factor, 1.0_rp / this%pr_turb)
       call field_add2(this%lambda, lambda_factor)
       call neko_scratch_registry%relinquish_field(index)
    end if

    ! Since cp is a fields and we use the %x(1,1,1,1) of the
    ! host array data to pass constant  material properties
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%cp%x, this%cp%x_d, this%cp%size(), &
            DEVICE_TO_HOST, sync=.false.)
    end if

  end subroutine adjoint_scalar_scheme_update_material_properties

  !> Set lamdba and cp.
  !! @param[inout] this The object.
  !! @param params_primal The case file configuration dictionary.
  !! @param user The user interface.
  subroutine adjoint_scalar_scheme_set_material_properties(this, &
       params_primal, user)
    class(adjoint_scalar_scheme_t), intent(inout) :: this
    type(json_file), intent(inout) :: params_primal
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties_intf), pointer :: dummy_mp_ptr
    real(kind=rp) :: const_cp, const_lambda
    ! Dummy time state set to 0
    type(time_state_t) :: time

    dummy_mp_ptr => dummy_user_material_properties

    ! Fill lambda field with the physical value
    call this%lambda%init(this%dm_Xh, "lambda")
    call this%cp%init(this%dm_Xh, "cp")

    call this%material_properties%init(2)
    call this%material_properties%assign_to_field(1, this%cp)
    call this%material_properties%assign_to_field(2, this%lambda)

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user " // &
            "file!"
       call neko_log%message(log_buf)
       this%user_material_properties => user%material_properties
       call user%material_properties(this%name, this%material_properties, time)
    else
       this%user_material_properties => dummy_user_material_properties
       if (params_primal%valid_path('Pe') .and. &
            (params_primal%valid_path('lambda') .or. &
            params_primal%valid_path('cp'))) then
          call neko_error("To set the material properties for the scalar, " // &
               "either provide Pe OR lambda and cp in the case file.")
          ! Non-dimensional case
       else if (params_primal%valid_path('Pe')) then
          write(log_buf, '(A)') 'Non-dimensional scalar material properties' //&
               ' input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Specific heat capacity will be set to 1,'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'conductivity to 1/Pe. Assumes density is 1.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Pe into lambda for further manipulation.
          call json_get(params_primal, 'Pe', const_lambda)
          write(log_buf, '(A,ES13.6)') 'Pe         :', const_lambda
          call neko_log%message(log_buf)

          ! Set cp and rho to 1 since the setup is non-dimensional.
          const_cp = 1.0_rp
          ! Invert the Pe to get conductivity
          const_lambda = 1.0_rp/const_lambda
          ! Dimensional case
       else
          call json_get(params_primal, 'lambda', const_lambda)
          call json_get(params_primal, 'cp', const_cp)
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

    ! Since cp is a field and we use the %x(1,1,1,1) of the
    ! host array data to pass constant material properties
    ! to some routines, we need to make sure that the host
    ! values are also filled
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%cp%x, this%cp%x_d, this%cp%size(), &
            DEVICE_TO_HOST, sync=.false.)
    end if
  end subroutine adjoint_scalar_scheme_set_material_properties

end module adjoint_scalar_scheme
