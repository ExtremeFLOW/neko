! Copyright (c) 2020-2024, The Neko Authors
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
module fluid_scheme
  use gather_scatter, only : gs_t
  use mean_sqr_flow, only : mean_sqr_flow_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use checkpoint, only : chkp_t
  use mean_flow, only : mean_flow_t
  use num_types, only : rp
  use comm
  use fluid_source_term, only: fluid_source_term_t
  use field, only : field_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t, krylov_solver_factory, krylov_solver_destroy, &
                     KSP_MAX_ITER
  use coefs, only: coef_t
  use wall, only : no_slip_wall_t
  use inflow, only : inflow_t
  use usr_inflow, only : usr_inflow_t, usr_inflow_eval
  use blasius, only : blasius_t
  use dirichlet, only : dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use symmetry, only : symmetry_t
  use field_dirichlet, only : field_dirichlet_t
  use field_dirichlet_vector, only: field_dirichlet_vector_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use fluid_stats, only : fluid_stats_t
  use bc, only : bc_t, bc_list_t, bc_list_init, bc_list_add, bc_list_free, &
       bc_list_apply_scalar, bc_list_apply_vector
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use math, only : cfill, add2s2
  use device_math, only : device_cfill, device_add2s2
  use time_scheme_controller, only : time_scheme_controller_t
  use operators, only : cfl
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use scratch_registry, only : scratch_registry_t
  use user_intf, only : user_t, dummy_user_material_properties, &
                        user_material_properties
  use utils, only : neko_error
  use field_series, only : field_series_t
  use time_step_controller, only : time_step_controller_t
  use field_math, only : field_cfill
  use gradient_jump_penalty, only : gradient_jump_penalty_t
  implicit none
  private

  !> Base type of all fluid formulations
  type, abstract :: fluid_scheme_t
     type(field_t), pointer :: u => null() !< x-component of Velocity
     type(field_t), pointer :: v => null() !< y-component of Velocity
     type(field_t), pointer :: w => null() !< z-component of Velocity
     type(field_t), pointer :: p => null() !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$
     !> The source term for the momentum equation.
     type(fluid_source_term_t) :: source_term
     !> X-component of the right-hand side.
     type(field_t), pointer :: f_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_z => null()
     class(ksp_t), allocatable :: ksp_vel     !< Krylov solver for velocity
     class(ksp_t), allocatable :: ksp_prs     !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel        !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs        !< Velocity Preconditioner
     integer :: vel_projection_dim         !< Size of the projection space for ksp_vel
     integer :: pr_projection_dim          !< Size of the projection space for ksp_pr
     integer :: vel_projection_activ_step  !< Steps to activate projection for ksp_vel
     integer :: pr_projection_activ_step   !< Steps to activate projection for ksp_pr
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     class(bc_t), allocatable :: bc_inflow !< Dirichlet inflow for velocity
     !> Gradient jump panelty
     logical :: if_gradient_jump_penalty
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_u
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_v
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_w

     ! Attributes for field dirichlet BCs
     type(field_dirichlet_vector_t) :: user_field_bc_vel   !< User-computed Dirichlet velocity condition
     type(field_dirichlet_t) :: user_field_bc_prs   !< User-computed Dirichlet pressure condition
     type(dirichlet_t) :: bc_prs               !< Dirichlet pressure condition
     type(dong_outflow_t) :: bc_dong           !< Dong outflow condition
     type(symmetry_t) :: bc_sym                !< Symmetry plane for velocity
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_prs              !< List of pressure conditions
     type(field_t) :: bdry                     !< Boundary markings
     type(json_file), pointer :: params        !< Parameters
     type(mesh_t), pointer :: msh => null()    !< Mesh
     type(chkp_t) :: chkp                      !< Checkpoint
     type(mean_flow_t) :: mean                 !< Mean flow field
     type(fluid_stats_t) :: stats              !< Fluid statistics
     type(mean_sqr_flow_t) :: mean_sqr         !< Mean squared flow field
     logical :: forced_flow_rate = .false.     !< Is the flow rate forced?
     logical :: freeze = .false.               !< Freeze velocity at initial condition?
     !> Dynamic viscosity
     real(kind=rp) :: mu
     !> The variable mu field
     type(field_t) :: mu_field
     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name
     !> Is mu varying in time? Currently only due to LES models.
     logical :: variable_material_properties = .false.
     !> Density
     real(kind=rp) :: rho
     !> The variable density field
     type(field_t) :: rho_field
     type(scratch_registry_t) :: scratch       !< Manager for temporary fields
     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)
   contains
     !> Constructor for the base type
     procedure, pass(this) :: fluid_scheme_init_all
     procedure, pass(this) :: fluid_scheme_init_common
     generic :: scheme_init => fluid_scheme_init_all, fluid_scheme_init_common
     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => fluid_scheme_validate
     !> Apply pressure boundary conditions
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     !> Apply velocity boundary conditions
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     !> Set the user inflow procedure
     procedure, pass(this) :: set_usr_inflow => fluid_scheme_set_usr_inflow
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => fluid_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: set_material_properties => &
          fluid_scheme_set_material_properties
     !> Constructor
     procedure(fluid_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor
     procedure(fluid_scheme_free_intrf), pass(this), deferred :: free
     !> Advance one step in time
     procedure(fluid_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint
     procedure(fluid_scheme_restart_intrf), pass(this), deferred :: restart
     procedure, private, pass(this) :: set_bc_type_output => &
       fluid_scheme_set_bc_type_output
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
       fluid_scheme_update_material_properties
  end type fluid_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_scheme_init_intrf(this, msh, lx, params, user, &
          time_scheme)
       import fluid_scheme_t
       import json_file
       import mesh_t
       import user_t
       import time_scheme_controller_t
       class(fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), intent(in) :: user
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
     end subroutine fluid_scheme_init_intrf
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine fluid_scheme_free_intrf(this)
       import fluid_scheme_t
       class(fluid_scheme_t), intent(inout) :: this
     end subroutine fluid_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine fluid_scheme_step_intrf(this, t, tstep, dt, ext_bdf, &
                                        dt_controller)
       import fluid_scheme_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(fluid_scheme_t), target, intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(inout) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine fluid_scheme_step_intrf
  end interface

  !> Abstract interface to restart a fluid scheme
  abstract interface
     subroutine fluid_scheme_restart_intrf(this, dtlag, tlag)
       import fluid_scheme_t
       import rp
       class(fluid_scheme_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)

     end subroutine fluid_scheme_restart_intrf
  end interface

  interface
     !> Initialise a fluid scheme
     module subroutine fluid_scheme_factory(object, type_name)
       class(fluid_scheme_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine fluid_scheme_factory
  end interface

  public :: fluid_scheme_t, fluid_scheme_factory

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_common(this, msh, lx, params, scheme, user, &
      kspv_init)
    implicit none
    class(fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical, intent(in) :: kspv_init
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), allocatable :: real_vec(:)
    real(kind=rp) :: real_val
    logical :: logical_val
    integer :: integer_val, ierr
    character(len=:), allocatable :: string_val1, string_val2

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
    this%scratch = scratch_registry_t(this%dm_Xh, 10, 2)

    ! Case parameters
    this%params => params


    !
    ! First section of fluid log
    !

    call neko_log%section('Fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)

    !
    ! Material properties
    !
    call this%set_material_properties(params, user)

    !
    ! Turbulence modelling and variable material properties
    !
    if (params%valid_path('case.fluid.nut_field')) then
       call json_get(params, 'case.fluid.nut_field', this%nut_field_name)
       this%variable_material_properties = .true.
    else
       this%nut_field_name = ""
    end if

    ! Fill mu and rho field with the physical value

    call this%mu_field%init(this%dm_Xh, "mu")
    call this%rho_field%init(this%dm_Xh, "mu")
    call field_cfill(this%mu_field, this%mu, this%mu_field%size())
    call field_cfill(this%rho_field, this%rho, this%mu_field%size())

    ! Since mu, rho is a field, and the none-stress simulation fetches
    ! data from the host arrays, we need to mirror the constant
    ! material properties on the host
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(this%mu_field%x, this%mu, this%mu_field%size())
       call cfill(this%rho_field%x, this%rho, this%rho_field%size())
    end if


    ! Projection spaces
    call json_get_or_default(params, &
                        'case.fluid.velocity_solver.projection_space_size', &
                        this%vel_projection_dim, 20)
    call json_get_or_default(params, &
                        'case.fluid.pressure_solver.projection_space_size', &
                        this%pr_projection_dim, 20)
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
    write(log_buf, '(A, I0)') 'DoFs       : ', this%dm_Xh%size()
    call neko_log%message(log_buf)

    write(log_buf, '(A,ES13.6)') 'rho        :',  this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'mu         :',  this%mu
    call neko_log%message(log_buf)

    call json_get(params, 'case.numerics.dealias', logical_val)
    write(log_buf, '(A, L1)') 'Dealias    : ',  logical_val
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.output_boundary', logical_val, &
                             .false.)
    write(log_buf, '(A, L1)') 'Save bdry  : ',  logical_val
    call neko_log%message(log_buf)
    if (logical_val) then
       write(log_buf, '(A)') 'bdry keys: '
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'No-slip wall             ''w'' = 1'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Inlet/velocity dirichlet ''v'' = 2'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Outlet                   ''o'' = 3'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Symmetry               ''sym'' = 4'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Outlet-normal           ''on'' = 5'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Periodic                     = 6'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Dirichlet on velocity components: '
       call neko_log%message(log_buf)
       write(log_buf, '(A)') ' ''d_vel_u, d_vel_v, d_vel_w'' = 7'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'Pressure dirichlet  ''d_pres'' = 8'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') '''d_vel_(u,v,w)'' and ''d_pres'' = 8'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') 'No boundary condition set    = 0'
       call neko_log%message(log_buf)
    end if


    !
    ! Setup velocity boundary conditions
    !
    allocate(this%bc_labels(NEKO_MSH_MAX_ZLBLS))
    this%bc_labels = "not"

    if (params%valid_path('case.fluid.boundary_types')) then
       call json_get(params, &
                     'case.fluid.boundary_types', &
                     this%bc_labels)
    end if

    call bc_list_init(this%bclst_vel)

    call this%bc_sym%init_base(this%c_Xh)
    call this%bc_sym%mark_zone(msh%sympln)
    call this%bc_sym%mark_zones_from_list(msh%labeled_zones, &
                        'sym', this%bc_labels)
    call this%bc_sym%finalize()
    call this%bc_sym%init(this%c_Xh)
    call bc_list_add(this%bclst_vel, this%bc_sym)

    !
    ! Inflow
    !
    if (params%valid_path('case.fluid.inflow_condition')) then
       call json_get(params, 'case.fluid.inflow_condition.type', string_val1)
       if (trim(string_val1) .eq. "uniform") then
          allocate(inflow_t::this%bc_inflow)
       else if (trim(string_val1) .eq. "blasius") then
          allocate(blasius_t::this%bc_inflow)
       else if (trim(string_val1) .eq. "user") then
          allocate(usr_inflow_t::this%bc_inflow)
       else
          call neko_error('Invalid inflow condition '//string_val1)
       end if

       call this%bc_inflow%init_base(this%c_Xh)
       call this%bc_inflow%mark_zone(msh%inlet)
       call this%bc_inflow%mark_zones_from_list(msh%labeled_zones, &
                        'v', this%bc_labels)
       call this%bc_inflow%finalize()
       call bc_list_add(this%bclst_vel, this%bc_inflow)

       if (trim(string_val1) .eq. "uniform") then
          call json_get(params, 'case.fluid.inflow_condition.value', real_vec)
          select type (bc_if => this%bc_inflow)
          type is (inflow_t)
             call bc_if%set_inflow(real_vec)
          end select
       else if (trim(string_val1) .eq. "blasius") then
          select type (bc_if => this%bc_inflow)
          type is (blasius_t)
             call json_get(params, 'case.fluid.blasius.delta', real_val)
             call json_get(params, 'case.fluid.blasius.approximation', &
                           string_val2)
             call json_get(params, 'case.fluid.blasius.freestream_velocity', &
                           real_vec)

             call bc_if%set_params(real_vec, real_val, string_val2)

          end select
       else if (trim(string_val1) .eq. "user") then
       end if
    end if

    call this%bc_wall%init_base(this%c_Xh)
    call this%bc_wall%mark_zone(msh%wall)
    call this%bc_wall%mark_zones_from_list(msh%labeled_zones, &
                        'w', this%bc_labels)
    call this%bc_wall%finalize()
    call bc_list_add(this%bclst_vel, this%bc_wall)

    ! Setup field dirichlet bc for u-velocity
    call this%user_field_bc_vel%bc_u%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_u%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_u', this%bc_labels)
    call this%user_field_bc_vel%bc_u%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_u%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0)  then
      call this%user_field_bc_vel%bc_u%init_field('d_vel_u')
    end if

    ! Setup field dirichlet bc for v-velocity
    call this%user_field_bc_vel%bc_v%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_v%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_v', this%bc_labels)
    call this%user_field_bc_vel%bc_v%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_v%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0)  then
      call this%user_field_bc_vel%bc_v%init_field('d_vel_v')
    end if

    ! Setup field dirichlet bc for w-velocity
    call this%user_field_bc_vel%bc_w%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_w%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_w', this%bc_labels)
    call this%user_field_bc_vel%bc_w%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_w%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0)  then
      call this%user_field_bc_vel%bc_w%init_field('d_vel_w')
    end if

    ! Setup our global field dirichlet bc
    call this%user_field_bc_vel%init_base(this%c_Xh)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_u', this%bc_labels)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_v', this%bc_labels)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones, &
                        'd_vel_w', this%bc_labels)
    call this%user_field_bc_vel%finalize()

    ! Add the field bc to velocity bcs
    call bc_list_add(this%bclst_vel, this%user_field_bc_vel)

    !
    ! Associate our field dirichlet update to the user one.
    !
    this%user_field_bc_vel%update => user%user_dirichlet_update

    !
    ! Initialize field list and bc list for user_dirichlet_update
    !

    ! Note, some of these are potentially not initialized !
    call this%user_field_bc_vel%field_list%init(4)
    call this%user_field_bc_vel%field_list%assign_to_field(1, &
            this%user_field_bc_vel%bc_u%field_bc)
    call this%user_field_bc_vel%field_list%assign_to_field(2, &
            this%user_field_bc_vel%bc_v%field_bc)
    call this%user_field_bc_vel%field_list%assign_to_field(3, &
            this%user_field_bc_vel%bc_w%field_bc)
    call this%user_field_bc_vel%field_list%assign_to_field(4, &
            this%user_field_bc_prs%field_bc)

    call bc_list_init(this%user_field_bc_vel%bc_list, size = 4)
    ! Note, bc_list_add only adds if the bc is not empty
    call bc_list_add(this%user_field_bc_vel%bc_list, &
                     this%user_field_bc_vel%bc_u)
    call bc_list_add(this%user_field_bc_vel%bc_list, &
                     this%user_field_bc_vel%bc_v)
    call bc_list_add(this%user_field_bc_vel%bc_list, &
                     this%user_field_bc_vel%bc_w)

    !
    ! Check if we need to output boundary types to a separate field
    call fluid_scheme_set_bc_type_output(this, params)

    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_x)
    allocate(this%f_y)
    allocate(this%f_z)
    call this%f_x%init(this%dm_Xh, fld_name = "fluid_rhs_x")
    call this%f_y%init(this%dm_Xh, fld_name = "fluid_rhs_y")
    call this%f_z%init(this%dm_Xh, fld_name = "fluid_rhs_z")

    ! Initialize the source term
    call this%source_term%init(this%f_x, this%f_y, this%f_z, this%c_Xh, user)
    call this%source_term%add(params, 'case.fluid.source_terms')
    
    ! If case.output_boundary is true, set the values for the bc types in the
    ! output of the field.
    call this%set_bc_type_output(params)

    ! Initialize velocity solver
    if (kspv_init) then
       call neko_log%section("Velocity solver")
       call json_get_or_default(params, &
                                'case.fluid.velocity_solver.max_iterations', &
                                integer_val, KSP_MAX_ITER)
       call json_get(params, 'case.fluid.velocity_solver.type', string_val1)
       call json_get(params, 'case.fluid.velocity_solver.preconditioner', &
                     string_val2)
       call json_get(params, 'case.fluid.velocity_solver.absolute_tolerance', &
                     real_val)
       call json_get_or_default(params, &
                                'case.fluid.velocity_solver.monitor', &
                                logical_val, .false.)
       
       call neko_log%message('Type       : ('// trim(string_val1) // &
           ', ' // trim(string_val2) // ')')

       write(log_buf, '(A,ES13.6)') 'Abs tol    :',  real_val
       call neko_log%message(log_buf)
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            string_val1, integer_val, real_val, logical_val)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, string_val2)
       call neko_log%end_section()
    end if

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, 'u')
    call neko_field_registry%add_field(this%dm_Xh, 'v')
    call neko_field_registry%add_field(this%dm_Xh, 'w')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    !! Initialize time-lag fields
    call this%ulag%init(this%u, 2)
    call this%vlag%init(this%v, 2)
    call this%wlag%init(this%w, 2)


  end subroutine fluid_scheme_init_common

  !> Initialize all components of the current scheme
  subroutine fluid_scheme_init_all(this, msh, lx, params, kspv_init, &
       kspp_init, scheme, user)
    implicit none
    class(fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical :: kspv_init
    logical :: kspp_init    
    character(len=*), intent(in) :: scheme
    real(kind=rp) :: abs_tol
    integer :: integer_val, ierr
    logical :: logical_val
    character(len=:), allocatable :: solver_type, precon_type
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: GJP_param_a, GJP_param_b

    call fluid_scheme_init_common(this, msh, lx, params, scheme, user, &
         kspv_init)

    call neko_field_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_field_registry%get_field('p')

    !
    ! Setup pressure boundary conditions
    !
    call bc_list_init(this%bclst_prs)
    call this%bc_prs%init_base(this%c_Xh)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones, &
                        'o', this%bc_labels)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones, &
                        'on', this%bc_labels)

    ! Field dirichlet pressure bc
    call this%user_field_bc_prs%init_base(this%c_Xh)
    call this%user_field_bc_prs%mark_zones_from_list(msh%labeled_zones, &
                        'd_pres', this%bc_labels)
    call this%user_field_bc_prs%finalize()
    call MPI_Allreduce(this%user_field_bc_prs%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    if (integer_val .gt. 0) call this%user_field_bc_prs%init_field('d_pres')
    call bc_list_add(this%bclst_prs, this%user_field_bc_prs)
    call bc_list_add(this%user_field_bc_vel%bc_list, this%user_field_bc_prs)

    if (msh%outlet%size .gt. 0) then
       call this%bc_prs%mark_zone(msh%outlet)
    end if
    if (msh%outlet_normal%size .gt. 0) then
       call this%bc_prs%mark_zone(msh%outlet_normal)
    end if

    call this%bc_prs%finalize()
    call this%bc_prs%set_g(0.0_rp)
    call bc_list_add(this%bclst_prs, this%bc_prs)
    call this%bc_dong%init_base(this%c_Xh)
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones, &
                        'o+dong', this%bc_labels)
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones, &
                        'on+dong', this%bc_labels)
    call this%bc_dong%finalize()


    call this%bc_dong%init(this%c_Xh, params)

    call bc_list_add(this%bclst_prs, this%bc_dong)

    ! Pressure solver
    if (kspp_init) then
       call neko_log%section("Pressure solver")

       call json_get_or_default(params, &
                               'case.fluid.pressure_solver.max_iterations', &
                               integer_val, KSP_MAX_ITER)
       call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
       call json_get(params, 'case.fluid.pressure_solver.preconditioner', &
                     precon_type)
       call json_get(params, 'case.fluid.pressure_solver.absolute_tolerance', &
                     abs_tol)
       call json_get_or_default(params, &
                                'case.fluid.pressure_solver.monitor', &
                                logical_val, .false.)
       call neko_log%message('Type       : ('// trim(solver_type) // &
             ', ' // trim(precon_type) // ')')
       write(log_buf, '(A,ES13.6)') 'Abs tol    :',  abs_tol
       call neko_log%message(log_buf)

       call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Xh%size(), &
            solver_type, integer_val, abs_tol, logical_val)
       call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs, precon_type)

       call neko_log%end_section()

    end if

    ! Initiate gradient jump penalty
    call json_get_or_default(params, &
                            'case.fluid.gradient_jump_penalty.enabled',&
                            this%if_gradient_jump_penalty, .false.)

    if (this%if_gradient_jump_penalty .eqv. .true.) then
       if ((this%dm_Xh%xh%lx - 1) .eq. 1) then
          call json_get_or_default(params, &
                            'case.fluid.gradient_jump_penalty.tau',&
                            GJP_param_a, 0.02_rp)
          GJP_param_b = 0.0_rp
       else
          call json_get_or_default(params, &
                        'case.fluid.gradient_jump_penalty.scaling_factor',&
                            GJP_param_a, 0.8_rp)
          call json_get_or_default(params, &
                        'case.fluid.gradient_jump_penalty.scaling_exponent',&
                            GJP_param_b, 4.0_rp)
       end if
       call this%gradient_jump_penalty_u%init(params, this%dm_Xh, this%c_Xh, &
                                              GJP_param_a, GJP_param_b)
       call this%gradient_jump_penalty_v%init(params, this%dm_Xh, this%c_Xh, &
                                              GJP_param_a, GJP_param_b)
       call this%gradient_jump_penalty_w%init(params, this%dm_Xh, this%c_Xh, &
                                              GJP_param_a, GJP_param_b)
    end if

    call neko_log%end_section()

  end subroutine fluid_scheme_init_all

  !> Deallocate a fluid formulation
  subroutine fluid_scheme_free(this)
    class(fluid_scheme_t), intent(inout) :: this

    call this%bdry%free()

    if (allocated(this%bc_inflow)) then
       call this%bc_inflow%free()
    end if

    call this%bc_wall%free()
    call this%bc_sym%free()

    !
    ! Free everything related to field_dirichlet BCs
    !
    call this%user_field_bc_prs%field_bc%free()
    call this%user_field_bc_prs%free()
    call this%user_field_bc_vel%bc_u%field_bc%free()
    call this%user_field_bc_vel%bc_v%field_bc%free()
    call this%user_field_bc_vel%bc_w%field_bc%free()
    call this%user_field_bc_vel%free()

    call this%Xh%free()

    if (allocated(this%ksp_vel)) then
       call krylov_solver_destroy(this%ksp_vel)
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call krylov_solver_destroy(this%ksp_prs)
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

    if (allocated(this%bc_labels)) then
       deallocate(this%bc_labels)
    end if

    call this%source_term%free()

    call this%gs_Xh%free()

    call this%c_Xh%free()

    call bc_list_free(this%bclst_vel)

    call this%scratch%free()

    nullify(this%params)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()


    if (associated(this%f_x)) then
       call this%f_x%free()
    end if

    if (associated(this%f_y)) then
       call this%f_y%free()
    end if

    if (associated(this%f_z)) then
       call this%f_z%free()
    end if

    nullify(this%f_x)
    nullify(this%f_y)
    nullify(this%f_z)

    call this%mu_field%free()

    ! Free gradient jump penalty
    if (this%if_gradient_jump_penalty .eqv. .true.) then
       call this%gradient_jump_penalty_u%free()
       call this%gradient_jump_penalty_v%free()
       call this%gradient_jump_penalty_w%free()
    end if


  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_t), target, intent(inout) :: this
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

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

    if (allocated(this%bc_inflow)) then
       select type (ip => this%bc_inflow)
       type is (usr_inflow_t)
          call ip%validate
       end select
    end if

    !
    ! Setup checkpoint structure (if everything is fine)
    !
    call this%chkp%init(this%u, this%v, this%w, this%p)

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! Here we perform additional gs operations to take care of
  !! shared points between elements that have different BCs, as done in Nek5000.
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this, t, tstep)
    class(fluid_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call bc_list_apply_vector(this%bclst_vel, &
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size(), t, tstep)

  end subroutine fluid_scheme_bc_apply_vel

  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this, t, tstep)
    class(fluid_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call bc_list_apply_scalar(this%bclst_prs, this%p%x, &
                              this%p%dof%size(), t, tstep)

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
  subroutine fluid_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, &
                                         pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst, &
                  trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)

  end subroutine fluid_scheme_precon_factory

  !> Initialize a user defined inflow condition
  subroutine fluid_scheme_set_usr_inflow(this, usr_eval)
    class(fluid_scheme_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval

    select type (bc_if => this%bc_inflow)
    type is (usr_inflow_t)
       call bc_if%set_eval(usr_eval)
    class default
       call neko_error("Not a user defined inflow condition")
    end select
  end subroutine fluid_scheme_set_usr_inflow

  !> Compute CFL
  function fluid_compute_cfl(this, dt) result(c)
    class(fluid_scheme_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u%x, this%v%x, this%w%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)

  end function fluid_compute_cfl

  !> Set boundary types for the diagnostic output.
  !! @param params The JSON case file.
  subroutine fluid_scheme_set_bc_type_output(this, params)
    class(fluid_scheme_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(dirichlet_t) :: bdry_mask
    logical :: found

    !
    ! Check if we need to output boundaries
    !
    call json_get_or_default(params, 'case.output_boundary', found, .false.)

    if (found) then
      call this%bdry%init(this%dm_Xh, 'bdry')
      this%bdry = 0.0_rp

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%wall)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'w', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(1.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%inlet)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'v', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(2.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%outlet)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'o', this%bc_labels)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'o+dong', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(3.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%sympln)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'sym', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(4.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%outlet_normal)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'on', this%bc_labels)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'on+dong', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(5.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zone(this%msh%periodic)
      call bdry_mask%finalize()
      call bdry_mask%set_g(6.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'd_vel_u', this%bc_labels)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'd_vel_v', this%bc_labels)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'd_vel_w', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(7.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

      call bdry_mask%init_base(this%c_Xh)
      call bdry_mask%mark_zones_from_list(this%msh%labeled_zones, &
                     'd_pres', this%bc_labels)
      call bdry_mask%finalize()
      call bdry_mask%set_g(8.0_rp)
      call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
      call bdry_mask%free()

    end if

  end subroutine fluid_scheme_set_bc_type_output

  !> Update the values of `mu_field` if necessary.
  subroutine fluid_scheme_update_material_properties(this)
    class(fluid_scheme_t), intent(inout) :: this
    type(field_t), pointer :: nut
    integer :: n

    if (this%variable_material_properties) then
      nut => neko_field_registry%get_field(this%nut_field_name)
      n = nut%size()

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(this%mu_field%x_d, this%mu, n)
         call device_add2s2(this%mu_field%x_d, nut%x_d, this%rho, n)
      else
         call cfill(this%mu_field%x, this%mu, n)
         call add2s2(this%mu_field%x, nut%x, this%rho, n)
      end if
    end if

  end subroutine fluid_scheme_update_material_properties

  !> Sets rho and mu
  !! @param params The case paramter file.
  !! @param user The user interface.
  subroutine fluid_scheme_set_material_properties(this, params, user)
    class(fluid_scheme_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties),  pointer :: dummy_mp_ptr
    logical :: nondimensional
    real(kind=rp) :: dummy_lambda, dummy_cp

    dummy_mp_ptr => dummy_user_material_properties

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user&
       & file!"
       call neko_log%message(log_buf)
       call user%material_properties(0.0_rp, 0, this%rho, this%mu, &
            dummy_cp, dummy_lambda, params)
    else
       ! Incorrect user input
       if (params%valid_path('case.fluid.Re') .and. &
           (params%valid_path('case.fluid.mu') .or. &
            params%valid_path('case.fluid.rho'))) then
          call neko_error("To set the material properties for the fluid,&
          & either provide Re OR mu and rho in the case file.")

          ! Non-dimensional case
       else if (params%valid_path('case.fluid.Re')) then

          write(log_buf, '(A)') 'Non-dimensional fluid material properties &
          & input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Density will be set to 1, dynamic viscosity to&
          & 1/Re.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Re into mu for further manipulation.
          call json_get(params, 'case.fluid.Re', this%mu)
          write(log_buf, '(A)') 'Read non-dimensional material properties'
          call neko_log%message(log_buf)
          write(log_buf, '(A,ES13.6)') 'Re         :',  this%mu
          call neko_log%message(log_buf)

          ! Set rho to 1 since the setup is non-dimensional.
          this%rho = 1.0_rp
          ! Invert the Re to get viscosity.
          this%mu = 1.0_rp/this%mu
          ! Dimensional case
       else
          call json_get(params, 'case.fluid.mu', this%mu)
          call json_get(params, 'case.fluid.rho', this%rho)
       end if

    end if
  end subroutine fluid_scheme_set_material_properties

end module fluid_scheme
