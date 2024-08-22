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
module adjoint_scheme
  use gather_scatter
  use mean_sqr_flow, only : mean_sqr_flow_t
  use neko_config
  use checkpoint, only : chkp_t
  use mean_flow, only : mean_flow_t
  use num_types
  use comm
  use fluid_user_source_term, only: fluid_user_source_term_t
  use fluid_source_term, only: fluid_source_term_t
  use field_list, only : field_list_t
  use field, only : field_t
  use space
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t
  use coefs, only: coef_t
  use wall, only : no_slip_wall_t
  use inflow, only : inflow_t
  use usr_inflow, only : usr_inflow_t, usr_inflow_eval
  use blasius, only : blasius_t
  use dirichlet, only : dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use symmetry, only : symmetry_t
  use non_normal, only : non_normal_t
  use field_dirichlet, only : field_dirichlet_t, field_dirichlet_update
  use field_dirichlet_vector, only: field_dirichlet_vector_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use precon, only : pc_t
  use krylov_fctry
  use precon_fctry
  use fluid_stats, only : fluid_stats_t
  use bc
  use mesh
  use math
  use device_math, only : device_cfill, device_add2s2
  use time_scheme_controller, only : time_scheme_controller_t
  use mathops
  use operators, only : cfl
  use logger
  use field_registry
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file, json_core, json_value
  use scratch_registry, only : scratch_registry_t
  use user_intf, only : user_t
  use utils, only : neko_warning, neko_error
  use material_properties, only : material_properties_t
  use field_series
  use time_step_controller
  use field_math, only : field_cfill
  implicit none

  !> Base type of all fluid formulations
  type, abstract :: adjoint_scheme_t
     type(field_t), pointer :: u_adj => null() !< x-component of Velocity
     type(field_t), pointer :: v_adj => null() !< y-component of Velocity
     type(field_t), pointer :: w_adj => null() !< z-component of Velocity
     type(field_t), pointer :: p_adj => null() !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)
     type(space_t) :: Xh !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh !< Coefficients associated with \f$ X_h \f$
     !> The source term for the momentum equation.
     type(fluid_source_term_t) :: source_term
     !> X-component of the right-hand side.
     type(field_t), pointer :: f_adj_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_adj_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_adj_z => null()
     class(ksp_t), allocatable :: ksp_vel !< Krylov solver for velocity
     class(ksp_t), allocatable :: ksp_prs !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs !< Velocity Preconditioner
     integer :: vel_projection_dim !< Size of the projection space for ksp_vel
     integer :: pr_projection_dim !< Size of the projection space for ksp_pr
     integer :: vel_projection_activ_step !< Steps to activate projection for ksp_vel
     integer :: pr_projection_activ_step !< Steps to activate projection for ksp_pr
     type(no_slip_wall_t) :: bc_wall !< No-slip wall for velocity
     class(bc_t), allocatable :: bc_inflow !< Dirichlet inflow for velocity

     ! Attributes for field dirichlet BCs
     type(field_dirichlet_vector_t) :: user_field_bc_vel !< User-computed Dirichlet velocity condition
     type(field_dirichlet_t) :: user_field_bc_prs !< User-computed Dirichlet pressure condition
     type(dirichlet_t) :: bc_prs !< Dirichlet pressure condition
     type(dong_outflow_t) :: bc_dong !< Dong outflow condition
     type(symmetry_t) :: bc_sym !< Symmetry plane for velocity
     type(bc_list_t) :: bclst_vel !< List of velocity conditions
     type(bc_list_t) :: bclst_prs !< List of pressure conditions
     type(field_t) :: bdry !< Boundary markings
     type(json_file), pointer :: params !< Parameters
     type(mesh_t), pointer :: msh => null() !< Mesh
     type(chkp_t) :: chkp !< Checkpoint
     type(mean_flow_t) :: mean !< Mean flow field
     type(fluid_stats_t) :: stats !< Fluid statistics
     type(mean_sqr_flow_t) :: mean_sqr !< Mean squared flow field
     logical :: forced_flow_rate = .false. !< Is the flow rate forced?
     logical :: freeze = .false. !< Freeze velocity at initial condition?
     !> Dynamic viscosity
     real(kind=rp), pointer :: mu => null()
     !> The variable mu field
     type(field_t) :: mu_field
     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name
     !> Is mu varying in time? Currently only due to LES models.
     logical :: variable_material_properties = .false.
     !> Density
     real(kind=rp), pointer :: rho => null()
     !> The variable density field
     type(field_t) :: rho_field
     type(scratch_registry_t) :: scratch !< Manager for temporary fields
     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)
   contains
     !> Constructor for the base type
     procedure, pass(this) :: adjoint_scheme_init_all
     procedure, pass(this) :: adjoint_scheme_init_common
     generic :: scheme_init => adjoint_scheme_init_all, adjoint_scheme_init_common
     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => adjoint_scheme_free
     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => adjoint_scheme_validate
     !> Apply pressure boundary conditions
     procedure, pass(this) :: bc_apply_vel => adjoint_scheme_bc_apply_vel
     !> Apply velocity boundary conditions
     procedure, pass(this) :: bc_apply_prs => adjoint_scheme_bc_apply_prs
     !> Set the user inflow procedure
     procedure, pass(this) :: set_usr_inflow => adjoint_scheme_set_usr_inflow
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => adjoint_compute_cfl
     !> Constructor
     procedure(adjoint_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor
     procedure(adjoint_scheme_free_intrf), pass(this), deferred :: free
     !> Advance one step in time
     procedure(adjoint_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint
     procedure(adjoint_scheme_restart_intrf), pass(this), deferred :: restart
     procedure, private, pass(this) :: set_bc_type_output => &
          adjoint_scheme_set_bc_type_output
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          adjoint_scheme_update_material_properties
  end type adjoint_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine adjoint_scheme_init_intrf(this, msh, lx, params, user, &
          material_properties)
       import adjoint_scheme_t
       import json_file
       import mesh_t
       import user_t
       import material_properties_t
       class(adjoint_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), intent(in) :: user
       type(material_properties_t), target, intent(inout) :: material_properties
     end subroutine adjoint_scheme_init_intrf
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine adjoint_scheme_free_intrf(this)
       import adjoint_scheme_t
       class(adjoint_scheme_t), intent(inout) :: this
     end subroutine adjoint_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine adjoint_scheme_step_intrf(this, t, tstep, dt, ext_bdf, &
          dt_controller)
       import adjoint_scheme_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(adjoint_scheme_t), target, intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(inout) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine adjoint_scheme_step_intrf
  end interface

  !> Abstract interface to restart a fluid scheme
  abstract interface
     subroutine adjoint_scheme_restart_intrf(this, dtlag, tlag)
       import adjoint_scheme_t
       import rp
       class(adjoint_scheme_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)

     end subroutine adjoint_scheme_restart_intrf
  end interface

contains

  !> Initialize common data for the current scheme
  subroutine adjoint_scheme_init_common(this, msh, lx, params, scheme, user, &
       material_properties, kspv_init)
    implicit none
    class(adjoint_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(material_properties_t), target, intent(inout) :: material_properties
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

    this%dm_Xh = dofmap_t(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Local scratch registry
    this%scratch = scratch_registry_t(this%dm_Xh, 10, 2)

    ! Case parameters
    this%params => params

    !
    ! Material properties
    !

    this%rho => material_properties%rho
    this%mu => material_properties%mu

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
         'case.fluid.velocity_solver.projection_space_size',&
         this%vel_projection_dim, 20)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_space_size',&
         this%pr_projection_dim, 20)
    call json_get_or_default(params, &
         'case.fluid.velocity_solver.projection_hold_steps',&
         this%vel_projection_activ_step, 5)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_hold_steps',&
         this%pr_projection_activ_step, 5)


    call json_get_or_default(params, 'case.fluid.freeze', this%freeze, .false.)

    if (params%valid_path("case.fluid.flow_rate_force")) then
       this%forced_flow_rate = .true.
    end if


    !
    ! First section of fluid log
    !

    call neko_log%section('Adjoint scheme')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
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

    write(log_buf, '(A,ES13.6)') 'rho        :', this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'mu         :', this%mu
    call neko_log%message(log_buf)

    call json_get(params, 'case.numerics.dealias', logical_val)
    write(log_buf, '(A, L1)') 'Dealias    : ', logical_val
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.output_boundary', logical_val, &
         .false.)
    write(log_buf, '(A, L1)') 'Save bdry  : ', logical_val
    call neko_log%message(log_buf)


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
    call this%bc_sym%mark_zones_from_list(msh%labeled_zones,&
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
       call this%bc_inflow%mark_zones_from_list(msh%labeled_zones,&
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
             call json_get(params, 'case.fluid.blasius.approximation',&
                  string_val2)
             call json_get(params, 'case.fluid.blasius.freestream_velocity',&
                  real_vec)

             call bc_if%set_params(real_vec, real_val, string_val2)

          end select
       else if (trim(string_val1) .eq. "user") then
       end if
    end if

    call this%bc_wall%init_base(this%c_Xh)
    call this%bc_wall%mark_zone(msh%wall)
    call this%bc_wall%mark_zones_from_list(msh%labeled_zones,&
         'w', this%bc_labels)
    call this%bc_wall%finalize()
    call bc_list_add(this%bclst_vel, this%bc_wall)

    ! Setup field dirichlet bc for u-velocity
    call this%user_field_bc_vel%bc_u%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_u%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_u', this%bc_labels)
    call this%user_field_bc_vel%bc_u%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_u%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0) then
       call this%user_field_bc_vel%bc_u%init_field('d_vel_u')
    end if

    ! Setup field dirichlet bc for v-velocity
    call this%user_field_bc_vel%bc_v%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_v%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_v', this%bc_labels)
    call this%user_field_bc_vel%bc_v%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_v%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0) then
       call this%user_field_bc_vel%bc_v%init_field('d_vel_v')
    end if

    ! Setup field dirichlet bc for w-velocity
    call this%user_field_bc_vel%bc_w%init_base(this%c_Xh)
    call this%user_field_bc_vel%bc_w%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_w', this%bc_labels)
    call this%user_field_bc_vel%bc_w%finalize()

    call MPI_Allreduce(this%user_field_bc_vel%bc_w%msk(0), integer_val, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    if (integer_val .gt. 0) then
       call this%user_field_bc_vel%bc_w%init_field('d_vel_w')
    end if

    ! Setup our global field dirichlet bc
    call this%user_field_bc_vel%init_base(this%c_Xh)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_u', this%bc_labels)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones,&
         'd_vel_v', this%bc_labels)
    call this%user_field_bc_vel%mark_zones_from_list(msh%labeled_zones,&
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
    call bc_list_add(this%user_field_bc_vel%bc_list,&
         this%user_field_bc_vel%bc_v)
    call bc_list_add(this%user_field_bc_vel%bc_list,&
         this%user_field_bc_vel%bc_w)

    !
    ! Check if we need to output boundary types to a separate field
    call adjoint_scheme_set_bc_type_output(this, params)

    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_adj_x)
    allocate(this%f_adj_y)
    allocate(this%f_adj_z)
    call this%f_adj_x%init(this%dm_Xh, fld_name = "adjoint_rhs_x")
    call this%f_adj_y%init(this%dm_Xh, fld_name = "adjoint_rhs_y")
    call this%f_adj_z%init(this%dm_Xh, fld_name = "adjoint_rhs_z")

    ! Initialize the source term
    call this%source_term%init(params, this%f_adj_x, this%f_adj_y, this%f_adj_z, this%c_Xh,&
         user)

    ! If case.output_boundary is true, set the values for the bc types in the
    ! output of the field.
    call this%set_bc_type_output(params)

    ! Initialize velocity solver
    if (kspv_init) then
       call neko_log%section("Velocity solver")
       call json_get_or_default(params, &
            'case.fluid.velocity_solver.max_iterations', &
            integer_val, 800)
       call json_get(params, 'case.fluid.velocity_solver.type', string_val1)
       call json_get(params, 'case.fluid.velocity_solver.preconditioner', &
            string_val2)
       call json_get(params, 'case.fluid.velocity_solver.absolute_tolerance', &
            real_val)

       call neko_log%message('Type       : ('// trim(string_val1) // &
            ', ' // trim(string_val2) // ')')

       write(log_buf, '(A,ES13.6)') 'Abs tol    :', real_val
       call neko_log%message(log_buf)
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            string_val1, integer_val, real_val)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, string_val2)
       call neko_log%end_section()
    end if

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, 'u_adj')
    call neko_field_registry%add_field(this%dm_Xh, 'v_adj')
    call neko_field_registry%add_field(this%dm_Xh, 'w_adj')
    this%u_adj => neko_field_registry%get_field('u_adj')
    this%v_adj => neko_field_registry%get_field('v_adj')
    this%w_adj => neko_field_registry%get_field('w_adj')

    !! Initialize time-lag fields
    call this%ulag%init(this%u_adj, 2)
    call this%vlag%init(this%v_adj, 2)
    call this%wlag%init(this%w_adj, 2)


  end subroutine adjoint_scheme_init_common

  !> Initialize all components of the current scheme
  subroutine adjoint_scheme_init_all(this, msh, lx, params, kspv_init, kspp_init,&
       scheme, user, material_properties)
    implicit none
    class(adjoint_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(material_properties_t), target, intent(inout) :: material_properties
    logical :: kspv_init
    logical :: kspp_init
    character(len=*), intent(in) :: scheme
    real(kind=rp) :: abs_tol
    integer :: integer_val, ierr
    character(len=:), allocatable :: solver_type, precon_type
    character(len=LOG_SIZE) :: log_buf

    call adjoint_scheme_init_common(this, msh, lx, params, scheme, user, &
         material_properties, kspv_init)

    call neko_field_registry%add_field(this%dm_Xh, 'p_adj')
    this%p_adj => neko_field_registry%get_field('p_adj')

    !
    ! Setup pressure boundary conditions
    !
    call bc_list_init(this%bclst_prs)
    call this%bc_prs%init_base(this%c_Xh)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
         'o', this%bc_labels)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
         'on', this%bc_labels)

    ! Field dirichlet pressure bc
    call this%user_field_bc_prs%init_base(this%c_Xh)
    call this%user_field_bc_prs%mark_zones_from_list(msh%labeled_zones,&
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
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
         'o+dong', this%bc_labels)
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
         'on+dong', this%bc_labels)
    call this%bc_dong%finalize()


    call this%bc_dong%init(this%c_Xh, params)

    call bc_list_add(this%bclst_prs, this%bc_dong)

    ! Pressure solver
    if (kspp_init) then
       call neko_log%section("Pressure solver")

       call json_get_or_default(params, &
            'case.fluid.pressure_solver.max_iterations', &
            integer_val, 800)
       call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
       call json_get(params, 'case.fluid.pressure_solver.preconditioner', &
            precon_type)
       call json_get(params, 'case.fluid.pressure_solver.absolute_tolerance', &
            abs_tol)
       call neko_log%message('Type       : ('// trim(solver_type) // &
            ', ' // trim(precon_type) // ')')
       write(log_buf, '(A,ES13.6)') 'Abs tol    :', abs_tol
       call neko_log%message(log_buf)

       call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Xh%size(), &
            solver_type, integer_val, abs_tol)
       call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs, precon_type)

       call neko_log%end_section()

    end if


    call neko_log%end_section()

  end subroutine adjoint_scheme_init_all

  !> Deallocate a fluid formulation
  subroutine adjoint_scheme_free(this)
    class(adjoint_scheme_t), intent(inout) :: this

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

    nullify(this%u_adj)
    nullify(this%v_adj)
    nullify(this%w_adj)
    nullify(this%p_adj)

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()


    if (associated(this%f_adj_x)) then
       call this%f_adj_x%free()
    end if

    if (associated(this%f_adj_y)) then
       call this%f_adj_y%free()
    end if

    if (associated(this%f_adj_z)) then
       call this%f_adj_z%free()
    end if

    nullify(this%f_adj_x)
    nullify(this%f_adj_y)
    nullify(this%f_adj_z)

    call this%mu_field%free()


  end subroutine adjoint_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine adjoint_scheme_validate(this)
    class(adjoint_scheme_t), target, intent(inout) :: this
    ! Variables for retrieving json parameters
    logical :: logical_val

    if ( (.not. associated(this%u_adj)) .or. &
         (.not. associated(this%v_adj)) .or. &
         (.not. associated(this%w_adj)) .or. &
         (.not. associated(this%p_adj))) then
       call neko_error('Fields are not registered')
    end if

    if ( (.not. allocated(this%u_adj%x)) .or. &
         (.not. allocated(this%v_adj%x)) .or. &
         (.not. allocated(this%w_adj%x)) .or. &
         (.not. allocated(this%p_adj%x))) then
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
    call this%chkp%init(this%u_adj, this%v_adj, this%w_adj, this%p_adj)

    !
    ! Setup mean flow fields if requested
    !
    if (this%params%valid_path('case.statistics')) then
       call json_get_or_default(this%params, 'case.statistics.enabled',&
            logical_val, .true.)
       if (logical_val) then
          call this%mean%init(this%u_adj, this%v_adj, this%w_adj, this%p_adj)
          call this%stats%init(this%c_Xh, this%mean%u, &
               this%mean%v, this%mean%w, this%mean%p)
       end if
    end if

  end subroutine adjoint_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! Here we perform additional gs operations to take care of
  !! shared points between elements that have different BCs, as done in Nek5000.
  !! @todo Why can't we call the interface here?
  subroutine adjoint_scheme_bc_apply_vel(this, t, tstep)
    class(adjoint_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call bc_list_apply_vector(this%bclst_vel,&
         this%u_adj%x, this%v_adj%x, this%w_adj%x, this%dm_Xh%size(), t, tstep)

  end subroutine adjoint_scheme_bc_apply_vel

  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine adjoint_scheme_bc_apply_prs(this, t, tstep)
    class(adjoint_scheme_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call bc_list_apply_scalar(this%bclst_prs, this%p_adj%x, &
         this%p_adj%dof%size(), t, tstep)

  end subroutine adjoint_scheme_bc_apply_prs

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver, max_iter, abstol)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*), intent(in) :: solver
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: abstol

    call krylov_solver_factory(ksp, n, solver, max_iter, abstol)

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
  subroutine adjoint_scheme_set_usr_inflow(this, usr_eval)
    class(adjoint_scheme_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval

    select type (bc_if => this%bc_inflow)
      type is (usr_inflow_t)
       call bc_if%set_eval(usr_eval)
      class default
       call neko_error("Not a user defined inflow condition")
    end select
  end subroutine adjoint_scheme_set_usr_inflow

  !> Compute CFL
  function adjoint_compute_cfl(this, dt) result(c)
    class(adjoint_scheme_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u_adj%x, this%v_adj%x, this%w_adj%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)

  end function adjoint_compute_cfl

  !> Set boundary types for the diagnostic output.
  !! @param params The JSON case file.
  subroutine adjoint_scheme_set_bc_type_output(this, params)
    class(adjoint_scheme_t), intent(inout) :: this
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
       call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
            'w', this%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(1.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init_base(this%c_Xh)
       call bdry_mask%mark_zone(this%msh%inlet)
       call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
            'v', this%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(2.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init_base(this%c_Xh)
       call bdry_mask%mark_zone(this%msh%outlet)
       call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
            'o', this%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(3.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init_base(this%c_Xh)
       call bdry_mask%mark_zone(this%msh%sympln)
       call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
            'sym', this%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(4.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init_base(this%c_Xh)
       call bdry_mask%mark_zone(this%msh%outlet_normal)
       call bdry_mask%mark_zones_from_list(this%msh%labeled_zones,&
            'on', this%bc_labels)
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
    end if

  end subroutine adjoint_scheme_set_bc_type_output

  !> Update the values of `mu_field` if necessary.
  subroutine adjoint_scheme_update_material_properties(this)
    class(adjoint_scheme_t), intent(inout) :: this
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

  end subroutine adjoint_scheme_update_material_properties

end module adjoint_scheme
