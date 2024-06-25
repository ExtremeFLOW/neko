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
!> Contains the scalar_scheme_t type.

! todo: module name
module scalar_scheme
  use gather_scatter, only : gs_t
  use checkpoint, only : chkp_t
  use num_types, only: rp
  use field, only : field_t
  use field_list, only: field_list_t
  use space, only : space_t
  use dofmap, only :  dofmap_t
  use krylov, only : ksp_t
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use krylov_fctry, only : krylov_solver_factory, krylov_solver_destroy
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use hsmg, only : hsmg_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use precon_fctry, only : precon_factory, precon_destroy
  use precon, only : pc_t
  use field_dirichlet, only: field_dirichlet_t, field_dirichlet_update
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use facet_zone, only : facet_zone_t
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, LOG_SIZE
  use field_registry, only : neko_field_registry
  use usr_scalar, only : usr_scalar_t, usr_scalar_bc_eval
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file, json_core, json_value
  use user_intf, only : user_t
  use material_properties, only : material_properties_t
  use utils, only : neko_error
  use comm, only: NEKO_COMM, MPI_INTEGER, MPI_SUM
  use scalar_source_term, only : scalar_source_term_t
  use field_series, only : field_series_t
  use time_step_controller, only : time_step_controller_t
  use bc_fctry, only : bc_factory
  implicit none

  !> Base type for a scalar advection-diffusion solver.
  type, abstract :: scalar_scheme_t
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
     type(coef_t), pointer  :: c_Xh
     !> Right-hand side.
     type(field_t), pointer :: f_Xh => null()
     !> The source term for equation.
     type(scalar_source_term_t) :: source_term
     !> Krylov solver.
     class(ksp_t), allocatable  :: ksp
     !> Max iterations in the Krylov solver.
     integer :: ksp_maxiter
     !> Projection space size.
     integer :: projection_dim
     !< Steps to activate projection for ksp
     integer :: projection_activ_step
     !> Preconditioner.
     class(pc_t), allocatable :: pc
     !> Field Dirichlet conditions.
     type(field_dirichlet_t) :: field_dir_bc
     !> List of BC objects to pass to user_dirichlet_update
     type(bc_list_t) :: field_dirichlet_bcs
     !> Number of strong  bcs.
     integer :: n_strong = 0
     !> List of boundary conditions, including the user one.
     type(bc_list_t) :: bcs
     !> User Dirichlet conditions.
     type(usr_scalar_t) :: user_bc
     !> Case paramters.
     type(json_file), pointer :: params
     !> Mesh.
     type(mesh_t), pointer :: msh => null()
     !> Checkpoint for restarts.
     type(chkp_t) :: chkp
     !> Thermal diffusivity.
     real(kind=rp), pointer :: lambda
     !> Density.
     real(kind=rp), pointer :: rho
     !> Specific heat capacity.
     real(kind=rp), pointer :: cp
   contains
     !> Constructor for the base type.
     procedure, pass(this) :: scheme_init => scalar_scheme_init
     !> Destructor for the base type.
     procedure, pass(this) :: scheme_free => scalar_scheme_free
     !> Validate successful initialization.
     procedure, pass(this) :: validate => scalar_scheme_validate
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
     subroutine scalar_scheme_init_intrf(this, msh, coef, gs, params, user,&
                                         material_properties)
       import scalar_scheme_t
       import json_file
       import coef_t
       import gs_t
       import mesh_t
       import user_t
       import material_properties_t
       class(scalar_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       type(coef_t), target, intent(inout) :: coef
       type(gs_t), target, intent(inout) :: gs
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       type(material_properties_t), intent(inout) :: material_properties
     end subroutine scalar_scheme_init_intrf
  end interface

  !> Abstract interface to restart a scalar formulation
  abstract interface
     subroutine scalar_scheme_restart_intrf(this, dtlag, tlag)
       import scalar_scheme_t
       import chkp_t
       import rp
       class(scalar_scheme_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)
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
     subroutine scalar_scheme_step_intrf(this, t, tstep, dt, ext_bdf, dt_controller)
       import scalar_scheme_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(scalar_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(inout) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine scalar_scheme_step_intrf
  end interface

contains

  !> Initialize boundary conditions
  !! @param user The user object binding the user-defined routines.
  subroutine scalar_scheme_setup_bcs(this, user)
    class(scalar_scheme_t), intent(inout) :: this
    type(user_t), target, intent(in) :: user
    integer :: i, j, n_bcs, ierr
    real(kind=rp) :: dir_value, flux_value
    logical :: bc_exists
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found

    call this%bcs%init()

    if (this%params%valid_path('case.scalar.boundary_conditions')) then
       call this%params%info('case.scalar.boundary_conditions', n_children=n_bcs)
       call this%params%get_core(core)
       call this%params%get('case.scalar.boundary_conditions', bc_object, found)
       call this%bcs%init(n_bcs)

       ! Gotta set this manually, becase we don't append to the list but
       ! rather directly allocate in the factory.
       this%bcs%size_ = n_bcs


       do i=1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)

          call bc_factory(this%bcs%items(i)%ptr, bc_subdict, &
                          this%c_Xh, user)

          if (this%bcs%strong(i)) then
             this%n_strong = this%n_strong + 1
          end if
       end do
    end if
  end subroutine scalar_scheme_setup_bcs

  !> Initialize all related components of the current scheme
  !! @param msh The mesh.
  !! @param c_Xh The coefficients.
  !! @param gs_Xh The gather-scatter.
  !! @param params The case parameter file in json.
  !! @param scheme The name of the scalar scheme.
  !! @param user Type with user-defined procedures.
  subroutine scalar_scheme_init(this, msh, c_Xh, gs_Xh, params, scheme, user, &
                                material_properties)
    class(scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    type(coef_t), target, intent(inout) :: c_Xh
    type(gs_t), target, intent(inout) :: gs_Xh
    type(json_file), target, intent(inout) :: params
    character(len=*), intent(in) :: scheme
    type(user_t), target, intent(in) :: user
    type(material_properties_t), target,  intent(inout) :: material_properties
    ! IO buffer for log output
    character(len=LOG_SIZE) :: log_buf
    ! Variables for retrieving json parameters
    logical :: logical_val
    real(kind=rp) :: real_val, solver_abstol
    integer :: integer_val, ierr
    character(len=:), allocatable :: solver_type, solver_precon

    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    call neko_log%section('Scalar')
    call json_get(params, 'case.fluid.velocity_solver.type', solver_type)
    call json_get(params, 'case.fluid.velocity_solver.preconditioner',&
                  solver_precon)
    call json_get(params, 'case.fluid.velocity_solver.absolute_tolerance',&
                  solver_abstol)

    !
    ! Material properties
    !
    this%rho => material_properties%rho
    this%lambda => material_properties%lambda
    this%cp => material_properties%cp

    write(log_buf, '(A,ES13.6)') 'rho        :',  this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'lambda     :',  this%lambda
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'cp         :',  this%cp
    call neko_log%message(log_buf)

    call json_get_or_default(params, &
                            'case.fluid.velocity_solver.projection_space_size',&
                            this%projection_dim, 20)
    call json_get_or_default(params, &
                            'case.fluid.velocity_solver.projection_hold_steps',&
                            this%projection_activ_step, 5)


    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    call neko_log%message('Ksp scalar : ('// trim(solver_type) // &
         ', ' // trim(solver_precon) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :',  solver_abstol
    call neko_log%message(log_buf)

    this%Xh => this%u%Xh
    this%dm_Xh => this%u%dof
    this%params => params
    this%msh => msh
    if (.not. neko_field_registry%field_exists('s')) then
       call neko_field_registry%add_field(this%dm_Xh, 's')
    end if
    this%s => neko_field_registry%get_field('s')

    call this%slag%init(this%s, 2)

    this%gs_Xh => gs_Xh
    this%c_Xh => c_Xh

    !
    ! Setup right-hand side field.
    !
    allocate(this%f_Xh)
    call this%f_Xh%init(this%dm_Xh, fld_name="scalar_rhs")

    ! Initialize the source term
    call this%source_term%init(params, this%f_Xh, this%c_Xh, user)

    ! Set up boundary conditions
    call scalar_scheme_setup_bcs(this, user)


!! COMMENTING USER STUFF
    ! Mark BC zones
!    call this%user_bc%mark_zone(msh%wall)
!    call this%user_bc%mark_zone(msh%inlet)
!    call this%user_bc%mark_zone(msh%outlet)
!    call this%user_bc%mark_zone(msh%outlet_normal)
!    call this%user_bc%mark_zone(msh%sympln)
!    call this%user_bc%finalize()
!    if (this%user_bc%msk(0) .gt. 0) call bc_list_add(this%bclst_dirichlet,&
!                                                     this%user_bc)

    ! Add field dirichlet BCs
!    call this%field_dir_bc%init_base(this%c_Xh)
!    call this%field_dir_bc%mark_zones_from_list(msh%labeled_zones, &
!         'd_s', this%bc_labels)
!    call this%field_dir_bc%finalize()
!    call MPI_Allreduce(this%field_dir_bc%msk(0), integer_val, 1, &
!         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
!    if (integer_val .gt. 0) call this%field_dir_bc%init_field('d_s')

!    call bc_list_add(this%bclst_dirichlet, this%field_dir_bc)

    !
    ! Associate our field dirichlet update to the user one.
    !
!    this%field_dir_bc%update => user%user_dirichlet_update

!    call bc_list_init(this%field_dirichlet_bcs, size=1)
!    call bc_list_add(this%field_dirichlet_bcs, this%field_dir_bc)

    ! todo parameter file ksp tol should be added
    call json_get_or_default(params, 'case.fluid.velocity_solver.max_iterations',&
                             integer_val, 800)
    call scalar_scheme_solver_factory(this%ksp, this%dm_Xh%size(), &
         solver_type, integer_val, solver_abstol)
    call scalar_scheme_precon_factory(this%pc, this%ksp, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs, solver_precon)

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
       call krylov_solver_destroy(this%ksp)
       deallocate(this%ksp)
    end if

    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if

    call this%source_term%free()

    call this%bcs%free()

    call this%slag%free()

    ! Free everything related to field dirichlet BCs
    call this%field_dirichlet_bcs%free()
    call this%field_dir_bc%free()

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

    !
    ! Setup checkpoint structure (if everything is fine)
    !
!    @todo no io for now
!    call this%chkp%init(this%u, this%v, this%w, this%p)

  end subroutine scalar_scheme_validate

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine scalar_scheme_solver_factory(ksp, n, solver, max_iter, abstol)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    integer, intent(in) :: max_iter
    character(len=*), intent(in) :: solver
    real(kind=rp) :: abstol

    call krylov_solver_factory(ksp, n, solver, max_iter, abstol)

  end subroutine scalar_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine scalar_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type(pcp => pc)
    type is(jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is(hsmg_t)
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

  end subroutine scalar_scheme_precon_factory

end module scalar_scheme
