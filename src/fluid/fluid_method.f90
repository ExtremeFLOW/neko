! Copyright (c) 2020-2022, The Neko Authors
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
module fluid_method
  use gather_scatter
  use mean_sqr_flow    
  use neko_config
  use parameters
  use checkpoint
  use mean_flow
  use num_types
  use source
  use field
  use space
  use dofmap
  use krylov
  use coefs
  use wall
  use inflow
  use usr_inflow
  use blasius
  use dirichlet
  use dong_outflow
  use symmetry
  use non_normal
  use krylov_fctry
  use precon_fctry
  use bc
  use mesh
  use math
  use abbdf
  use mathops
  use operators
  use logger
  use field_registry
  implicit none
  
  !> Base type of all fluid formulations
  type, abstract :: fluid_scheme_t
     type(field_t), pointer :: u => null() !< x-component of Velocity
     type(field_t), pointer :: v => null() !< y-component of Velocity
     type(field_t), pointer :: w => null() !< z-component of Velocity
     type(field_t), pointer :: p => null() !< Pressure
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$
     type(source_t) :: f_Xh     !< Source term associated with \f$ X_h \f$
     class(ksp_t), allocatable  :: ksp_vel     !< Krylov solver for velocity
     class(ksp_t), allocatable  :: ksp_prs     !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel        !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs        !< Velocity Preconditioner
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     class(inflow_t), allocatable :: bc_inflow !< Dirichlet inflow for velocity
     type(dirichlet_t) :: bc_prs               !< Dirichlet pressure condition
     type(dong_outflow_t) :: bc_dong           !< Dong outflow condition
     type(symmetry_t) :: bc_sym                !< Symmetry plane for velocity
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_prs              !< List of pressure conditions
     type(field_t) :: bdry                     !< Boundary markings     
     type(param_t), pointer :: params          !< Parameters          
     type(mesh_t), pointer :: msh => null()    !< Mesh
     type(chkp_t) :: chkp                      !< Checkpoint
     type(mean_flow_t) :: mean                 !< Mean flow field
     type(mean_sqr_flow_t) :: mean_sqr         !< Mean squared flow field
   contains
     procedure, pass(this) :: fluid_scheme_init_all
     procedure, pass(this) :: fluid_scheme_init_uvw
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     procedure, pass(this) :: validate => fluid_scheme_validate
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     procedure, pass(this) :: set_source => fluid_scheme_set_source
     procedure, pass(this) :: set_usr_inflow => fluid_scheme_set_usr_inflow
     procedure, pass(this) :: compute_cfl => fluid_compute_cfl
     procedure(fluid_method_init), pass(this), deferred :: init
     procedure(fluid_method_free), pass(this), deferred :: free
     procedure(fluid_method_step), pass(this), deferred :: step
     generic :: scheme_init => fluid_scheme_init_all, fluid_scheme_init_uvw
  end type fluid_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_method_init(this, msh, lx, param)
       import fluid_scheme_t
       import param_t
       import mesh_t
       class(fluid_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh       
       integer, intent(inout) :: lx
       type(param_t), target, intent(inout) :: param              
     end subroutine fluid_method_init
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine fluid_method_free(this)
       import fluid_scheme_t
       class(fluid_scheme_t), intent(inout) :: this
     end subroutine fluid_method_free
  end interface
  
  !> Abstract interface to compute a time-step
  abstract interface
     subroutine fluid_method_step(this, t, tstep, ab_bdf)
       import fluid_scheme_t
       import abbdf_t
       import rp
       class(fluid_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       type(abbdf_t), intent(inout) :: ab_bdf
     end subroutine fluid_method_step
  end interface

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_common(this, msh, lx, params, scheme)
    class(fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(param_t), target, intent(inout) :: params
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    
    call neko_log%section('Fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'lx         : ', lx
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'lx         : ', lx
    else
       write(log_buf, '(A, I3)') 'lx         : ', lx
    end if
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'Re         :',  params%Re
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'rho        :',  params%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'mu         :',  params%mu
    call neko_log%message(log_buf)
    call neko_log%message('Ksp vel.   : ('// trim(params%ksp_vel) // &
         ', ' // trim(params%pc_vel) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :',  params%abstol_vel
    call neko_log%message(log_buf)
    call neko_log%message('Ksp prs.   : ('// trim(params%ksp_prs) // &
         ', ' // trim(params%pc_prs) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :',  params%abstol_prs
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'Dealias    : ',  params%dealias
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'Save bdry  : ',  params%output_bdry
    call neko_log%message(log_buf)


    if (msh%gdim .eq. 2) then
       call space_init(this%Xh, GLL, lx, lx)
    else
       call space_init(this%Xh, GLL, lx, lx, lx)
    end if

    this%dm_Xh = dofmap_t(msh, this%Xh)

    this%params => params

    this%msh => msh

    call gs_init(this%gs_Xh, this%dm_Xh)

    call coef_init(this%c_Xh, this%gs_Xh)

    call source_init(this%f_Xh, this%dm_Xh)

    !
    ! Setup velocity boundary conditions
    !
    call bc_list_init(this%bclst_vel)

    call this%bc_sym%init(this%dm_Xh)
    call this%bc_sym%mark_zone(msh%sympln)
    call this%bc_sym%mark_zones_from_list(msh%labeled_zones,&
                        'sym', this%params%bc_labels)
    call this%bc_sym%finalize()
    call this%bc_sym%init_msk(this%c_Xh)    
    call bc_list_add(this%bclst_vel, this%bc_sym)

    if (trim(params%fluid_inflow) .eq. "default") then
       allocate(inflow_t::this%bc_inflow)
    else if (trim(params%fluid_inflow) .eq. "blasius") then
       allocate(blasius_t::this%bc_inflow)
    else if (trim(params%fluid_inflow) .eq. "user") then
       allocate(usr_inflow_t::this%bc_inflow)
    else
       call neko_error('Invalid Inflow condition')
    end if
    
    call this%bc_inflow%init(this%dm_Xh)
    call this%bc_inflow%mark_zone(msh%inlet)
    call this%bc_inflow%mark_zones_from_list(msh%labeled_zones,&
                        'v', this%params%bc_labels)
    call this%bc_inflow%finalize()
    call this%bc_inflow%set_inflow(params%uinf)
    call bc_list_add(this%bclst_vel, this%bc_inflow)

    if (trim(params%fluid_inflow) .eq. "blasius") then
       select type(bc_if => this%bc_inflow)
       type is(blasius_t)
          call bc_if%set_coef(this%C_Xh)
          call bc_if%set_params(params%delta, params%blasius_approx)
       end select
    end if
    
    if (trim(params%fluid_inflow) .eq. "user") then
       select type(bc_if => this%bc_inflow)
       type is(usr_inflow_t)
          call bc_if%set_coef(this%C_Xh)
       end select
    end if
    
    call this%bc_wall%init(this%dm_Xh)
    call this%bc_wall%mark_zone(msh%wall)
    call this%bc_wall%mark_zones_from_list(msh%labeled_zones,&
                        'w', this%params%bc_labels)
    call this%bc_wall%finalize()
    call bc_list_add(this%bclst_vel, this%bc_wall)
       
    if (params%output_bdry) then       
       call field_init(this%bdry, this%dm_Xh, 'bdry')
       this%bdry = 0.0_rp
       
       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%wall)
       call bdry_mask%mark_zones_from_list(msh%labeled_zones,&
                      'w', this%params%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(1.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%inlet)
       call bdry_mask%mark_zones_from_list(msh%labeled_zones,&
                      'v', this%params%bc_labels)

       call bdry_mask%finalize()
       call bdry_mask%set_g(2.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%outlet)
       call bdry_mask%mark_zones_from_list(msh%labeled_zones,&
                      'o', this%params%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(3.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%sympln)
       call bdry_mask%mark_zones_from_list(msh%labeled_zones,&
                      'sym', this%params%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(4.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%periodic)
       call bdry_mask%finalize()
       call bdry_mask%set_g(5.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%outlet_normal)
       call bdry_mask%mark_zones_from_list(msh%labeled_zones,&
                      'on', this%params%bc_labels)
       call bdry_mask%finalize()
       call bdry_mask%set_g(6.0_rp)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%size())
       call bdry_mask%free()

    end if
    
  end subroutine fluid_scheme_init_common

  !> Initialize all velocity related components of the current scheme
  subroutine fluid_scheme_init_uvw(this, msh, lx, params, kspv_init, scheme)
    class(fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), target, intent(inout) :: params
    logical :: kspv_init
    character(len=*), intent(in) :: scheme

    call fluid_scheme_init_common(this, msh, lx, params, scheme)
    
    call neko_field_registry%add_field(this%dm_Xh, 'u')
    call neko_field_registry%add_field(this%dm_Xh, 'v')
    call neko_field_registry%add_field(this%dm_Xh, 'w')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    if (kspv_init) then
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            params%ksp_vel, params%abstol_vel)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, params%pc_vel)
    end if

    call neko_log%end_section()
  end subroutine fluid_scheme_init_uvw

  !> Initialize all components of the current scheme
  subroutine fluid_scheme_init_all(this, msh, lx, params, &
                                   kspv_init, kspp_init, scheme)
    class(fluid_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), target, intent(inout) :: params
    logical :: kspv_init
    logical :: kspp_init
    character(len=*), intent(in) :: scheme

    call fluid_scheme_init_common(this, msh, lx, params, scheme)

    call neko_field_registry%add_field(this%dm_Xh, 'u')
    call neko_field_registry%add_field(this%dm_Xh, 'v')
    call neko_field_registry%add_field(this%dm_Xh, 'w')
    call neko_field_registry%add_field(this%dm_Xh, 'p')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')
    this%p => neko_field_registry%get_field('p')

    !
    ! Setup pressure boundary conditions
    !
    call bc_list_init(this%bclst_prs)
    call this%bc_prs%init(this%dm_Xh)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
                        'o', this%params%bc_labels)
    call this%bc_prs%mark_zones_from_list(msh%labeled_zones,&
                        'on', this%params%bc_labels)

    if (msh%outlet%size .gt. 0) then
       call this%bc_prs%mark_zone(msh%outlet)
    end if
    if (msh%outlet_normal%size .gt. 0) then
       call this%bc_prs%mark_zone(msh%outlet_normal)
    end if

    call this%bc_prs%finalize()
    call this%bc_prs%set_g(0.0_rp)
    call bc_list_add(this%bclst_prs, this%bc_prs)
    call this%bc_dong%init(this%dm_Xh)
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
                        'o+dong', this%params%bc_labels)
    call this%bc_dong%mark_zones_from_list(msh%labeled_zones,&
                        'on+dong', this%params%bc_labels)
    call this%bc_dong%finalize()
    call this%bc_dong%set_vars(this%c_Xh, this%u, this%v, this%w,&
         params%dong_uchar, params%dong_delta)

    call bc_list_add(this%bclst_prs, this%bc_dong)

    if (kspv_init) then
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            params%ksp_vel, params%abstol_vel)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, params%pc_vel)
    end if

    if (kspp_init) then
       call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Xh%size(), &
            params%ksp_prs, params%abstol_prs)
       call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs, params%pc_prs)
    end if


    call neko_log%end_section()
    
  end subroutine fluid_scheme_init_all

  !> Deallocate a fluid formulation
  subroutine fluid_scheme_free(this)
    class(fluid_scheme_t), intent(inout) :: this

    call field_free(this%bdry)

    if (allocated(this%bc_inflow)) then
       call this%bc_inflow%free()
    end if

    call this%bc_wall%free()
    call this%bc_sym%free()

    call space_free(this%Xh)    

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

    call gs_free(this%gs_Xh)

    call coef_free(this%c_Xh)

    call source_free(this%f_Xh)

    call bc_list_free(this%bclst_vel)

    nullify(this%params)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    
    
  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_t), target, intent(inout) :: this

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

    if (.not. associated(this%f_Xh%eval)) then
       call neko_error('No source term defined')
    end if

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

    select type(ip => this%bc_inflow)
    type is(usr_inflow_t)
       call ip%validate
    end select

    !
    ! Setup checkpoint structure (if everything is fine)
    !
    call this%chkp%init(this%u, this%v, this%w, this%p)

    !
    ! Setup mean flow fields if requested
    !
    if (this%params%stats_mean_flow) then
       call this%mean%init(this%u, this%v, this%w, this%p)
    end if

    if (this%params%stats_mean_sqr_flow) then
       call this%mean_sqr%init(this%u, this%v, this%w, this%p)
    end if

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_vector(this%bclst_vel,&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size())
  end subroutine fluid_scheme_bc_apply_vel
  
  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_scalar(this%bclst_prs, this%p%x, this%p%dof%size())
  end subroutine fluid_scheme_bc_apply_prs
  
  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver, abstol)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=20), intent(inout) :: solver
    real(kind=rp) :: abstol

    call krylov_solver_factory(ksp, n, solver, abstol)
    
  end subroutine fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine fluid_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=20) :: pctype
    
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
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, &
                  bclst, trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)
    
  end subroutine fluid_scheme_precon_factory

  !> Initialize source term
  subroutine fluid_scheme_set_source(this, source_term_type, usr_f, usr_f_vec)
    class(fluid_scheme_t), intent(inout) :: this
    character(len=*) :: source_term_type
    procedure(source_term_pw), optional :: usr_f
    procedure(source_term), optional :: usr_f_vec

    if (trim(source_term_type) .eq. 'noforce') then
       call source_set_type(this%f_Xh, source_eval_noforce)
    else if (trim(source_term_type) .eq. 'user' .and. present(usr_f)) then
       call source_set_pw_type(this%f_Xh, usr_f)
    else if (trim(source_term_type) .eq. 'user_vector' .and. present(usr_f_vec)) then
       call source_set_type(this%f_Xh, usr_f_vec)
    else if (trim(source_term_type) .eq. '') then
       if (pe_rank .eq. 0) then
          call neko_warning('No source term defined, using default (noforce)')
       end if
       call source_set_type(this%f_Xh, source_eval_noforce)
    else
       call neko_error('Invalid source term')
    end if

  end subroutine fluid_scheme_set_source

  !> Initialize a user defined inflow condition
  subroutine fluid_scheme_set_usr_inflow(this, usr_eval)
    class(fluid_scheme_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval

    select type(bc_if => this%bc_inflow)
    type is(usr_inflow_t)
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
     
end module fluid_method
