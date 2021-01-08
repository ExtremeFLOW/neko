!> Fluid formulations
module fluid_method
  use gather_scatter
  use parameters
  use num_types
  use source
  use field
  use space
  use dofmap
  use krylov
  use coefs
  use wall
  use inflow
  use dirichlet
  use cg
  use bc
  use jacobi
  use gmres
  use mesh
  use math
  use abbdf
  use mathops
  use operators
  use hsmg
  implicit none
  
  !> Base type of all fluid formulations
  type, abstract :: fluid_scheme_t
     type(field_t) :: u         !< x-component of Velocity
     type(field_t) :: v         !< y-component of Velocity
     type(field_t) :: w         !< z-component of Velocity
     type(field_t) :: p         !< Pressure
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
     type(inflow_t) :: bc_inflow               !< Dirichlet inflow for velocity
     type(dirichlet_t) :: bc_prs               !< Dirichlet pressure condition
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_prs              !< List of pressure conditions
     type(field_t) :: bdry                     !< Boundary markings
     type(param_t), pointer :: params          !< Parameters          
     type(mesh_t), pointer :: msh => null()    !< Mesh
   contains
     procedure, pass(this) :: fluid_scheme_init_all
     procedure, pass(this) :: fluid_scheme_init_uvw
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     procedure, pass(this) :: validate => fluid_scheme_validate
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     procedure(fluid_method_init), pass(this), deferred :: init
     procedure(fluid_method_free), pass(this), deferred :: free
     procedure(fluid_method_step), pass(this), deferred :: step
     generic :: scheme_init => fluid_scheme_init_all, fluid_scheme_init_uvw
  end type fluid_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_method_init(this, msh, lx, param, vel, prs)
       import fluid_scheme_t
       import param_t
       import mesh_t
       class(fluid_scheme_t), intent(inout) :: this
       type(mesh_t), intent(inout) :: msh       
       integer, intent(inout) :: lx
       type(param_t), intent(inout) :: param              
       character(len=80), intent(inout) :: vel
       character(len=80), intent(inout) :: prs
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
       import dp
       class(fluid_scheme_t), intent(inout) :: this
       real(kind=dp), intent(inout) :: t
       integer, intent(inout) :: tstep
       type(abbdf_t), intent(inout) :: ab_bdf
     end subroutine fluid_method_step
  end interface

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_common(this, msh, lx, params)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout), target :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout), target :: params
    type(dirichlet_t) :: bdry_mask
    
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
    call this%bc_wall%init(this%dm_Xh)
    call this%bc_wall%mark_zone(msh%wall)
    call this%bc_wall%finalize()

    call this%bc_inflow%init(this%dm_Xh)
    call this%bc_inflow%mark_zone(msh%inlet)
    call this%bc_inflow%finalize()
    call this%bc_inflow%set_inflow(params%uinf)

    call bc_list_init(this%bclst_vel)
    call bc_list_add(this%bclst_vel, this%bc_inflow)
    call bc_list_add(this%bclst_vel, this%bc_wall)
    
    if (params%output_bdry) then

       if (pe_rank .eq. 0) then
          write(*,*) 'Saving boundary markings'
       end if
       
       call field_init(this%bdry, this%dm_Xh, 'bdry')
       this%bdry = 0d0
       
       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%wall)
       call bdry_mask%finalize()
       call bdry_mask%set_g(1d0)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%inlet)
       call bdry_mask%finalize()
       call bdry_mask%set_g(2d0)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%outlet)
       call bdry_mask%finalize()
       call bdry_mask%set_g(3d0)
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()
    end if

  end subroutine fluid_scheme_init_common

  !> Initialize all velocity related components of the current scheme
  subroutine fluid_scheme_init_uvw(this, msh, lx, params, solver_vel)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: params
    character(len=80), intent(inout) :: solver_vel

    call fluid_scheme_init_common(this, msh, lx, params)
    
    call field_init(this%u, this%dm_Xh, 'u')
    call field_init(this%v, this%dm_Xh, 'v')
    call field_init(this%w, this%dm_Xh, 'w')

    call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), solver_vel)

  end subroutine fluid_scheme_init_uvw

  !> Initialize all components of the current scheme
  subroutine fluid_scheme_init_all(this, msh, lx, params, solver_vel, solver_prs)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: params      
    character(len=80), intent(inout) :: solver_vel
    character(len=80), intent(inout) :: solver_prs

    call fluid_scheme_init_common(this, msh, lx, params)
    
    call field_init(this%u, this%dm_Xh, 'u')
    call field_init(this%v, this%dm_Xh, 'v')
    call field_init(this%w, this%dm_Xh, 'w')
    call field_init(this%p, this%dm_Xh, 'p')

    !
    ! Setup pressure boundary conditions
    !
    call this%bc_prs%init(this%dm_Xh)
    call this%bc_prs%mark_zone(msh%outlet)
    call this%bc_prs%finalize()
    call this%bc_prs%set_g(0d0)
    call bc_list_init(this%bclst_prs)
    call bc_list_add(this%bclst_prs, this%bc_prs)

    call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), solver_vel)
    call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
         this%c_Xh, this%dm_Xh, this%gs_Xh)

    call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Xh%size(), solver_prs)
    !call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
    !     this%c_Xh, this%dm_Xh, this%gs_Xh)
    call fluid_scheme_hsmg_factory(this%pc_prs, this%ksp_prs, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs)

  end subroutine fluid_scheme_init_all

  !> Deallocate a fluid formulation
  subroutine fluid_scheme_free(this)
    class(fluid_scheme_t), intent(inout) :: this

    call field_free(this%u)
    call field_free(this%v)
    call field_free(this%w)
    call field_free(this%p)
    call field_free(this%bdry)

    call this%bc_inflow%free()
    call this%bc_wall%free()

    call space_free(this%Xh)

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if

    call gs_free(this%gs_Xh)

    call coef_free(this%c_Xh)

    call source_free(this%f_Xh)

    call bc_list_free(this%bclst_vel)

    nullify(this%params)
    
  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_t), intent(inout) :: this

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

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_vector(this%bclst_vel,&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%n_dofs)
  end subroutine fluid_scheme_bc_apply_vel
  
  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_scalar(this%bclst_prs, this%p%x, this%p%dof%n_dofs)
  end subroutine fluid_scheme_bc_apply_prs
  
  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver)
    class(ksp_t), allocatable, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=80), intent(inout) :: solver

    if (trim(solver) .eq. 'cg') then
       allocate(cg_t::ksp)
    else if (trim(solver) .eq. 'gmres') then
       allocate(gmres_t::ksp)
    else
       call neko_error('Unknown linear solver')
    end if

    select type(kp => ksp)
    type is(cg_t)
       call kp%init(n)
    type is(gmres_t)
       call kp%init(n)
    end select
    
  end subroutine fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  !! @note Currently hardcoded to jacobi
  subroutine fluid_scheme_precon_factory(pc, ksp, coef, dof, gs)
    class(pc_t), allocatable, intent(inout), target :: pc
    class(ksp_t), allocatable, intent(inout) :: ksp
    type(coef_t), intent(inout) :: coef
    type(dofmap_t), intent(inout) :: dof
    type(gs_t), intent(inout) :: gs

    allocate(jacobi_t::pc)

    select type(pcp => pc)
    type is(jacobi_t)
       call pcp%init(coef, dof, gs)
    end select

    ksp%M => pc
    
  end subroutine fluid_scheme_precon_factory
 
  !> Initialize a Krylov preconditioner
  !! @note Currently hardcoded to hsmg
  subroutine fluid_scheme_hsmg_factory(pc, ksp, coef, dof, gs, bclst)
    class(pc_t), allocatable, intent(inout), target :: pc
    class(ksp_t), allocatable, intent(inout) :: ksp
    type(coef_t), intent(inout) :: coef
    type(dofmap_t), intent(inout) :: dof
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bclst

    allocate(hsmg_t::pc)

    select type(pcp => pc)
    type is(hsmg_t)
        call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
    end select

    ksp%M => pc
    
  end subroutine fluid_scheme_hsmg_factory
    
end module fluid_method
