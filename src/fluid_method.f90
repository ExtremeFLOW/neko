!> Fluid formulations
module fluid_method
  use field
  use space
  use dofmap
  use krylov
  use mesh
  implicit none
  
  !> Base type of all fluid formulations
  type, abstract :: fluid_scheme_t
     type(field_t) :: u         !< x-component of Velocity
     type(field_t) :: v         !< y-component of Velocity
     type(field_t) :: w         !< z-component of Velocity
     type(field_t) :: p         !< Pressure
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dof      !< Dofmap assoicated with \f$ X_h \f$
     class(ksp_t), allocatable  :: ksp_vel     !< Krylov solver for velocity
     class(ksp_t), allocatable  :: ksp_prs     !< Krylov solver for pressure
   contains
     procedure, pass(this) :: scheme_init => fluid_scheme_init
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     procedure(fluid_method_init), pass(this), deferred :: init
     procedure(fluid_method_free), pass(this), deferred :: free
     procedure(fluid_method_step), pass(this), deferred :: step
  end type fluid_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_method_init(this, msh, lx)
       import fluid_scheme_t
       import mesh_t
       class(fluid_scheme_t), intent(inout) :: this
       type(mesh_t), intent(inout) :: msh
       integer, intent(inout) :: lx
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
     subroutine fluid_method_step(this)
       import fluid_scheme_t
       class(fluid_scheme_t), intent(inout) :: this
     end subroutine fluid_method_step
  end interface

contains

  subroutine fluid_scheme_init(this, msh, lx)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx

    if (msh%gdim .eq. 2) then
       call space_init(this%Xh, GLL, lx, lx)
    else
       call space_init(this%Xh, GLL, lx, lx, lx)
    end if

    this%dof = dofmap_t(msh, this%Xh)

    call field_init(this%u, this%dof)
    call field_init(this%v, this%dof)
    call field_init(this%w, this%dof)
    
  end subroutine fluid_scheme_init

  subroutine fluid_scheme_free(this)
    class(fluid_scheme_t), intent(inout) :: this

    call field_free(this%u)
    call field_free(this%v)
    call field_free(this%w)
    call field_free(this%p)

    call space_free(this%Xh)

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if
    
  end subroutine fluid_scheme_free
  
end module fluid_method
