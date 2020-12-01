!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4
  use fluid_method
  implicit none

  type, extends(fluid_scheme_t) :: fluid_plan4_t
     !>@todo Add plan4 related data, ax, precon ect
   contains
     procedure, pass(this) :: init => fluid_plan4_init
     procedure, pass(this) :: free => fluid_plan4_free
     procedure, pass(this) :: step => fluid_plan4_step
  end type fluid_plan4_t

contains

  subroutine fluid_plan4_init(this, msh, lx, vel, prs)
    class(fluid_plan4_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=80), intent(inout) :: vel
    character(len=80), intent(inout) :: prs
    
    !> Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, solver_vel=vel, solver_prs=prs)

  end subroutine fluid_plan4_init

  subroutine fluid_plan4_free(this)
    class(fluid_plan4_t), intent(inout) :: this

    ! Deallocate velocity and pressure fields
    call this%scheme_free()
    
  end subroutine fluid_plan4_free
  
  subroutine fluid_plan4_step(this)
    class(fluid_plan4_t), intent(inout) :: this

  end subroutine fluid_plan4_step
  
end module fluid_plan4
