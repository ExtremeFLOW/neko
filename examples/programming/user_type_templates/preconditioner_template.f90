!> Template for a user-defined Krylov preconditioner.
module preconditioner_template
  use num_types, only : rp
  use precon, only : pc_t, precon_allocate, register_precon
  implicit none
  private

  type, extends(pc_t) :: preconditioner_template_t
   contains
     procedure :: solve => preconditioner_template_solve
     procedure :: update => preconditioner_template_update
  end type preconditioner_template_t

  public :: preconditioner_template_register_types

contains

  subroutine preconditioner_template_solve(this, z, r, n)
    class(preconditioner_template_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: z(n), r(n)
    ! TODO: Replace the identity operation with M z = r.
    z = r
  end subroutine preconditioner_template_solve

  subroutine preconditioner_template_update(this)
    class(preconditioner_template_t), intent(inout) :: this
    ! TODO: Update state after a change in geometry or operator coefficients.
  end subroutine preconditioner_template_update

  subroutine preconditioner_template_register_types()
    procedure(precon_allocate), pointer :: allocator
    allocator => preconditioner_template_allocate
    call register_precon('preconditioner_template', allocator)
  end subroutine preconditioner_template_register_types

  subroutine preconditioner_template_allocate(obj)
    class(pc_t), allocatable, intent(inout) :: obj
    allocate(preconditioner_template_t :: obj)
  end subroutine preconditioner_template_allocate
end module preconditioner_template
