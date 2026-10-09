!> Template for a user-defined Krylov preconditioner.
module preconditioner_template
  use num_types, only : rp
  use precon, only : pc_t, precon_allocate, register_precon
  use coefs, only : coef_t
  use bc_list, only : bc_list_t
  use json_module, only : json_file
  implicit none
  private

  type, extends(pc_t) :: preconditioner_template_t
   contains
     procedure :: init => preconditioner_template_init
     procedure :: solve => preconditioner_template_solve
     procedure :: update => preconditioner_template_update
     procedure :: free => preconditioner_template_free
  end type preconditioner_template_t

  public :: preconditioner_template_register_types

contains

  subroutine preconditioner_template_init(this, coef, bclst, json)
    class(preconditioner_template_t), intent(inout), target :: this
    type(coef_t), intent(in), target :: coef
    type(bc_list_t), intent(inout), target :: bclst
    type(json_file), intent(inout) :: json
    ! TODO: Set up the preconditioner from the SEM coefficients, the boundary
    ! conditions of the system and the preconditioner's case-file dictionary.
  end subroutine preconditioner_template_init

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

  subroutine preconditioner_template_free(this)
    class(preconditioner_template_t), intent(inout) :: this
    ! TODO: Release any resources held by the preconditioner.
  end subroutine preconditioner_template_free

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
