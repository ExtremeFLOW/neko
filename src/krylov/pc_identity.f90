!> Krylov preconditioner
module identity
  use math
  use utils
  use precon
  use ax_product
  use num_types
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type, public, extends(pc_t) :: ident_t
  contains
     procedure, pass(this) :: solve => ident_solve
     procedure, pass(this) :: update => ident_update
  end type ident_t

contains
  !> The (default) naive preconditioner \f$ I z = r \f$
  subroutine ident_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(ident_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    call copy(z, r, n)    
  end subroutine ident_solve
  !> Mandatory update routine
  subroutine ident_update(this)
    class(ident_t), intent(inout) :: this
  end subroutine ident_update
end module identity
