!> Identity Krylov preconditioner for accelerators
module device_identity
  use utils
  use precon
  use device
  use device_math
  use num_types
  implicit none
  
  !> Defines a canonical Krylov preconditioner for accelerators
  type, public, extends(pc_t) :: device_ident_t
  contains
     procedure, pass(this) :: solve => device_ident_solve
     procedure, pass(this) :: update => device_ident_update
  end type device_ident_t

contains

  !> The (default) naive preconditioner \f$ I z = r \f$
  subroutine device_ident_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(device_ident_t), intent(inout) :: this
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    
    z_d = device_get_ptr(z, n)
    r_d = device_get_ptr(r, n)
    
    call device_copy(z_d, r_d, n)
    
  end subroutine device_ident_solve
  
  !> Mandatory update routine (NOP)
  subroutine device_ident_update(this)
    class(device_ident_t), intent(inout) :: this
  end subroutine device_ident_update
  
end module device_identity
