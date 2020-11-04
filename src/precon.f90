!> Krylov preconditioner
module precon
  use math
  use utils
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type :: pc_t
     procedure(pc_solve), nopass, pointer :: solve => pc_ident
  end type pc_t

  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
  abstract interface
     subroutine pc_solve(z, r, n)
       import dp
       implicit none
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), dimension(n), intent(inout) :: r
       integer, intent(inout) :: n
     end subroutine pc_solve
  end interface

contains

  !> The (default) naive preconditioner \f$ I z = r \f$
  subroutine pc_ident(z, r, n)
    real(kind=dp), dimension(n), intent(inout) :: z
    real(kind=dp), dimension(n), intent(inout) :: r
    integer, intent(inout) :: n    
    call copy(z, r, n)    
  end subroutine pc_ident

end module precon
