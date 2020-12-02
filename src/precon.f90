!> Krylov preconditioner
module precon
  use math
  use utils
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type, public, abstract :: pc_t
   contains
     procedure(pc_solve), pass(this), deferred :: solve
  end type pc_t

  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
  abstract interface
     subroutine pc_solve(this, z, r, n)
       import dp
       import :: pc_t
       implicit none
       integer, intent(inout) :: n
       class(pc_t), intent(inout) :: this
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), dimension(n), intent(inout) :: r
     end subroutine pc_solve
  end interface
end module precon
