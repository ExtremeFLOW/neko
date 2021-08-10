!> Krylov preconditioner
module precon
  use math
  use utils
  implicit none
  
  !> Defines a canonical Krylov preconditioner
  type, public, abstract :: pc_t
   contains
     procedure(pc_solve), pass(this), deferred :: solve
     procedure(pc_update), pass(this), deferred :: update
  end type pc_t

  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
  abstract interface
     subroutine pc_solve(this, z, r, n)
       import rp
       import :: pc_t
       implicit none
       integer, intent(inout) :: n
       class(pc_t), intent(inout) :: this
       real(kind=rp), dimension(n), intent(inout) :: z
       real(kind=rp), dimension(n), intent(inout) :: r
     end subroutine pc_solve
     subroutine pc_update(this)
       import rp
       import :: pc_t
       implicit none
       class(pc_t), intent(inout) :: this
     end subroutine pc_update
  end interface
end module precon
