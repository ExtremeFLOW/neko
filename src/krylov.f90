!> Krylov solver
module krylov
  use ax_product
  use num_types
  use precon
  use mesh
  use field
  use utils
  implicit none

  integer, public, parameter :: KSP_MAX_ITER = 1e4 !< Maximum number of iters.

  !> Base type for a canonical Krylov method, solving \f$ Ax = f \f$
  type, abstract :: ksp_t
     type(pc_t) :: M            !< Preconditioner
   contains
     procedure, pass(this) :: ksp => krylov_init
     procedure, pass(this) :: ksp_free => krylov_free
     procedure(ksp_method), pass(this), deferred :: solve
  end type ksp_t
  
  !> Abstract interface for a Krylov method's solve routine
  !!
  !! @param x field to solve for
  !! @param f right hand side 
  !! @param n integer, size of vectors
  !! @param niter iteration trip count
  abstract interface
     subroutine ksp_method(this, Ax, x, f, n, niter)
       import :: field_t
       import :: ksp_t
       import :: ax_t
       import dp
       implicit none
       class(ksp_t), intent(inout) :: this
       class(ax_t), intent(inout) :: Ax
       type(field_t), intent(inout) :: x
       real(kind=dp), dimension(n), intent(inout) :: f
       integer, intent(inout) :: n
       integer, optional, intent(in) :: niter
     end subroutine ksp_method
  end interface
  
contains

  !> Create a krylov solver with @a nvec of size @a n
  subroutine krylov_init(this)    
    class(ksp_t), intent(inout) :: this
    integer :: i
    
    call krylov_free(this)

  end subroutine krylov_init
  
  !> Deallocate a Krylov solver
  subroutine krylov_free(this)
    class(ksp_t), intent(inout) :: this

  end subroutine krylov_free
  
end module krylov
