!> Krylov solver
module krylov
  use gather_scatter
  use num_types
  use precon
  use mesh
  use field
  use space
  use utils
  implicit none

  !> Base type for a canonical Krylov method, solving \f$ Ax = f \f$
  type, abstract :: ksp_t
     type(pc_t) :: M            !< Preconditioner
     procedure(ksp_ax), nopass, pointer :: Ax     
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
  !! @param niter iteartion trip count
  abstract interface
     subroutine ksp_method(this, x, f, n, niter)
       import :: field_t
       import :: ksp_t
       import dp
       implicit none
       class(ksp_t), intent(inout) :: this
       type(field_t), intent(inout) :: x
       real(kind=dp), dimension(n), intent(inout) :: f
       integer, intent(inout) :: n
       integer, intent(in) :: niter
     end subroutine ksp_method
  end interface

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of length @a n
  !! @param z vector of length @a n
  !! @param g geometric factors
  !! @param gs_h gather-scatter handle
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  !! @param n integer, size of vectors
  abstract interface
     subroutine ksp_ax(w, z, g, gs_h, msh, Xh, n)
       import space_t
       import mesh_t
       import gs_t
       import dp
       implicit none
       type(gs_t), intent(inout) :: gs_h
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       integer, intent(inout) :: n
       real(kind=dp), dimension(n), intent(inout) :: w
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), intent(inout) :: g(6, Xh%lx, Xh%ly, Xh%lz)
     end subroutine ksp_ax
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
