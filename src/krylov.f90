!> Krylov solver
module krylov
  use gather_scatter
  use ax_product
  use num_types
  use precon
  use mesh
  use field
  use utils
  use bc
  implicit none

  integer, public, parameter :: KSP_MAX_ITER = 1e4 !< Maximum number of iters.
  real(kind=dp), public, parameter :: KSP_ABS_TOL = 1d-20 !< Absolut tolerance
  real(kind=dp), public, parameter :: KSP_REL_TOL = 1d-10 !< Relative tolerance

  !> Base type for a canonical Krylov method, solving \f$ Ax = f \f$
  type, public, abstract :: ksp_t
     type(pc_t) :: M            !< Preconditioner
     real(kind=dp) :: rel_tol   !< Relative tolerance
     real(kind=dp) :: abs_tol   !< Absolute tolerance
   contains
     procedure, pass(this) :: ksp_init => krylov_init
     procedure, pass(this) :: ksp_free => krylov_free
     procedure(ksp_method), pass(this), deferred :: solve
  end type ksp_t
  
  !> Abstract interface for a Krylov method's solve routine
  !!
  !! @param x field to solve for
  !! @param f right hand side 
  !! @param n integer, size of vectors
  !! @param blst list of  boundary conditions
  !! @param gs_h Gather-scatter handle
  !! @param niter iteration trip count
  abstract interface
     function ksp_method(this, Ax, x, f, n, blst, gs_h, niter) result(iter)
       import :: bc_list_t       
       import :: field_t
       import :: ksp_t
       import :: gs_t
       import :: ax_t
       import dp
       implicit none
       class(ksp_t), intent(inout) :: this
       class(ax_t), intent(inout) :: Ax
       type(field_t), intent(inout) :: x
       integer, intent(inout) :: n
       real(kind=dp), dimension(n), intent(inout) :: f
       type(bc_list_t), intent(inout) :: blst
       type(gs_t), intent(inout) :: gs_h              
       integer, optional, intent(in) :: niter       
       integer :: iter
     end function ksp_method
  end interface
  
contains

  !> Create a krylov solver
  subroutine krylov_init(this, rel_tol, abs_tol)    
    class(ksp_t), intent(inout) :: this
    real(kind=dp), optional, intent(in) :: rel_tol
    real(kind=dp), optional, intent(in) :: abs_tol
    integer :: i
    
    call krylov_free(this)

    if (present(rel_tol)) then
       this%rel_tol = rel_tol
    else
       this%rel_tol = KSP_REL_TOL
    end if

    if (present(abs_tol)) then
       this%abs_tol = abs_tol
    else
       this%abs_tol = KSP_ABS_TOL
    end if

    if (.not. associated(this%M%solve)) then
       this%M%solve => pc_ident
    end if

  end subroutine krylov_init
  
  !> Deallocate a Krylov solver
  subroutine krylov_free(this)
    class(ksp_t), intent(inout) :: this

    !> @todo add calls to destroy precon. if necessary

  end subroutine krylov_free
  
end module krylov
