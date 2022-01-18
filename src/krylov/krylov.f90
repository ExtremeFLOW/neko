!> Krylov solver
module krylov
  use gather_scatter
  use ax_product
  use num_types
  use precon
  use coefs    
  use mesh
  use field
  use utils
  use bc
  use identity
  use device_identity
  use neko_config
  implicit none

  integer, public, parameter :: KSP_MAX_ITER = 1e4 !< Maximum number of iters.
  real(kind=rp), public, parameter :: KSP_ABS_TOL = 1d-9 !< Absolut tolerance
  real(kind=rp), public, parameter :: KSP_REL_TOL = 1d-9 !< Relative tolerance

  type, public :: ksp_monitor_t
    integer :: iter
    real(kind=rp) :: res_start
    real(kind=rp) :: res_final
  end type ksp_monitor_t

  !> Base type for a canonical Krylov method, solving \f$ Ax = f \f$
  type, public, abstract :: ksp_t
     class(pc_t), pointer :: M            !< Preconditioner
     real(kind=rp) :: rel_tol   !< Relative tolerance
     real(kind=rp) :: abs_tol   !< Absolute tolerance
   contains
     procedure, pass(this) :: ksp_init => krylov_init
     procedure, pass(this) :: ksp_free => krylov_free
     procedure, pass(this) :: set_pc => krylov_set_pc
     procedure(ksp_method), pass(this), deferred :: solve
     procedure(ksp_t_free), pass(this), deferred :: free
  end type ksp_t

  
  !> Abstract interface for a Krylov method's solve routine
  !!
  !! @param x field to solve for
  !! @param f right hand side 
  !! @param n integer, size of vectors
  !! @param coef Coefficients
  !! @param blst list of  boundary conditions
  !! @param gs_h Gather-scatter handle
  !! @param niter iteration trip count
  abstract interface
     function ksp_method(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
       import :: bc_list_t       
       import :: field_t
       import :: ksp_t
       import :: coef_t
       import :: gs_t
       import :: ax_t
       import :: ksp_monitor_t
       import rp
       implicit none
       class(ksp_t), intent(inout) :: this
       class(ax_t), intent(inout) :: Ax
       type(field_t), intent(inout) :: x
       integer, intent(inout) :: n
       real(kind=rp), dimension(n), intent(inout) :: f
       type(coef_t), intent(inout) :: coef
       type(bc_list_t), intent(inout) :: blst
       type(gs_t), intent(inout) :: gs_h              
       integer, optional, intent(in) :: niter       
       type(ksp_monitor_t) :: ksp_results
     end function ksp_method
  end interface

  !> Abstract interface for deallocating a Krylov method
  abstract interface
     subroutine ksp_t_free(this)
       import :: ksp_t
       class(ksp_t), intent(inout) :: this
     end subroutine ksp_t_free
  end interface
  
contains

  !> Create a krylov solver
  subroutine krylov_init(this, rel_tol, abs_tol, M)    
    class(ksp_t), intent(inout) :: this
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    class(pc_t), optional, target, intent(in) :: M
    class(pc_t), allocatable, target :: M_ident
    
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

    if (present(M)) then
       this%M => M
    else
       if (.not. associated(this%M)) then
          if (NEKO_BCKND_HIP .eq. 1 .or. NEKO_BCKND_CUDA .eq. 1) then
             allocate(device_ident_t::M_ident)
          else
             allocate(ident_t::M_ident)
          end if
          this%M => M_ident
       end if
    end if

  end subroutine krylov_init
  
  !> Deallocate a Krylov solver
  subroutine krylov_free(this)
    class(ksp_t), intent(inout) :: this

    !> @todo add calls to destroy precon. if necessary

  end subroutine krylov_free

  !> Setup a Krylov solvers preconditioners
  subroutine krylov_set_pc(this, M)
    class(ksp_t), intent(inout) :: this
    class(pc_t), optional, target, intent(in) :: M

    if (associated(this%M)) then
       select type(pc => this%M)
       type is (ident_t)
       class default
          call neko_error('Preconditioner already defined')
       end select
    end if
    
    this%M => M
    
  end subroutine krylov_set_pc
  
end module krylov
