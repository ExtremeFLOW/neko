!> Defines various Conjugate Gradient methods
module cg
  use krylov
  use math
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_t
     real(kind=dp), allocatable :: w(:)
     real(kind=dp), allocatable :: r(:)
     real(kind=dp), allocatable :: p(:)
     real(kind=dp), allocatable :: z(:)
   contains
     procedure, pass(this) :: init => cg_init
     procedure, pass(this) :: free => cg_free
     procedure, pass(this) :: solve => cg_solve
  end type cg_t

contains

  !> Initialise a standard PCG solver
  subroutine cg_init(this, n, M, rel_tol, abs_tol)
    class(cg_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=dp), optional, intent(inout) :: rel_tol
    real(kind=dp), optional, intent(inout) :: abs_tol

        
    call this%free()
    
    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%p(n))
    allocate(this%z(n))
    if (present(M)) then 
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(rel_tol=rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(abs_tol=abs_tol)
    else
       call this%ksp_init()
    end if
          
  end subroutine cg_init

  !> Deallocate a standard PCG solver
  subroutine cg_free(this)
    class(cg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if
    
    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    nullify(this%M)


  end subroutine cg_free
  
  !> Standard PCG solve
  function cg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(inout) :: n
    real(kind=dp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=dp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=dp) :: beta, pap, alpha, alphm, eps, norm_fac
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1./sqrt(coef%volume)

    rtz1 = 1d0
    call rzero(x%x, n)
    call copy(this%r, f, n)

    rtr = sqrt(glsc3(this%r, coef%mult, this%r, n))
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    if(rnorm .eq. 0d0) return
    do iter = 1, max_iter
       call this%M%solve(this%z, this%r, n)
       rtz2 = rtz1
       rtz1 = glsc3(this%r, coef%mult, this%z, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = 0d0
       call add2s1(this%p, this%z, beta, n)
       
       call Ax%compute(this%w, this%p, coef, x%msh, x%Xh)
       call gs_op(gs_h, this%w, n, GS_OP_ADD)
       call bc_list_apply(blst, this%w, n)

       pap = glsc3(this%w, coef%mult, this%p, n)

       alpha = rtz1 / pap
       alphm = -alpha
       call add2s2(x%x, this%p, alpha, n)
       call add2s2(this%r, this%w, alphm, n)

       rtr = glsc3(this%r, coef%mult, this%r, n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr)*norm_fac
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function cg_solve

end module cg
  

