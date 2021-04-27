!> Defines various Conjugate Gradient methods
module cg
  use krylov
  use math
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_t
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: z(:)
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
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol

        
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
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    real(kind=rp), parameter :: one = 1.0
    real(kind=rp), parameter :: zero = 0.0
    integer :: iter, max_iter
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=rp) :: beta, pap, alpha, alphm, eps, norm_fac
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = one/sqrt(coef%volume)

    rtz1 = one
    call rzero(x%x, n)
    call copy(this%r, f, n)

    rtr = glsc3(this%r, coef%mult, this%r, n)
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    if(rnorm .eq. zero) return
    do iter = 1, max_iter
       call this%M%solve(this%z, this%r, n)
       rtz2 = rtz1
       rtz1 = glsc3(this%r, coef%mult, this%z, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = zero
       call add2s1(this%p, this%z, beta, n)
       
       call Ax%compute(this%w, this%p, coef, x%msh, x%Xh)
       call gs_op(gs_h, this%w, n, GS_OP_ADD)
       call bc_list_apply(blst, this%w, n)

       pap = glsc3(this%w, coef%mult, this%p, n)

       alpha = rtz1 / pap
       alphm = -alpha
       call add2s2(x%x, this%p, alpha, n)
       call add2s2(this%r, this%w, alphm, n)
       !compute residuals
       if (x%Xh%lx .eq. 2) then
          rnorm = sqrt(glsc3(this%r, coef%mult, this%r,n))*norm_fac
       else
          if ( pe_rank .eq. 0) write (*,*) 'Iteration:', iter
          rnorm = calc_norms(this%r, f, coef%mult, coef%Binv, norm_fac, n)
          call Ax%compute(this%w, x%x, coef, x%msh, x%Xh)
          call gs_op(gs_h, this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          call sub2(this%w,f,n)
          rnorm = calc_norms(this%w, f, coef%mult, coef%Binv, norm_fac, n)
       end if
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function cg_solve

  function calc_norms(r, f, mult, Binv, norm_fac, n) result(rnorm)
    integer :: n
    real(kind=rp) :: r(n), mult(n), Binv(n), f(n)
    real(kind=rp) :: norm_fac
    real(kind=rp) :: alg_norm, vel_norm
    real(kind=rp) :: l2_norm, max_norm, rel_l2_norm, max_rel_norm, rnorm
    
    l2_norm = sqrt(glsc3(r, mult, r, n))
    
    alg_norm = l2_norm*norm_fac
    vel_norm = sqrt(glsc4(r, mult, Binv, r, n))*norm_fac
    
    max_norm = max(glmax(r, n), -glmin(r, n))
    max_rel_norm = max_norm/max(glmax(f,n), -glmin(f,n))
    rel_l2_norm = l2_norm/sqrt(glsc3(f,mult,f, n))
    rnorm = alg_norm
    if( pe_rank .eq. 0) write (*,*) alg_norm, vel_norm, l2_norm, rel_l2_norm, max_norm, max_rel_norm

  end function calc_norms

end module cg
  

