!> Defines various Conjugate Gradient methods
module pipecg
  use krylov
  use math
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: pipecg_t
     real(kind=dp), allocatable :: p(:)
     real(kind=dp), allocatable :: q(:)
     real(kind=dp), allocatable :: r(:)
     real(kind=dp), allocatable :: s(:)
     real(kind=dp), allocatable :: u(:)
     real(kind=dp), allocatable :: w(:)
     real(kind=dp), allocatable :: z(:)
     real(kind=dp), allocatable :: mi(:)
     real(kind=dp), allocatable :: ni(:)
   contains
     procedure, pass(this) :: init => pipecg_init
     procedure, pass(this) :: free => pipecg_free
     procedure, pass(this) :: solve => pipecg_solve
  end type pipecg_t

contains

  !> Initialise a standard PCG solver
  subroutine pipecg_init(this, n, M, rel_tol, abs_tol)
    class(pipecg_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=dp), optional, intent(inout) :: rel_tol
    real(kind=dp), optional, intent(inout) :: abs_tol

        
    call this%free()
    
    allocate(this%p(n))
    allocate(this%q(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%u(n))
    allocate(this%w(n))
    allocate(this%z(n))
    allocate(this%mi(n))
    allocate(this%ni(n))
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
          
  end subroutine pipecg_init

  !> Deallocate a standard PCG solver
  subroutine pipecg_free(this)
    class(pipecg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%p)) then
       deallocate(this%p)
    end if
    if (allocated(this%q)) then
       deallocate(this%q)
    end if
    if (allocated(this%r)) then
       deallocate(this%r)
    end if
    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    if (allocated(this%u)) then
       deallocate(this%u)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%z)) then
       deallocate(this%z)
    end if
    if (allocated(this%mi)) then
       deallocate(this%mi)
    end if
    if (allocated(this%ni)) then
       deallocate(this%ni)
    end if

    nullify(this%M)


  end subroutine pipecg_free
  
  !> Standard PCG solve
  function pipecg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(pipecg_t), intent(inout) :: this
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
    real(kind=dp) :: rnorm, rtr, reduction(3), norm_fac 
    real(kind=dp) :: alpha, beta, gamma1, gamma2, delta
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1./sqrt(coef%volume)

    call rzero(x%x, n)
    call copy(this%r, f, n)
    call this%M%solve(this%u, this%r, n)
    call Ax%compute(this%w, this%u, coef, x%msh, x%Xh)
    call gs_op(gs_h, this%w, n, GS_OP_ADD)
    call bc_list_apply(blst, this%w, n)
    
    rtr = glsc4(this%r, coef%Binv,coef%mult, this%r, n)
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    if(rnorm .eq. 0d0) return
    do iter = 1, max_iter
       reduction = pipecg_reduction(this%r,this%w, this%u,coef%Binv, coef%mult, n)
       gamma2 = gamma1
       gamma1 = reduction(1)
       delta = reduction(2)
       rtr = reduction(3)
       rnorm = sqrt(rtr)*norm_fac
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
       call this%M%solve(this%mi, this%w, n)
       call Ax%compute(this%ni, this%mi, coef, x%msh, x%Xh)
       call gs_op(gs_h, this%ni, n, GS_OP_ADD)
       call bc_list_apply(blst, this%ni, n)
 
       if (iter .gt. 1) then
          beta = gamma1 / gamma2
          alpha = gamma1 / (delta - beta * gamma1/alpha)
       else 
          beta = 0d0 
          alpha = gamma1/delta
       end if
 
       call add2s1(this%z, this%ni, beta, n)
       call add2s1(this%q, this%mi, beta, n)
       call add2s1(this%s, this%w, beta, n)
       call add2s1(this%p, this%u, beta, n)
 
       
       call add2s2(x%x, this%p, alpha, n)
       call add2s2(this%r, this%s, -alpha, n)
       call add2s2(this%u, this%q, -alpha, n)
       call add2s2(this%w, this%z, -alpha, n)
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function pipecg_solve
   
  function pipecg_reduction(r, w, u, binv, mult, n) result(reduction)
    integer, intent(in) :: n
    real(kind=dp), dimension(n), intent(in) :: r
    real(kind=dp), dimension(n), intent(in) :: w
    real(kind=dp), dimension(n), intent(in) :: u
    real(kind=dp), dimension(n), intent(in) :: binv
    real(kind=dp), dimension(n), intent(in) :: mult
    real(kind=dp) :: tmp1, tmp2, tmp3, temp(3), reduction(3)
    integer :: i, ierr

    tmp1 = 0d0
    tmp2 = 0d0
    tmp3 = 0d0
    do i = 1, n
       tmp1 = tmp1 + r(i) * mult(i) * u(i)
       tmp2 = tmp2 + w(i) * mult(i) * u(i)
       ! This way of calculating the norm is really questionable
       tmp3 = tmp3 + r(i)**2 * mult(i) * binv(i)
    end do
    temp(1) = tmp1
    temp(2) = tmp2
    temp(3) = tmp3
    call MPI_Allreduce(temp, reduction, 3, &
         MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM, ierr)
  end function pipecg_reduction

end module pipecg
  

