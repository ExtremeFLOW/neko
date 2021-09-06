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
     real(kind=rp), allocatable :: p(:,:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: alpha(:)
     integer :: p_space
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
    this%p_space = 50
    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%p(n,this%p_space))
    allocate(this%z(n))
    allocate(this%alpha(this%p_space))
    
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
    integer, parameter :: BLOCK_SIZE = 50000
    integer :: iter, max_iter, x_update, i, j, p_cur, p_prev, k
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1, x_plus(BLOCK_SIZE)
    real(kind=rp) :: beta, pap, alpha, alphm, eps, norm_fac
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = one / sqrt(coef%volume)

    rtz1 = one
    call rzero(x%x, n)
    call rzero(this%p, n)
    call copy(this%r, f, n)

    rtr = glsc3(this%r, coef%mult, this%r, n)
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    p_prev = this%p_space
    p_cur = 1
    if(rnorm .eq. zero) return
    do iter = 1, max_iter
       call this%M%solve(this%z, this%r, n)
       rtz2 = rtz1
       rtz1 = glsc3(this%r, coef%mult, this%z, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = zero
       !call add2s1(this%p, this%z, beta, n)
       do i = 1, n
          this%p(i,p_cur) = this%z(i) + beta*this%p(i,p_prev)
       end do
       
       call Ax%compute(this%w, this%p(1,p_cur), coef, x%msh, x%Xh)
       call gs_op(gs_h, this%w, n, GS_OP_ADD)
       call bc_list_apply(blst, this%w, n)

       pap = glsc3(this%w, coef%mult, this%p(1,p_cur), n)

       this%alpha(p_cur) = rtz1 / pap
       call second_cg_part(rtr, this%r, coef%mult, this%w, this%alpha(p_cur), n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr) * norm_fac
       if (p_cur .eq. this%p_space .or. rnorm .lt. this%abs_tol .or. iter .eq. max_iter) then
           do i = 0,n,BLOCK_SIZE
              if (i + BLOCK_SIZE .le. n) then
                 do k = 1, BLOCK_SIZE
                    x_plus(k) = 0.0
                 end do
                 do j = 1, p_cur
                    do k = 1, BLOCK_SIZE
                       x_plus(k) = x_plus(k) + this%alpha(j)*this%p(i+k,j)
                    end do
                 end do
                 do k = 1, BLOCK_SIZE
                    x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                 end do
              else 
                 do k = 1, n-i
                    x_plus(1) = 0.0
                    do j = 1, p_cur
                       x_plus(1) = x_plus(1) + this%alpha(j)*this%p(i+k,j)
                    end do
                    x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                 end do
              end if
          end do 
          p_prev = p_cur
          p_cur = 1
          if (rnorm .lt. this%abs_tol) exit
       else
          p_prev = p_cur
          p_cur = p_cur + 1
       end if
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function cg_solve
  subroutine second_cg_part(rtr, r, mult, w, alpha, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n), rtr
    real(kind=rp), intent(in) ::mult(n), w(n), alpha 
    real(kind=rp) :: tmp
    integer :: i, ierr
    tmp = 0.0
    do i = 1,n
       r(i) =r(i) - alpha*w(i)
       tmp = tmp + r(i) * r(i) * mult(i)
    end do
    call MPI_Allreduce(tmp, rtr, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
  end subroutine second_cg_part 

end module cg
  

