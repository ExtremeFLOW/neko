!> Defines a pipelined Conjugate Gradient methods
module pipecg
  use krylov
  use math
  use num_types
  implicit none

  !> Pipelined preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: pipecg_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: q(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: u(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: mi(:)
     real(kind=rp), allocatable :: ni(:)
   contains
     procedure, pass(this) :: init => pipecg_init
     procedure, pass(this) :: free => pipecg_free
     procedure, pass(this) :: solve => pipecg_solve
  end type pipecg_t

contains

  !> Initialise a pipelined PCG solver
  subroutine pipecg_init(this, n, M, rel_tol, abs_tol)
    class(pipecg_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
        
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

  !> Deallocate a pipelined PCG solver
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
  
  !> Pipelined PCG solve
  function pipecg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(pipecg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(inout) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, ierr
    real(kind=rp) :: rnorm, rtr, reduction(3), norm_fac 
    real(kind=rp) :: alpha, beta, gamma1, gamma2, delta
    real(kind=rp) :: tmp1, tmp2, tmp3
    type(MPI_Request) :: request
    type(MPI_Status) :: status
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(p => this%p, q => this%q, r => this%r, s => this%s, &
         u => this%u, w => this%w, z => this%z, mi => this%mi, ni => this%ni)
      
      call rzero(x%x, n)
      call rzero(z, n)
      call rzero(q, n)
      call rzero(p, n)
      call rzero(s, n)
      call copy(r, f, n)
      call this%M%solve(u, r, n)
      call Ax%compute(w, u, coef, x%msh, x%Xh)
      call gs_op(gs_h, w, n, GS_OP_ADD)
      call bc_list_apply(blst, w, n)
    
      rtr = glsc3(r, coef%mult, r, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(rnorm .eq. 0.0_rp) return

      gamma1 = 0.0_rp
      
      do iter = 1, max_iter

         tmp1 = 0.0_rp
         tmp2 = 0.0_rp
         tmp3 = 0.0_rp
         do i = 1, n
            tmp1 = tmp1 + r(i) * coef%mult(i,1,1,1) * u(i)
            tmp2 = tmp2 + w(i) * coef%mult(i,1,1,1) * u(i)
            tmp3 = tmp3 + r(i) * coef%mult(i,1,1,1) * r(i)
         end do
         reduction(1) = tmp1
         reduction(2) = tmp2
         reduction(3) = tmp3

         call MPI_Iallreduce(MPI_IN_PLACE, reduction, 3, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, request, ierr)
         
         call this%M%solve(mi, w, n)
         call Ax%compute(ni, mi, coef, x%msh, x%Xh)
         call gs_op(gs_h, ni, n, GS_OP_ADD)
         call bc_list_apply(blst, ni, n)

         call MPI_Wait(request, status, ierr)
         gamma2 = gamma1       
         gamma1 = reduction(1)
         delta = reduction(2)
         rtr = reduction(3)

         rnorm = sqrt(rtr)*norm_fac
         if (rnorm .lt. this%abs_tol) then
            exit
         end if
         
         if (iter .gt. 1) then
            beta = gamma1 / gamma2
            alpha = gamma1 / (delta - (beta * gamma1/alpha))
         else 
            beta = 0.0_rp
            alpha = gamma1/delta
         end if
         
         call add2s1(z, ni, beta, n)
         call add2s1(q, mi, beta, n)
         call add2s1(s, w, beta, n)
         call add2s1(p, u, beta, n)
 
         call add2s2(x%x, p, alpha, n)
         call add2s2(r, s, -alpha, n)
         call add2s2(u, q, -alpha, n)
         call add2s2(w, z, -alpha, n)
         
      end do
      
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
      
    end associate
    
  end function pipecg_solve
   
end module pipecg
  

