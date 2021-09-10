!> Defines various Conjugate Gradient methods
module cg_gs
  use krylov
  use math
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_gs_t
     real(kind=rp), allocatable :: w(:), w_s(:)
     real(kind=rp), allocatable :: r(:), r_s(:)
     real(kind=rp), allocatable :: p(:,:), p_s(:)
     real(kind=rp), allocatable :: x(:)
     real(kind=rp), allocatable :: alpha(:)
     integer :: p_space
   contains
     procedure, pass(this) :: init => cg_gs_init
     procedure, pass(this) :: free => cg_gs_free
     procedure, pass(this) :: solve => cg_gs_solve
  end type cg_gs_t

contains

  !> Initialise a standard PCG solver
  subroutine cg_gs_init(this, n, M, rel_tol, abs_tol)
    class(cg_gs_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
        
    call this%free()
    this%p_space = 50
    allocate(this%w(n), this%w_s(n))
    allocate(this%r(n), this%r_s(n))
    allocate(this%p(n,this%p_space), this%p_s(n))
    allocate(this%x(n))
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
          
  end subroutine cg_gs_init

  !> Deallocate a standard PCG solver
  subroutine cg_gs_free(this)
    class(cg_gs_t), intent(inout) :: this

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
    
    
    if (allocated(this%alpha)) then
       deallocate(this%alpha)
    end if

    nullify(this%M)

  end subroutine cg_gs_free
  
  !> Standard PCG solve
  function cg_gs_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cg_gs_t), intent(inout) :: this
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
    integer :: iter, max_iter, i, j, k, p_cur, p_prev, n_dofs, it
    real(kind=rp) :: rnorm, rtr, rtz2, rtz1, x_plus(NEKO_BLK_SIZE)
    real(kind=rp) :: beta, pap, norm_fac
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = one / sqrt(coef%volume)

    associate(w => this%w, r => this%r, p => this%p, alpha => this%alpha)
      n_dofs = n/(x%Xh%lx**3)*((x%Xh%lx-2)**3)
      it = 0
      do i = 1,gs_h%nlocal_blks
         it = it + 1
      end do
      do i = gs_h%local_facet_offset, gs_h%nlocal,2
         it = it + 1 
      end do

      if (pe_size .gt. 1) then
         k = 0
         do i = 1,gs_h%nshared_blks
            it = it + 1
         end do
         do i = gs_h%shared_facet_offset, gs_h%nshared
            it = it + 1
         end do
      end if
      n_dofs = n_dofs + it
      print *, n, n_dofs, it

      rtz1 = one
      call rzero(x%x, n)
      call rzero(p, n)
      call copy(this%r_s, f, n)
      call col2(this%r_s, coef%mult,n)
      call gather(r, this%r_s, x%Xh, gs_h, n) 
      rtr = glsc2(r, r, n_dofs)
      rnorm = sqrt(rtr) * norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      p_prev = this%p_space
      p_cur = 1
      if(rnorm .eq. zero) return
      do iter = 1, max_iter
         !call this%M%solve(z, r, n)
         rtz2 = rtz1
         rtz1 = glsc2(r, r, n_dofs)
      
         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = zero
         do i = 1, n_dofs
            p(i,p_cur) = r(i) + beta * p(i,p_prev)
         end do
   
         call scatter(p(1,p_cur), this%p_s, x%Xh, gs_h, n) 
       
         call Ax%compute(this%w_s, this%p_s, coef, x%msh, x%Xh)
         call bc_list_apply(blst, this%w_s, n)
         call gather(w, this%w_s, x%Xh, gs_h, n) 
         pap = glsc2(w, p(1,p_cur), n_dofs)
         
         alpha(p_cur) = rtz1 / pap
         call second_cg_gs_part(rtr, r, coef%mult, w, alpha(p_cur), n_dofs)
         rnorm = sqrt(rtr) * norm_fac

         if ((p_cur .eq. this%p_space) .or. &
             (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
            do i = 0, n_dofs, NEKO_BLK_SIZE
               if (i + NEKO_BLK_SIZE .le. n) then
                  do k = 1, NEKO_BLK_SIZE
                     x_plus(k) = 0.0
                  end do
                  do j = 1, p_cur
                     do k = 1, NEKO_BLK_SIZE
                        x_plus(k) = x_plus(k) + alpha(j) * p(i+k,j)
                     end do
                  end do
                  do k = 1, NEKO_BLK_SIZE
                     this%x(i+k) = this%x(i+k) + x_plus(k)
                  end do
               else 
                  do k = 1, n_dofs-i
                     x_plus(1) = 0.0
                     do j = 1, p_cur
                        x_plus(1) = x_plus(1) + alpha(j) * p(i+k,j)
                     end do
                     this%x(i+k) = this%x(i+k) + x_plus(1)
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
    end associate
    call scatter(this%x,x%x, x%Xh, gs_h, n)
    !WHY?
    call cmult(x%x, 0.5_rp,n)
    ksp_results%res_final = rnorm
    ksp_results%iter = iter

  end function cg_gs_solve

  subroutine second_cg_gs_part(rtr, r, mult, w, alpha, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n), rtr
    real(kind=rp), intent(in) ::mult(n), w(n), alpha 
    integer :: i, ierr

    rtr = 0.0
    do i = 1, n
       r(i) = r(i) - alpha*w(i)
       rtr = rtr + r(i) * r(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, rtr, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end subroutine second_cg_gs_part 

  subroutine scatter(u,u_s,Xh,gs,n)
    integer, intent(inout) :: n
    type(gs_t), intent(inout) :: gs
    type(space_t), intent(inout) :: Xh
    real(kind=rp), intent(inout) :: u(n), u_s(Xh%lx,Xh%ly,Xh%lz,n/(Xh%lxyz))
    integer :: it, e, i, j, k
    it = 0
    do e = 1,n/Xh%lxyz
       do k = 2, Xh%lx-1
          do j = 2, Xh%lx-1
             do i = 2, Xh%lx-1
                it = it + 1
                u_s(i,j,k,e) = u(it)
             end do
          end do
       end do
    end do
    k = 0
    do i = 1,gs%nlocal_blks
       it = it + 1
       gs%local_gs(gs%local_dof_gs(k+1)) = u(it)
       k = k + gs%local_blk_len(i)
    end do
    do i = gs%local_facet_offset, gs%nlocal,2
       it = it + 1 
       gs%local_gs(gs%local_dof_gs(i)) = u(it)
    end do

    if (pe_size .gt. 1) then
       k = 0
       do i = 1,gs%nshared_blks
          it = it + 1
          gs%shared_gs(gs%shared_dof_gs(k+1)) = u(it)
          k = k + gs%shared_blk_len(i)
       end do
       do i = gs%shared_facet_offset, gs%nshared
          it = it + 1
          gs%shared_gs(gs%shared_dof_gs(i)) = u(it)
       end do
    end if
    call gs_scatter_vector(gs, u_s, n, GS_OP_ADD)
  end subroutine scatter

  subroutine gather(u,u_s,Xh,gs,n)
    integer, intent(inout) :: n
    type(gs_t), intent(inout) :: gs
    type(space_t), intent(inout) :: Xh
    real(kind=rp), intent(inout) :: u(n), u_s(Xh%lx,Xh%ly,Xh%lz,n/(Xh%lxyz))
    integer :: it, e, i, j, k
    it = 0
    call gs_gather_vector(gs, u_s, n, GS_OP_ADD)
    do e = 1,n/Xh%lxyz
       do k = 2, Xh%lx-1
          do j = 2, Xh%lx-1
             do i = 2, Xh%lx-1
                it = it + 1
                u(it) = u_s(i,j,k,e) 
             end do
          end do
       end do
    end do
    k = 0
    do i = 1,gs%nlocal_blks
       it = it + 1
       u(it) = gs%local_gs(gs%local_dof_gs(k+1))
       k = k + gs%local_blk_len(i)
    end do
    do i = gs%local_facet_offset, gs%nlocal,2
       it = it + 1 
       u(it) = gs%local_gs(gs%local_dof_gs(i))
    end do

    if (pe_size .gt. 1) then
       k = 0
       do i = 1,gs%nshared_blks
          it = it + 1
          u(it) = gs%shared_gs(gs%shared_dof_gs(k+1))
          k = k + gs%shared_blk_len(i)
       end do
       do i = gs%shared_facet_offset, gs%nshared
          it = it + 1
          u(it) = gs%shared_gs(gs%shared_dof_gs(i))
       end do
    end if
  end subroutine gather

end module cg_gs
  

