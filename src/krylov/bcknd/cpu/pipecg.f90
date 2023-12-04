! Copyright (c) 2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Defines a pipelined Conjugate Gradient methods
module pipecg
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, copy
  use comm
  implicit none
  private

  integer, parameter :: PIPECG_P_SPACE = 7

  !> Pipelined preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: pipecg_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: q(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: u(:,:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: mi(:)
     real(kind=rp), allocatable :: ni(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => pipecg_init
     !> Destructor.
     procedure, pass(this) :: free => pipecg_free
     !> Solve the linear system.
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
    allocate(this%u(n,PIPECG_P_SPACE+1))
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
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, j, k, ierr, p_cur, p_prev, u_prev
    real(kind=rp) :: rnorm, rtr, reduction(3), norm_fac
    real(kind=rp) :: alpha(PIPECG_P_SPACE), beta(PIPECG_P_SPACE)
    real(kind=rp) :: gamma1, gamma2, delta
    real(kind=rp) :: tmp1, tmp2, tmp3, x_plus(NEKO_BLK_SIZE)
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
      
      p_prev = PIPECG_P_SPACE
      u_prev = PIPECG_P_SPACE+1
      p_cur = 1
      call rzero(x%x, n)
      call rzero(z, n)
      call rzero(q, n)
      call rzero(p, n)
      call rzero(s, n)
      call copy(r, f, n)
      call this%M%solve(u(1,u_prev), r, n)
      call Ax%compute(w, u(1,u_prev), coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD)
      call bc_list_apply(blst, w, n)
      
      rtr = glsc3(r, coef%mult, r, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(rnorm .eq. 0.0_rp) return
      
      gamma1 = 0.0_rp
      tmp1 = 0.0_rp
      tmp2 = 0.0_rp
      tmp3 = 0.0_rp
      do i = 1, n
         tmp1 = tmp1 + r(i) * coef%mult(i,1,1,1) * u(i,u_prev)
         tmp2 = tmp2 + w(i) * coef%mult(i,1,1,1) * u(i,u_prev)
         tmp3 = tmp3 + r(i) * coef%mult(i,1,1,1) * r(i)
      end do
      reduction(1) = tmp1
      reduction(2) = tmp2
      reduction(3) = tmp3
      
      do iter = 1, max_iter
         call MPI_Iallreduce(MPI_IN_PLACE, reduction, 3, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, request, ierr)
         
         call this%M%solve(mi, w, n)
         call Ax%compute(ni, mi, coef, x%msh, x%Xh)
         call gs_h%op(ni, n, GS_OP_ADD)
         call bc_list_apply(blst, ni, n)
         
         call MPI_Wait(request, status, ierr)
         gamma2 = gamma1       
         gamma1 = reduction(1)
         delta = reduction(2)
         rtr = reduction(3)
         
         rnorm = sqrt(rtr)*norm_fac
         if (rnorm .lt. this%abs_tol) exit

         if (iter .gt. 1) then
            beta(p_cur) = gamma1 / gamma2
            alpha(p_cur) = gamma1 / (delta - (beta(p_cur) * gamma1/alpha(p_prev)))
         else 
            beta(p_cur) = 0.0_rp
            alpha(p_cur) = gamma1/delta
         end if
         
         tmp1 = 0.0_rp
         tmp2 = 0.0_rp
         tmp3 = 0.0_rp
         do i = 0, n, NEKO_BLK_SIZE
            if (i + NEKO_BLK_SIZE .le. n) then
               do k = 1, NEKO_BLK_SIZE
                  z(i+k) = beta(p_cur) * z(i+k) + ni(i+k)
                  q(i+k) = beta(p_cur) * q(i+k) + mi(i+k)
                  s(i+k) = beta(p_cur) * s(i+k) + w(i+k)
                  r(i+k) =  r(i+k) - alpha(p_cur) * s(i+k)
                  u(i+k,p_cur) =  u(i+k,u_prev) - alpha(p_cur) * q(i+k)
                  w(i+k) =  w(i+k) - alpha(p_cur) * z(i+k)
                  tmp1 = tmp1 + r(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp2 = tmp2 + w(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp3 = tmp3 + r(i+k) * coef%mult(i+k,1,1,1) * r(i+k)
               end do
            else
               do k = 1, n-i
                  z(i+k) = beta(p_cur) * z(i+k) + ni(i+k)
                  q(i+k) = beta(p_cur) * q(i+k) + mi(i+k)
                  s(i+k) = beta(p_cur) * s(i+k) + w(i+k)
                  r(i+k) =  r(i+k) - alpha(p_cur) * s(i+k)
                  u(i+k,p_cur) =  u(i+k,u_prev) - alpha(p_cur) * q(i+k)
                  w(i+k) =  w(i+k) - alpha(p_cur) * z(i+k)
                  tmp1 = tmp1 + r(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp2 = tmp2 + w(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp3 = tmp3 + r(i+k) * coef%mult(i+k,1,1,1) * r(i+k)
               end do
            end if
         end do
         
         reduction(1) = tmp1
         reduction(2) = tmp2
         reduction(3) = tmp3
         
         if (p_cur .eq. PIPECG_P_SPACE) then
            do i = 0, n, NEKO_BLK_SIZE
               if (i + NEKO_BLK_SIZE .le. n) then
                  do k = 1, NEKO_BLK_SIZE
                     x_plus(k) = 0.0_rp
                  end do
                  p_prev = PIPECG_P_SPACE+1
                  do j = 1, p_cur
                     do k = 1, NEKO_BLK_SIZE
                        p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                        x_plus(k) = x_plus(k) + alpha(j) * p(i+k)
                     end do
                     p_prev = j
                  end do
                  do k = 1, NEKO_BLK_SIZE
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                     u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
                  end do
               else 
                  do k = 1, n-i
                     x_plus(1) = 0.0_rp
                     p_prev = PIPECG_P_SPACE + 1
                     do j = 1, p_cur
                        p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                        x_plus(1) = x_plus(1) + alpha(j) * p(i+k)
                        p_prev = j
                     end do
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                     u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
                  end do
               end if
            end do
            p_prev = p_cur
            u_prev = PIPECG_P_SPACE+1
            alpha(1) = alpha(p_cur) 
            beta(1) = beta(p_cur)
            p_cur = 1
         else
            u_prev = p_cur
            p_prev = p_cur
            p_cur = p_cur + 1
         end if
      end do
      
      if ( p_cur .ne. 1) then
         do i = 0, n, NEKO_BLK_SIZE
            if (i + NEKO_BLK_SIZE .le. n) then
               do k = 1, NEKO_BLK_SIZE
                  x_plus(k) = 0.0_rp
               end do
               p_prev = PIPECG_P_SPACE+1
               do j = 1, p_cur
                  do k = 1, NEKO_BLK_SIZE
                     p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                     x_plus(k) = x_plus(k) + alpha(j) * p(i+k)
                  end do
                  p_prev = j
               end do
               do k = 1, NEKO_BLK_SIZE
                  x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                  u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
               end do
            else 
               do k = 1, n-i
                  x_plus(1) = 0.0_rp
                  p_prev = PIPECG_P_SPACE + 1
                  do j = 1, p_cur
                     p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                     x_plus(1) = x_plus(1) + alpha(j) * p(i+k)
                     p_prev = j
                  end do
                  x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                  u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
               end do
            end if
         end do
      end if
      
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
      
    end associate
    
  end function pipecg_solve
   
end module pipecg
  

