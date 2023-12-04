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
!> Defines a pipelined Conjugate Gradient methods SX-Aurora backend
module pipecg_sx
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
  
  !> Pipelined preconditioned conjugate gradient method for SX-Aurora
  type, public, extends(ksp_t) :: sx_pipecg_t
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
     procedure, pass(this) :: init => sx_pipecg_init
     procedure, pass(this) :: free => sx_pipecg_free
     procedure, pass(this) :: solve => sx_pipecg_solve
  end type sx_pipecg_t

contains

  !> Initialise a pipelined PCG solver
  subroutine sx_pipecg_init(this, n, M, rel_tol, abs_tol)
    class(sx_pipecg_t), intent(inout) :: this
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
          
  end subroutine sx_pipecg_init

  !> Deallocate a pipelined PCG solver
  subroutine sx_pipecg_free(this)
    class(sx_pipecg_t), intent(inout) :: this

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


  end subroutine sx_pipecg_free
  
  !> Pipelined PCG solve
  function sx_pipecg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(sx_pipecg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
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

    do i = 1, n
       x%x(i,1,1,1) = 0.0_rp
       this%z(i) = 0.0_rp
       this%q(i) = 0.0_rp
       this%p(i) = 0.0_rp
       this%s(i) = 0.0_rp
       this%r(i) = f(i)
    end do

    call this%M%solve(this%u, this%r, n)
    call Ax%compute(this%w, this%u, coef, x%msh, x%Xh)
    call gs_h%op(this%w, n, GS_OP_ADD)
    call bc_list_apply(blst, this%w, n)
    
    rtr = glsc3(this%r, coef%mult, this%r, n)
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
          tmp1 = tmp1 + this%r(i) * coef%mult(i,1,1,1) * this%u(i)
          tmp2 = tmp2 + this%w(i) * coef%mult(i,1,1,1) * this%u(i)
          tmp3 = tmp3 + this%r(i) * coef%mult(i,1,1,1) * this%r(i)
       end do
       reduction(1) = tmp1
       reduction(2) = tmp2
       reduction(3) = tmp3
       
       call MPI_Iallreduce(MPI_IN_PLACE, reduction, 3, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, request, ierr)
       
       call this%M%solve(this%mi, this%w, n)
       call Ax%compute(this%ni, this%mi, coef, x%msh, x%Xh)
       call gs_h%op(this%ni, n, GS_OP_ADD)
       call bc_list_apply(blst, this%ni, n)

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
       
       do i = 1, n
          this%z(i) = beta * this%z(i) + this%ni(i)
          this%q(i) = beta * this%q(i) + this%mi(i)
          this%s(i) = beta * this%s(i) + this%w(i)
          this%p(i) = beta * this%p(i) + this%u(i)
       end do

       do i = 1, n
          x%x(i,1,1,1) = x%x(i,1,1,1) + alpha * this%p(i)
          this%r(i) = this%r(i) - alpha * this%s(i)
          this%u(i) = this%u(i) - alpha * this%q(i)
          this%w(i) = this%w(i) - alpha * this%z(i)
       end do
       
    end do
    
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    
  end function sx_pipecg_solve
   
end module pipecg_sx
  

