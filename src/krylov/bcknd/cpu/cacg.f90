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
!> Defines a communication avoiding Conjugate Gradient method
module cacg
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply, bc_list_apply_scalar
  use math, only : glsc3, rzero, copy, x_update
  use utils, only : neko_warning
  use comm
  use mxm_wrapper
  implicit none
  private
  
  !> S-step communication avoiding  preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cacg_t
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: PR(:,:)
     integer :: s
   contains
     procedure, pass(this) :: init => cacg_init
     procedure, pass(this) :: free => cacg_free
     procedure, pass(this) :: solve => cacg_solve
  end type cacg_t

contains

  !> Initialise a s-step CA  PCG solver
  subroutine cacg_init(this, n, M, s, rel_tol, abs_tol)
    class(cacg_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    integer, optional, intent(inout) :: s
    call this%free()

    if (present(s)) then
       this%s = s
    else 
       this%s = 4
    end if
    if (pe_rank .eq. 0) then
       call neko_warning("Communication Avoiding CG chosen, be aware of potential instabilities")
    end if
    
    allocate(this%r(n))
    allocate(this%p(n))
    allocate(this%PR(n,4*this%s+1))
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
          
  end subroutine cacg_init

  !> Deallocate a s-step CA PCG solver
  subroutine cacg_free(this)
    class(cacg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%PR)) then
       deallocate(this%PR)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if
    
    nullify(this%M)


  end subroutine cacg_free
  
  !> S-step CA PCG solve
  function cacg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cacg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: i, j, k, l, iter, max_iter, s, ierr, it
    real(kind=rp) :: rnorm, rtr, rtz1, tmp
    real(kind=rp) :: beta(this%s+1), alpha(this%s+1), alpha1, alpha2, norm_fac
    real(kind=rp), dimension(4*this%s+1,4*this%s+1) :: Tt, G, GTt, temp, temp2
    real(kind=rp) :: p_c(4*this%s+1,this%s+1)
    real(kind=rp) :: r_c(4*this%s+1,this%s+1)
    real(kind=rp) :: z_c(4*this%s+1,this%s+1)
    real(kind=rp) :: x_c(4*this%s+1,this%s+1)
    
    associate(PR => this%PR, r => this%r, p => this%p)
      s = this%s
      if (present(niter)) then
         max_iter = niter
      else
         max_iter = KSP_MAX_ITER
      end if
      norm_fac = 1.0_rp / sqrt(coef%volume)
      
      rtz1 = 1.0_rp
      call rzero(x%x, n)
      call copy(r, f, n)
      call this%M%solve(p, r, n)
      
      rtr = glsc3(r, coef%mult, r, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      iter = 0
      if(rnorm .eq. 0.0_rp) return
      do while (iter < max_iter)

         call copy(PR,p, n)
         call copy(PR(1,2*s+2), r, n)

         !Here we have hardcoded a monomial basis atm. 
         do i = 2, 2*s + 1
            if (mod(i,2) .eq. 0) then
               call Ax%compute(PR(1,i), PR(1,i-1), coef, x%msh, x%Xh)
               call gs_h%gs_op_vector(PR(1,i), n, GS_OP_ADD)
               call bc_list_apply_scalar(blst, PR(1,i), n)
            else
               call this%M%solve(PR(1,i), PR(1,i-1), n)
            end if
         end do

         do i = 2*s+2, 4*s
            if (mod(i,2) == 0) then
               call this%M%solve(PR(1,i+1), PR(1,i), n)
            else
               call Ax%compute(PR(1,i+1), PR(1,i), coef, x%msh, x%Xh)
               call gs_h%gs_op_vector(PR(1,i+1), n, GS_OP_ADD)
               call bc_list_apply_scalar(blst, PR(1,1+i), n)
            end if
         end do

         call construct_basis_matrix(Tt, s)
         call rzero(p_c, (4*s+1) * (s+1))
         p_c(1,1) = 1.0_rp 
         call rzero(r_c, (4*s+1) * (s+1))
         r_c(2*s+2,1) = 1.0_rp
         call mxm(Tt, 4*s+1, r_c, 4*s+1, z_c,s+1)
         call rzero(x_c, (4*s+1) * (s+1))
         call rzero(temp, (4*s+1)**2)

         do i = 0, n, NEKO_BLK_SIZE
            it = 0
            if (i + NEKO_BLK_SIZE .le. n) then
               do j = 1, 4*s+1
                  do l = 1, j
                     it = it + 1
                     do k = 1, NEKO_BLK_SIZE
                        temp(it,1) = temp(it,1) &
                                   + PR(i+k,j) * PR(i+k,l) * coef%mult(i+k,1,1,1)
                     end do
                  end do
               end do
            else
               do j = 1, 4*s+1
                  do l = 1, j
                     it = it + 1
                     do k = 1, n-i
                        temp(it,1) = temp(it,1) &
                                   + PR(i+k,j) * PR(i+k,l) * coef%mult(i+k,1,1,1)
                     end do
                  end do
               end do
            end if
         end do

         call MPI_Allreduce(temp, temp2, it, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
         it = 0
         do j = 1, 4*s+1
            do k = 1, j
               it = it + 1
               G(j,k) = temp2(it,1)
               G(k,j) = temp2(it,1)
            end do
         end do
         
         call mxm(G,4*s+1, Tt, 4*s+1,GTt,4*s+1)
         
         do j = 1, s
            iter = iter + 1
          
            call mxm(G, 4*s+1, r_c(1,j), 4*s+1,temp, 1)
            call mxm(GTt, 4*s+1, p_c(1,j), 4*s+1,temp2, 1)
            alpha1 = 0.0_rp
            alpha2 = 0.0_rp
            do i = 1,4*s+1
               alpha1 = alpha1 + temp(i,1) * z_c(i,j)
               alpha2 = alpha2 + temp2(i,1) * p_c(i,j)
            end do
            alpha(j) = alpha1/alpha2

            do i = 1, 4*s+1
               x_c(i,j+1) = x_c(i,j) + alpha(j) * p_c(i,j)
               tmp = 0.0_rp
               do k = 1, 4*s+1
                  tmp = tmp + Tt(i,k) * p_c(k,j)
               end do
               r_c(i,j+1) =  r_c(i,j) - alpha(j)*tmp
               tmp = 0.0_rp
               do k = 1, 4*s+1
                  tmp = tmp + Tt(i,k)*r_c(k,j+1)
               end do
               z_c(i,j+1) = tmp
            end do

            call mxm(G,4*s+1,r_c(1,j+1),4*s+1,temp2,1)
            alpha2 = 0.0_rp
            do i = 1,4*s+1
               alpha2 = alpha2 + temp2(i,1)*z_c(i,j+1)
            end do
            beta(j) = alpha2 / alpha1
            do i = 1,4*s+1
               p_c(i,j+1) = z_c(i,j+1) +  beta(j)*p_c(i,j)
            end do
         end do

         call rzero(p, n)
         call rzero(r, n)
         rtr = 0.0_rp
         do i = 0, n, NEKO_BLK_SIZE
            if (i + NEKO_BLK_SIZE .le. n) then
               do j = 1, 4*s + 1
                  do k = 1, NEKO_BLK_SIZE
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + PR(i+k,j) * x_c(j,s+1)
                     p(i+k) = p(i+k) + PR(i+k,j) * p_c(j,s+1)
                     tmp = PR(i+k,j) * r_c(j,s+1)
                     r(i+k) = r(i+k) + tmp
                  end do
               end do
               do k = 1, NEKO_BLK_SIZE
                  rtr = rtr + r(i+k)**2 * coef%mult(i+k,1,1,1)
               end do
            else 
               do j = 1,4*s+1
                  do k = 1, n-i
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + PR(i+k,j) * x_c(j,s+1)
                     p(i+k) = p(i+k) + PR(i+k,j) * p_c(j,s+1)
                     tmp = PR(i+k,j) * r_c(j,s+1)
                     r(i+k) = r(i+k) + tmp
                  end do
               end do
               do k = 1, n-i
                  rtr = rtr + r(i+k)**2 * coef%mult(i+k,1,1,1)
               end do
            end if
         end do

         call MPI_Allreduce(rtr, tmp, 1, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
         rnorm = norm_fac*sqrt(tmp)
         if( rnorm <= this%abs_tol) exit 
      end do

      ksp_results%res_final = rnorm
      ksp_results%iter = iter

    end associate

  end function cacg_solve

  !> Monomial matrix constuction, not sparse
  subroutine construct_basis_matrix(Tt, s)
     integer, intent(in) :: s
     real(kind=rp), intent(inout) :: Tt(4*s+1,4*s+1)
     integer :: mlen, i
     mlen = (4*s+1)*(4*s+1)
     call rzero(Tt,mlen)
     do i = 1, 2*s 
        Tt(i+1,i) = 1.0_rp
     end do
     do i = 1, (2*s-1)
        Tt(2*s+2+i,2*s+1+i) = 1.0_rp
     end do
  end subroutine construct_basis_matrix

end module cacg
  

