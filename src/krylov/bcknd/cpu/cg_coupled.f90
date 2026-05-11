! Copyright (c) 2024, The Neko Authors
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
!> Defines a coupled Conjugate Gradient methods
module cg_cpld
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use scratch_registry, only : neko_scratch_registry
  use field_math, only : field_copy
  use math, only : glsc2, abscmp
  use utils, only : neko_error
  use operators, only : rotate_cyc
  implicit none
  private

  !> Coupled preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_cpld_t
     real(kind=rp), allocatable :: w1(:)
     real(kind=rp), allocatable :: w2(:)
     real(kind=rp), allocatable :: w3(:)
     real(kind=rp), allocatable :: r1(:)
     real(kind=rp), allocatable :: r2(:)
     real(kind=rp), allocatable :: r3(:)
     real(kind=rp), allocatable :: p1(:)
     real(kind=rp), allocatable :: p2(:)
     real(kind=rp), allocatable :: p3(:)
     real(kind=rp), allocatable :: z1(:)
     real(kind=rp), allocatable :: z2(:)
     real(kind=rp), allocatable :: z3(:)
     real(kind=rp), allocatable :: tmp(:)
   contains
     procedure, pass(this) :: init => cg_cpld_init
     procedure, pass(this) :: free => cg_cpld_free
     procedure, pass(this) :: solve => cg_cpld_nop
     procedure, pass(this) :: solve_coupled => cg_cpld_solve
  end type cg_cpld_t

contains

  !> Initialise a coupled PCG solver
  subroutine cg_cpld_init(this, n, max_iter, M, abs_tol, rel_tol, monitor)
    class(cg_cpld_t), target, intent(inout) :: this
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: abs_tol
    real(kind=rp), intent(in) :: rel_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%w1(n))
    allocate(this%w2(n))
    allocate(this%w3(n))
    allocate(this%r1(n))
    allocate(this%r2(n))
    allocate(this%r3(n))
    allocate(this%p1(n))
    allocate(this%p2(n))
    allocate(this%p3(n))
    allocate(this%z1(n))
    allocate(this%z2(n))
    allocate(this%z3(n))
    allocate(this%tmp(n))

    if (present(M)) then
       this%M => M
    end if

    if (present(monitor)) then
       call this%ksp_init(max_iter, abs_tol, rel_tol, monitor = monitor)
    else
       call this%ksp_init(max_iter, abs_tol, rel_tol)
    end if

  end subroutine cg_cpld_init

  !> Deallocate a coupled PCG solver
  subroutine cg_cpld_free(this)
    class(cg_cpld_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w1)) then
       deallocate(this%w1)
    end if

    if (allocated(this%w2)) then
       deallocate(this%w2)
    end if

    if (allocated(this%w3)) then
       deallocate(this%w3)
    end if

    if (allocated(this%r1)) then
       deallocate(this%r1)
    end if

    if (allocated(this%r2)) then
       deallocate(this%r2)
    end if

    if (allocated(this%r3)) then
       deallocate(this%r3)
    end if

    if (allocated(this%p1)) then
       deallocate(this%p1)
    end if

    if (allocated(this%p2)) then
       deallocate(this%p2)
    end if

    if (allocated(this%p3)) then
       deallocate(this%p3)
    end if

    if (allocated(this%z1)) then
       deallocate(this%z1)
    end if

    if (allocated(this%z2)) then
       deallocate(this%z2)
    end if

    if (allocated(this%z3)) then
       deallocate(this%z3)
    end if

    if (allocated(this%tmp)) then
       deallocate(this%tmp)
    end if

    nullify(this%M)

  end subroutine cg_cpld_free

  function cg_cpld_nop(this, Ax, x, f, n, coef, blst, gs_h, niter, ref) &
       result(ksp_results)
    class(cg_cpld_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    type(field_t), optional, intent(in) :: ref

    ! Throw and error
    call neko_error('The cpldcg solver is only defined for coupled solves')

    ksp_results%res_final = 0.0
    ksp_results%iter = 0
  end function cg_cpld_nop

  !> Coupled PCG solve
  function cg_cpld_solve(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter, refx, refy, refz) &
       result(ksp_results)
    class(cg_cpld_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter
    type(field_t), optional, intent(in) :: refx
    type(field_t), optional, intent(in) :: refy
    type(field_t), optional, intent(in) :: refz
    integer :: i, iter, max_iter, weight_idx, refx_idx, refy_idx, refz_idx
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=rp) :: beta, pap, alpha, alphm, norm_scale
    type(field_t), pointer :: weight, refx_, refy_, refz_

     if (present(niter)) then
        max_iter = niter
     else
        max_iter = this%max_iter
     end if
     call neko_scratch_registry%request_field(weight, weight_idx, .false.)
     call neko_scratch_registry%request_field(refx_, refx_idx, .false.)
     call neko_scratch_registry%request_field(refy_, refy_idx, .false.)
     call neko_scratch_registry%request_field(refz_, refz_idx, .false.)
     call this%residual%compute_weight(weight, coef, n)
     if (present(refx)) then
        call field_copy(refx_, refx, n)
     else
        call field_copy(refx_, x, n)
     end if
     if (present(refy)) then
        call field_copy(refy_, refy, n)
     else
        call field_copy(refy_, y, n)
     end if
     if (present(refz)) then
        call field_copy(refz_, refz, n)
     else
        call field_copy(refz_, z, n)
     end if
     select case (trim(this%residual%normalization_type))
     case ('none')
        norm_scale = 1.0_rp
     case ('volume')
        norm_scale = sqrt(coef%volume)
     case default
        norm_scale = sqrt( &
             this%residual%compute_normalization(Ax, refx_, fx, coef, gs_h, &
             blstx, weight, n)**2 + &
             this%residual%compute_normalization(Ax, refy_, fy, coef, gs_h, &
             blsty, weight, n)**2 + &
             this%residual%compute_normalization(Ax, refz_, fz, coef, gs_h, &
             blstz, weight, n)**2)
     end select

     associate (p1 => this%p1, p2 => this%p2, p3 => this%p3, z1 => this%z1, &
          z2 => this%z2, z3 => this%z3, r1 => this%r1, r2 => this%r2, &
         r3 => this%r3, tmp => this%tmp, w1 => this%w1, w2 => this%w2, &
         w3 => this%w3)

      rtz1 = 1.0_rp
      do concurrent (i = 1:n)
         x%x(i,1,1,1) = 0.0_rp
         y%x(i,1,1,1) = 0.0_rp
         z%x(i,1,1,1) = 0.0_rp
         p1(i) = 0.0_rp
         p2(i) = 0.0_rp
         p3(i) = 0.0_rp
         z1(i) = 0.0_rp
         z2(i) = 0.0_rp
         z3(i) = 0.0_rp
         r1(i) = fx(i)
         r2(i) = fy(i)
         r3(i) = fz(i)
         tmp(i) = r1(i)**2 + r2(i)**2 + r3(i)**2
      end do

       rtr = glsc2(tmp, weight%x, n)
       rnorm = sqrt(rtr) / norm_scale
       ksp_results%res_start = rnorm
       ksp_results%res_final = rnorm
       ksp_results%iter = 0
       if (abscmp(rnorm, 0.0_rp)) then
          ksp_results%converged = .true.
          call neko_scratch_registry%relinquish_field(refz_idx)
          call neko_scratch_registry%relinquish_field(refy_idx)
          call neko_scratch_registry%relinquish_field(refx_idx)
          call neko_scratch_registry%relinquish_field(weight_idx)
          return
       end if

      call this%monitor_start('cpldCG')
      do iter = 1, max_iter
         call this%M%solve(z1, this%r1, n)
         call this%M%solve(z2, this%r2, n)
         call this%M%solve(z3, this%r3, n)
         rtz2 = rtz1

         do concurrent (i = 1:n)
            this%tmp(i) = z1(i) * r1(i) &
                 + z2(i) * r2(i) &
                 + z3(i) * r3(i)
         end do

          rtz1 = glsc2(tmp, weight%x, n)

         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         do concurrent (i = 1:n)
            p1(i) = p1(i) * beta + z1(i)
            p2(i) = p2(i) * beta + z2(i)
            p3(i) = p3(i) * beta + z3(i)
         end do

         call Ax%compute_vector(w1, w2, w3, p1, p2, p3, coef, x%msh, x%Xh)

         call rotate_cyc(w1, w2, w3, 1, coef)
         call gs_h%op(w1, n, GS_OP_ADD)
         call gs_h%op(w2, n, GS_OP_ADD)
         call gs_h%op(w3, n, GS_OP_ADD)
         call rotate_cyc(w1, w2, w3, 0, coef)

         call blstx%apply_scalar(w1, n)
         call blsty%apply_scalar(w2, n)
         call blstz%apply_scalar(w3, n)

         do concurrent (i = 1:n)
            tmp(i) = w1(i) * p1(i) &
                 + w2(i) * p2(i) &
                 + w3(i) * p3(i)
         end do

          pap = glsc2(tmp, weight%x, n)

         alpha = rtz1 / pap
         alphm = -alpha
         do concurrent (i = 1:n)
            x%x(i,1,1,1) = x%x(i,1,1,1) + alpha * p1(i)
            y%x(i,1,1,1) = y%x(i,1,1,1) + alpha * p2(i)
            z%x(i,1,1,1) = z%x(i,1,1,1) + alpha * p3(i)
            r1(i) = r1(i) + alphm * w1(i)
            r2(i) = r2(i) + alphm * w2(i)
            r3(i) = r3(i) + alphm * w3(i)
            tmp(i) = r1(i)**2 + r2(i)**2 + r3(i)**2
         end do

          rtr = glsc2(tmp, weight%x, n)
          if (iter .eq. 1) rtr0 = rtr
          rnorm = sqrt(rtr) / norm_scale
          call this%monitor_iter(iter, rnorm)
          if (rnorm .lt. this%abs_tol) then
             exit
          end if
       end do
     end associate
     call neko_scratch_registry%relinquish_field(refz_idx)
     call neko_scratch_registry%relinquish_field(refy_idx)
     call neko_scratch_registry%relinquish_field(refx_idx)
     call neko_scratch_registry%relinquish_field(weight_idx)
     call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)
  end function cg_cpld_solve

end module cg_cpld
