! Copyright (c) 2024-2025, The Neko Authors
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
!> Implements smoothers for use with TreeAMG matrix vector product
module tree_amg_smoother
  use tree_amg, only : tamg_hierarchy_t
  use tree_amg_utils, only : tamg_sample_matrix_val
  use num_types, only : rp
  use math, only : col2, add2, add2s2, glsc2, glsc3, sub2, cmult, &
       cmult2, copy, add3s2
  use device_math, only : device_glsc2, device_glsc3, device_rzero, &
       device_cmult2, device_sub2, device_add2, device_add3s2, &
       device_copy
  use krylov, only : ksp_monitor_t
  use bc_list, only: bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use logger, only : neko_log, LOG_SIZE
  use device, only: device_map, device_free, device_memcpy, HOST_TO_DEVICE, &
       device_deassociate
  use device_tree_amg_smoother, only : amg_device_cheby_solve_part1, &
       amg_device_cheby_solve_part2
  use neko_config, only: NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Type for Chebyshev iteration using TreeAMG matvec
  type, public :: amg_jacobi_t
     real(kind=rp), allocatable :: d(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp) :: omega
     integer :: lvl
     integer :: n
     integer :: max_iter = 10
     logical :: recompute_diag = .true.
   contains
     procedure, pass(this) :: init => amg_jacobi_init
     procedure, pass(this) :: solve => amg_jacobi_solve
     procedure, pass(this) :: comp_diag => amg_jacobi_diag
     procedure, pass(this) :: free => amg_jacobi_free
  end type amg_jacobi_t

  !> Type for Chebyshev iteration using TreeAMG matvec
  type, public :: amg_cheby_t
     real(kind=rp), allocatable :: d(:)
     type(c_ptr) :: d_d = C_NULL_PTR
     real(kind=rp), allocatable :: w(:)
     type(c_ptr) :: w_d = C_NULL_PTR
     real(kind=rp), allocatable :: r(:)
     type(c_ptr) :: r_d = C_NULL_PTR
     real(kind=rp) :: tha, dlt
     integer :: lvl
     integer :: n
     integer :: power_its = 250
     integer :: max_iter = 10
     logical :: recompute_eigs = .true.
   contains
     procedure, pass(this) :: init => amg_cheby_init
     procedure, pass(this) :: solve => amg_cheby_solve
     procedure, pass(this) :: comp_eig => amg_cheby_power
     procedure, pass(this) :: device_solve => amg_device_cheby_solve
     procedure, pass(this) :: device_comp_eig => amg_device_cheby_power
     procedure, pass(this) :: free => amg_cheby_free
  end type amg_cheby_t

contains

  !> Initialization of chebyshev
  !! @param n Number of dofs
  !! @param lvl The tamg hierarchy level on which the iterations are to be applied
  !! @param max_iter The number of iterations (chebyshev degree)
  subroutine amg_cheby_init(this, n, lvl, max_iter)
    class(amg_cheby_t), intent(inout), target :: this
    integer, intent(in) :: n
    integer, intent(in) :: lvl
    integer, intent(in) :: max_iter

    allocate(this%d(n))
    allocate(this%w(n))
    allocate(this%r(n))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%d, this%d_d, n)
       call device_map(this%w, this%w_d, n)
       call device_map(this%r, this%r_d, n)
    end if
    this%n = n
    this%lvl = lvl
    this%max_iter = max_iter
    this%recompute_eigs = .true.

    call amg_smoo_monitor(lvl, this)

  end subroutine amg_cheby_init

  !> free cheby data
  subroutine amg_cheby_free(this)
    class(amg_cheby_t), intent(inout), target :: this
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_deassociate(this%d)
       call device_deassociate(this%w)
       call device_deassociate(this%r)
    end if
    if (allocated(this%d)) then
       deallocate(this%d)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%r)) then
       deallocate(this%r)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_free(this%d_d)
       call device_free(this%w_d)
       call device_free(this%r_d)
    end if
  end subroutine amg_cheby_free


  !> Power method to approximate largest eigenvalue
  !! @param amg TreeAMG object
  !! @param n Number of dofs
  subroutine amg_cheby_power(this, amg, n)
    class(amg_cheby_t), intent(inout) :: this
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: n
    real(kind=rp) :: lam, b, a, rn
    real(kind=rp), parameter :: boost = 1.1_rp
    real(kind=rp), parameter :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer :: i
    associate(w => this%w, d => this%d, coef => amg%coef, gs_h => amg%gs_h, &
         msh => amg%msh, Xh => amg%Xh, blst => amg%blst)

      do i = 1, n
         !call random_number(rn)
         !d(i) = rn + 10.0_rp
         d(i) = sin(real(i))
      end do
      if (this%lvl .eq. 0) then
         call gs_h%op(d, n, GS_OP_ADD)!TODO
         call blst%apply(d, n)
      end if
      !Power method to get lamba max
      do i = 1, this%power_its
         call amg%matvec(w, d, this%lvl)

         if (this%lvl .eq. 0) then
            wtw = glsc3(w, coef%mult, w, n)
         else
            wtw = glsc2(w, w, n)
         end if

         call cmult2(d, w, 1.0_rp/sqrt(wtw), n)
      end do

      call amg%matvec(w, d, this%lvl)

      if (this%lvl .eq. 0) then
         dtw = glsc3(d, coef%mult, w, n)
         dtd = glsc3(d, coef%mult, d, n)
      else
         dtw = glsc2(d, w, n)
         dtd = glsc2(d, d, n)
      end if
      lam = dtw / dtd
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0_rp
      this%dlt = (b-a)/2.0_rp

      this%recompute_eigs = .false.
      call amg_cheby_monitor(this%lvl, lam)
    end associate
  end subroutine amg_cheby_power

  !> Chebyshev smoother
  !> From Saad's iterative methods textbook
  !! @param x The solution to be returned
  !! @param f The right-hand side
  !! @param n Number of dofs
  !! @param amg The TreeAMG object
  subroutine amg_cheby_solve(this, x, f, n, amg, zero_init)
    class(amg_cheby_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(inout) :: f
    class(tamg_hierarchy_t), intent(inout) :: amg
    type(ksp_monitor_t) :: ksp_results
    logical, optional, intent(in) :: zero_init
    integer :: iter, max_iter, i
    real(kind=rp) :: rtr, rnorm
    real(kind=rp) :: rhok, rhokp1, s1, thet, delt, tmp1, tmp2
    logical :: zero_initial_guess

    if (this%recompute_eigs) then
       call this%comp_eig(amg, n)
    end if
    if (present(zero_init)) then
       zero_initial_guess = zero_init
    else
       zero_initial_guess = .false.
    end if
    max_iter = this%max_iter

    associate( w => this%w, r => this%r, d => this%d, blst => amg%blst)
      call copy(r, f, n)
      if (.not. zero_initial_guess) then
         call amg%matvec(w, x, this%lvl)
         call sub2(r, w, n)
      end if

      thet = this%tha
      delt = this%dlt
      s1 = thet / delt
      rhok = 1.0_rp / s1

      ! First iteration
      do concurrent (i = 1:n)
         d(i) = 1.0_rp/thet * r(i)
         x(i) = x(i) + d(i)
      end do

      ! Rest of iterations
      do iter = 2, max_iter
         call amg%matvec(w, d, this%lvl)

         rhokp1 = 1.0_rp / (2.0_rp * s1 - rhok)
         tmp1 = rhokp1 * rhok
         tmp2 = 2.0_rp * rhokp1 / delt
         rhok = rhokp1

         do concurrent (i = 1:n)
            r(i) = r(i) - w(i)
            d(i) = tmp1 * d(i) + tmp2 * r(i)
            x(i) = x(i) + d(i)
         end do

      end do
    end associate
  end subroutine amg_cheby_solve

  !> Power method to approximate largest eigenvalue
  !! @param amg TreeAMG object
  !! @param n Number of dofs
  subroutine amg_device_cheby_power(this, amg, n)
    class(amg_cheby_t), intent(inout) :: this
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: n
    real(kind=rp) :: lam, b, a, rn
    real(kind=rp), parameter :: boost = 1.1_rp
    real(kind=rp), parameter :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer :: i
    associate(w => this%w, d => this%d, coef => amg%coef, gs_h => amg%gs_h, &
         msh => amg%msh, Xh => amg%Xh, blst => amg%blst)
      do i = 1, n
         !TODO: replace with a better way to initialize power method
         d(i) = sin(real(i))
      end do
      call device_memcpy(this%d, this%d_d, n, HOST_TO_DEVICE, .true.)
      if (this%lvl .eq. 0) then
         call gs_h%op(d, n, GS_OP_ADD)!TODO
         call blst%apply(d, n)
      end if
      do i = 1, this%power_its
         call amg%device_matvec(w, d, this%w_d, this%d_d, this%lvl)

         if (this%lvl .eq. 0) then
            wtw = device_glsc3(this%w_d, coef%mult_d, this%w_d, n)
         else
            wtw = device_glsc2(this%w_d, this%w_d, n)
         end if

         call device_cmult2(this%d_d, this%w_d, 1.0_rp/sqrt(wtw), n)
      end do

      call amg%device_matvec(w, d, this%w_d, this%d_d, this%lvl)

      if (this%lvl .eq. 0) then
         dtw = device_glsc3(this%d_d, coef%mult_d, this%w_d, n)
         dtd = device_glsc3(this%d_d, coef%mult_d, this%d_d, n)
      else
         dtw = device_glsc2(this%d_d, this%w_d, n)
         dtd = device_glsc2(this%d_d, this%d_d, n)
      end if
      lam = dtw / dtd
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0_rp
      this%dlt = (b-a)/2.0_rp

      this%recompute_eigs = .false.
      call amg_cheby_monitor(this%lvl, lam)
    end associate
  end subroutine amg_device_cheby_power

  !> Chebyshev smoother
  !> From Saad's iterative methods textbook
  !! @param x The solution to be returned
  !! @param f The right-hand side
  !! @param n Number of dofs
  !! @param amg The TreeAMG object
  subroutine amg_device_cheby_solve(this, x, f, x_d, f_d, n, amg, zero_init)
    class(amg_cheby_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(inout) :: f
    type(c_ptr) :: x_d
    type(c_ptr) :: f_d
    class(tamg_hierarchy_t), intent(inout) :: amg
    type(ksp_monitor_t) :: ksp_results
    logical, optional, intent(in) :: zero_init
    integer :: iter, max_iter
    real(kind=rp) :: rtr, rnorm
    real(kind=rp) :: rhok, rhokp1, s1, thet, delt, tmp1, tmp2
    logical :: zero_initial_guess

    if (this%recompute_eigs) then
       call this%device_comp_eig(amg, n)
    end if
    if (present(zero_init)) then
       zero_initial_guess = zero_init
    else
       zero_initial_guess = .false.
    end if
    max_iter = this%max_iter

    associate( w_d => this%w_d, r_d => this%r_d, d_d => this%d_d, &
         blst => amg%blst)

      if (.not. zero_initial_guess) then
         call amg%device_matvec(this%w, x, w_d, x_d, this%lvl)
      end if

      thet = this%tha
      delt = this%dlt
      s1 = thet / delt
      rhok = 1.0_rp / s1

      ! First iteration
      tmp1 = 1.0_rp / thet
      call amg_device_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
           tmp1, n, zero_initial_guess)
      ! Rest of iterations
      do iter = 2, max_iter
         call amg%device_matvec(this%w, this%d, w_d, d_d, this%lvl)

         rhokp1 = 1.0_rp / (2.0_rp * s1 - rhok)
         tmp1 = rhokp1 * rhok
         tmp2 = 2.0_rp * rhokp1 / delt
         rhok = rhokp1

         call amg_device_cheby_solve_part2(r_d, w_d, d_d, x_d, tmp1, tmp2, n)

      end do
    end associate
  end subroutine amg_device_cheby_solve

  !> Initialization of Jacobi (this is expensive...)
  !! @param n Number of dofs
  !! @param lvl The tamg hierarchy level on which the iterations are to be applied
  !! @param max_iter The number of iterations
  subroutine amg_jacobi_init(this, n, lvl, max_iter)
    class(amg_jacobi_t), intent(inout), target :: this
    integer, intent(in) :: n
    integer, intent(in) :: lvl
    integer, intent(in) :: max_iter

    allocate(this%d(n))
    allocate(this%w(n))
    allocate(this%r(n))
    this%n = n
    this%lvl = lvl
    this%max_iter = max_iter
    this%omega = 0.7_rp

  end subroutine amg_jacobi_init

  !> free jacobi data
  subroutine amg_jacobi_free(this)
    class(amg_jacobi_t), intent(inout), target :: this
    if (allocated(this%d)) then
       deallocate(this%d)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%r)) then
       deallocate(this%r)
    end if
  end subroutine amg_jacobi_free

  !> SAMPLE MATRIX DIAGONAL VALUES (DO NOT USE, EXPENSIVE)
  !! @param amg TreeAMG object
  !! @param n Number of dofs
  subroutine amg_jacobi_diag(this, amg, n)
    class(amg_jacobi_t), intent(inout) :: this
    type(tamg_hierarchy_t), intent(inout) :: amg
    integer, intent(in) :: n
    real(kind=rp) :: val
    integer :: i
    do i = 1, n
       call tamg_sample_matrix_val(val, amg, this%lvl, i, i)
       this%d(i) = 1.0_rp / val
    end do
    this%recompute_diag = .false.
  end subroutine amg_jacobi_diag

  !> Jacobi smoother
  !! @param x The solution to be returned
  !! @param f The right-hand side
  !! @param n Number of dofs
  !! @param amg The TreeAMG object
  subroutine amg_jacobi_solve(this, x, f, n, amg, niter)
    class(amg_jacobi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(inout) :: f
    class(tamg_hierarchy_t), intent(inout) :: amg
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: rtr, rnorm
    integer :: i

    if (this%recompute_diag) then
       call this%comp_diag(amg, n)
    end if

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if

    ! x = x + omega * Dinv( f - Ax )
    associate( w => this%w, r => this%r, d => this%d)
      do iter = 1, max_iter
         w = 0.0_rp
         !> w = A x
         call amg%matvec(w, x, this%lvl)
         !> r = f - Ax
         call copy(r, f, n)
         call sub2(r, w, n)
         !> r = Dinv * (f - Ax)
         call col2(r, d, n)
         !> x = x + omega * Dinv * (f - Ax)
         call add2s2(x, r, this%omega, n)
      end do
    end associate
  end subroutine amg_jacobi_solve

  subroutine amg_smoo_monitor(lvl, smoo)
    integer, intent(in) :: lvl
    class(amg_cheby_t), intent(in) :: smoo
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A8,I2,A28)') '-- level', lvl, '-- init smoother: Chebyshev'
    call neko_log%message(log_buf)
    write(log_buf, '(A22,I6)') 'Iterations:', smoo%max_iter
    call neko_log%message(log_buf)
  end subroutine amg_smoo_monitor

  subroutine amg_cheby_monitor(lvl, lam)
    integer, intent(in) :: lvl
    real(kind=rp), intent(in) :: lam
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A12,I2,A29,F12.3)') '-- AMG level', lvl, &
         '-- Chebyshev approx. max eig', lam
    call neko_log%message(log_buf)
  end subroutine amg_cheby_monitor

end module tree_amg_smoother
