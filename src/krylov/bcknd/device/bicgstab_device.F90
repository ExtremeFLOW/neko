! Copyright (c) 2026, The Neko Authors
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
!> Provides a device implementation of the BiCGStab method.
module bicgstab_device
  use num_types, only : rp, xp, c_rp, c_xp
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t, &
       vector_bc_projector_components
  use math, only : NEKO_EPS
  use device, only : device_map, device_unmap, device_get_ptr, &
       device_event_create, device_event_destroy, device_event_sync, &
       glb_cmd_queue
  use device_math, only : device_rzero, device_copy, device_glsc3, &
       device_add2s2
  use utils, only : neko_error
  use comm, only : pe_size, NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated, &
       c_int
  implicit none
  private

  !> Device implementation of the right-preconditioned BiCGStab method.
  !!
  !! The recurrence and the breakdown criteria mirror the CPU implementation
  !! in [bicgstab_t](#bicgstab::bicgstab_t). The solve starts from a zero
  !! initial guess and uses the preconditioner provided by
  !! [ksp_t](#krylov::ksp_t). The coupled interface solves the three
  !! components independently and does not apply a coupled operator.
  !!
  !! Each vector update is fused with the reduction that follows it, so an
  !! iteration issues four kernels and four global reductions rather than
  !! thirteen and seven.
  type, public, extends(ksp_t) :: bicgstab_device_t
     !> Search direction \f$p\f$.
     real(kind=rp), allocatable :: p(:)
     !> Preconditioned search direction \f$\hat{p} = M^{-1}p\f$.
     real(kind=rp), allocatable :: p_hat(:)
     !> Residual \f$r\f$.
     real(kind=rp), allocatable :: r(:)
     !> Intermediate residual \f$s = r - \alpha v\f$.
     real(kind=rp), allocatable :: s(:)
     !> Preconditioned intermediate residual \f$\hat{s} = M^{-1}s\f$.
     real(kind=rp), allocatable :: s_hat(:)
     !> Operator action \f$t = A\hat{s}\f$.
     real(kind=rp), allocatable :: t(:)
     !> Operator action \f$v = A\hat{p}\f$.
     real(kind=rp), allocatable :: v(:)
     type(c_ptr) :: p_d = C_NULL_PTR
     type(c_ptr) :: p_hat_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: s_d = C_NULL_PTR
     type(c_ptr) :: s_hat_d = C_NULL_PTR
     type(c_ptr) :: t_d = C_NULL_PTR
     type(c_ptr) :: v_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     !> Initialise a device BiCGStab solver.
     procedure, pass(this) :: init => bicgstab_device_init
     !> Free a device BiCGStab solver.
     procedure, pass(this) :: free => bicgstab_device_free
     !> Solve a linear system with the device BiCGStab method.
     procedure, pass(this) :: solve => bicgstab_device_solve
     !> Solve three independent systems with the device BiCGStab method.
     procedure, pass(this) :: solve_coupled => bicgstab_device_solve_coupled
  end type bicgstab_device_t

#if HAVE_HIP
  interface
     subroutine hip_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n) &
          bind(c, name = 'hip_bicgstab_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, r_d, v_d
       real(c_rp) :: beta, omega
       integer(c_int) :: n
     end subroutine hip_bicgstab_update_p
  end interface

  interface
     subroutine hip_bicgstab_product_and_norm(a_d, b_d, mult_d, res, n) &
          bind(c, name = 'hip_bicgstab_product_and_norm')
       use, intrinsic :: iso_c_binding
       import c_xp
       implicit none
       type(c_ptr), value :: a_d, b_d, mult_d
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine hip_bicgstab_product_and_norm
  end interface

  interface
     real(c_xp) function hip_bicgstab_part1(s_d, r_d, v_d, mult_d, &
          alpha, n) bind(c, name = 'hip_bicgstab_part1')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: s_d, r_d, v_d, mult_d
       real(c_rp) :: alpha
       integer(c_int) :: n
     end function hip_bicgstab_part1
  end interface

  interface
     subroutine hip_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, &
          t_d, f_d, mult_d, alpha, omega, res, n) &
          bind(c, name = 'hip_bicgstab_part2')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d
       type(c_ptr), value :: mult_d
       real(c_rp) :: alpha, omega
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine hip_bicgstab_part2
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n) &
          bind(c, name = 'cuda_bicgstab_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, r_d, v_d
       real(c_rp) :: beta, omega
       integer(c_int) :: n
     end subroutine cuda_bicgstab_update_p
  end interface

  interface
     subroutine cuda_bicgstab_product_and_norm(a_d, b_d, mult_d, res, n) &
          bind(c, name = 'cuda_bicgstab_product_and_norm')
       use, intrinsic :: iso_c_binding
       import c_xp
       implicit none
       type(c_ptr), value :: a_d, b_d, mult_d
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine cuda_bicgstab_product_and_norm
  end interface

  interface
     real(c_xp) function cuda_bicgstab_part1(s_d, r_d, v_d, mult_d, &
          alpha, n) bind(c, name = 'cuda_bicgstab_part1')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: s_d, r_d, v_d, mult_d
       real(c_rp) :: alpha
       integer(c_int) :: n
     end function cuda_bicgstab_part1
  end interface

  interface
     subroutine cuda_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, &
          t_d, f_d, mult_d, alpha, omega, res, n) &
          bind(c, name = 'cuda_bicgstab_part2')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d
       type(c_ptr), value :: mult_d
       real(c_rp) :: alpha, omega
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine cuda_bicgstab_part2
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n, &
          strm) bind(c, name = 'opencl_bicgstab_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, r_d, v_d, strm
       real(c_rp) :: beta, omega
       integer(c_int) :: n
     end subroutine opencl_bicgstab_update_p
  end interface

  interface
     subroutine opencl_bicgstab_product_and_norm(a_d, b_d, mult_d, res, &
          n, strm) bind(c, name = 'opencl_bicgstab_product_and_norm')
       use, intrinsic :: iso_c_binding
       import c_xp
       implicit none
       type(c_ptr), value :: a_d, b_d, mult_d, strm
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine opencl_bicgstab_product_and_norm
  end interface

  interface
     real(c_xp) function opencl_bicgstab_part1(s_d, r_d, v_d, mult_d, &
          alpha, n, strm) bind(c, name = 'opencl_bicgstab_part1')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: s_d, r_d, v_d, mult_d, strm
       real(c_rp) :: alpha
       integer(c_int) :: n
     end function opencl_bicgstab_part1
  end interface

  interface
     subroutine opencl_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, &
          t_d, f_d, mult_d, alpha, omega, res, n, strm) &
          bind(c, name = 'opencl_bicgstab_part2')
       use, intrinsic :: iso_c_binding
       import c_rp, c_xp
       implicit none
       type(c_ptr), value :: x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d
       type(c_ptr), value :: mult_d, strm
       real(c_rp) :: alpha, omega
       real(c_xp) :: res(2)
       integer(c_int) :: n
     end subroutine opencl_bicgstab_part2
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n, &
          strm) bind(c, name = 'metal_bicgstab_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, r_d, v_d, strm
       real(c_rp) :: beta, omega
       integer(c_int) :: n
     end subroutine metal_bicgstab_update_p
  end interface

  interface
     subroutine metal_bicgstab_product_and_norm(a_d, b_d, mult_d, res, n, &
          strm) bind(c, name = 'metal_bicgstab_product_and_norm')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, mult_d, strm
       real(c_rp) :: res(2)
       integer(c_int) :: n
     end subroutine metal_bicgstab_product_and_norm
  end interface

  interface
     real(c_rp) function metal_bicgstab_part1(s_d, r_d, v_d, mult_d, &
          alpha, n, strm) bind(c, name = 'metal_bicgstab_part1')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: s_d, r_d, v_d, mult_d, strm
       real(c_rp) :: alpha
       integer(c_int) :: n
     end function metal_bicgstab_part1
  end interface

  interface
     subroutine metal_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, &
          t_d, f_d, mult_d, alpha, omega, res, n, strm) &
          bind(c, name = 'metal_bicgstab_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d
       type(c_ptr), value :: mult_d, strm
       real(c_rp) :: alpha, omega, res(2)
       integer(c_int) :: n
     end subroutine metal_bicgstab_part2
  end interface
#endif

contains

  !> Search direction update \f$p = r + \beta (p - \omega v)\f$.
  !! @param p_d Search direction, updated in place.
  !! @param r_d Residual.
  !! @param v_d Operator action \f$A\hat{p}\f$ of the previous iteration.
  !! @param beta BiCGStab \f$\beta\f$.
  !! @param omega BiCGStab \f$\omega\f$ of the previous iteration.
  !! @param n Number of vector entries.
  subroutine device_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n)
    type(c_ptr) :: p_d, r_d, v_d
    real(kind=rp) :: beta, omega
    integer :: n

#if HAVE_HIP
    call hip_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n)
#elif HAVE_CUDA
    call cuda_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n)
#elif HAVE_OPENCL
    call opencl_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n, &
         glb_cmd_queue)
#elif HAVE_METAL
    call metal_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n, &
         glb_cmd_queue)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_bicgstab_update_p

  !> Weighted inner product and squared norm in one reduction.
  !!
  !! Combining the two values avoids a global synchronisation solely for the
  !! scale used by the breakdown checks.
  !! @param product Weighted inner product of `a` and `b`.
  !! @param norm_squared Weighted squared norm of `b`.
  !! @param a_d First vector in the inner product.
  !! @param b_d Second vector in the inner product and vector whose norm is
  !! taken.
  !! @param mult_d Spectral element multiplicity weights.
  !! @param n Number of vector entries.
  subroutine device_bicgstab_product_and_norm(product, norm_squared, &
       a_d, b_d, mult_d, n)
    real(kind=rp), intent(out) :: product
    real(kind=rp), intent(out) :: norm_squared
    type(c_ptr) :: a_d, b_d, mult_d
    integer :: n
    real(kind=xp) :: res_xp(2)
    integer :: ierr
#if HAVE_METAL
    real(kind=rp) :: res(2)
#endif

    res_xp = 0.0_xp
#if HAVE_HIP
    call hip_bicgstab_product_and_norm(a_d, b_d, mult_d, res_xp, n)
#elif HAVE_CUDA
    call cuda_bicgstab_product_and_norm(a_d, b_d, mult_d, res_xp, n)
#elif HAVE_OPENCL
    call opencl_bicgstab_product_and_norm(a_d, b_d, mult_d, res_xp, n, &
         glb_cmd_queue)
#elif HAVE_METAL
    call metal_bicgstab_product_and_norm(a_d, b_d, mult_d, res, n, &
         glb_cmd_queue)
    res_xp = real(res, kind=xp)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res_xp, 2, &
            MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
    product = real(res_xp(1), kind=rp)
    norm_squared = real(res_xp(2), kind=rp)

  end subroutine device_bicgstab_product_and_norm

  !> BiCGStab part 1, \f$s = r - \alpha v\f$.
  !! @param s_d Intermediate residual, written.
  !! @param r_d Residual.
  !! @param v_d Operator action \f$A\hat{p}\f$.
  !! @param mult_d Spectral element multiplicity weights.
  !! @param alpha BiCGStab \f$\alpha\f$.
  !! @param n Number of vector entries.
  !! @return The weighted squared norm \f$s^T M s\f$.
  function device_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n) &
       result(res)
    type(c_ptr) :: s_d, r_d, v_d, mult_d
    real(kind=rp) :: alpha
    integer :: n
    real(kind=rp) :: res
    real(kind=xp) :: res_xp
    integer :: ierr

    res_xp = 0.0_xp
#if HAVE_HIP
    res_xp = hip_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n)
#elif HAVE_CUDA
    res_xp = cuda_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n)
#elif HAVE_OPENCL
    res_xp = opencl_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n, &
         glb_cmd_queue)
#elif HAVE_METAL
    res = metal_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n, &
         glb_cmd_queue)
    res_xp = real(res, kind=xp)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res_xp, 1, &
            MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
    res = real(res_xp, kind=rp)

  end function device_bicgstab_part1

  !> BiCGStab part 2, \f$x = x + \alpha \hat{p} + \omega \hat{s}\f$ and
  !! \f$r = s - \omega t\f$.
  !!
  !! The rho inner product of the next iteration is reduced together with
  !! the residual norm, since both are taken against the residual this
  !! routine writes.
  !! @param rtr Weighted squared norm \f$r^T M r\f$ of the new residual.
  !! @param rho Weighted inner product \f$f^T M r\f$ of the new residual.
  !! @param x_d Solution, updated in place.
  !! @param r_d Residual, written.
  !! @param p_hat_d Preconditioned search direction.
  !! @param s_hat_d Preconditioned intermediate residual.
  !! @param s_d Intermediate residual.
  !! @param t_d Operator action \f$A\hat{s}\f$.
  !! @param f_d Right-hand side, which is also the shadow residual.
  !! @param mult_d Spectral element multiplicity weights.
  !! @param alpha BiCGStab \f$\alpha\f$.
  !! @param omega BiCGStab \f$\omega\f$.
  !! @param n Number of vector entries.
  subroutine device_bicgstab_part2(rtr, rho, x_d, r_d, p_hat_d, s_hat_d, &
       s_d, t_d, f_d, mult_d, alpha, omega, n)
    real(kind=rp), intent(out) :: rtr
    real(kind=rp), intent(out) :: rho
    type(c_ptr) :: x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d, mult_d
    real(kind=rp) :: alpha, omega
    integer :: n
    real(kind=xp) :: res_xp(2)
    integer :: ierr
#if HAVE_METAL
    real(kind=rp) :: res(2)
#endif

    res_xp = 0.0_xp
#if HAVE_HIP
    call hip_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d, &
         mult_d, alpha, omega, res_xp, n)
#elif HAVE_CUDA
    call cuda_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d, &
         mult_d, alpha, omega, res_xp, n)
#elif HAVE_OPENCL
    call opencl_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d, &
         mult_d, alpha, omega, res_xp, n, glb_cmd_queue)
#elif HAVE_METAL
    call metal_bicgstab_part2(x_d, r_d, p_hat_d, s_hat_d, s_d, t_d, f_d, &
         mult_d, alpha, omega, res, n, glb_cmd_queue)
    res_xp = real(res, kind=xp)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res_xp, 2, &
            MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif
    rtr = real(res_xp(1), kind=rp)
    rho = real(res_xp(2), kind=rp)

  end subroutine device_bicgstab_part2

  !> Initialise a device BiCGStab solver.
  !! @param n Number of degrees of freedom.
  !! @param max_iter Maximum number of iterations.
  !! @param M Optional preconditioner. An identity preconditioner is used if
  !! absent.
  !! @param rel_tol Optional relative convergence tolerance.
  !! @param abs_tol Optional absolute convergence tolerance.
  !! @param monitor Optional switch for logging the residual at each iteration.
  subroutine bicgstab_device_init(this, n, max_iter, M, rel_tol, abs_tol, &
       monitor)
    class(bicgstab_device_t), target, intent(inout) :: this
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%p(n))
    allocate(this%p_hat(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%s_hat(n))
    allocate(this%t(n))
    allocate(this%v(n))

    call device_map(this%p, this%p_d, n)
    call device_map(this%p_hat, this%p_hat_d, n)
    call device_map(this%r, this%r_d, n)
    call device_map(this%s, this%s_d, n)
    call device_map(this%s_hat, this%s_hat_d, n)
    call device_map(this%t, this%t_d, n)
    call device_map(this%v, this%v_d, n)

    if (present(M)) then
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(monitor) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, monitor = monitor)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else if (present(monitor)) then
       call this%ksp_init(max_iter, monitor = monitor)
    else
       call this%ksp_init(max_iter)
    end if

    call device_event_create(this%gs_event, 2)

  end subroutine bicgstab_device_init

  !> Free a device BiCGStab solver.
  subroutine bicgstab_device_free(this)
    class(bicgstab_device_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%v)) then
       if (c_associated(this%v_d)) then
          call device_unmap(this%v, this%v_d)
       end if
       deallocate(this%v)
    end if

    if (allocated(this%r)) then
       if (c_associated(this%r_d)) then
          call device_unmap(this%r, this%r_d)
       end if
       deallocate(this%r)
    end if

    if (allocated(this%t)) then
       if (c_associated(this%t_d)) then
          call device_unmap(this%t, this%t_d)
       end if
       deallocate(this%t)
    end if

    if (allocated(this%p)) then
       if (c_associated(this%p_d)) then
          call device_unmap(this%p, this%p_d)
       end if
       deallocate(this%p)
    end if

    if (allocated(this%p_hat)) then
       if (c_associated(this%p_hat_d)) then
          call device_unmap(this%p_hat, this%p_hat_d)
       end if
       deallocate(this%p_hat)
    end if

    if (allocated(this%s)) then
       if (c_associated(this%s_d)) then
          call device_unmap(this%s, this%s_d)
       end if
       deallocate(this%s)
    end if

    if (allocated(this%s_hat)) then
       if (c_associated(this%s_hat_d)) then
          call device_unmap(this%s_hat, this%s_hat_d)
       end if
       deallocate(this%s_hat)
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine bicgstab_device_free

  !> Solve a linear system with the device BiCGStab method.
  !!
  !! The initial guess in `x` is discarded and the iteration starts from zero.
  !! @param Ax Linear operator.
  !! @param x Solution field.
  !! @param f Right-hand side.
  !! @param n Number of degrees of freedom.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param bc_projector Projector for dirichlet boundary nodes.
  !! @param gs_h Gather-scatter handle used to assemble the operator result.
  !! @param niter Optional maximum number of iterations, overriding the
  !! configured value.
  !! @return Convergence information for the solve.
  function bicgstab_device_solve(this, Ax, x, f, n, coef, bc_projector, &
       gs_h, niter) result(ksp_results)
    class(bicgstab_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    class(scalar_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: rnorm, rtr, norm_fac, gamma
    real(kind=rp) :: r_norm, s_norm, shadow_norm, t_norm, v_norm
    ! s^T s, f^T v, v^T v, s^T t, t^T t
    real(kind=rp) :: sts, ftv, vtv, stt, ttt
    real(kind=rp) :: beta, alpha, omega, rho_1, rho_2, rho_next
    type(c_ptr) :: f_d

    f_d = device_get_ptr(f)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(r_d => this%r_d, t_d => this%t_d, s_d => this%s_d, &
         v_d => this%v_d, p_d => this%p_d, s_hat_d => this%s_hat_d, &
         p_hat_d => this%p_hat_d, mult_d => coef%mult_d)

      call device_rzero(x%x_d, n)
      call device_copy(r_d, f_d, n)
      rtr = device_glsc3(r_d, mult_d, r_d, n)

      ! The implementation deliberately starts from x = 0, so f is both the
      ! initial residual and the fixed shadow residual used by BiCGStab.
      r_norm = bicgstab_device_sqrt(rtr, 'initial residual')
      shadow_norm = r_norm
      rnorm = r_norm * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      ! Avoid entering the recurrence if the zero initial guess already meets
      ! either stopping criterion. Apart from saving work, this prevents a
      ! small right-hand side from being mistaken for a rho breakdown.
      if (r_norm .le. 0.0_rp .or. rnorm .lt. this%abs_tol .or. &
           rnorm .lt. gamma) then
         ksp_results%converged = .true.
         return
      end if

      ! Every later rho comes out of part 2 together with the residual norm.
      ! On the first iteration r is still f, so rho is the norm just reduced.
      rho_1 = rtr

      call this%monitor_start('BiCGStab')
      do iter = 1, max_iter

         ! BiCGStab breaks down if the residual becomes orthogonal to the
         ! shadow residual. Test the inner product relative to the norms of
         ! both vectors so that uniformly scaling the system does not change
         ! the decision.
         call bicgstab_device_check_inner_product(rho_1, shadow_norm, &
              r_norm, 'rho inner product')

         if (iter .eq. 1) then
            call device_copy(p_d, r_d, n)
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            if (.not. ieee_is_finite(beta)) then
               call neko_error('BiCGStab failure: non-finite beta')
            end if
            call device_bicgstab_update_p(p_d, r_d, v_d, beta, omega, n)
         end if

         call this%M%solve(this%p_hat, this%p, n)
         call Ax%compute(this%v, this%p_hat, coef, x%msh, x%Xh)
         call gs_h%op(this%v, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call bc_projector%apply(this%v, n)

         ! The alpha denominator is another BiCG breakdown point. Reducing it
         ! together with ||v|| permits a scale-aware orthogonality check
         ! without an additional global synchronisation.
         call device_bicgstab_product_and_norm(ftv, vtv, f_d, v_d, &
              mult_d, n)
         v_norm = bicgstab_device_sqrt(vtv, 'operator result v')
         call bicgstab_device_check_inner_product(ftv, shadow_norm, &
              v_norm, 'alpha denominator')
         alpha = rho_1 / ftv
         if (.not. ieee_is_finite(alpha)) then
            call neko_error('BiCGStab failure: non-finite alpha')
         end if

         sts = device_bicgstab_part1(s_d, r_d, v_d, mult_d, alpha, n)

         s_norm = bicgstab_device_sqrt(sts, 'intermediate residual')
         rnorm = s_norm * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call device_add2s2(x%x_d, p_hat_d, alpha, n)
            call this%monitor_iter(iter, rnorm)
            exit
         end if

         call this%M%solve(this%s_hat, this%s, n)
         call Ax%compute(this%t, this%s_hat, coef, x%msh, x%Xh)
         call gs_h%op(this%t, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call bc_projector%apply(this%t, n)

         call device_bicgstab_product_and_norm(stt, ttt, s_d, t_d, &
              mult_d, n)
         t_norm = bicgstab_device_sqrt(ttt, 'operator result t')
         if (t_norm .le. 0.0_rp) then
            call neko_error('BiCGStab breakdown: zero omega denominator')
         end if
         if (.not. ieee_is_finite(stt)) then
            call neko_error('BiCGStab failure: non-finite omega numerator')
         end if
         omega = stt / ttt
         if (.not. ieee_is_finite(omega)) then
            call neko_error('BiCGStab failure: non-finite omega')
         end if

         call device_bicgstab_part2(rtr, rho_next, x%x_d, r_d, p_hat_d, &
              s_hat_d, s_d, t_d, f_d, mult_d, alpha, omega, n)

         r_norm = bicgstab_device_sqrt(rtr, 'recursive residual')
         rnorm = r_norm * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         ! A negative omega is valid. Breakdown occurs only when the local
         ! minimal-residual step stagnates, i.e. when t and s are numerically
         ! orthogonal.
         call bicgstab_device_check_inner_product(stt, t_norm, s_norm, &
              'omega numerator')
         rho_2 = rho_1
         rho_1 = rho_next

      end do
    end associate
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)

  end function bicgstab_device_solve

  !> Check an inner product for a BiCGStab breakdown.
  !!
  !! The comparison is normalised by both vector norms. Dividing by the larger
  !! norm first avoids overflow in their product and preserves scale invariance.
  !! @param inner_product Weighted inner product of the two vectors.
  !! @param norm_a Norm of the first vector.
  !! @param norm_b Norm of the second vector.
  !! @param quantity Name of the algorithmic quantity for error reporting.
  subroutine bicgstab_device_check_inner_product(inner_product, norm_a, &
       norm_b, quantity)
    real(kind=rp), intent(in) :: inner_product
    real(kind=rp), intent(in) :: norm_a
    real(kind=rp), intent(in) :: norm_b
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: large_norm, small_norm

    if (.not. ieee_is_finite(inner_product)) then
       call neko_error('BiCGStab failure: non-finite ' // trim(quantity))
    end if

    large_norm = max(norm_a, norm_b)
    small_norm = min(norm_a, norm_b)
    if (large_norm .le. 0.0_rp .or. &
         abs(inner_product) / large_norm .le. NEKO_EPS * small_norm) then
       call neko_error('BiCGStab breakdown: near-zero ' // trim(quantity))
    end if

  end subroutine bicgstab_device_check_inner_product

  !> Return the square root of a valid squared norm.
  !! @param value Weighted squared norm.
  !! @param quantity Name of the vector for error reporting.
  !! @return The non-negative square root.
  function bicgstab_device_sqrt(value, quantity) result(root)
    real(kind=rp), intent(in) :: value
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: root

    if (.not. ieee_is_finite(value) .or. value .lt. 0.0_rp) then
       call neko_error('BiCGStab failure: invalid ' // trim(quantity) // &
            ' norm')
    end if
    root = sqrt(value)

  end function bicgstab_device_sqrt

  !> Solve three independent systems with the device BiCGStab method.
  !!
  !! This routine sequentially invokes the scalar solver for each component.
  !! @param Ax Linear operator.
  !! @param x Solution field for the first component.
  !! @param y Solution field for the second component.
  !! @param z Solution field for the third component.
  !! @param fx Right-hand side for the first component.
  !! @param fy Right-hand side for the second component.
  !! @param fz Right-hand side for the third component.
  !! @param n Number of degrees of freedom per component.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param bc_projector Boundary conditions for the three components.
  !! @param gs_h Gather-scatter handle used to assemble operator results.
  !! @param niter Optional maximum number of iterations, overriding the
  !! configured value.
  !! @return Convergence information for each component.
  function bicgstab_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, bc_projector, gs_h, niter) result(ksp_results)
    class(bicgstab_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    class(vector_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter
    type(scalar_bc_projector_t), pointer :: bc_x, bc_y, bc_z

    call vector_bc_projector_components(bc_projector, bc_x, bc_y, bc_z)
    ksp_results(1) = this%solve(Ax, x, fx, n, coef, bc_x, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, bc_y, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, bc_z, gs_h, niter)

  end function bicgstab_device_solve_coupled

end module bicgstab_device
