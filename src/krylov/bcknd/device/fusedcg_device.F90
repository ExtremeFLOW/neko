! Copyright (c) 2021-2024, The Neko Authors
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
!> Defines a fused Conjugate Gradient method for accelerators
module fusedcg_device
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp, c_rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, copy, abscmp
  use device_math, only : device_rzero, device_copy, device_glsc3
  use device
  use comm
  implicit none
  private

  integer, parameter :: DEVICE_FUSEDCG_P_SPACE = 10

  !> Fused preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: fusedcg_device_t
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: p(:,:)
     real(kind=rp), allocatable :: alpha(:)
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: z_d = C_NULL_PTR
     type(c_ptr) :: alpha_d = C_NULL_PTR
     type(c_ptr) :: p_d_d = C_NULL_PTR
     type(c_ptr), allocatable :: p_d(:)
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => fusedcg_device_init
     procedure, pass(this) :: free => fusedcg_device_free
     procedure, pass(this) :: solve => fusedcg_device_solve
     procedure, pass(this) :: solve_coupled => fusedcg_device_solve_coupled
  end type fusedcg_device_t

#ifdef HAVE_CUDA
  interface
     subroutine cuda_fusedcg_update_p(p_d, z_d, po_d, beta, n) &
          bind(c, name='cuda_fusedcg_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, z_d, po_d
       real(c_rp) :: beta
       integer(c_int) :: n
     end subroutine cuda_fusedcg_update_p
  end interface

  interface
     subroutine cuda_fusedcg_update_x(x_d, p_d, alpha, p_cur, n) &
          bind(c, name='cuda_fusedcg_update_x')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x_d, p_d, alpha
       integer(c_int) :: p_cur, n
     end subroutine cuda_fusedcg_update_x
  end interface

  interface
     real(c_rp) function cuda_fusedcg_part2(a_d, b_d, c_d, alpha_d, alpha, &
          p_cur, n) bind(c, name='cuda_fusedcg_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, alpha_d
       real(c_rp) :: alpha
       integer(c_int) :: n, p_cur
     end function cuda_fusedcg_part2
  end interface
#elif HAVE_HIP
  interface
     subroutine hip_fusedcg_update_p(p_d, z_d, po_d, beta, n) &
          bind(c, name='hip_fusedcg_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, z_d, po_d
       real(c_rp) :: beta
       integer(c_int) :: n
     end subroutine hip_fusedcg_update_p
  end interface

  interface
     subroutine hip_fusedcg_update_x(x_d, p_d, alpha, p_cur, n) &
          bind(c, name='hip_fusedcg_update_x')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x_d, p_d, alpha
       integer(c_int) :: p_cur, n
     end subroutine hip_fusedcg_update_x
  end interface

  interface
     real(c_rp) function hip_fusedcg_part2(a_d, b_d, c_d, alpha_d, alpha, &
          p_cur, n) bind(c, name='hip_fusedcg_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d, alpha_d
       real(c_rp) :: alpha
       integer(c_int) :: n, p_cur
     end function hip_fusedcg_part2
  end interface
#endif

contains

  subroutine device_fusedcg_update_p(p_d, z_d, po_d, beta, n)
    type(c_ptr), value :: p_d, z_d, po_d
    real(c_rp) :: beta
    integer(c_int) :: n
#ifdef HAVE_HIP
    call hip_fusedcg_update_p(p_d, z_d, po_d, beta, n)
#elif HAVE_CUDA
    call cuda_fusedcg_update_p(p_d, z_d, po_d, beta, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_fusedcg_update_p

  subroutine device_fusedcg_update_x(x_d, p_d, alpha, p_cur, n)
    type(c_ptr), value :: x_d, p_d, alpha
    integer(c_int) :: p_cur, n
#ifdef HAVE_HIP
    call hip_fusedcg_update_x(x_d, p_d, alpha, p_cur, n)
#elif HAVE_CUDA
    call cuda_fusedcg_update_x(x_d, p_d, alpha, p_cur, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_fusedcg_update_x

  function device_fusedcg_part2(a_d, b_d, c_d, alpha_d, alpha, &
                                p_cur, n) result(res)
    type(c_ptr), value :: a_d, b_d, c_d, alpha_d
    real(c_rp) :: alpha
    integer :: n, p_cur
    real(kind=rp) :: res
    integer :: ierr
#ifdef HAVE_HIP
    res = hip_fusedcg_part2(a_d, b_d, c_d, alpha_d, alpha, p_cur, n)
#elif HAVE_CUDA
    res = cuda_fusedcg_part2(a_d, b_d, c_d, alpha_d, alpha, p_cur, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif

  end function device_fusedcg_part2

  !> Initialise a fused PCG solver
  subroutine fusedcg_device_init(this, n, max_iter, M, rel_tol, abs_tol, &
                                 monitor)
    class(fusedcg_device_t), target, intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    logical, optional, intent(in) :: monitor
    type(c_ptr) :: ptr
    integer(c_size_t) :: p_size
    integer :: i

    call this%free()

    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%z(n))
    allocate(this%p(n, DEVICE_FUSEDCG_P_SPACE))
    allocate(this%p_d(DEVICE_FUSEDCG_P_SPACE))
    allocate(this%alpha(DEVICE_FUSEDCG_P_SPACE))

    if (present(M)) then
       this%M => M
    end if

    call device_map(this%w, this%w_d, n)
    call device_map(this%r, this%r_d, n)
    call device_map(this%z, this%z_d, n)
    call device_map(this%alpha, this%alpha_d, DEVICE_FUSEDCG_P_SPACE)
    do i = 1, DEVICE_FUSEDCG_P_SPACE+1
       this%p_d(i) = C_NULL_PTR
       call device_map(this%p(:,i), this%p_d(i), n)
    end do

    p_size = c_sizeof(C_NULL_PTR) * (DEVICE_FUSEDCG_P_SPACE)
    call device_alloc(this%p_d_d, p_size)
    ptr = c_loc(this%p_d)
    call device_memcpy(ptr, this%p_d_d, p_size, &
                       HOST_TO_DEVICE, sync=.false.)
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

  end subroutine fusedcg_device_init

  !> Deallocate a pipelined PCG solver
  subroutine fusedcg_device_free(this)
    class(fusedcg_device_t), intent(inout) :: this
    integer :: i

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if


    if (allocated(this%z)) then
       deallocate(this%z)
    end if


    if (allocated(this%alpha)) then
       deallocate(this%alpha)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if

    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if

    if (c_associated(this%r_d)) then
       call device_free(this%r_d)
    end if

    if (c_associated(this%z_d)) then
       call device_free(this%z_d)
    end if

    if (c_associated(this%alpha_d)) then
       call device_free(this%alpha_d)
    end if

    if (allocated(this%p_d)) then
       do i = 1, DEVICE_FUSEDCG_P_SPACE
          if (c_associated(this%p_d(i))) then
             call device_free(this%p_d(i))
          end if
       end do
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine fusedcg_device_free

  !> Pipelined PCG solve
  function fusedcg_device_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(fusedcg_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, ierr, i, p_cur, p_prev
    real(kind=rp) :: rnorm, rtr, norm_fac,  rtz1, rtz2
    real(kind=rp) :: pap, beta
    type(c_ptr) :: f_d
    f_d = device_get_ptr(f)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(w => this%w, r => this%r, p => this%p, z => this%z, &
         alpha => this%alpha, alpha_d => this%alpha_d, &
         w_d => this%w_d, r_d => this%r_d, z_d => this%z_d, &
         p_d => this%p_d, p_d_d => this%p_d_d)

      rtz1 = 1.0_rp
      p_prev = DEVICE_FUSEDCG_P_SPACE
      p_cur = 1
      call device_rzero(x%x_d, n)
      call device_rzero(p_d(1), n)
      call device_copy(r_d, f_d, n)

      rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(abscmp(rnorm, 0.0_rp)) return

      call this%monitor_start('FusedCG')
      do iter = 1, max_iter
         call this%M%solve(z, r, n)
         rtz2 = rtz1
         rtz1 = device_glsc3(r_d, coef%mult_d, z_d, n)
         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         call device_fusedcg_update_p(p_d(p_cur), z_d, p_d(p_prev), beta, n)

         call Ax%compute(w, p(1, p_cur), coef, x%msh, x%Xh)
         call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call bc_list_apply(blst, w, n)

         pap = device_glsc3(w_d, coef%mult_d, this%p_d(p_cur), n)

         alpha(p_cur) = rtz1 / pap
         rtr = device_fusedcg_part2(r_d, coef%mult_d, w_d, &
                                    alpha_d, alpha(p_cur), p_cur, n)
         rnorm = sqrt(rtr)*norm_fac
         call this%monitor_iter(iter, rnorm)
         if ((p_cur .eq. DEVICE_FUSEDCG_P_SPACE) .or. &
              (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
            call device_fusedcg_update_x(x%x_d, p_d_d, alpha_d, p_cur, n)
            p_prev = p_cur
            p_cur = 1
            if (rnorm .lt. this%abs_tol) exit
         else
            p_prev = p_cur
            p_cur = p_cur + 1
         end if
      end do
      call this%monitor_stop()
      ksp_results%res_final = rnorm
      ksp_results%iter = iter

    end associate

  end function fusedcg_device_solve

  !> Pipelined PCG solve coupled solve
  function fusedcg_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(fusedcg_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: fx
    real(kind=rp), dimension(n), intent(inout) :: fy
    real(kind=rp), dimension(n), intent(inout) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter

    ksp_results(1) =  this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) =  this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) =  this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function fusedcg_device_solve_coupled
  
end module fusedcg_device


