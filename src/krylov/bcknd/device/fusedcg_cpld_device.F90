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
module fusedcg_cpld_device
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp, c_rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, copy, abscmp
  use device_math, only : device_rzero, device_copy, device_glsc3, device_glsc2
  use device
  use comm
  implicit none
  private

  integer, parameter :: DEVICE_FUSEDCG_CPLD_P_SPACE = 10

  !> Fused preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: fusedcg_cpld_device_t
     real(kind=rp), allocatable :: w1(:)
     real(kind=rp), allocatable :: w2(:)
     real(kind=rp), allocatable :: w3(:)
     real(kind=rp), allocatable :: r1(:)
     real(kind=rp), allocatable :: r2(:)
     real(kind=rp), allocatable :: r3(:)
     real(kind=rp), allocatable :: z1(:)
     real(kind=rp), allocatable :: z2(:)
     real(kind=rp), allocatable :: z3(:)
     real(kind=rp), allocatable :: tmp(:)
     real(kind=rp), allocatable :: p1(:,:)
     real(kind=rp), allocatable :: p2(:,:)
     real(kind=rp), allocatable :: p3(:,:)
     real(kind=rp), allocatable :: alpha(:)
     type(c_ptr) :: w1_d = C_NULL_PTR
     type(c_ptr) :: w2_d = C_NULL_PTR
     type(c_ptr) :: w3_d = C_NULL_PTR     
     type(c_ptr) :: r1_d = C_NULL_PTR
     type(c_ptr) :: r2_d = C_NULL_PTR
     type(c_ptr) :: r3_d = C_NULL_PTR
     type(c_ptr) :: z1_d = C_NULL_PTR
     type(c_ptr) :: z2_d = C_NULL_PTR
     type(c_ptr) :: z3_d = C_NULL_PTR
     type(c_ptr) :: alpha_d = C_NULL_PTR
     type(c_ptr) :: p1_d_d = C_NULL_PTR
     type(c_ptr) :: p2_d_d = C_NULL_PTR
     type(c_ptr) :: p3_d_d = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
     type(c_ptr), allocatable :: p1_d(:)
     type(c_ptr), allocatable :: p2_d(:)
     type(c_ptr), allocatable :: p3_d(:)
     type(c_ptr) :: gs_event1 = C_NULL_PTR
     type(c_ptr) :: gs_event2 = C_NULL_PTR
     type(c_ptr) :: gs_event3 = C_NULL_PTR
   contains
     procedure, pass(this) :: init => fusedcg_cpld_device_init
     procedure, pass(this) :: free => fusedcg_cpld_device_free
     procedure, pass(this) :: solve => fusedcg_cpld_device_solve
     procedure, pass(this) :: solve_coupled => fusedcg_cpld_device_solve_coupled
  end type fusedcg_cpld_device_t

#ifdef HAVE_CUDA
  interface
     subroutine cuda_fusedcg_cpld_part1(a1_d, a2_d, a3_d, &
          b1_d, b2_d, b3_d, tmp_d, n) bind(c, name='cuda_fusedcg_cpld_part1')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, tmp_d
       integer(c_int) :: n
     end subroutine cuda_fusedcg_cpld_part1
  end interface

  interface
     subroutine cuda_fusedcg_cpld_update_p(p1_d, p2_d, p3_d, z1_d, z2_d, z3_d, &
          po1_d, po2_d, po3_d, beta, n) bind(c, name='cuda_fusedcg_cpld_update_p')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p1_d, p2_d, p3_d, z1_d, z2_d, z3_d
       type(c_ptr), value :: po1_d, po2_d, po3_d
       real(c_rp) :: beta
       integer(c_int) :: n
     end subroutine cuda_fusedcg_cpld_update_p
  end interface

  interface
     subroutine cuda_fusedcg_cpld_update_x(x1_d, x2_d, x3_d, p1_d, p2_d, p3_d, &
          alpha, p_cur, n) bind(c, name='cuda_fusedcg_cpld_update_x')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x1_d, x2_d, x3_d, p1_d, p2_d, p3_d, alpha
       integer(c_int) :: p_cur, n
     end subroutine cuda_fusedcg_cpld_update_x
  end interface

  interface
     real(c_rp) function cuda_fusedcg_cpld_part2(a1_d, a2_d, a3_d, b_d, &
          c1_d, c2_d, c3_d, alpha_d, alpha,  p_cur, n) &
          bind(c, name='cuda_fusedcg_cpld_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a1_d, a2_d, a3_d, b_d
       type(c_ptr), value :: c1_d, c2_d, c3_d, alpha_d
       real(c_rp) :: alpha
       integer(c_int) :: n, p_cur
     end function cuda_fusedcg_cpld_part2
  end interface
#elif HAVE_HIP
#endif

contains

  subroutine device_fusedcg_cpld_part1(a1_d, a2_d, a3_d, &
                                       b1_d, b2_d, b3_d, tmp_d, n)
    type(c_ptr), value :: a1_d, a2_d, a3_d, b1_d, b2_d, b3_d
    type(c_ptr), value :: tmp_d
    integer(c_int) :: n
#ifdef HAVE_HIP
#elif HAVE_CUDA
    call cuda_fusedcg_cpld_part1(a1_d, a2_d, a3_d, b1_d, b2_d, b3_d, tmp_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_fusedcg_cpld_part1

  subroutine device_fusedcg_cpld_update_p(p1_d, p2_d, p3_d, z1_d, z2_d, z3_d, &
                                          po1_d, po2_d, po3_d, beta, n)
    type(c_ptr), value :: p1_d, p2_d, p3_d, z1_d, z2_d, z3_d
    type(c_ptr), value :: po1_d, po2_d, po3_d
    real(c_rp) :: beta
    integer(c_int) :: n
#ifdef HAVE_HIP
#elif HAVE_CUDA
    call cuda_fusedcg_cpld_update_p(p1_d, p2_d, p3_d, z1_d, z2_d, z3_d, &
                                    po1_d, po2_d, po3_d, beta, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_fusedcg_cpld_update_p

  subroutine device_fusedcg_cpld_update_x(x1_d, x2_d, x3_d, &
                                          p1_d, p2_d, p3_d, alpha, p_cur, n)
    type(c_ptr), value :: x1_d, x2_d, x3_d, p1_d, p2_d, p3_d, alpha
    integer(c_int) :: p_cur, n
#ifdef HAVE_HIP
#elif HAVE_CUDA
    call cuda_fusedcg_cpld_update_x(x1_d, x2_d, x3_d, &
                                    p1_d, p2_d, p3_d, alpha, p_cur, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_fusedcg_cpld_update_x

  function device_fusedcg_cpld_part2(a1_d, a2_d, a3_d, b_d, &
       c1_d, c2_d, c3_d, alpha_d, alpha, p_cur, n) result(res)
    type(c_ptr), value :: a1_d, a2_d, a3_d, b_d
    type(c_ptr), value :: c1_d, c2_d, c3_d, alpha_d
    real(c_rp) :: alpha
    integer :: n, p_cur
    real(kind=rp) :: res
    integer :: ierr
#ifdef HAVE_HIP
#elif HAVE_CUDA
    res = cuda_fusedcg_cpld_part2(a1_d, a2_d, a3_d, b_d, &
         c1_d, c2_d, c3_d, alpha_d, alpha, p_cur, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif

  end function device_fusedcg_cpld_part2

  !> Initialise a fused PCG solver
  subroutine fusedcg_cpld_device_init(this, n, max_iter, M, &
                                      rel_tol, abs_tol, monitor)
    class(fusedcg_cpld_device_t), target, intent(inout) :: this
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

    allocate(this%w1(n))
    allocate(this%w2(n))
    allocate(this%w3(n))
    allocate(this%r1(n))
    allocate(this%r2(n))
    allocate(this%r3(n))
    allocate(this%z1(n))
    allocate(this%z2(n))
    allocate(this%z3(n))
    allocate(this%tmp(n))
    allocate(this%p1(n, DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%p2(n, DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%p3(n, DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%p1_d(DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%p2_d(DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%p3_d(DEVICE_FUSEDCG_CPLD_P_SPACE))
    allocate(this%alpha(DEVICE_FUSEDCG_CPLD_P_SPACE))

    if (present(M)) then
       this%M => M
    end if

    call device_map(this%w1, this%w1_d, n)
    call device_map(this%w2, this%w2_d, n)
    call device_map(this%w3, this%w3_d, n)
    call device_map(this%r1, this%r1_d, n)
    call device_map(this%r2, this%r2_d, n)
    call device_map(this%r3, this%r3_d, n)
    call device_map(this%z1, this%z1_d, n)
    call device_map(this%z2, this%z2_d, n)
    call device_map(this%z3, this%z3_d, n)
    call device_map(this%tmp, this%tmp_d, n)
    call device_map(this%alpha, this%alpha_d, DEVICE_FUSEDCG_CPLD_P_SPACE)
    do i = 1, DEVICE_FUSEDCG_CPLD_P_SPACE+1
       this%p1_d(i) = C_NULL_PTR
       call device_map(this%p1(:,i), this%p1_d(i), n)
       
       this%p2_d(i) = C_NULL_PTR
       call device_map(this%p2(:,i), this%p2_d(i), n)
       
       this%p3_d(i) = C_NULL_PTR
       call device_map(this%p3(:,i), this%p3_d(i), n)
    end do

    p_size = c_sizeof(C_NULL_PTR) * (DEVICE_FUSEDCG_CPLD_P_SPACE)
    call device_alloc(this%p1_d_d, p_size)
    call device_alloc(this%p2_d_d, p_size)
    call device_alloc(this%p3_d_d, p_size)
    ptr = c_loc(this%p1_d)
    call device_memcpy(ptr, this%p1_d_d, p_size, &
                       HOST_TO_DEVICE, sync=.false.)
    ptr = c_loc(this%p2_d)
    call device_memcpy(ptr, this%p2_d_d, p_size, &
                      HOST_TO_DEVICE, sync=.false.)
    ptr = c_loc(this%p3_d)
    call device_memcpy(ptr, this%p3_d_d, p_size, &
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

    call device_event_create(this%gs_event1, 2)
    call device_event_create(this%gs_event2, 2)
    call device_event_create(this%gs_event3, 2)

  end subroutine fusedcg_cpld_device_init

  !> Deallocate a pipelined PCG solver
  subroutine fusedcg_cpld_device_free(this)
    class(fusedcg_cpld_device_t), intent(inout) :: this
    integer :: i

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

    if (allocated(this%alpha)) then
       deallocate(this%alpha)
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

    if (c_associated(this%w1_d)) then
       call device_free(this%w1_d)
    end if

    if (c_associated(this%w2_d)) then
       call device_free(this%w2_d)
    end if

    if (c_associated(this%w3_d)) then
       call device_free(this%w3_d)
    end if

    if (c_associated(this%r1_d)) then
       call device_free(this%r1_d)
    end if

    if (c_associated(this%r2_d)) then
       call device_free(this%r2_d)
    end if

    if (c_associated(this%r3_d)) then
       call device_free(this%r3_d)
    end if

    if (c_associated(this%z1_d)) then
       call device_free(this%z1_d)
    end if

    if (c_associated(this%z2_d)) then
       call device_free(this%z2_d)
    end if

    if (c_associated(this%z3_d)) then
       call device_free(this%z3_d)
    end if

    if (c_associated(this%alpha_d)) then
       call device_free(this%alpha_d)
    end if

    if (c_associated(this%tmp_d)) then
       call device_free(this%tmp_d)
    end if

    if (allocated(this%p1_d)) then
       do i = 1, DEVICE_FUSEDCG_CPLD_P_SPACE
          if (c_associated(this%p1_d(i))) then
             call device_free(this%p1_d(i))
          end if
       end do
    end if

    if (allocated(this%p2_d)) then
       do i = 1, DEVICE_FUSEDCG_CPLD_P_SPACE
          if (c_associated(this%p2_d(i))) then
             call device_free(this%p2_d(i))
          end if
       end do
    end if

    if (allocated(this%p3_d)) then
       do i = 1, DEVICE_FUSEDCG_CPLD_P_SPACE
          if (c_associated(this%p3_d(i))) then
             call device_free(this%p3_d(i))
          end if
       end do
    end if

    nullify(this%M)

    if (c_associated(this%gs_event1)) then
       call device_event_destroy(this%gs_event1)
    end if

    if (c_associated(this%gs_event2)) then
       call device_event_destroy(this%gs_event2)
    end if
    
    if (c_associated(this%gs_event3)) then
       call device_event_destroy(this%gs_event3)
    end if

  end subroutine fusedcg_cpld_device_free

  !> Pipelined PCG solve coupled solve
  function fusedcg_cpld_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(fusedcg_cpld_device_t), intent(inout) :: this
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
    integer :: iter, max_iter, ierr, i, p_cur, p_prev
    real(kind=rp) :: rnorm, rtr, norm_fac,  rtz1, rtz2
    real(kind=rp) :: pap, beta
    type(c_ptr) :: fx_d
    type(c_ptr) :: fy_d
    type(c_ptr) :: fz_d

    fx_d = device_get_ptr(fx)
    fy_d = device_get_ptr(fy)
    fz_d = device_get_ptr(fz)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(w1 => this%w1, w2 => this%w2, w3 => this%w3, r1 => this%r1, &
         r2 => this%r2, r3 => this%r3, p1 => this%p1, p2 => this%p2, &
         p3 => this%p3, z1 => this%z1, z2 => this%z2, z3 => this%z3,  &
         tmp_d => this%tmp_d, alpha => this%alpha, alpha_d => this%alpha_d, &
         w1_d => this%w1_d, w2_d => this%w2_d, w3_d => this%w3_d, &
         r1_d => this%r1_d, r2_d => this%r2_d, r3_d => this%r3_d, &
         z1_d => this%z1_d, z2_d => this%z2_d, z3_d => this%z3_d, &
         p1_d => this%p1_d, p2_d => this%p2_d, p3_d => this%p3_d, &
         p1_d_d => this%p1_d_d, p2_d_d => this%p2_d_d, p3_d_d => this%p3_d_d)

      rtz1 = 1.0_rp
      p_prev = DEVICE_FUSEDCG_CPLD_P_SPACE
      p_cur = 1
          
 
      call device_rzero(x%x_d, n)
      call device_rzero(y%x_d, n)
      call device_rzero(z%x_d, n)
      call device_rzero(p1_d(1), n)
      call device_rzero(p2_d(1), n)
      call device_rzero(p3_d(1), n)
      call device_copy(r1_d, fx_d, n)
      call device_copy(r2_d, fy_d, n)
      call device_copy(r3_d, fz_d, n)

      call device_fusedcg_cpld_part1(r1_d, r2_d, r3_d, r1_d, &
                                     r2_d, r3_d, tmp_d, n)

      rtr = device_glsc3(tmp_d, coef%mult_d, coef%binv_d, n)
      
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results(1)%iter = 0
      ksp_results(2:3)%iter = -1
      if(abscmp(rnorm, 0.0_rp)) return
      call this%monitor_start('fcpldCG')
      do iter = 1, max_iter
         call this%M%solve(z1, r1, n)
         call this%M%solve(z2, r2, n)
         call this%M%solve(z3, r3, n)
         rtz2 = rtz1

         call device_fusedcg_cpld_part1(z1_d, z2_d, z3_d, &
                                        r1_d, r2_d, r3_d, tmp_d, n)
         rtz1 = device_glsc2(tmp_d, coef%mult_d, n)

         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         
         call device_fusedcg_cpld_update_p(p1_d(p_cur), p2_d(p_cur), p3_d(p_cur), &
              z1_d, z2_d, z3_d, p1_d(p_prev), p2_d(p_prev), p3_d(p_prev), beta, n)

         call Ax%compute_vector(w1, w2, w3, &
              p1(1, p_cur), p2(1, p_cur), p3(1, p_cur), coef, x%msh, x%Xh)
         call gs_h%op(w1, n, GS_OP_ADD, this%gs_event1)
         call gs_h%op(w2, n, GS_OP_ADD, this%gs_event2)
         call gs_h%op(w3, n, GS_OP_ADD, this%gs_event3)
         call device_event_sync(this%gs_event1)
         call device_event_sync(this%gs_event2)
         call device_event_sync(this%gs_event3)
         call bc_list_apply(blstx, w1, n)
         call bc_list_apply(blsty, w2, n)
         call bc_list_apply(blstz, w3, n)

         call device_fusedcg_cpld_part1(w1_d, w2_d, w3_d,  p1_d(p_cur), &
                                        p2_d(p_cur), p3_d(p_cur), tmp_d, n)

         pap = device_glsc2(tmp_d, coef%mult_d, n)
                  
         alpha(p_cur) = rtz1 / pap
         rtr = device_fusedcg_cpld_part2(r1_d, r2_d, r3_d, coef%mult_d, &
              w1_d, w2_d, w3_d, alpha_d, alpha(p_cur), p_cur, n)
         rnorm = sqrt(rtr)*norm_fac
         call this%monitor_iter(iter, rnorm)
         if ((p_cur .eq. DEVICE_FUSEDCG_CPLD_P_SPACE) .or. &
              (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
            call device_fusedcg_cpld_update_x(x%x_d, y%x_d, z%x_d, &
                 p1_d_d, p2_d_d, p3_d_d, alpha_d, p_cur, n)
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

  end function fusedcg_cpld_device_solve_coupled

  !> Pipelined PCG solve
  function fusedcg_cpld_device_solve(this, Ax, x, f, n, coef, blst, &
       gs_h, niter)  result(ksp_results)
    class(fusedcg_cpld_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter

    ! Throw and error
    call neko_error('Only defined for coupled solves')

    ksp_results%res_final = 0.0
    ksp_results%iter = 0
    
  end function fusedcg_cpld_device_solve
  
end module fusedcg_cpld_device


