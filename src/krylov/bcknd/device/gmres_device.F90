! Copyright (c) 2022-2024, The Neko Authors
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
!> Defines various GMRES methods
module gmres_device
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp, c_rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use device_identity, only : device_ident_t
  use math, only : rone, rzero, abscmp
  use device_math, only : device_rzero, device_copy, device_glsc3, &
                          device_add2s2, device_add2s1, device_rone, &
                          device_cmult2, device_add2s2_many, device_glsc3_many,&
                          device_sub2
  use device
  use comm
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Standard preconditioned generalized minimal residual method
  type, public, extends(ksp_t) :: gmres_device_t
     integer :: m_restart
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: c(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: z(:,:)
     real(kind=rp), allocatable :: h(:,:)
     real(kind=rp), allocatable :: v(:,:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: gam(:)
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: c_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: s_d = C_NULL_PTR
     type(c_ptr) :: gam_d = C_NULL_PTR
     type(c_ptr), allocatable :: z_d(:), h_d(:), v_d(:)
     type(c_ptr) :: z_d_d = C_NULL_PTR
     type(c_ptr) :: h_d_d = C_NULL_PTR
     type(c_ptr) :: v_d_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gmres_device_init
     procedure, pass(this) :: free => gmres_device_free
     procedure, pass(this) :: solve => gmres_device_solve
     procedure, pass(this) :: solve_coupled => gmres_device_solve_coupled
  end type gmres_device_t

#ifdef HAVE_HIP
  interface
     real(c_rp) function hip_gmres_part2(w_d, v_d_d, h_d, mult_d, j, n) &
          bind(c, name = 'hip_gmres_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: h_d, w_d, v_d_d, mult_d
       integer(c_int) :: j, n
     end function hip_gmres_part2
  end interface
#elif HAVE_CUDA

  interface
     real(c_rp) function cuda_gmres_part2(w_d, v_d_d, h_d, mult_d, j, n) &
          bind(c, name = 'cuda_gmres_part2')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: h_d, w_d, v_d_d, mult_d
       integer(c_int) :: j, n
     end function cuda_gmres_part2
  end interface
#endif

contains

  function device_gmres_part2(w_d, v_d_d, h_d, mult_d, j, n) result(alpha)
    type(c_ptr), value :: h_d, w_d, v_d_d, mult_d
    integer(c_int) :: j, n
    real(c_rp) :: alpha
    integer :: ierr
#ifdef HAVE_HIP
    alpha = hip_gmres_part2(w_d, v_d_d, h_d, mult_d, j, n)
#elif HAVE_CUDA
    alpha = cuda_gmres_part2(w_d, v_d_d, h_d, mult_d, j, n)
#else
    call neko_error('No device backend configured')
#endif

#ifndef HAVE_DEVICE_MPI
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, alpha, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
#endif

  end function device_gmres_part2

  !> Initialise a standard GMRES solver
  subroutine gmres_device_init(this, n, max_iter, M, m_restart, &
       rel_tol, abs_tol)
    class(gmres_device_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(inout), target :: M
    integer, optional, intent(inout) :: m_restart
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    type(device_ident_t), target :: M_ident
    type(c_ptr) :: ptr
    integer(c_size_t) :: z_size
    integer :: i

    if (present(m_restart)) then
       this%m_restart = m_restart
    else
       this%m_restart = 30
    end if


    call this%free()

    if (present(M)) then
       this%M => M
    else
       this%M => M_ident
    end if

    allocate(this%w(n))
    allocate(this%r(n))
    call device_map(this%w, this%w_d, n)
    call device_map(this%r, this%r_d, n)

    allocate(this%c(this%m_restart))
    allocate(this%s(this%m_restart))
    allocate(this%gam(this%m_restart + 1))
    call device_map(this%c, this%c_d, this%m_restart)
    call device_map(this%s, this%s_d, this%m_restart)
    call device_map(this%gam, this%gam_d, this%m_restart+1)

    allocate(this%z(n, this%m_restart))
    allocate(this%v(n, this%m_restart))
    allocate(this%h(this%m_restart, this%m_restart))
    allocate(this%z_d(this%m_restart))
    allocate(this%v_d(this%m_restart))
    allocate(this%h_d(this%m_restart))
    do i = 1, this%m_restart
       this%z_d(i) = c_null_ptr
       call device_map(this%z(:,i), this%z_d(i), n)

       this%v_d(i) = c_null_ptr
       call device_map(this%v(:,i), this%v_d(i), n)

       this%h_d(i) = c_null_ptr
       call device_map(this%h(:,i), this%h_d(i), this%m_restart)
    end do

    z_size = c_sizeof(C_NULL_PTR) * (this%m_restart)
    call device_alloc(this%z_d_d, z_size)
    call device_alloc(this%v_d_d, z_size)
    call device_alloc(this%h_d_d, z_size)
    ptr = c_loc(this%z_d)
    call device_memcpy(ptr, this%z_d_d, z_size, &
                       HOST_TO_DEVICE, sync = .false.)
    ptr = c_loc(this%v_d)
    call device_memcpy(ptr, this%v_d_d, z_size, &
                       HOST_TO_DEVICE, sync = .false.)
    ptr = c_loc(this%h_d)
    call device_memcpy(ptr, this%h_d_d, z_size, &
                       HOST_TO_DEVICE, sync = .false.)


    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else
       call this%ksp_init(max_iter)
    end if

    call device_event_create(this%gs_event, 2)

  end subroutine gmres_device_init

  !> Deallocate a standard GMRES solver
  subroutine gmres_device_free(this)
    class(gmres_device_t), intent(inout) :: this
    integer :: i

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%c)) then
       deallocate(this%c)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    if (allocated(this%h)) then
       deallocate(this%h)
    end if

    if (allocated(this%v)) then
       deallocate(this%v)
    end if

    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    if (allocated(this%gam)) then
       deallocate(this%gam)
    end if

    if (allocated(this%v_d)) then
       do i = 1, this%m_restart
          if (c_associated(this%v_d(i))) then
             call device_free(this%v_d(i))
          end if
       end do
    end if

    if (allocated(this%z_d)) then
       do i = 1, this%m_restart
          if (c_associated(this%z_d(i))) then
             call device_free(this%z_d(i))
          end if
       end do
    end if
    if (allocated(this%h_d)) then
       do i = 1, this%m_restart
          if (c_associated(this%h_d(i))) then
             call device_free(this%h_d(i))
          end if
       end do
    end if



    if (c_associated(this%gam_d)) then
       call device_free(this%gam_d)
    end if
    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if
    if (c_associated(this%c_d)) then
       call device_free(this%c_d)
    end if
    if (c_associated(this%r_d)) then
       call device_free(this%r_d)
    end if
    if (c_associated(this%s_d)) then
       call device_free(this%s_d)
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine gmres_device_free

  !> Standard GMRES solve
  function gmres_device_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) &
       result(ksp_results)
    class(gmres_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    integer :: i, j, k
    real(kind=rp) :: rnorm, alpha, temp, lr, alpha2, norm_fac
    logical :: conv
    type(c_ptr) :: f_d

    f_d = device_get_ptr(f)

    conv = .false.
    iter = 0

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if

    associate(w => this%w, c => this%c, r => this%r, z => this%z, h => this%h, &
         v => this%v, s => this%s, gam => this%gam, v_d => this%v_d, &
         w_d => this%w_d, r_d => this%r_d, h_d => this%h_d, &
         v_d_d => this%v_d_d, x_d => x%x_d, z_d_d => this%z_d_d, &
         c_d => this%c_d)

      norm_fac = 1.0_rp / sqrt(coef%volume)
      call rzero(gam, this%m_restart + 1)
      call rone(s, this%m_restart)
      call rone(c, this%m_restart)
      call rzero(h, this%m_restart * this%m_restart)
      call device_rzero(x%x_d, n)
      call device_rzero(this%gam_d, this%m_restart + 1)
      call device_rone(this%s_d, this%m_restart)
      call device_rone(this%c_d, this%m_restart)

      call rzero(this%h, this%m_restart**2)
!       do j = 1, this%m_restart
!          call device_rzero(h_d(j), this%m_restart)
!       end do
      do while (.not. conv .and. iter .lt. max_iter)

         if (iter .eq. 0) then
            call device_copy(r_d, f_d, n)
         else
            call device_copy(r_d, f_d, n)
            call Ax%compute(w, x%x, coef, x%msh, x%Xh)
            call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
            call device_event_sync(this%gs_event)
            call bc_list_apply(blst, w, n)
            call device_sub2(r_d, w_d, n)
         end if

         gam(1) = sqrt(device_glsc3(r_d, r_d, coef%mult_d, n))
         if (iter .eq. 0) then
            ksp_results%res_start = gam(1) * norm_fac
         end if

         if (abscmp(gam(1), 0.0_rp)) return

         rnorm = 0.0_rp
         temp = 1.0_rp / gam(1)
         call device_cmult2(v_d(1), r_d, temp, n)
         do j = 1, this%m_restart
            iter = iter+1

            call this%M%solve(z(1,j), v(1,j), n)

            call Ax%compute(w, z(1,j), coef, x%msh, x%Xh)
            call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
            call device_event_sync(this%gs_event)
            call bc_list_apply(blst, w, n)

            if (NEKO_BCKND_OPENCL .eq. 1) then
               do i = 1, j
                  h(i,j) = device_glsc3(w_d, v_d(i), coef%mult_d, n)

                  call device_add2s2(w_d, v_d(i), -h(i,j), n)

                  alpha2 = device_glsc3(w_d, w_d, coef%mult_d, n)
               end do
            else
               call device_glsc3_many(h(1,j), w_d, v_d_d, coef%mult_d, j, n)

               call device_memcpy(h(:,j), h_d(j), j, &
                                   HOST_TO_DEVICE, sync = .false.)

               alpha2 = device_gmres_part2(w_d, v_d_d, h_d(j), &
                    coef%mult_d, j, n)

            end if

            alpha = sqrt(alpha2)
            do i = 1, j-1
               temp = h(i,j)
               h(i,j) = c(i)*temp + s(i) * h(i+1,j)
               h(i+1,j) = -s(i)*temp + c(i) * h(i+1,j)
            end do

            rnorm = 0.0_rp
            if (abscmp(alpha, 0.0_rp)) then
               conv = .true.
               exit
            end if

            lr = sqrt(h(j,j) * h(j,j) + alpha2)
            temp = 1.0_rp / lr
            c(j) = h(j,j) * temp
            s(j) = alpha * temp
            h(j,j) = lr
            call device_memcpy(h(:,j), h_d(j), j, &
                                HOST_TO_DEVICE, sync = .false.)
            gam(j+1) = -s(j) * gam(j)
            gam(j) = c(j) * gam(j)

            rnorm = abs(gam(j+1)) * norm_fac
            if (rnorm .lt. this%abs_tol) then
               conv = .true.
               exit
            end if

            if (iter + 1 .gt. max_iter) exit

            if (j .lt. this%m_restart) then
               temp = 1.0_rp / alpha
               call device_cmult2(v_d(j+1), w_d, temp, n)
            end if

         end do

         j = min(j, this%m_restart)
         do k = j, 1, -1
            temp = gam(k)
            do i = j, k+1, -1
               temp = temp - h(k,i) * c(i)
            end do
            c(k) = temp / h(k,k)
         end do

         if (NEKO_BCKND_OPENCL .eq. 1) then
            do i = 1, j
               call device_add2s2(x_d, this%z_d(i), c(i), n)
            end do
         else
            call device_memcpy(c, c_d, j, HOST_TO_DEVICE, sync = .false.)
            call device_add2s2_many(x_d, z_d_d, c_d, j, n)
         end if
      end do

    end associate

    ksp_results%res_final = rnorm
    ksp_results%iter = iter

  end function gmres_device_solve

  !> Standard GMRES coupled solve
  function gmres_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(gmres_device_t), intent(inout) :: this
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

    ksp_results(1) = this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function gmres_device_solve_coupled

end module gmres_device


