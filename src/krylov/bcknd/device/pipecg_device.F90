! Copyright (c) 2021-2023, The Neko Authors
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
module pipecg_device
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp, c_rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, copy
  use device_math, only : device_rzero, device_copy, device_glsc3
  use device
  use comm
  implicit none
  private

  integer, parameter :: DEVICE_PIPECG_P_SPACE = 10

  !> Pipelined preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: pipecg_device_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: q(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: u(:,:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: mi(:)
     real(kind=rp), allocatable :: ni(:)
     real(kind=rp), allocatable :: alpha(:)
     real(kind=rp), allocatable :: beta(:)
     type(c_ptr) :: p_d = C_NULL_PTR
     type(c_ptr) :: q_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: s_d = C_NULL_PTR
     type(c_ptr) :: u_d_d = C_NULL_PTR
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: z_d = C_NULL_PTR
     type(c_ptr) :: mi_d = C_NULL_PTR
     type(c_ptr) :: ni_d = C_NULL_PTR
     type(c_ptr) :: alpha_d = C_NULL_PTR
     type(c_ptr) :: beta_d = C_NULL_PTR
     type(c_ptr), allocatable :: u_d(:)
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => pipecg_device_init
     procedure, pass(this) :: free => pipecg_device_free
     procedure, pass(this) :: solve => pipecg_device_solve
  end type pipecg_device_t
  
#ifdef HAVE_CUDA
  interface
     subroutine cuda_pipecg_vecops(p_d, q_d, r_d, s_d, u_d1, u_d2, &
          w_d, z_d, ni_d, mi_d, alpha, beta, mult_d, reduction,n) &
          bind(c, name='cuda_pipecg_vecops')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, q_d, r_d, s_d, u_d1, u_d2
       type(c_ptr), value :: w_d, ni_d, mi_d, z_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: alpha, beta, reduction(3)
     end subroutine cuda_pipecg_vecops
  end interface
  
  interface
     subroutine cuda_cg_update_xp(x_d, p_d, u_d_d, alpha, beta, &
                                  p_cur, p_space, n) &
          bind(c, name='cuda_cg_update_xp')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x_d, p_d, u_d_d, alpha, beta
       integer(c_int) :: p_cur, n, p_space
     end subroutine cuda_cg_update_xp
  end interface
#elif HAVE_HIP
  interface
     subroutine hip_pipecg_vecops(p_d, q_d, r_d, s_d, u_d1, u_d2, &
          w_d, z_d, ni_d, mi_d, alpha, beta, mult_d, reduction,n) &
          bind(c, name='hip_pipecg_vecops')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_d, q_d, r_d, s_d, u_d1, u_d2
       type(c_ptr), value :: w_d, ni_d, mi_d, z_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: alpha, beta, reduction(3)
     end subroutine hip_pipecg_vecops
  end interface

  interface
     subroutine hip_cg_update_xp(x_d, p_d, u_d_d, alpha, beta, &
                                 p_cur, p_space, n) &
          bind(c, name='hip_cg_update_xp')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x_d, p_d, u_d_d, alpha, beta
       integer(c_int) :: p_cur, n, p_space
     end subroutine hip_cg_update_xp
  end interface
#endif
  
contains
  
  subroutine device_pipecg_vecops(p_d, q_d, r_d, s_d, u_d1, u_d2, &
       w_d, z_d, ni_d, mi_d, alpha, beta, mult_d, reduction,n)
    type(c_ptr), value :: p_d, q_d, r_d, s_d, u_d1, u_d2
    type(c_ptr), value :: w_d, ni_d, mi_d, z_d, mult_d
    integer(c_int) :: n
    real(c_rp) :: alpha, beta, reduction(3)
#ifdef HAVE_HIP
    call hip_pipecg_vecops(p_d, q_d, r_d,&
       s_d, u_d1, u_d2, w_d, z_d, ni_d, mi_d, alpha, beta, mult_d, reduction,n) 
#elif HAVE_CUDA
    call cuda_pipecg_vecops(p_d, q_d, r_d,&
       s_d, u_d1, u_d2, w_d, z_d, ni_d, mi_d, alpha, beta, mult_d, reduction,n) 
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_pipecg_vecops
  
  subroutine device_cg_update_xp(x_d, p_d, u_d_d, alpha, beta, p_cur, p_space, n)
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: x_d, p_d, u_d_d, alpha, beta
    integer(c_int) :: p_cur, n, p_space
#ifdef HAVE_HIP
    call hip_cg_update_xp(x_d, p_d, u_d_d, alpha, beta, p_cur, p_space, n)
#elif HAVE_CUDA
    call cuda_cg_update_xp(x_d, p_d, u_d_d, alpha, beta, p_cur, p_space, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_cg_update_xp
  
  !> Initialise a pipelined PCG solver
  subroutine pipecg_device_init(this, n, M, rel_tol, abs_tol)
    class(pipecg_device_t), target, intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    type(c_ptr) :: ptr
    integer(c_size_t) :: u_size
    integer :: i
        
    call this%free()

    allocate(this%p(n))
    allocate(this%q(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%u(n, DEVICE_PIPECG_P_SPACE+1))
    allocate(this%u_d(DEVICE_PIPECG_P_SPACE+1))
    allocate(this%w(n))
    allocate(this%z(n))
    allocate(this%mi(n))
    allocate(this%ni(n))
    allocate(this%alpha(DEVICE_PIPECG_P_SPACE))
    allocate(this%beta(DEVICE_PIPECG_P_SPACE))
    
    if (present(M)) then 
       this%M => M
    end if

    call device_map(this%p, this%p_d, n)
    call device_map(this%q, this%q_d, n)
    call device_map(this%r, this%r_d, n)
    call device_map(this%s, this%s_d, n)
    call device_map(this%w, this%w_d, n)
    call device_map(this%z, this%z_d, n)
    call device_map(this%mi, this%mi_d, n)
    call device_map(this%ni, this%ni_d, n)
    call device_map(this%alpha, this%alpha_d, DEVICE_PIPECG_P_SPACE)
    call device_map(this%beta, this%beta_d, DEVICE_PIPECG_P_SPACE)
    do i = 1, DEVICE_PIPECG_P_SPACE+1
       this%u_d(i) = C_NULL_PTR
       call device_map(this%u(:,i), this%u_d(i), n)
    end do
    !Did not work with 4 for some reason...
    u_size = 8*(DEVICE_PIPECG_P_SPACE+1)
    call device_alloc(this%u_d_d, u_size)
    ptr = c_loc(this%u_d)
    call device_memcpy(ptr,this%u_d_d, u_size, HOST_TO_DEVICE)

    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(rel_tol=rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(abs_tol=abs_tol)
    else
       call this%ksp_init()
    end if

    call device_event_create(this%gs_event, 2)
          
  end subroutine pipecg_device_init

  !> Deallocate a pipelined PCG solver
  subroutine pipecg_device_free(this)
    class(pipecg_device_t), intent(inout) :: this
    integer :: i

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
    if (allocated(this%alpha)) then
       deallocate(this%alpha)
    end if
    if (allocated(this%beta)) then
       deallocate(this%beta)
    end if
    

    if (c_associated(this%p_d)) then
       call device_free(this%p_d)
    end if
    if (c_associated(this%q_d)) then
       call device_free(this%q_d)
    end if
    if (c_associated(this%r_d)) then
       call device_free(this%r_d)
    end if
    if (c_associated(this%s_d)) then
       call device_free(this%s_d)
    end if
    if (c_associated(this%u_d_d)) then
       call device_free(this%u_d_d)
    end if
    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if
    if (c_associated(this%z_d)) then
       call device_free(this%z_d)
    end if
    if (c_associated(this%mi_d)) then
       call device_free(this%mi_d)
    end if
    if (c_associated(this%ni_d)) then
       call device_free(this%ni_d)
    end if
    if (c_associated(this%alpha_d)) then
       call device_free(this%alpha_d)
    end if
    if (c_associated(this%beta_d)) then
       call device_free(this%beta_d)
    end if
    if (allocated(this%u_d)) then
       do i = 1, DEVICE_PIPECG_P_SPACE
          if (c_associated(this%u_d(i))) then
             call device_free(this%u_d(i))
          end if
       end do
    end if

    nullify(this%M)

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine pipecg_device_free
  
  !> Pipelined PCG solve
  function pipecg_device_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(pipecg_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, ierr, p_cur, p_prev, u_prev
    real(kind=rp) :: rnorm, rtr, reduction(3), norm_fac
    real(kind=rp) :: gamma1, gamma2, delta
    real(kind=rp) :: tmp1, tmp2, tmp3
    type(MPI_Request) :: request
    type(MPI_Status) :: status
    type(c_ptr) :: f_d
    f_d = device_get_ptr(f)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(p => this%p, q => this%q, r => this%r, s => this%s, &
         u => this%u, w => this%w, z => this%z, mi => this%mi, ni => this%ni, &
         alpha => this%alpha, beta => this%beta,  &
         alpha_d => this%alpha_d, beta_d => this%beta_d, &
         p_d => this%p_d, q_d => this%q_d, r_d => this%r_d, &
         s_d => this%s_d, u_d => this%u_d, u_d_d => this%u_d_d, &
         w_d => this%w_d, z_d => this%z_d, mi_d => this%mi_d, ni_d => this%ni_d)
      
      p_prev = DEVICE_PIPECG_P_SPACE !this%p_space
      u_prev = DEVICE_PIPECG_P_SPACE + 1 !this%p_space+1
      p_cur = 1
      call device_rzero(x%x_d, n)
      call device_rzero(z_d, n)
      call device_rzero(q_d, n)
      call device_rzero(p_d, n)
      call device_rzero(s_d, n)
      call device_copy(r_d, f_d, n)
      !apply u=M^-1r
      !call device_copy(u_d(u_prev), r_d, n)
      call this%M%solve(u(1,u_prev), r, n)
      call Ax%compute(w, u(1,u_prev), coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
      call device_event_sync(this%gs_event)
      call bc_list_apply(blst, w, n)
    
      rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(rnorm .eq. 0.0_rp) return

      gamma1 = 0.0_rp
      tmp1 = 0.0_rp
      tmp2 = 0.0_rp
      tmp3 = 0.0_rp
      tmp1 = device_glsc3(r_d,coef%mult_d,u_d(u_prev),n)
      tmp2 = device_glsc3(w_d,coef%mult_d,u_d(u_prev),n)
      tmp3 = device_glsc3(r_d,coef%mult_d,r_d,n)
      reduction(1) = tmp1
      reduction(2) = tmp2
      reduction(3) = tmp3

      do iter = 1, max_iter
         call MPI_Iallreduce(MPI_IN_PLACE, reduction, 3, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, request, ierr)
         
         call this%M%solve(mi, w, n)
         call Ax%compute(ni, mi, coef, x%msh, x%Xh)
         call gs_h%op(ni, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
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
                  
         call device_pipecg_vecops(p_d, q_d, r_d,&
                                 s_d, u_d(u_prev), u_d(p_cur),&
                                 w_d, z_d, ni_d,&
                                 mi_d, alpha(p_cur), beta(p_cur),&
                                 coef%mult_d, reduction,n) 
         if (p_cur .eq. DEVICE_PIPECG_P_SPACE) then
            call device_memcpy(alpha, alpha_d, p_cur, HOST_TO_DEVICE)
            call device_memcpy(beta, beta_d, p_cur, HOST_TO_DEVICE)
            call device_cg_update_xp(x%x_d, p_d, u_d_d, alpha_d, beta_d, p_cur, &
                                     DEVICE_PIPECG_P_SPACE, n)
            p_prev = p_cur
            u_prev = DEVICE_PIPECG_P_SPACE + 1
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
         call device_memcpy(alpha, alpha_d, p_cur, HOST_TO_DEVICE)
         call device_memcpy(beta, beta_d, p_cur, HOST_TO_DEVICE)
         call device_cg_update_xp(x%x_d, p_d, u_d_d, alpha_d, beta_d, p_cur, &
                                  DEVICE_PIPECG_P_SPACE, n)
      end if

      ksp_results%res_final = rnorm
      ksp_results%iter = iter
      
    end associate
    
  end function pipecg_device_solve
   
end module pipecg_device
  

