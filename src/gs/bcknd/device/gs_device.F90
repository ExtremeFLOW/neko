! Copyright (c) 2021-2022, The Neko Authors
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
!> Generic Gather-scatter backend for accelerators
module gs_device
  use num_types
  use gs_bcknd
  use device    
  use gs_ops
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Gather-scatter backend for offloading devices
  type, public, extends(gs_bcknd_t) :: gs_device_t
     integer, allocatable :: local_blk_off(:)    !< Local block offset
     integer, allocatable :: shared_blk_off(:)   !< Shared block offset
     type(c_ptr) :: local1_gs_d = C_NULL_PTR     !< Dev. ptr local gs-ops
     type(c_ptr) :: local2_gs_d = C_NULL_PTR     !< Dev. ptr local gs-ops
     type(c_ptr) :: local3_gs_d = C_NULL_PTR     !< Dev. ptr local gs-ops
     type(c_ptr) :: local_dof_gs_d = C_NULL_PTR  !< Dev. ptr local dof to gs map.
     type(c_ptr) :: local_gs_dof_d = C_NULL_PTR  !< Dev. ptr local gs to dof map.
     type(c_ptr) :: shared1_gs_d = C_NULL_PTR    !< Dev. ptr shared gs-ops
     type(c_ptr) :: shared2_gs_d = C_NULL_PTR    !< Dev. ptr shared gs-ops
     type(c_ptr) :: shared3_gs_d = C_NULL_PTR    !< Dev. ptr shared gs-ops
     type(c_ptr) :: shared_dof_gs_d = C_NULL_PTR !< Dev. ptr shrd dof to gs map.
     type(c_ptr) :: shared_gs_dof_d = C_NULL_PTR !< Dev. ptr shrd gs to dof map.
     type(c_ptr) :: local_blk_len_d = C_NULL_PTR !< Dev. ptr local n-f blocks
     type(c_ptr) :: shared_blk_len_d = C_NULL_PTR!< Dev. ptr shared n-f blocks
     type(c_ptr) :: local_blk_off_d = C_NULL_PTR !< Dev. ptr local blk offset
     type(c_ptr) :: shared_blk_off_d = C_NULL_PTR!< Dev. ptr shared blk offset
     integer :: nlocal              
     integer :: nshared
   contains
     procedure, pass(this) :: init => gs_device_init
     procedure, pass(this) :: free => gs_device_free
     procedure, pass(this) :: gather => gs_gather_device
     procedure, pass(this) :: scatter => gs_scatter_device
     procedure, pass(this) :: gather_many => gs_gather_many_device
     procedure, pass(this) :: scatter_many => gs_scatter_many_device
  end type gs_device_t

#ifdef HAVE_HIP
  interface
     subroutine hip_gather_kernel(v, m, o, dg, u, n, gd, nb, b, bo, op) &
          bind(c, name='hip_gather_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb, o, op
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine hip_gather_kernel
  end interface

  interface
     subroutine hip_scatter_kernel(v, m, dg, u, n, gd, nb, b, bo) &
          bind(c, name='hip_scatter_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine hip_scatter_kernel
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_gather_kernel(v, m, o, dg, u, n, gd, nb, b, bo, op) &
          bind(c, name='cuda_gather_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb, o, op
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine cuda_gather_kernel
  end interface

  interface
     subroutine cuda_scatter_kernel(v, m, dg, u, n, gd, nb, b, bo) &
          bind(c, name='cuda_scatter_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine cuda_scatter_kernel
  end interface

  interface
     subroutine cuda_gather_many_kernel(v1, v2, v3, m, o, dg, &
                                          u1, u2, u3, n, gd, nb, b, bo, op) &
          bind(c, name='cuda_gather_many_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb, o, op
       type(c_ptr), value :: v1, v2, v3, u1, u2, u3, dg, gd, b, bo
     end subroutine cuda_gather_many_kernel
  end interface

  interface
     subroutine cuda_scatter_many_kernel(v1, v2, v3, m, dg, &
                                           u1, u2, u3, n, gd, nb, b, bo) &
          bind(c, name='cuda_scatter_many_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb
       type(c_ptr), value :: v1, v2, v3, u1, u2, u3, dg, gd, b, bo
     end subroutine cuda_scatter_many_kernel
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_gather_kernel(v, m, o, dg, u, n, gd, nb, b, bo, op) &
          bind(c, name='opencl_gather_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb, o, op
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine opencl_gather_kernel
  end interface

  interface
     subroutine opencl_scatter_kernel(v, m, dg, u, n, gd, nb, b, bo) &
          bind(c, name='opencl_scatter_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb
       type(c_ptr), value :: v, u, dg, gd, b, bo
     end subroutine opencl_scatter_kernel
  end interface

  interface
     subroutine opencl_gather_many_kernel(v1, v2, v3, m, o, dg, &
                                          u1, u2, u3, n, gd, nb, b, bo, op) &
          bind(c, name='opencl_gather_many_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb, o, op
       type(c_ptr), value :: v1, v2, v3, u1, u2, u3, dg, gd, b, bo
     end subroutine opencl_gather_many_kernel
  end interface

  interface
     subroutine opencl_scatter_many_kernel(v1, v2, v3, m, dg, &
                                           u1, u2, u3, n, gd, nb, b, bo) &
          bind(c, name='opencl_scatter_many_kernel')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, nb
       type(c_ptr), value :: v1, v2, v3, u1, u2, u3, dg, gd, b, bo
     end subroutine opencl_scatter_many_kernel
  end interface
#endif

contains
  
  !> Accelerator backend initialisation
  subroutine gs_device_init(this, nlocal, nshared, nlcl_blks, nshrd_blks)
    class(gs_device_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared
    integer, intent(in) :: nlcl_blks
    integer, intent(in) :: nshrd_blks

    call this%free()

    this%nlocal = nlocal
    this%nshared = nshared

    allocate(this%local_blk_off(nlcl_blks))
    allocate(this%shared_blk_off(nshrd_blks))

    this%local1_gs_d = C_NULL_PTR
    this%local2_gs_d = C_NULL_PTR
    this%local3_gs_d = C_NULL_PTR
    this%local_dof_gs_d = C_NULL_PTR
    this%local_gs_dof_d = C_NULL_PTR
    this%local_blk_len_d = C_NULL_PTR
    this%local_blk_off_d = C_NULL_PTR
    this%shared1_gs_d = C_NULL_PTR
    this%shared2_gs_d = C_NULL_PTR
    this%shared3_gs_d = C_NULL_PTR
    this%shared_dof_gs_d = C_NULL_PTR
    this%shared_gs_dof_d = C_NULL_PTR
    this%shared_blk_len_d = C_NULL_PTR
    this%shared_blk_off_d = C_NULL_PTR
      
  end subroutine gs_device_init

  !> Dummy backend deallocation
  subroutine gs_device_free(this)
    class(gs_device_t), intent(inout) :: this

    if (allocated(this%local_blk_off)) then
       deallocate(this%local_blk_off)
    end if

    if (allocated(this%shared_blk_off)) then
       deallocate(this%shared_blk_off)
    end if

    if (c_associated(this%local1_gs_d)) then
       call device_free(this%local1_gs_d)
    end if

    if (c_associated(this%local2_gs_d)) then
       call device_free(this%local2_gs_d)
    end if

    if (c_associated(this%local3_gs_d)) then
       call device_free(this%local3_gs_d)
    end if

    if (c_associated(this%shared1_gs_d)) then
       call device_free(this%shared1_gs_d)
    end if

    if (c_associated(this%shared2_gs_d)) then
       call device_free(this%shared2_gs_d)
    end if

    if (c_associated(this%shared3_gs_d)) then
       call device_free(this%shared3_gs_d)
    end if

    if (c_associated(this%local_dof_gs_d)) then
       call device_free(this%local_dof_gs_d)
    end if

    if (c_associated(this%local_gs_dof_d)) then
       call device_free(this%local_gs_dof_d)
    end if

    if (c_associated(this%local_blk_len_d)) then
       call device_free(this%local_blk_len_d)
    end if
    
    if (c_associated(this%shared_blk_len_d)) then
       call device_free(this%shared_blk_len_d)
    end if

    if (c_associated(this%local_blk_off_d)) then
       call device_free(this%local_blk_off_d)
    end if

    if (c_associated(this%shared_blk_off_d)) then
       call device_free(this%shared_blk_off_d)
    end if

    this%nlocal = 0
    this%nshared = 0
    
  end subroutine gs_device_free

  !> Gather kernel
  subroutine gs_gather_device(this, v, m, o, dg, u, n, gd, nb, b, op)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer :: op, i
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)
        
    if (this%nlocal .eq. m) then       
       associate(v_d=>this%local1_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, b_d=>this%local_blk_len_d, &
            bo=>this%local_blk_off, bo_d=>this%local_blk_off_d)

         if (.not. c_associated(v_d)) then
            call device_map(v, v_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(b_d)) then
            call device_map(b, b_d, nb)
            call device_memcpy(b, b_d, nb, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(bo_d)) then
            call device_map(bo, bo_d, nb)
            bo(1) = 0
            do  i = 2, nb
               bo(i) = bo(i - 1) + b(i - 1)
            end do
            call device_memcpy(bo, bo_d, nb, HOST_TO_DEVICE)
         end if
         
#ifdef HAVE_HIP
         call hip_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
                                nb, b_d, bo_d, op)
#elif HAVE_CUDA
         call cuda_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
              nb, b_d, bo_d, op)
#elif HAVE_OPENCL
         call opencl_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
                                   nb, b_d, bo_d, op)
#else
         call neko_error('No device backend configured')
#endif
         
       end associate
    else if (this%nshared .eq. m) then
       associate(v_d=>this%shared1_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, b_d=>this%shared_blk_len_d, &
            bo=>this%shared_blk_off, bo_d=>this%shared_blk_off_d)

         if (.not. c_associated(v_d)) then
            call device_map(v, v_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(b_d)) then
            call device_map(b, b_d, nb)
            call device_memcpy(b, b_d, nb, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(bo_d)) then
            call device_map(bo, bo_d, nb)
            bo(1) = 0
            do  i = 2, nb
               bo(i) = bo(i - 1) + b(i - 1)
            end do
            call device_memcpy(bo, bo_d, nb, HOST_TO_DEVICE)
         end if

         
#ifdef HAVE_HIP   
         call hip_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
                                nb, b_d, bo_d, op)
#elif HAVE_CUDA
         call cuda_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
              nb, b_d, bo_d, op)
#elif HAVE_OPENCL
         call opencl_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
                                   nb, b_d, bo_d, op)
#else
         call neko_error('No device backend configured')
#endif
         if (this%nshared .eq. m) then
            call device_memcpy(v, v_d, m, DEVICE_TO_HOST)
         end if
       end associate
    end if

  end subroutine gs_gather_device
 
  !> Scatter kernel
  subroutine gs_scatter_device(this, v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)

    if (this%nlocal .eq. m) then
       associate(v_d=>this%local1_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, b_d=>this%local_blk_len_d, &
            bo_d=>this%local_blk_off_d)
#ifdef HAVE_HIP
         call hip_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_CUDA
         call cuda_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_OPENCL
         call opencl_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    else if (this%nshared .eq. m) then
       associate(v_d=>this%shared1_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, b_d=>this%shared_blk_len_d, &
            bo_d=>this%shared_blk_off_d)

         call device_memcpy(v, v_d, m, HOST_TO_DEVICE)
         
#ifdef HAVE_HIP
         call hip_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_CUDA
         call cuda_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_OPENCL
         call opencl_scatter_kernel(v_d, m, dg_d, u_d, n, gd_d, nb, b_d, bo_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    end if

  end subroutine gs_scatter_device

  !> Gather kernel
  subroutine gs_gather_many_device(this, v1, v2, v3, m, o, dg, &
                                         u1, u2, u3, n, gd, nb, b, op)
    integer, intent(inout) :: m
    integer, intent(inout) :: n
    integer, intent(inout) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(inout) :: o
    integer :: op, i
    type(c_ptr) :: u1_d, u2_d, u3_d

    u1_d = device_get_ptr(u1, n)
    u2_d = device_get_ptr(u2, n)
    u3_d = device_get_ptr(u3, n)
        
    if (this%nlocal .eq. m) then       
       associate(v1_d=>this%local1_gs_d, v2_d=>this%local2_gs_d, &
            v3_d=>this%local3_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, b_d=>this%local_blk_len_d, &
            bo=>this%local_blk_off, bo_d=>this%local_blk_off_d)

         if (.not. c_associated(v1_d)) then
            call device_map(v1, v1_d, m)
            call device_map(v2, v2_d, m)
            call device_map(v3, v3_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(b_d)) then
            call device_map(b, b_d, nb)
            call device_memcpy(b, b_d, nb, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(bo_d)) then
            call device_map(bo, bo_d, nb)
            bo(1) = 0
            do  i = 2, nb
               bo(i) = bo(i - 1) + b(i - 1)
            end do
            call device_memcpy(bo, bo_d, nb, HOST_TO_DEVICE)
         end if
         
#ifdef HAVE_HIP
         call hip_gather_kernel(v1_d, m, o, dg_d, u1_d, n, gd_d, &
                                nb, b_d, bo_d, op)
         call hip_gather_kernel(v2_d, m, o, dg_d, u2_d, n, gd_d, &
                                nb, b_d, bo_d, op)
         call hip_gather_kernel(v3_d, m, o, dg_d, u3_d, n, gd_d, &
                                nb, b_d, bo_d, op)
#elif HAVE_CUDA
         call cuda_gather_many_kernel(v1_d, v2_d, v3_d, m, o, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d, op)
#elif HAVE_OPENCL
         call opencl_gather_many_kernel(v1_d, v2_d, v3_d, m, o, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d, op)
#else
         call neko_error('No device backend configured')
#endif
         
       end associate
    else if (this%nshared .eq. m) then
       associate(v1_d=>this%shared1_gs_d, v2_d=>this%shared2_gs_d, &
            v3_d=>this%shared3_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, b_d=>this%shared_blk_len_d, &
            bo=>this%shared_blk_off, bo_d=>this%shared_blk_off_d)

         if (.not. c_associated(v1_d)) then
            call device_map(v1, v1_d, m)
            call device_map(v2, v2_d, m)
            call device_map(v3, v3_d, m)
         end if

         if (.not. c_associated(dg_d)) then
            call device_map(dg, dg_d, m)
            call device_memcpy(dg, dg_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(gd_d)) then
            call device_map(gd, gd_d, m)
            call device_memcpy(gd, gd_d, m, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(b_d)) then
            call device_map(b, b_d, nb)
            call device_memcpy(b, b_d, nb, HOST_TO_DEVICE)
         end if

         if (.not. c_associated(bo_d)) then
            call device_map(bo, bo_d, nb)
            bo(1) = 0
            do  i = 2, nb
               bo(i) = bo(i - 1) + b(i - 1)
            end do
            call device_memcpy(bo, bo_d, nb, HOST_TO_DEVICE)
         end if

         
#ifdef HAVE_HIP   
         call hip_gather_kernel(v_d, m, o, dg_d, u_d, n, gd_d, &
                                nb, b_d, bo_d, op)
#elif HAVE_CUDA
         call cuda_gather_many_kernel(v1_d, v2_d, v3_d, m, o, &
              dg_d, u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d, op)
#elif HAVE_OPENCL
         call opencl_gather_many_kernel(v1_d, v2_d, v3_d, m, o, &
              dg_d, u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d, op)
#else
         call neko_error('No device backend configured')
#endif
         if (this%nshared .eq. m) then
            call device_memcpy(v1, v1_d, m, DEVICE_TO_HOST)
            call device_memcpy(v2, v2_d, m, DEVICE_TO_HOST)
            call device_memcpy(v3, v3_d, m, DEVICE_TO_HOST)
         end if
       end associate
    end if

  end subroutine gs_gather_many_device

  !> Scatter many kernel
  subroutine gs_scatter_many_device(this, v1, v2, v3, m, dg, &
                                          u1, u2, u3, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_device_t), intent(inout) :: this
    real(kind=rp), dimension(m), intent(inout) :: v1
    real(kind=rp), dimension(m), intent(inout) :: v2
    real(kind=rp), dimension(m), intent(inout) :: v3
    integer, dimension(m), intent(inout) :: dg
    real(kind=rp), dimension(n), intent(inout) :: u1
    real(kind=rp), dimension(n), intent(inout) :: u2
    real(kind=rp), dimension(n), intent(inout) :: u3
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    type(c_ptr) :: u1_d, u2_d, u3_d

    u1_d = device_get_ptr(u1, n)
    u2_d = device_get_ptr(u2, n)
    u3_d = device_get_ptr(u3, n)

    if (this%nlocal .eq. m) then
       associate(v1_d=>this%local1_gs_d, v2_d=>this%local2_gs_d, &
            v3_d=>this%local3_gs_d, dg_d=>this%local_dof_gs_d, &
            gd_d=>this%local_gs_dof_d, b_d=>this%local_blk_len_d, &
            bo_d=>this%local_blk_off_d)
#ifdef HAVE_HIP
         call hip_scatter_kernel(v1_d, m, dg_d, u1_d, n, gd_d, nb, b_d, bo_d)
         call hip_scatter_kernel(v2_d, m, dg_d, u2_d, n, gd_d, nb, b_d, bo_d)
         call hip_scatter_kernel(v3_d, m, dg_d, u3_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_CUDA
         call cuda_scatter_many_kernel(v1_d, v2_d, v3_d, m, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_OPENCL
         call opencl_scatter_many_kernel(v1_d, v2_d, v3_d, m, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    else if (this%nshared .eq. m) then
       associate(v1_d=>this%shared1_gs_d, v2_d=>this%shared2_gs_d, &
            v3_d=>this%shared3_gs_d, dg_d=>this%shared_dof_gs_d, &
            gd_d=>this%shared_gs_dof_d, b_d=>this%shared_blk_len_d, &
            bo_d=>this%shared_blk_off_d)

         call device_memcpy(v1, v1_d, m, HOST_TO_DEVICE)
         call device_memcpy(v2, v2_d, m, HOST_TO_DEVICE)
         call device_memcpy(v3, v3_d, m, HOST_TO_DEVICE)
         
#ifdef HAVE_HIP
         call hip_scatter_kernel(v1_d, m, dg_d, u1_d, n, gd_d, nb, b_d, bo_d)
         call hip_scatter_kernel(v2_d, m, dg_d, u2_d, n, gd_d, nb, b_d, bo_d)
         call hip_scatter_kernel(v3_d, m, dg_d, u3_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_CUDA
         call cuda_scatter_many_kernel(v1_d, v2_d, v3_d, m, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d)
#elif HAVE_OPENCL
         call opencl_scatter_many_kernel(v1_d, v2_d, v3_d, m, dg_d, &
              u1_d, u2_d, u3_d, n, gd_d, nb, b_d, bo_d)
#else
         call neko_error('No device backend configured')
#endif
       end associate
    end if

  end subroutine gs_scatter_many_device

end module gs_device
