module device_schwarz
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx, nelv) &
       bind(c, name='hip_schwarz_extrude')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: arr1_d, arr2_d
       integer(c_int) :: l1, l2, nx, nelv
       real(c_rp) :: f1, f2
     end subroutine hip_schwarz_extrude
     subroutine hip_schwarz_toext3d(a_d,b_d,nx, nelv) &
       bind(c, name='hip_schwarz_toext3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nc, nelv
     end subroutine hip_schwarz_toext3d
     subroutine hip_schwarz_toreg3d(b_d,a_d,nx, nelv) &
       bind(c, name='hip_schwarz_toreg3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nc, nelv
     end subroutine hip_schwarz_toreg3d
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx, nelv) &
       bind(c, name='cuda_schwarz_extrude')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: arr1_d, arr2_d
       integer(c_int) :: l1, l2, nx, nelv
       real(c_rp) :: f1, f2
     end subroutine cuda_schwarz_extrude
     subroutine cuda_schwarz_toext3d(a_d,b_d,nx, nelv) &
       bind(c, name='cuda_schwarz_toext3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine cuda_schwarz_toext3d
     subroutine cuda_schwarz_toreg3d(b_d,a_d,nx, nelv) &
       bind(c, name='cuda_schwarz_toreg3d')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: a_d, b_d 
       integer(c_int) :: nx, nelv
     end subroutine cuda_schwarz_toreg3d
  end interface

#endif
contains
  subroutine device_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,ny,nz, nelv)
    integer, intent(in) :: l1,l2,nx,ny,nz, nelv
    type(c_ptr), intent(inout) :: arr1_d,arr2_d
    real(kind=rp), intent(in) :: f1,f2
#ifdef HAVE_HIP
    call hip_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_extrude(arr1_d,l1,f1,arr2_d,l2,f2,nx,nelv)
#elif HAVE_OPENCL
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_schwarz_extrude

  subroutine device_schwarz_toext3d(a_d,b_d,nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: a_d, b_d
#ifdef HAVE_HIP
    call hip_schwarz_toext3d(a_d,b_d,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_toext3d(a_d,b_d,nx,nelv)
#elif HAVE_OPENCL
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_schwarz_toext3d

  subroutine device_schwarz_toreg3d(b_d,a_d,nx, nelv)
    integer, intent(in) :: nx, nelv
    type(c_ptr) :: a_d, b_d
#ifdef HAVE_HIP
    call hip_schwarz_toreg3d(b_d,a_d,nx,nelv)
#elif HAVE_CUDA
    call cuda_schwarz_toreg3d(b_d,a_d,nx,nelv)
#elif HAVE_OPENCL
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_schwarz_toreg3d
end module device_schwarz
