module device_wall
  use num_types
  use utils
  use, intrinsic :: iso_c_binding

#ifdef HAVE_HIP
    interface
     subroutine hip_no_slip_wall_apply_scalar(msk, x, m) &
          bind(c, name='hip_no_slip_wall_apply_scalar')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x
     end subroutine hip_no_slip_wall_apply_scalar
  end interface
  
  interface
     subroutine hip_no_slip_wall_apply_vector(msk, x, y, z, m) &
          bind(c, name='hip_no_slip_wall_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z
     end subroutine hip_no_slip_wall_apply_vector
  end interface
#elif HAVE_CUDA
    interface
     subroutine cuda_no_slip_wall_apply_scalar(msk, x, m) &
          bind(c, name='cuda_no_slip_wall_apply_scalar')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x
     end subroutine cuda_no_slip_wall_apply_scalar
  end interface
  
  interface
     subroutine cuda_no_slip_wall_apply_vector(msk, x, y, z, m) &
          bind(c, name='cuda_no_slip_wall_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z
     end subroutine cuda_no_slip_wall_apply_vector
  end interface
#endif
  
contains

  subroutine device_no_slip_wall_apply_scalar(msk, x, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x

#ifdef HAVE_HIP
    call hip_no_slip_wall_apply_scalar(msk, x, m)
#elif HAVE_CUDA
    call cuda_no_slip_wall_apply_scalar(msk, x, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_no_slip_wall_apply_scalar

  subroutine device_no_slip_wall_apply_vector(msk, x, y, z, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x, y, z

#ifdef HAVE_HIP
    call hip_no_slip_wall_apply_vector(msk, x, y, z, m)
#elif HAVE_CUDA
    call cuda_no_slip_wall_apply_vector(msk, x, y, z, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_no_slip_wall_apply_vector
  
end module device_wall
