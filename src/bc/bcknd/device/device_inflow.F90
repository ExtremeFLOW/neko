module device_inflow
  use num_types
  use utils
  use, intrinsic :: iso_c_binding

#ifdef HAVE_HIP

  interface
     subroutine hip_inflow_apply_vector(msk, x, y, z, g, m) &
          bind(c, name='hip_inflow_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, g
     end subroutine hip_inflow_apply_vector
  end interface

#elif HAVE_CUDA

  interface
     subroutine cuda_inflow_apply_vector(msk, x, y, z, g, m) &
          bind(c, name='cuda_inflow_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, g
     end subroutine cuda_inflow_apply_vector
  end interface

#endif

contains

    subroutine device_inflow_apply_vector(msk, x, y, z, g, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x, y, z, g

#ifdef HAVE_HIP
    call hip_inflow_apply_vector(msk, x, y, z, g, m)
#elif HAVE_CUDA
    call cuda_inflow_apply_vector(msk, x, y, z, g, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_inflow_apply_vector
  
end module device_inflow

