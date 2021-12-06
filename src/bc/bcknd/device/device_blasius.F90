module device_blasius
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_blasius_apply_vector(msk, x, y, z, bla_x, bla_y, bla_z, m) &
          bind(c, name='hip_blasius_apply_vector')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, bla_x, bla_y, bla_z
     end subroutine hip_blasius_apply_vector
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_blasius_apply_vector(msk, x, y, z, bla_x, bla_y, bla_z, m) &
          bind(c, name='cuda_blasius_apply_vector')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z, bla_x, bla_y, bla_z
     end subroutine cuda_blasius_apply_vector
  end interface
#elif HAVE_OPENCL
#endif

contains

  subroutine device_blasius_apply_vector(msk, x, y, z, bla_x, bla_y, bla_z, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x, y, z, bla_x, bla_y, bla_z

#ifdef HAVE_HIP
    call hip_blasius_apply_vector(msk, x, y, z, bla_x, bla_y, bla_z, m)
#elif HAVE_CUDA
    call cuda_blasius_apply_vector(msk, x, y, z, bla_x, bla_y, bla_z, m)
#elif HAVE_OPENCL
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_blasius_apply_vector
  
end module device_blasius
