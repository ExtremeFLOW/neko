module device_symmetry
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l) &
          bind(c, name='hip_symmetry_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, l
       type(c_ptr), value :: xmsk, ymsk, zmsk, x, y, z
     end subroutine hip_symmetry_apply_vector
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l) &
          bind(c, name='cuda_symmetry_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: m, n, l
       type(c_ptr), value :: xmsk, ymsk, zmsk, x, y, z
     end subroutine cuda_symmetry_apply_vector
  end interface
#endif
  
contains

  subroutine device_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
    integer, intent(in) :: m, n, l
    type(c_ptr) :: xmsk, ymsk, zmsk, x, y, z

#ifdef HAVE_HIP
    call hip_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
#elif HAVE_CUDA
    call cuda_symmetry_apply_vector(xmsk, ymsk, zmsk, x, y, z, m, n, l)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_symmetry_apply_vector
  
end module device_symmetry
