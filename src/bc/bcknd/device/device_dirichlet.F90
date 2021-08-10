module device_dirichlet
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_dirichlet_apply_scalar(msk, x, g, m) &
          bind(c, name='hip_dirichlet_apply_scalar')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double) :: g
       integer(c_int) :: m
       type(c_ptr), value :: msk, x
     end subroutine hip_dirichlet_apply_scalar
  end interface
  
  interface
     subroutine hip_dirichlet_apply_vector(msk, x, y, z, g, m) &
          bind(c, name='hip_dirichlet_apply_vector')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double) :: g
       integer(c_int) :: m
       type(c_ptr), value :: msk, x, y, z
     end subroutine hip_dirichlet_apply_vector
  end interface
#endif
  
contains

  subroutine device_dirichlet_apply_scalar(msk, x, g, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x
    real(kind=rp), intent(in) :: g

#ifdef HAVE_HIP
    call hip_dirichlet_apply_scalar(msk, x, g, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_dirichlet_apply_scalar

  subroutine device_dirichlet_apply_vector(msk, x, y, z, g, m)
    integer, intent(in) :: m
    type(c_ptr) :: msk, x, y, z
    real(kind=rp), intent(in) :: g

#ifdef HAVE_HIP
    call hip_dirichlet_apply_vector(msk, x, y, z, g, m)
#else
    call neko_error('No device backend configured')
#endif
    
  end subroutine device_dirichlet_apply_vector
  
end module device_dirichlet
