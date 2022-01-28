module fdm_device
  use num_types
  use utils
  use device
  use, intrinsic :: iso_c_binding
  implicit none
#ifdef HAVE_HIP
  interface
     subroutine hip_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv) &
          bind(c, name='hip_fdm_do_fast')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: e_d, r_d, s_d, d_d
       integer(c_int) :: nl, nelv
     end subroutine hip_fdm_do_fast
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv) &
          bind(c, name='cuda_fdm_do_fast')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: e_d, r_d, s_d, d_d
       integer(c_int) :: nl, nelv
     end subroutine cuda_fdm_do_fast
  end interface
#endif
contains

  subroutine fdm_do_fast_device(e, r, s, d, nl, ldim, nelv)
    integer, intent(in) :: nl, nelv, ldim
    real(kind=rp), intent(inout) :: e(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: r(nl**ldim, nelv)
    real(kind=rp), intent(inout) :: s(nl*nl,2,ldim, nelv)
    real(kind=rp), intent(inout) :: d(nl**ldim, nelv)    
    integer ::  ie, nn, i
    type(c_ptr) :: e_d, r_d, s_d, d_d

    e_d = device_get_ptr(e, nelv)
    r_d = device_get_ptr(r, nelv)
    s_d = device_get_ptr(s, nelv)
    d_d = device_get_ptr(d, nelv)
    if (ldim .ne. 3) call neko_error('fdm dim not supported')

#ifdef HAVE_HIP
    call hip_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv)
#elif HAVE_CUDA
    call cuda_fdm_do_fast(e_d, r_d, s_d, d_d, nl, nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine fdm_do_fast_device

end module fdm_device
