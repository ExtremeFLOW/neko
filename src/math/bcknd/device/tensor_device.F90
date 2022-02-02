module tensor_device
  use num_types
  use utils
  use comm
  use, intrinsic :: iso_c_binding
  implicit none
#ifdef HAVE_HIP
  interface
     subroutine hip_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv) &
          bind(c, name='hip_tnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d
       integer(c_int) :: nu, nv, nelv
     end subroutine hip_tnsr3d
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv) &
          bind(c, name='cuda_tnsr3d')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: v_d, u_d, A_d, Bt_d, Ct_d
       integer(c_int) :: nu, nv, nelv
     end subroutine cuda_tnsr3d
  end interface
#endif
contains

  subroutine tnsr3d_device(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
    type(c_ptr) :: v_d, u_d, A_d, Bt_d, Ct_d
    integer(c_int) :: nu, nv, nelv
#ifdef HAVE_HIP
    call hip_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
#elif HAVE_CUDA
    call cuda_tnsr3d(v_d, nv, u_d, nu, A_d, Bt_d, Ct_d, nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine tnsr3d_device

end module tensor_device
