!> Implements the device kernel for the `rough_log_law_t` type.
module rough_log_law_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_rough_log_law_compute(u_d, v_d, w_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, B, z0, tstep) &
          bind(c, name = 'hip_rough_log_law_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       real(c_rp) :: kappa, B, z0
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine hip_rough_log_law_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_rough_log_law_compute(u_d, v_d, w_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, B, z0, tstep) &
          bind(c, name = 'cuda_rough_log_law_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       real(c_rp) :: kappa, B, z0
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine cuda_rough_log_law_compute
  end interface
#elif HAVE_OPENCL
#endif
  public :: rough_log_law_compute_device

contains
  !> Compute the wall shear stress on device using the rough log-law model.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine rough_log_law_compute_device(u_d, v_d, w_d, &
       ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
       n_x_d, n_y_d, n_z_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
       n_nodes, lx, kappa, B, z0, tstep)
    integer, intent(in) :: n_nodes, lx, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d
    type(c_ptr), intent(in) :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, B, z0

#if HAVE_HIP
    call hip_rough_log_law_compute(u_d, v_d, w_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, z0, tstep)
#elif HAVE_CUDA
    call cuda_rough_log_law_compute(u_d, v_d, w_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, z0, tstep)
#elif HAVE_OPENCL
    call neko_error("OPENCL is not implemented for the rough log-law model")
#else
    call neko_error('No device backend configured')
#endif

  end subroutine rough_log_law_compute_device
end module rough_log_law_device
