!> Device dispatch for `caisagaut_model_ii_t`.
module caisagaut_model_ii_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, &
          ind_s_d, ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, &
          h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s) &
          bind(c, name = 'hip_caisagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes, lx
     end subroutine hip_caisagaut_model_ii_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, &
          ind_s_d, ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, &
          h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s) &
          bind(c, name = 'cuda_caisagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes, lx
     end subroutine cuda_caisagaut_model_ii_compute
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, &
          ind_s_d, ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, &
          h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s) &
          bind(c, name = 'opencl_caisagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes, lx
     end subroutine opencl_caisagaut_model_ii_compute
  end interface
#endif
  public :: caisagaut_model_ii_compute_device

contains
  !> Evaluate the device wall-model kernel for Model-II.
  !! @param u_d The sampled x-velocity field on the device.
  !! @param v_d The sampled y-velocity field on the device.
  !! @param w_d The sampled z-velocity field on the device.
  !! @param ind_r_d The r-index array for sampled GLL points.
  !! @param ind_s_d The s-index array for sampled GLL points.
  !! @param ind_t_d The t-index array for sampled GLL points.
  !! @param ind_e_d The element-index array for sampled GLL points.
  !! @param n_x_d The x-component of the wall normals.
  !! @param n_y_d The y-component of the wall normals.
  !! @param n_z_d The z-component of the wall normals.
  !! @param nu_d The sampled kinematic viscosity at wall points.
  !! @param rho_w_d The sampled density at wall points.
  !! @param h_d The wall-model sampling distances.
  !! @param tau_x_d The x-component of the wall shear stress.
  !! @param tau_y_d The y-component of the wall shear stress.
  !! @param tau_z_d The z-component of the wall shear stress.
  !! @param n_nodes The number of wall points.
  !! @param lx The one-dimensional polynomial order.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param p The blending exponent.
  !! @param s The blending scale.
  subroutine caisagaut_model_ii_compute_device(u_d, v_d, w_d, ind_r_d, &
       ind_s_d, ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
       tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s)
    integer, intent(in) :: n_nodes, lx
    type(c_ptr), intent(in) :: u_d, v_d, w_d, rho_w_d
    type(c_ptr), intent(in) :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d, nu_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, B, p, s

#if HAVE_HIP
    call hip_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, ind_s_d, &
         ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s)
#elif HAVE_CUDA
    call cuda_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, ind_s_d, &
         ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s)
#elif HAVE_OPENCL
    call opencl_caisagaut_model_ii_compute(u_d, v_d, w_d, ind_r_d, ind_s_d, &
         ind_t_d, ind_e_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, p, s)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine caisagaut_model_ii_compute_device
end module caisagaut_model_ii_device
