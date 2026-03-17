!> Implements the device kernel for the `most_t` type.

module most_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr, c_char, c_null_char
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_most_compute(u_d, v_d, w_d, temp_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, z0, z0h_in, bc_type, bc_value, tstep) &
          bind(c, name = 'hip_most_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, temp_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       type(c_ptr), value :: bc_type ! pointer to first char of the string
       real(c_rp) :: kappa, z0, z0h_in, bc_value
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine hip_most_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_most_compute(u_d, v_d, w_d, temp_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, z0, z0h_in, bc_type, bc_value, tstep) &
          bind(c, name = 'cuda_most_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, temp_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       type(c_ptr) :: bc_type
       real(c_rp) :: kappa, z0, z0h_in, bc_value
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine cuda_most_compute
  end interface
#elif HAVE_OPENCL
#endif
  public :: most_compute_device

contains
  !> Compute the wall shear stress on device using the rough log-law model.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine most_compute_device(u_d, v_d, w_d, temp_d, &
       ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
       n_x_d, n_y_d, n_z_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
       n_nodes, lx, kappa, z0, z0h_in, bc_type, bc_value, tstep)
    use, intrinsic :: iso_c_binding, only : c_loc ! needed for C-string handling
    integer, intent(in) :: n_nodes, lx, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d, temp_d
    type(c_ptr), intent(in) :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, z0, z0h_in, bc_value
    character(len=*), intent(in) :: bc_type ! passed as a normal Fortran string

    ! bc_type must be C-compatible string
    character(kind=c_char, len=len_trim(bc_type)+1), target :: bc_type_c
    ! add null terminator for C/C++ compatibility
    bc_type_c = trim(bc_type) // c_null_char

#if HAVE_HIP
    call hip_most_compute(u_d, v_d, w_d,temp_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         lx, kappa, z0, z0h_in, c_loc(bc_type_c), bc_value, tstep)
#elif HAVE_CUDA
    call cuda_most_compute(u_d, v_d, w_d,temp_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         lx, kappa, z0, z0h_in, c_loc(bc_type_c), bc_value, tstep)
#elif HAVE_OPENCL
    call neko_error("OPENCL is not implemented for the MOST wall model")
#else
    call neko_error('No device backend configured')
#endif

  end subroutine most_compute_device
end module most_device
