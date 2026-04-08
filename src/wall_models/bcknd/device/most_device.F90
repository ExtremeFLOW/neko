!> Implements the device kernel for the `most_t` type.

module most_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_most_compute(u_d, v_d, w_d, temp_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, mu, rho, g, z0, z0h_in, bc_type_int, bc_value, tstep, &
          Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn, ti_diagn, q_diagn) &
          bind(c, name = 'hip_most_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, temp_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       real(c_rp) :: kappa, mu, rho, z0, z0h_in, bc_value
       real(c_rp) :: g(3)
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep, bc_type_int
       type(c_ptr), value :: Ri_b_diagn, L_ob_diagn
       type(c_ptr), value :: utau_diagn, magu_diagn
       type(c_ptr), value :: ti_diagn, q_diagn
     end subroutine hip_most_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_most_compute(u_d, v_d, w_d, temp_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, mu, rho, g, z0, z0h_in, bc_type_int, bc_value, tstep, &
          Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn, ti_diagn, q_diagn) &
          bind(c, name = 'cuda_most_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, temp_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d
       real(c_rp) :: kappa, mu, rho, z0, z0h_in, bc_value
       real(c_rp) :: g(3)
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep, bc_type_int
       type(c_ptr), value :: Ri_b_diagn, L_ob_diagn
       type(c_ptr), value :: utau_diagn, magu_diagn
       type(c_ptr), value :: ti_diagn, q_diagn
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
       n_nodes, lx, kappa, mu, rho, g, z0, z0h_in, bc_type, bc_value, tstep, &
       Ri_b_diagn, L_ob_diagn, utau_diagn, magu_diagn, ti_diagn, q_diagn)
    integer, intent(in) :: n_nodes, lx, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d, temp_d
    type(c_ptr), intent(in) :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, mu, rho, z0, z0h_in, bc_value
    real(kind=rp) :: g(3)
    character(len=*), intent(in) :: bc_type ! passed as a normal Fortran string
    integer :: bc_type_int
    type(c_ptr), value :: Ri_b_diagn, L_ob_diagn
    type(c_ptr), value :: utau_diagn, magu_diagn
    type(c_ptr), value :: ti_diagn, q_diagn

    ! convert bc_type to integer to avoid cross-language passing of strings
    select case (trim(adjustl(bc_type))) ! (trimmed, lowercase-consistent)
    case ("neumann")
       bc_type_int = 0
    case ("dirichlet")
       bc_type_int = 1
    case default
       call neko_error("Neumann/Dirichlet bc not specified correctly (most_device)")
    end select

#if HAVE_HIP
    call hip_most_compute(u_d, v_d, w_d,temp_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         lx, kappa, mu, rho, g, z0, z0h_in, &
         bc_type_int, bc_value, tstep, Ri_b_diagn, &
         L_ob_diagn, utau_diagn, magu_diagn, ti_diagn, q_diagn)
#elif HAVE_CUDA
    call cuda_most_compute(u_d, v_d, w_d,temp_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         lx, kappa, mu, rho, g, z0, z0h_in, &
         bc_type_int, bc_value, tstep, Ri_b_diagn, &
         L_ob_diagn, utau_diagn, magu_diagn, ti_diagn, q_diagn)
#elif HAVE_OPENCL
    call neko_error("OPENCL is not implemented for the MOST wall model")
#else
    call neko_error('No device backend configured')
#endif

  end subroutine most_compute_device
end module most_device
