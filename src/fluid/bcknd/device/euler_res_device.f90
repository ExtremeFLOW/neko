module euler_res_device
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use num_types, only : rp, c_rp
  use scratch_registry, only: neko_scratch_registry
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int

  type, public, extends(euler_rhs_t) :: euler_res_device_t
   contains
     procedure, nopass :: compute => euler_res_device_compute
  end type euler_res_device_t

#ifdef HAVE_HIP
  interface
    subroutine euler_res_part1_hip(rhs_rho_field_d, Binv_d, lap_rho_d, c_avisc, n) &
        bind(c, name = 'euler_res_part1_hip')
      use, intrinsic :: iso_c_binding
      import c_rp
      implicit none
      type(c_ptr), value :: rhs_rho_field_d, Binv_d, lap_rho_d
      real(c_rp) :: c_avisc
      integer(c_int) :: n
    end subroutine euler_res_part1_hip
  end interface
#endif

contains
  subroutine euler_res_device_compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                c_Xh, gs_Xh)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    integer :: n
    real(kind=rp) :: h, c_avisc
    type(field_t), pointer :: temp, f_x, f_y, f_z
    integer :: temp_indices(4)

    h = 0.01_rp / 1.0_rp ! grid size / polynomial degreedm
    c_avisc = 0.5_rp*h
    n = c_Xh%dof%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1))
    call neko_scratch_registry%request_field(f_x, temp_indices(2))
    call neko_scratch_registry%request_field(f_y, temp_indices(3))
    call neko_scratch_registry%request_field(f_z, temp_indices(4))

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, c_Xh)
    ! artificial diffusion for rho
    call Ax%compute(temp%x, rho_field%x, c_Xh, p%msh, p%Xh)
    call gs_Xh%op(temp, GS_OP_ADD)

#ifdef HAVE_HIP
    call euler_res_part1_hip(rhs_rho_field%x_d, c_Xh%Binv_d, c_avisc, temp%x_d, n)
#elif HAVE_CUDA
    call neko_error("CUDA not supported")
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine euler_res_device_compute

end module euler_res_device