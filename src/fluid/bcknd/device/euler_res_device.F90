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
  use operators, only : div

  type, public, extends(euler_rhs_t) :: euler_res_device_t
   contains
     procedure, nopass :: compute => euler_res_device_compute
  end type euler_res_device_t

#ifdef HAVE_HIP
  interface
    subroutine euler_res_part_visc_hip(rhs_rho_field_d, Binv_d, lap_rho_d, c_avisc_low, n) &
        bind(c, name = 'euler_res_part_visc_hip')
      use, intrinsic :: iso_c_binding
      import c_rp
      implicit none
      type(c_ptr), value :: rhs_rho_field_d, Binv_d, lap_rho_d
      real(c_rp) :: c_avisc_low
      integer(c_int) :: n
    end subroutine euler_res_part_visc_hip
  end interface

  interface
    subroutine euler_res_part_mx_flux_hip(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
        bind(c, name = 'euler_res_part_mx_flux_hip')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
      integer(c_int) :: n
    end subroutine euler_res_part_mx_flux_hip
  end interface

  interface
    subroutine euler_res_part_my_flux_hip(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
        bind(c, name = 'euler_res_part_my_flux_hip')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
      integer(c_int) :: n
    end subroutine euler_res_part_my_flux_hip
  end interface

  interface
    subroutine euler_res_part_mz_flux_hip(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
        bind(c, name = 'euler_res_part_mz_flux_hip')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
      integer(c_int) :: n
    end subroutine euler_res_part_mz_flux_hip
  end interface

  interface
    subroutine euler_res_part_E_flux_hip(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E, n) &
        bind(c, name = 'euler_res_part_E_flux_hip')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
      integer(c_int) :: n
    end subroutine euler_res_part_E_flux_hip
  end interface

  interface
    subroutine euler_res_part_coef_mult_hip(rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
                                            rhs_E_d, mult_d, n) &
        bind(c, name = 'euler_res_part_coef_mult_hip')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
                            rhs_E_d, mult_d
      integer(c_int) :: n
    end subroutine euler_res_part_coef_mult_hip
  end interface
#elif HAVE_CUDA
  interface
  subroutine euler_res_part_visc_cuda(rhs_rho_field_d, Binv_d, lap_rho_d, c_avisc_low, n) &
      bind(c, name = 'euler_res_part_visc_cuda')
    use, intrinsic :: iso_c_binding
    import c_rp
    implicit none
    type(c_ptr), value :: rhs_rho_field_d, Binv_d, lap_rho_d
    real(c_rp) :: c_avisc_low
    integer(c_int) :: n
  end subroutine euler_res_part_visc_cuda
  end interface

  interface
  subroutine euler_res_part_mx_flux_cuda(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
      bind(c, name = 'euler_res_part_mx_flux_cuda')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
    integer(c_int) :: n
  end subroutine euler_res_part_mx_flux_cuda
  end interface

  interface
  subroutine euler_res_part_my_flux_cuda(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
      bind(c, name = 'euler_res_part_my_flux_cuda')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
    integer(c_int) :: n
  end subroutine euler_res_part_my_flux_cuda
  end interface

  interface
  subroutine euler_res_part_mz_flux_cuda(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n) &
      bind(c, name = 'euler_res_part_mz_flux_cuda')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
    integer(c_int) :: n
  end subroutine euler_res_part_mz_flux_cuda
  end interface

  interface
  subroutine euler_res_part_E_flux_cuda(f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E, n) &
      bind(c, name = 'euler_res_part_E_flux_cuda')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
    integer(c_int) :: n
  end subroutine euler_res_part_E_flux_cuda
  end interface

  interface
    subroutine euler_res_part_coef_mult_cuda(rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
                                            rhs_E_d, mult_d, n) &
        bind(c, name = 'euler_res_part_coef_mult_cuda')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
                            rhs_E_d, mult_d
      integer(c_int) :: n
    end subroutine euler_res_part_coef_mult_cuda
  end interface
#endif

contains
  subroutine euler_res_device_compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
                rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
                coef, gs, h, c_avisc_low)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n
    real(kind=rp) :: c_avisc_low
    type(field_t), pointer :: temp, f_x, f_y, f_z, &
                              visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: temp_indices(9)

    n = coef%dof%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1))
    call neko_scratch_registry%request_field(f_x, temp_indices(2))
    call neko_scratch_registry%request_field(f_y, temp_indices(3))
    call neko_scratch_registry%request_field(f_z, temp_indices(4))

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

    !> m = m - dt * div(rho * u * u^T + p*I)
    ! m_x
#ifdef HAVE_HIP
    call euler_res_part_mx_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call euler_res_part_mx_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_y
#ifdef HAVE_HIP
    call euler_res_part_my_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
                                    m_x%x_d, m_y%x_d, m_z%x_d, &
                                    rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call euler_res_part_my_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
                                    m_x%x_d, m_y%x_d, m_z%x_d, &
                                    rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_z
#ifdef HAVE_HIP
    call euler_res_part_mz_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
                                   m_x%x_d, m_y%x_d, m_z%x_d, &
                                   rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call euler_res_part_mz_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
                                  m_x%x_d, m_y%x_d, m_z%x_d, &
                                  rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

    !> E = E - dt * div(u * (E + p))
#ifdef HAVE_HIP
    call euler_res_part_E_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
                                 m_x%x_d, m_y%x_d, m_z%x_d, &
                                 rho_field%x_d, p%x_d, E%x_d, n)
#elif HAVE_CUDA
    call euler_res_part_E_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
                                m_x%x_d, m_y%x_d, m_z%x_d, &
                                rho_field%x_d, p%x_d, E%x_d, n)
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call gs%op(rhs_m_x, GS_OP_ADD)
    call gs%op(rhs_m_y, GS_OP_ADD)
    call gs%op(rhs_m_z, GS_OP_ADD)
    call gs%op(rhs_E, GS_OP_ADD)

    call neko_scratch_registry%request_field(visc_rho, temp_indices(5))
    call neko_scratch_registry%request_field(visc_m_x, temp_indices(6))
    call neko_scratch_registry%request_field(visc_m_y, temp_indices(7))
    call neko_scratch_registry%request_field(visc_m_z, temp_indices(8))
    call neko_scratch_registry%request_field(visc_E, temp_indices(9))

    ! Calculate artificial diffusion
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    call gs%op(visc_rho, GS_OP_ADD)
    call gs%op(visc_m_x, GS_OP_ADD)
    call gs%op(visc_m_y, GS_OP_ADD)
    call gs%op(visc_m_z, GS_OP_ADD)
    call gs%op(visc_E, GS_OP_ADD)

#ifdef HAVE_HIP
    call euler_res_part_coef_mult_hip(rhs_rho_field%x_d, rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, &
                                      rhs_E%x_d, coef%mult_d, n)
    call euler_res_part_visc_hip(rhs_rho_field%x_d, coef%Binv_d, visc_rho%x_d, c_avisc_low, n)
    call euler_res_part_visc_hip(rhs_m_x%x_d, coef%Binv_d, visc_m_x%x_d, c_avisc_low, n)
    call euler_res_part_visc_hip(rhs_m_y%x_d, coef%Binv_d, visc_m_y%x_d, c_avisc_low, n)
    call euler_res_part_visc_hip(rhs_m_z%x_d, coef%Binv_d, visc_m_z%x_d, c_avisc_low, n)
    call euler_res_part_visc_hip(rhs_E%x_d, coef%Binv_d, visc_E%x_d, c_avisc_low, n)
#elif HAVE_CUDA
    call euler_res_part_coef_mult_cuda(rhs_rho_field%x_d, rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, &
                                      rhs_E%x_d, coef%mult_d, n)
    call euler_res_part_visc_cuda(rhs_rho_field%x_d, coef%Binv_d, visc_rho%x_d, c_avisc_low, n)
    call euler_res_part_visc_cuda(rhs_m_x%x_d, coef%Binv_d, visc_m_x%x_d, c_avisc_low, n)
    call euler_res_part_visc_cuda(rhs_m_y%x_d, coef%Binv_d, visc_m_y%x_d, c_avisc_low, n)
    call euler_res_part_visc_cuda(rhs_m_z%x_d, coef%Binv_d, visc_m_z%x_d, c_avisc_low, n)
    call euler_res_part_visc_cuda(rhs_E%x_d, coef%Binv_d, visc_E%x_d, c_avisc_low, n)
#elif HAVE_OPENCL
    call neko_error("OpenCL not supported")
#endif

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine euler_res_device_compute

end module euler_res_device