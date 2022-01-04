module ax_helm_device
  use device_math
  use ax_product
  use num_types
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public, extends(ax_t) :: ax_helm_device_t
   contains
     procedure, nopass :: compute => ax_helm_device_compute
  end type ax_helm_device_t

#ifdef HAVE_HIP
    interface
     subroutine hip_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, dxt_d, dyt_d, dzt_d, &
          h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d, nelv, lx) &
          bind(c, name='hip_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine hip_ax_helm
  end interface
#elif HAVE_CUDA
    interface
     subroutine cuda_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, h1_d, g11_d, g22_d, g33_d, &
          g12_d, g13_d, g23_d, nelv, lx) bind(c, name='cuda_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine cuda_ax_helm
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_ax_helm(w_d, u_d, &
          dx_d, dy_d, dz_d, dxt_d, dyt_d, dzt_d, &
          h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d, nelv, lx) &
          bind(c, name='opencl_ax_helm')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d
       type(c_ptr), value :: h1_d, g11_d, g22_d, g33_d, g12_d, g13_d, g23_d
       integer(c_int) :: nel, lx
     end subroutine opencl_ax_helm
  end interface
#endif

contains

  subroutine ax_helm_device_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(c_ptr) :: u_d, w_d

    u_d = device_get_ptr(u, size(u))
    w_d = device_get_ptr(w, size(u))

#ifdef HAVE_HIP
    call hip_ax_helm(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_CUDA
    call cuda_ax_helm(w_d, u_d, &
         Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call opencl_ax_helm(w_d, u_d, Xh%dx_d, Xh%dy_d, Xh%dz_d, &
         Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%h1_d, &
         coef%G11_d, coef%G22_d, coef%G33_d, &
         coef%G12_d, coef%G13_d, coef%G23_d, &
         msh%nelv, Xh%lx)
#endif

    if (coef%ifh2) then
       call device_addcol4(w_d ,coef%h2_d, coef%B_d, u_d, coef%dof%n_dofs)
    end if
    
  end subroutine ax_helm_device_compute
  
end module ax_helm_device


