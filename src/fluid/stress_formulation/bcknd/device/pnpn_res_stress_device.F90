!> Residuals in the Pn-Pn formulation (device version)
module pnpn_res_stress_device
  use gather_scatter, only : gs_t, GS_OP_ADD
  use utils, only : neko_error
  use operators, only : dudxyz, cdtp, curl, opgrad, strain_rate
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t
  use scratch_registry, only: neko_scratch_registry
  use mesh, only : mesh_t
  use num_types, only : rp, c_rp
  use space, only : space_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use device_mathops, only : device_opcolv
  use device_math, only : device_rzero, device_vdot3, device_cmult, &
                          device_sub2, device_col2, device_copy, device_invcol1
  implicit none
  private

  !> Device implementation of the pressure residual for the PnPn fluid with
  !! full viscous stress formulation.
  type, public, extends(pnpn_prs_res_t) :: pnpn_prs_res_stress_device_t
   contains
     procedure, nopass :: compute => pnpn_prs_res_stress_device_compute
  end type pnpn_prs_res_stress_device_t

  !> Device implementation of the velocity residual for the PnPn fluid with
  !! full viscous stress formulation.
  type, public, extends(pnpn_vel_res_t) :: pnpn_vel_res_stress_device_t
   contains
     procedure, nopass :: compute => pnpn_vel_res_stress_device_compute
  end type pnpn_vel_res_stress_device_t

#ifdef HAVE_HIP
  interface
     subroutine pnpn_prs_res_part1_hip(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, mu, rho, n) &
          bind(c, name = 'pnpn_prs_res_part1_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d
       real(c_rp) :: mu, rho
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part1_hip
  end interface

  interface
     subroutine pnpn_prs_res_part2_hip(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name = 'pnpn_prs_res_part2_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_hip
  end interface

  interface
     subroutine pnpn_prs_res_part3_hip(p_res_d, ta1_d, ta2_d, ta3_d, dtbd, n) &
          bind(c, name = 'pnpn_prs_res_part3_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part3_hip
  end interface

  interface
     subroutine pnpn_vel_res_update_hip(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name = 'pnpn_vel_res_update_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_hip
  end interface
#elif HAVE_CUDA
  interface
     subroutine pnpn_prs_stress_res_part1_cuda(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, rho_d, n) &
          bind(c, name = 'pnpn_prs_stress_res_part1_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d, rho_d
       integer(c_int) :: n
     end subroutine pnpn_prs_stress_res_part1_cuda
  end interface

  interface
     subroutine pnpn_prs_res_part2_cuda(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name = 'pnpn_prs_res_part2_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_cuda
  end interface

  interface
     subroutine pnpn_prs_stress_res_part3_cuda(p_res_d, ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, dtbd, n) &
          bind(c, name = 'pnpn_prs_stress_res_part3_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_stress_res_part3_cuda
  end interface

  interface
     subroutine pnpn_vel_res_update_cuda(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name = 'pnpn_vel_res_update_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_cuda
  end interface
#elif HAVE_OPENCL
  interface
     subroutine pnpn_prs_res_part1_opencl(ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, f_u_d, f_v_d, f_w_d, &
          B_d, h1_d, mu, rho, n) &
          bind(c, name = 'pnpn_prs_res_part1_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d, h1_d
       real(c_rp) :: mu, rho
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part1_opencl
  end interface

  interface
     subroutine pnpn_prs_res_part2_opencl(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name = 'pnpn_prs_res_part2_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part2_opencl
  end interface

  interface
     subroutine pnpn_prs_res_part3_opencl(p_res_d, ta1_d, ta2_d, ta3_d, &
          wa1_d, wa2_d, wa3_d, dtbd, n) &
          bind(c, name = 'pnpn_prs_res_part3_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: p_res_d, ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: wa1_d, wa2_d, wa3_d
       real(c_rp) :: dtbd
       integer(c_int) :: n
     end subroutine pnpn_prs_res_part3_opencl
  end interface

  interface
     subroutine pnpn_vel_res_update_opencl(u_res_d, v_res_d, w_res_d, &
          ta1_d, ta2_d, ta3_d, f_u_d, f_v_d, f_w_d, n) &
          bind(c, name = 'pnpn_vel_res_update_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine pnpn_vel_res_update_opencl
  end interface
#endif

contains

  subroutine pnpn_prs_res_stress_device_compute(p, p_res, u, v, w, u_e, v_e,&
       w_e, f_x, f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd,&
       dt, mu, rho)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: p_res
    type(field_t), intent(inout) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(inout) :: bc_prs_surface
    type(facet_normal_t), intent(inout) :: bc_sym_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=rp), intent(inout) :: bd
    real(kind=rp), intent(in) :: dt
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp) :: dtbd
    integer :: n, nelv, lxyz, gdim
    integer :: i, e
    ! Work arrays
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, work3
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    integer :: temp_indices(15)

    ! Work arrays
    call neko_scratch_registry%request_field(ta1, temp_indices(1))
    call neko_scratch_registry%request_field(ta2, temp_indices(2))
    call neko_scratch_registry%request_field(ta3, temp_indices(3))
    call neko_scratch_registry%request_field(wa1, temp_indices(4))
    call neko_scratch_registry%request_field(wa2, temp_indices(5))
    call neko_scratch_registry%request_field(wa3, temp_indices(6))
    call neko_scratch_registry%request_field(work1, temp_indices(7))
    call neko_scratch_registry%request_field(work2, temp_indices(8))
    call neko_scratch_registry%request_field(work3, temp_indices(9))

   ! Stress tensor
    call neko_scratch_registry%request_field(s11, temp_indices(10))
    call neko_scratch_registry%request_field(s22, temp_indices(11))
    call neko_scratch_registry%request_field(s33, temp_indices(12))
    call neko_scratch_registry%request_field(s12, temp_indices(13))
    call neko_scratch_registry%request_field(s13, temp_indices(14))
    call neko_scratch_registry%request_field(s23, temp_indices(15))

    n = c_Xh%dof%size()
    lxyz = c_Xh%Xh%lxyz
    nelv = c_Xh%msh%nelv
    gdim = c_Xh%msh%gdim

    call device_copy(c_Xh%h1_d, rho%x_d, n)
    call device_invcol1(c_Xh%h1_d, n)
    call device_rzero(c_Xh%h2_d, n)
    c_Xh%ifh2 = .false.

    ! mu times the double curl of the velocity
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

    call device_col2(wa1%x_d, mu%x_d, n)
    call device_col2(wa2%x_d, mu%x_d, n)
    call device_col2(wa3%x_d, mu%x_d, n)


    ! The strain rate tensor
    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, &
                     u_e, v_e, w_e, c_Xh)


    ! Gradient of viscosity * 2
    !call opgrad(ta1%x, ta2%x, ta3%x, mu%x, c_Xh)
    call dudxyz(ta1%x, mu%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call dudxyz(ta2%x, mu%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call dudxyz(ta3%x, mu%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    call device_cmult(ta1%x_d, 2.0_rp, n)
    call device_cmult(ta2%x_d, 2.0_rp, n)
    call device_cmult(ta3%x_d, 2.0_rp, n)

    ! S^T grad \mu
    call device_vdot3 (work1%x_d, ta1%x_d, ta2%x_d, ta3%x_d, &
                       s11%x_d, s12%x_d, s13%x_d, n)

    call device_vdot3 (work2%x_d, ta1%x_d, ta2%x_d, ta3%x_d, &
                       s12%x_d, s22%x_d, s23%x_d, lxyz)

    call device_vdot3 (work3%x_d, ta1%x_d, ta2%x_d, ta3%x_d, &
                       s13%x_d, s23%x_d, s33%x_d, lxyz)

    ! Subtract the two terms of the viscous stress to get
    ! \nabla x \nabla u - S^T \nabla \mu
    ! The sign is consitent with the fact that we subtract the term
    ! below.
    call device_sub2(wa1%x_d, work1%x_d, n)
    call device_sub2(wa2%x_d, work2%x_d, n)
    call device_sub2(wa3%x_d, work3%x_d, n)

#if HAVE_CUDA
    call pnpn_prs_stress_res_part1_cuda(ta1%x_d, ta2%x_d, ta3%x_d, &
         wa1%x_d, wa2%x_d, wa3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, &
         c_Xh%B_d, c_Xh%h1_d, rho%x_d, n)
#else
    call neko_error('No device backend configured')
#endif

    call gs_Xh%op(ta1, GS_OP_ADD)
    call gs_Xh%op(ta2, GS_OP_ADD)
    call gs_Xh%op(ta3, GS_OP_ADD)

    call device_opcolv(ta1%x_d, ta2%x_d, ta3%x_d, c_Xh%Binv_d, gdim, n)

    ! Compute the components of the divergence of the rhs
    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    ! The laplacian of the pressure
    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

#ifdef HAVE_HIP
    call pnpn_prs_res_part2_hip(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n)
#elif HAVE_CUDA
    call pnpn_prs_res_part2_cuda(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n)
#elif HAVE_OPENCL
    call pnpn_prs_res_part2_opencl(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n)
#endif

    !
    ! Surface velocity terms
    !
    call device_rzero(wa1%x_d, n)
    call device_rzero(wa2%x_d, n)
    call device_rzero(wa3%x_d, n)

    call bc_sym_surface%apply_surfvec_dev(wa1%x_d, wa2%x_d, wa3%x_d, &
                                          ta1%x_d , ta2%x_d, ta3%x_d)

    dtbd = bd / dt
    call device_rzerO(ta1%x_d, n)
    call device_rzerO(ta2%x_d, n)
    call device_rzerO(ta3%x_d, n)

    call bc_prs_surface%apply_surfvec_dev(ta1%x_d, ta2%x_d, ta3%x_d, &
                                          u%x_D, v%x_d, w%x_d)

#if HAVE_CUDA
    call pnpn_prs_stress_res_part3_cuda(p_res%x_d, ta1%x_d, ta2%x_d, ta3%x_d, &
                                        wa1%x_d, wa2%x_d, wa3%x_d, dtbd, n)
#else
    call neko_error('No device backend configured')
#endif

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine pnpn_prs_res_stress_device_compute

  subroutine pnpn_vel_res_stress_device_compute(Ax, u, v, w, u_res, v_res, &
       w_res, p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(inout) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: bddt
    integer :: temp_indices(3)
    type(field_t), pointer :: ta1, ta2, ta3
    integer, intent(in) :: n
    integer :: i

    call device_copy(c_Xh%h1_d, mu%x_d, n)
    call device_copy(c_Xh%h2_d, rho%x_d, n)

    bddt = bd / dt
    call device_cmult(c_Xh%h2_d, bddt, n)

    c_Xh%ifh2 = .true.

    ! Viscous stresses
    call Ax%compute_vector(u_res%x, v_res%x, w_res%x, u%x, v%x, w%x, c_Xh,&
                                msh, Xh)

    call neko_scratch_registry%request_field(ta1, temp_indices(1))
    call neko_scratch_registry%request_field(ta2, temp_indices(2))
    call neko_scratch_registry%request_field(ta3, temp_indices(3))

    ! Pressure gradient
    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    ! Sum all the terms
#ifdef HAVE_HIP
    call pnpn_vel_res_update_hip(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_CUDA
    call pnpn_vel_res_update_cuda(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_OPENCL
    call pnpn_vel_res_update_opencl(u_res%x_d, v_res%x_d, w_res%x_d, &
         ta1%x_d, ta2%x_d, ta3%x_d, f_x%x_d, f_y%x_d, f_z%x_d, n)
#endif

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine pnpn_vel_res_stress_device_compute

end module pnpn_res_stress_device
