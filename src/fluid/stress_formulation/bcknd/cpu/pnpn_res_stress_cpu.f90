!> Residuals in the Pn-Pn formulation (CPU version)
module pnpn_res_stress_cpu
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators, only : dudxyz, cdtp, curl, opgrad, strain_rate, rotate_cyc
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use pnpn_residual, only : pnpn_prs_res_t, pnpn_vel_res_t
  use scratch_registry, only : neko_scratch_registry
  use mesh, only : mesh_t
  use num_types, only : rp
  use space, only : space_t
  use math, only : rzero, col2, copy, invers2, cmult2
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> CPU implementation of the pressure residual for the PnPn fluid with
  !! full viscous stress formulation.
  type, public, extends(pnpn_prs_res_t) :: pnpn_prs_res_stress_cpu_t
   contains
     procedure, nopass :: compute => pnpn_prs_res_stress_cpu_compute
  end type pnpn_prs_res_stress_cpu_t

  !> CPU implementation of the velocity residual for the PnPn fluid with
  !! full viscous stress formulation.
  type, public, extends(pnpn_vel_res_t) :: pnpn_vel_res_stress_cpu_t
   contains
     procedure, nopass :: compute => pnpn_vel_res_stress_cpu_compute
  end type pnpn_vel_res_stress_cpu_t

contains

  subroutine pnpn_prs_res_stress_cpu_compute(p, p_res, u, v, w, u_e, v_e, w_e,&
       f_x, f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt,&
       mu, rho, event)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(in) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: p_res
    type(field_t), intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(in) :: bc_prs_surface
    type(facet_normal_t), intent(in) :: bc_sym_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    type(c_ptr), intent(inout) :: event
    real(kind=rp) :: dtbd
    real(kind=rp) :: w1, w2, w3
    integer :: n
    integer :: i
    ! Work arrays
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3, work1, work2
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    integer :: temp_indices(14)

    ! Work arrays
    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(work1, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(work2, temp_indices(8), .false.)

    ! Stress tensor
    call neko_scratch_registry%request_field(s11, temp_indices(9), .false.)
    call neko_scratch_registry%request_field(s22, temp_indices(10), .false.)
    call neko_scratch_registry%request_field(s33, temp_indices(11), .false.)
    call neko_scratch_registry%request_field(s12, temp_indices(12), .false.)
    call neko_scratch_registry%request_field(s13, temp_indices(13), .false.)
    call neko_scratch_registry%request_field(s23, temp_indices(14), .false.)

    n = c_Xh%dof%size()

    call invers2(c_Xh%h1, rho%x, n)
    call rzero(c_Xh%h2, n)
    c_Xh%ifh2 = .false.

    ! mu times the double curl of the velocity
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

    call col2(wa1%x, mu%x, n)
    call col2(wa2%x, mu%x, n)
    call col2(wa3%x, mu%x, n)


    ! The strain rate tensor
    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, &
         u_e%x, v_e%x, w_e%x, c_Xh)


    ! Gradient of viscosity
    call dudxyz(ta1%x, mu%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call dudxyz(ta2%x, mu%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call dudxyz(ta3%x, mu%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    ! Subtract the two terms of the viscous stress to get
    ! \nabla x \nabla u - S^T \nabla \mu, and form
    ! ta = f / rho - (wa / rho) * B from it. The sign is consistent with
    ! the fact that we subtract the term below.
    !
    ! This mirrors prs_stress_res_part1 on the device backends: the whole
    ! block is a single pass rather than a chain of cmult/vdot3/sub2 calls.
    ! wa1..wa3 are dead from here until cdtp overwrites them below, so the
    ! viscous stress is kept in registers instead of being written back.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do private(i, w1, w2, w3)
    do i = 1, n
       w1 = wa1%x(i,1,1,1) - 2.0_rp * (ta1%x(i,1,1,1) * s11%x(i,1,1,1) &
            + ta2%x(i,1,1,1) * s12%x(i,1,1,1) &
            + ta3%x(i,1,1,1) * s13%x(i,1,1,1))
       w2 = wa2%x(i,1,1,1) - 2.0_rp * (ta1%x(i,1,1,1) * s12%x(i,1,1,1) &
            + ta2%x(i,1,1,1) * s22%x(i,1,1,1) &
            + ta3%x(i,1,1,1) * s23%x(i,1,1,1))
       w3 = wa3%x(i,1,1,1) - 2.0_rp * (ta1%x(i,1,1,1) * s13%x(i,1,1,1) &
            + ta2%x(i,1,1,1) * s23%x(i,1,1,1) &
            + ta3%x(i,1,1,1) * s33%x(i,1,1,1))

       ta1%x(i,1,1,1) = f_x%x(i,1,1,1) / rho%x(i,1,1,1) &
            - ((w1 / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1))
       ta2%x(i,1,1,1) = f_y%x(i,1,1,1) / rho%x(i,1,1,1) &
            - ((w2 / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1))
       ta3%x(i,1,1,1) = f_z%x(i,1,1,1) / rho%x(i,1,1,1) &
            - ((w3 / rho%x(i,1,1,1)) * c_Xh%B(i,1,1,1))
    end do
    !$omp end parallel do

    call rotate_cyc(ta1%x, ta2%x, ta3%x, 1, c_Xh)
    call gs_Xh%op(ta1%x, ta2%x, ta3%x, n, GS_OP_ADD)
    call rotate_cyc(ta1%x, ta2%x, ta3%x, 0, c_Xh)

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do
    do i = 1, n
       ta1%x(i,1,1,1) = ta1%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta2%x(i,1,1,1) = ta2%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta3%x(i,1,1,1) = ta3%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
    end do
    !$omp end parallel do

    ! Compute the components of the divergence of the rhs
    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    ! The laplacian of the pressure
    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

    !$omp parallel private (i)
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp do
    do i = 1, n
       p_res%x(i,1,1,1) = (-p_res%x(i,1,1,1)) &
            + wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1)
    end do
    !$omp end do

    !
    ! Surface velocity terms
    !
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp do
    do i = 1, n
       wa1%x(i,1,1,1) = 0.0_rp
       wa2%x(i,1,1,1) = 0.0_rp
       wa3%x(i,1,1,1) = 0.0_rp
    end do
    !$omp end do
    !$omp end parallel

    call bc_sym_surface%apply_surfvec(wa1%x, wa2%x, wa3%x, ta1%x, ta2%x, ta3%x,&
         n)

    dtbd = bd / dt
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do
    do i = 1, n
       ta1%x(i,1,1,1) = 0.0_rp
       ta2%x(i,1,1,1) = 0.0_rp
       ta3%x(i,1,1,1) = 0.0_rp
    end do
    !$omp end parallel do

    call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do
    do i = 1, n
       p_res%x(i,1,1,1) = p_res%x(i,1,1,1) &
            - (dtbd * (ta1%x(i,1,1,1) + ta2%x(i,1,1,1) + ta3%x(i,1,1,1)))&
            - (wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1))
    end do
    !$omp end parallel do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine pnpn_prs_res_stress_cpu_compute

  subroutine pnpn_vel_res_stress_cpu_compute(Ax, u, v, w, u_res, v_res, w_res, &
       p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer :: temp_indices(3)
    type(field_t), pointer :: ta1, ta2, ta3
    integer, intent(in) :: n
    integer :: i

    call copy(c_Xh%h1, mu%x, n)
    call cmult2(c_Xh%h2, rho%x, bd / dt, n)
    c_Xh%ifh2 = .true.

    ! Viscous stresses
    call Ax%compute_vector(u_res%x, v_res%x, w_res%x, u%x, v%x, w%x, c_Xh,&
         msh, Xh)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)

    ! Pressure gradient
    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    ! Sum all the terms
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do
    do i = 1, n
       u_res%x(i,1,1,1) = (-u_res%x(i,1,1,1)) - ta1%x(i,1,1,1) + f_x%x(i,1,1,1)
       v_res%x(i,1,1,1) = (-v_res%x(i,1,1,1)) - ta2%x(i,1,1,1) + f_y%x(i,1,1,1)
       w_res%x(i,1,1,1) = (-w_res%x(i,1,1,1)) - ta3%x(i,1,1,1) + f_z%x(i,1,1,1)
    end do
    !$omp end parallel do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine pnpn_vel_res_stress_cpu_compute

end module pnpn_res_stress_cpu
