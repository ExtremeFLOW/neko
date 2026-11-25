! Copyright (c) 2025-2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
module compressible_res_device
  use compressible_residual, only : compressible_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use num_types, only : rp, c_rp
  use scratch_registry, only : neko_scratch_registry
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use operators, only : grad, opgrad, rotate_cyc
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  use compressible_ops_device, only : compressible_ops_device_update_uvw, &
       compressible_ops_device_update_mxyz_p_ruvw, &
       compressible_ops_device_ns_flux_prepare, &
       compressible_ops_device_ns_flux_finalize, &
       compressible_ops_device_ns_flux_temperature
  use device_math, only : device_copy, device_add2, device_add3s2
  use device, only : glb_cmd_event
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t

  type, public, extends(compressible_rhs_t) :: compressible_res_device_t
   contains
     procedure, nopass :: step => advance_primitive_variables_device
     procedure, nopass :: evaluate_rhs => evaluate_rhs_device
  end type compressible_res_device_t

  !> Whether physical Navier-Stokes fluxes are active for the current step.
  !! Set by `compressible_rhs_set_physical_flux`, see `compressible_residual`.
  logical, public :: compressible_res_device_add_physical_flux = .false.
  !> Whether physical viscous stress is active for the current step.
  !! Set by `compressible_rhs_set_physical_flux`, see `compressible_residual`.
  logical, public :: compressible_res_device_add_physical_stress = .false.
  !> Module variable to store thermodynamic parameter set by factory.
  real(kind=rp), public :: compressible_res_device_gamma = 1.4_rp

#ifdef HAVE_HIP
  interface
     subroutine compressible_res_part_visc_hip(rhs_rho_d, rhs_m_x_d, &
          rhs_m_y_d, rhs_m_z_d, rhs_E_d, visc_rho_d, visc_m_x_d, &
          visc_m_y_d, visc_m_z_d, visc_E_d, mult_d, Binv_d, h1_d, n) &
          bind(c, name = 'compressible_res_part_visc_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, visc_rho_d, visc_m_x_d, visc_m_y_d, visc_m_z_d, &
            visc_E_d, mult_d, Binv_d, h1_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_hip
  end interface

  interface
     subroutine compressible_res_inviscid_div_hip(rhs_rho, rhs_m_x, &
          rhs_m_y, rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, &
          drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
          dx, dy, dz, jacinv, nel, lx) &
          bind(c, name = 'compressible_res_inviscid_div_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, p, u, v, w
       type(c_ptr), value :: drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
            drdz, dsdz, dtdz, dx, dy, dz, jacinv
       integer(c_int) :: nel, lx
     end subroutine compressible_res_inviscid_div_hip
  end interface

  interface
     subroutine compressible_res_part_rk_sum_hip(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, b_i, n) &
          bind(c, name = 'compressible_res_part_rk_sum_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, b_i
       integer(c_int) :: n
     end subroutine compressible_res_part_rk_sum_hip
  end interface
#elif HAVE_CUDA
  interface
     subroutine compressible_res_part_visc_cuda(rhs_rho_d, rhs_m_x_d, &
          rhs_m_y_d, rhs_m_z_d, rhs_E_d, visc_rho_d, visc_m_x_d, &
          visc_m_y_d, visc_m_z_d, visc_E_d, mult_d, Binv_d, h1_d, n) &
          bind(c, name = 'compressible_res_part_visc_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, visc_rho_d, visc_m_x_d, visc_m_y_d, visc_m_z_d, &
            visc_E_d, mult_d, Binv_d, h1_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_cuda
  end interface

  interface
     subroutine compressible_res_inviscid_div_cuda(rhs_rho, rhs_m_x, &
          rhs_m_y, rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, &
          drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
          dx, dy, dz, jacinv, nel, lx) &
          bind(c, name = 'compressible_res_inviscid_div_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, p, u, v, w
       type(c_ptr), value :: drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
            drdz, dsdz, dtdz, dx, dy, dz, jacinv
       integer(c_int) :: nel, lx
     end subroutine compressible_res_inviscid_div_cuda
  end interface

  interface
     subroutine compressible_res_part_rk_sum_cuda(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, c, n) &
          bind(c, name = 'compressible_res_part_rk_sum_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, c
       integer(c_int) :: n
     end subroutine compressible_res_part_rk_sum_cuda
  end interface
#elif HAVE_OPENCL
  interface
     subroutine compressible_res_part_visc_opencl(rhs_rho_d, rhs_m_x_d, &
          rhs_m_y_d, rhs_m_z_d, rhs_E_d, visc_rho_d, visc_m_x_d, &
          visc_m_y_d, visc_m_z_d, visc_E_d, mult_d, Binv_d, h1_d, n) &
          bind(c, name = 'compressible_res_part_visc_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, visc_rho_d, visc_m_x_d, visc_m_y_d, visc_m_z_d, &
            visc_E_d, mult_d, Binv_d, h1_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_opencl
  end interface

  interface
     subroutine compressible_res_inviscid_div_opencl(rhs_rho, rhs_m_x, &
          rhs_m_y, rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, &
          drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
          dx, dy, dz, jacinv, nel, lx) &
          bind(c, name = 'compressible_res_inviscid_div_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, p, u, v, w
       type(c_ptr), value :: drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
            drdz, dsdz, dtdz, dx, dy, dz, jacinv
       integer(c_int) :: nel, lx
     end subroutine compressible_res_inviscid_div_opencl
  end interface

  interface
     subroutine compressible_res_part_rk_sum_opencl(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, c, n) &
          bind(c, name = 'compressible_res_part_rk_sum_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, c
       integer(c_int) :: n
     end subroutine compressible_res_part_rk_sum_opencl
  end interface
#elif HAVE_METAL
  interface
     subroutine compressible_res_part_visc_metal(rhs_rho_d, rhs_m_x_d, &
          rhs_m_y_d, rhs_m_z_d, rhs_E_d, visc_rho_d, visc_m_x_d, &
          visc_m_y_d, visc_m_z_d, visc_E_d, mult_d, Binv_d, h1_d, n) &
          bind(c, name = 'compressible_res_part_visc_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, visc_rho_d, visc_m_x_d, visc_m_y_d, visc_m_z_d, &
            visc_E_d, mult_d, Binv_d, h1_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_metal
  end interface

  interface
     subroutine compressible_res_inviscid_div_metal(rhs_rho, rhs_m_x, &
          rhs_m_y, rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, &
          drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
          dx, dy, dz, jacinv, nel, lx) &
          bind(c, name = 'compressible_res_inviscid_div_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, p, u, v, w
       type(c_ptr), value :: drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
            drdz, dsdz, dtdz, dx, dy, dz, jacinv
       integer(c_int) :: nel, lx
     end subroutine compressible_res_inviscid_div_metal
  end interface

  interface
     subroutine compressible_res_part_rk_sum_metal(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, c, n) &
          bind(c, name = 'compressible_res_part_rk_sum_metal')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, c
       integer(c_int) :: n
     end subroutine compressible_res_part_rk_sum_metal
  end interface
#endif

contains
  subroutine advance_primitive_variables_device(rho_field, &
       m_x, m_y, m_z, E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, artificial_visc, mu, kappa, bcs_vel, time, &
       rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bcs_vel
    type(time_state_t), intent(in) :: time
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, ntmp
    real(kind=rp) :: c
    logical :: first
    type(field_t), pointer :: k_rho_i, k_m_x_i, k_m_y_i, k_m_z_i, k_E_i, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
         temp_p, temp_u, temp_v, temp_w, temp_ruvw
    integer :: temp_indices(30)
    type(field_list_t) :: k_rho, k_m_x, k_m_y, k_m_z, k_E

    n = p%dof%size()
    s = rk_scheme%order

    ! One stage derivative per conserved variable and stage. They are
    ! fully written by evaluate_rhs before they are read, so the scratch
    ! fields are requested without clearing them (no zero-fill kernels).
    call k_rho%init(s)
    call k_m_x%init(s)
    call k_m_y%init(s)
    call k_m_z%init(s)
    call k_E%init(s)
    do i = 1, s
       call neko_scratch_registry%request_field(k_rho_i, temp_indices(5*i-4), &
            .false.)
       call neko_scratch_registry%request_field(k_m_x_i, temp_indices(5*i-3), &
            .false.)
       call neko_scratch_registry%request_field(k_m_y_i, temp_indices(5*i-2), &
            .false.)
       call neko_scratch_registry%request_field(k_m_z_i, temp_indices(5*i-1), &
            .false.)
       call neko_scratch_registry%request_field(k_E_i, temp_indices(5*i), &
            .false.)
       call k_rho%assign(i, k_rho_i)
       call k_m_x%assign(i, k_m_x_i)
       call k_m_y%assign(i, k_m_y_i)
       call k_m_z%assign(i, k_m_z_i)
       call k_E%assign(i, k_E_i)
    end do
    ntmp = 5 * s
    call neko_scratch_registry%request_field(temp_rho, temp_indices(ntmp+1), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_x, temp_indices(ntmp+2), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_y, temp_indices(ntmp+3), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_z, temp_indices(ntmp+4), &
         .false.)
    call neko_scratch_registry%request_field(temp_E, temp_indices(ntmp+5), &
         .false.)
    call neko_scratch_registry%request_field(temp_p, temp_indices(ntmp+6), &
         .false.)
    call neko_scratch_registry%request_field(temp_u, temp_indices(ntmp+7), &
         .false.)
    call neko_scratch_registry%request_field(temp_v, temp_indices(ntmp+8), &
         .false.)
    call neko_scratch_registry%request_field(temp_w, temp_indices(ntmp+9), &
         .false.)
    call neko_scratch_registry%request_field(temp_ruvw, temp_indices(ntmp+10), &
         .false.)
    ntmp = ntmp + 10

    ! Runge-Kutta stages: temp = U + dt * sum_j A(i,j) k_j. Zero entries of
    ! the Butcher tableau are skipped (exact test, the tableau holds literal
    ! zeros; for the classical RK4 only the sub-diagonal is non-zero), and
    ! the first non-zero term is fused with the copy of U into one add3s2
    ! per field, so U is read once instead of twice.
    do i = 1, s
       first = .true.
       do j = 1, i-1
          if (rk_scheme%coeffs_A(i, j) .eq. 0.0_rp) cycle
          if (first) then
             c = dt * rk_scheme%coeffs_A(i, j)
             call device_add3s2(temp_rho%x_d, rho_field%x_d, &
                  k_rho%items(j)%ptr%x_d, 1.0_rp, c, n)
             call device_add3s2(temp_m_x%x_d, m_x%x_d, &
                  k_m_x%items(j)%ptr%x_d, 1.0_rp, c, n)
             call device_add3s2(temp_m_y%x_d, m_y%x_d, &
                  k_m_y%items(j)%ptr%x_d, 1.0_rp, c, n)
             call device_add3s2(temp_m_z%x_d, m_z%x_d, &
                  k_m_z%items(j)%ptr%x_d, 1.0_rp, c, n)
             call device_add3s2(temp_E%x_d, E%x_d, &
                  k_E%items(j)%ptr%x_d, 1.0_rp, c, n)
             first = .false.
          else
#ifdef HAVE_HIP
             call compressible_res_part_rk_sum_hip(temp_rho%x_d, &
                  temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
                  k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
                  k_m_y%items(j)%ptr%x_d, k_m_z%items(j)%ptr%x_d, &
                  k_E%items(j)%ptr%x_d, dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_CUDA
             call compressible_res_part_rk_sum_cuda(temp_rho%x_d, &
                  temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
                  k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
                  k_m_y%items(j)%ptr%x_d, k_m_z%items(j)%ptr%x_d, &
                  k_E%items(j)%ptr%x_d, dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_OPENCL
             call compressible_res_part_rk_sum_opencl(temp_rho%x_d, &
                  temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
                  k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
                  k_m_y%items(j)%ptr%x_d, k_m_z%items(j)%ptr%x_d, &
                  k_E%items(j)%ptr%x_d, dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_METAL
             call compressible_res_part_rk_sum_metal(temp_rho%x_d, &
                  temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
                  k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
                  k_m_y%items(j)%ptr%x_d, k_m_z%items(j)%ptr%x_d, &
                  k_E%items(j)%ptr%x_d, dt, rk_scheme%coeffs_A(i, j), n)
#endif
          end if
       end do

       ! No stage contribution at all (always the case for the first stage):
       ! plain copy of the solution state
       if (first) then
          call device_copy(temp_rho%x_d, rho_field%x_d, n)
          call device_copy(temp_m_x%x_d, m_x%x_d, n)
          call device_copy(temp_m_y%x_d, m_y%x_d, n)
          call device_copy(temp_m_z%x_d, m_z%x_d, n)
          call device_copy(temp_E%x_d, E%x_d, n)
       end if

       ! Compute f(U) = rhs(U) with primitive variables from the RK stage state.
       call compressible_ops_device_update_uvw(temp_u%x_d, temp_v%x_d, &
            temp_w%x_d, temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, &
            temp_rho%x_d, n)
       if (compressible_res_device_add_physical_stress) then
          call bcs_vel%apply_vector(temp_u%x, temp_v%x, temp_w%x, n, time, &
               strong = .true.)
       end if
       call compressible_ops_device_update_mxyz_p_ruvw(temp_m_x%x_d, &
            temp_m_y%x_d, temp_m_z%x_d, temp_p%x_d, temp_ruvw%x_d, &
            temp_u%x_d, temp_v%x_d, temp_w%x_d, temp_E%x_d, &
            temp_rho%x_d, compressible_res_device_gamma, n)

       call evaluate_rhs_device(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            temp_p, temp_u, temp_v, temp_w, Ax, &
            Ax_stress, coef, gs, artificial_visc, mu, kappa)
    end do

    ! Update the solution
    do i = 1, s
#ifdef HAVE_HIP
       call compressible_res_part_rk_sum_hip(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, &
            k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#elif HAVE_CUDA
       call compressible_res_part_rk_sum_cuda(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, &
            k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#elif HAVE_OPENCL
       call compressible_res_part_rk_sum_opencl(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, &
            k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#elif HAVE_METAL
       call compressible_res_part_rk_sum_metal(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
           k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, &
           k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#endif
    end do

    call neko_scratch_registry%relinquish_field(temp_indices(1:ntmp))
  end subroutine advance_primitive_variables_device

  subroutine evaluate_rhs_device(rhs_rho_field, rhs_m_x, rhs_m_y, &
       rhs_m_z, rhs_E, rho_field, &
       m_x, m_y, m_z, E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, artificial_visc, mu, kappa)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, &
         rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n, i
    type(field_t), pointer :: visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: temp_indices(5)

    n = coef%dof%size()

    !> Inviscid part, rhs = div(F) for all five conserved variables in
    !! one element-local pass: the mass flux m, the momentum flux
    !! rho u u^T + p I from m and rho, and the energy flux (E + p) u from
    !! the stage velocity. See the kernel for how the 15 flux components
    !! and the geometric factors stay on chip.
#ifdef HAVE_HIP
    call compressible_res_inviscid_div_hip(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, rho_field%x_d, m_x%x_d, &
         m_y%x_d, m_z%x_d, E%x_d, p%x_d, u%x_d, v%x_d, w%x_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, coef%jacinv_d, &
         coef%msh%nelv, coef%Xh%lx)
#elif HAVE_CUDA
    call compressible_res_inviscid_div_cuda(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, rho_field%x_d, m_x%x_d, &
         m_y%x_d, m_z%x_d, E%x_d, p%x_d, u%x_d, v%x_d, w%x_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, coef%jacinv_d, &
         coef%msh%nelv, coef%Xh%lx)
#elif HAVE_OPENCL
    call compressible_res_inviscid_div_opencl(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, rho_field%x_d, m_x%x_d, &
         m_y%x_d, m_z%x_d, E%x_d, p%x_d, u%x_d, v%x_d, w%x_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, coef%jacinv_d, &
         coef%msh%nelv, coef%Xh%lx)
#elif HAVE_METAL
    call compressible_res_inviscid_div_metal(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, rho_field%x_d, m_x%x_d, &
         m_y%x_d, m_z%x_d, E%x_d, p%x_d, u%x_d, v%x_d, w%x_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, coef%jacinv_d, &
         coef%msh%nelv, coef%Xh%lx)
#endif

    call neko_scratch_registry%request_field(visc_rho, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(visc_m_x, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(visc_m_y, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(visc_m_z, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_E, temp_indices(5), .false.)

    coef%ifh2 = .false.

    ! Density, momentum and energy are all stabilized by the same
    ! artificial viscosity. The momentum components share the operator,
    ! so they use the fused vector apply (one kernel launch, geometric
    ! factors and h1 are loaded once for all three components).
    call device_copy(coef%h1_d, artificial_visc%x_d, n)
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute_vector(visc_m_x%x, visc_m_y%x, visc_m_z%x, &
         m_x%x, m_y%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    if (compressible_res_device_add_physical_flux) then
       call add_navier_stokes_flux_device(visc_m_x, visc_m_y, visc_m_z, &
            visc_E, rho_field, p, u, v, w, mu, kappa, Ax, Ax_stress, coef)
    end if

    ! Combine the inviscid divergence and the viscous Laplacian before a
    ! single gather-scatter per field, see evaluate_rhs_cpu: both weights
    ! are continuous, so gs(a) * mult + Binv * gs(b) = gs(a * mult + Binv * b).
    ! rhs = -(rhs * mult + Binv * visc), and reset h1 to 1.
#ifdef HAVE_HIP
    call compressible_res_part_visc_hip(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, visc_rho%x_d, visc_m_x%x_d, &
         visc_m_y%x_d, visc_m_z%x_d, visc_E%x_d, coef%mult_d, coef%Binv_d, &
         coef%h1_d, n)
#elif HAVE_CUDA
    call compressible_res_part_visc_cuda(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, visc_rho%x_d, visc_m_x%x_d, &
         visc_m_y%x_d, visc_m_z%x_d, visc_E%x_d, coef%mult_d, coef%Binv_d, &
         coef%h1_d, n)
#elif HAVE_OPENCL
    call compressible_res_part_visc_opencl(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, visc_rho%x_d, visc_m_x%x_d, &
         visc_m_y%x_d, visc_m_z%x_d, visc_E%x_d, coef%mult_d, coef%Binv_d, &
         coef%h1_d, n)
#elif HAVE_METAL
    call compressible_res_part_visc_metal(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, rhs_E%x_d, visc_rho%x_d, visc_m_x%x_d, &
         visc_m_y%x_d, visc_m_z%x_d, visc_E%x_d, coef%mult_d, coef%Binv_d, &
         coef%h1_d, n)
#endif

    call gs%op(rhs_rho_field, GS_OP_ADD, glb_cmd_event)
    call rotate_cyc(rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, 1, coef)
    call gs%op(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, n, GS_OP_ADD, glb_cmd_event)
    call rotate_cyc(rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD, glb_cmd_event)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine evaluate_rhs_device

  !> Add the physical Navier-Stokes flux contribution to the viscous residual.
  !! @param visc_m_x Viscous residual for x-momentum.
  !! @param visc_m_y Viscous residual for y-momentum.
  !! @param visc_m_z Viscous residual for z-momentum.
  !! @param visc_E Viscous residual for total energy.
  !! @param rho_field Density field.
  !! @param p Pressure field.
  !! @param u X-velocity field.
  !! @param v Y-velocity field.
  !! @param w Z-velocity field.
  !! @param mu Dynamic viscosity field.
  !! @param kappa Thermal conductivity field.
  !! @param Ax Matrix-vector product operator for scalar diffusion.
  !! @param Ax_stress Matrix-vector product operator for full stress.
  !! @param coef Spatial discretization coefficients.
  subroutine add_navier_stokes_flux_device(visc_m_x, visc_m_y, visc_m_z, &
       visc_E, rho_field, p, u, v, w, mu, kappa, Ax, Ax_stress, coef)
    type(field_t), intent(inout) :: visc_m_x, visc_m_y, visc_m_z, visc_E
    type(field_t), intent(in) :: rho_field
    type(field_t), intent(in) :: p, u, v, w, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(field_t), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, &
         dwdx, dwdy, dwdz, f_x, f_y, f_z, div_flux, dissipation
    integer :: temp_indices(14)
    integer :: n

    n = coef%dof%size()

    call neko_scratch_registry%request_field(dudx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dudy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dudz, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(dvdx, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(dvdy, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(dvdz, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(dwdx, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(dwdy, temp_indices(8), .false.)
    call neko_scratch_registry%request_field(dwdz, temp_indices(9), .false.)
    call neko_scratch_registry%request_field(f_x, temp_indices(10), .false.)
    call neko_scratch_registry%request_field(f_y, temp_indices(11), .false.)
    call neko_scratch_registry%request_field(f_z, temp_indices(12), .false.)
    call neko_scratch_registry%request_field(div_flux, temp_indices(13), &
         .false.)
    call neko_scratch_registry%request_field(dissipation, temp_indices(14), &
         .false.)

    call grad(dudx%x_d, dudy%x_d, dudz%x_d, u%x_d, coef)
    call grad(dvdx%x_d, dvdy%x_d, dvdz%x_d, v%x_d, coef)
    call grad(dwdx%x_d, dwdy%x_d, dwdz%x_d, w%x_d, coef)

    call compressible_ops_device_ns_flux_prepare(div_flux%x_d, &
         dissipation%x_d, coef%h1_d, dudx%x_d, dudy%x_d, dudz%x_d, &
         dvdx%x_d, dvdy%x_d, dvdz%x_d, dwdx%x_d, dwdy%x_d, dwdz%x_d, &
         mu%x_d, n)

    call Ax_stress%compute_vector(f_x%x, f_y%x, f_z%x, u%x, v%x, w%x, &
         coef, p%msh, p%Xh)
    call opgrad(dudx%x, dudy%x, dudz%x, div_flux%x, coef)
    call compressible_ops_device_ns_flux_finalize(visc_m_x%x_d, &
         visc_m_y%x_d, visc_m_z%x_d, visc_E%x_d, f_x%x_d, f_y%x_d, &
         f_z%x_d, dudx%x_d, dudy%x_d, dudz%x_d, u%x_d, v%x_d, w%x_d, &
         coef%B_d, dissipation%x_d, n)

    call compressible_ops_device_ns_flux_temperature(div_flux%x_d, coef%h1_d, &
         p%x_d, rho_field%x_d, kappa%x_d, compressible_res_device_gamma, n)
    call Ax%compute(dudx%x, div_flux%x, coef, p%msh, p%Xh)
    call device_add2(visc_E%x_d, dudx%x_d, n)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine add_navier_stokes_flux_device

end module compressible_res_device
