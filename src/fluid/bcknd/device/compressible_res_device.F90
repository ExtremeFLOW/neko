! Copyright (c) 2025, The Neko Authors
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
  use operators, only : div, grad, opgrad, rotate_cyc
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  use compressible_ops_device, only : compressible_ops_device_update_uvw, &
       compressible_ops_device_update_mxyz_p_ruvw, &
       compressible_ops_device_ns_flux_prepare, &
       compressible_ops_device_ns_flux_finalize, &
       compressible_ops_device_ns_flux_temperature
  use device_math, only : device_copy, device_rone, device_col2, &
       device_cmult, device_sub2, device_add2
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t

  type, public, extends(compressible_rhs_t) :: compressible_res_device_t
   contains
     procedure, nopass :: step => advance_primitive_variables_device
     procedure, nopass :: evaluate_rhs => evaluate_rhs_device
  end type compressible_res_device_t

  !> Whether physical Navier-Stokes fluxes are active for the current step.
  logical :: compressible_res_device_add_physical_flux = .false.
  !> Whether physical viscous stress is active for the current step.
  logical :: compressible_res_device_add_physical_stress = .false.
  !> Module variable to store thermodynamic parameter set by factory.
  real(kind=rp), public :: compressible_res_device_gamma = 1.4_rp

#ifdef HAVE_HIP
  interface
     subroutine compressible_res_part_visc_hip(rhs_field_d, Binv_d, field_d, &
          artificial_visc_d, n) &
          bind(c, name = 'compressible_res_part_visc_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, artificial_visc_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_hip
  end interface

  interface
     subroutine inviscid_res_part_mx_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mx_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mx_flux_hip
  end interface

  interface
     subroutine inviscid_res_part_my_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_my_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_my_flux_hip
  end interface

  interface
     subroutine inviscid_res_part_mz_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mz_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mz_flux_hip
  end interface

  interface
     subroutine inviscid_res_part_E_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'inviscid_res_part_E_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine inviscid_res_part_E_flux_hip
  end interface

  interface
     subroutine compressible_res_part_coef_mult_hip(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'compressible_res_part_coef_mult_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine compressible_res_part_coef_mult_hip
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
     subroutine compressible_res_part_visc_cuda(rhs_field_d, Binv_d, field_d, &
          artificial_visc_d, n) &
          bind(c, name = 'compressible_res_part_visc_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, artificial_visc_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_cuda
  end interface

  interface
     subroutine inviscid_res_part_mx_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mx_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mx_flux_cuda
  end interface

  interface
     subroutine inviscid_res_part_my_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_my_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_my_flux_cuda
  end interface

  interface
     subroutine inviscid_res_part_mz_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mz_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mz_flux_cuda
  end interface

  interface
     subroutine inviscid_res_part_E_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'inviscid_res_part_E_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine inviscid_res_part_E_flux_cuda
  end interface

  interface
     subroutine compressible_res_part_coef_mult_cuda(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'compressible_res_part_coef_mult_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine compressible_res_part_coef_mult_cuda
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
     subroutine compressible_res_part_visc_opencl(rhs_field_d, Binv_d, &
          field_d, artificial_visc_d, n) &
          bind(c, name = 'compressible_res_part_visc_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, artificial_visc_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_opencl
  end interface

  interface
     subroutine inviscid_res_part_mx_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mx_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mx_flux_opencl
  end interface

  interface
     subroutine inviscid_res_part_my_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_my_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_my_flux_opencl
  end interface

  interface
     subroutine inviscid_res_part_mz_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mz_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mz_flux_opencl
  end interface

  interface
     subroutine inviscid_res_part_E_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'inviscid_res_part_E_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine inviscid_res_part_E_flux_opencl
  end interface

  interface
     subroutine compressible_res_part_coef_mult_opencl(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'compressible_res_part_coef_mult_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine compressible_res_part_coef_mult_opencl
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
     subroutine compressible_res_part_visc_metal(rhs_field_d, Binv_d, field_d, &
          effective_visc_d, n) &
          bind(c, name = 'compressible_res_part_visc_metal')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, effective_visc_d
       integer(c_int) :: n
     end subroutine compressible_res_part_visc_metal
  end interface

  interface
     subroutine inviscid_res_part_mx_flux_metal(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mx_flux_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mx_flux_metal
  end interface

  interface
     subroutine inviscid_res_part_my_flux_metal(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_my_flux_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_my_flux_metal
  end interface

  interface
     subroutine inviscid_res_part_mz_flux_metal(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'inviscid_res_part_mz_flux_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine inviscid_res_part_mz_flux_metal
  end interface

  interface
     subroutine inviscid_res_part_E_flux_metal(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'inviscid_res_part_E_flux_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine inviscid_res_part_E_flux_metal
  end interface

  interface
     subroutine compressible_res_part_coef_mult_metal(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'compressible_res_part_coef_mult_metal')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine compressible_res_part_coef_mult_metal
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
       Ax_stress, coef, gs, h, artificial_visc, mu, kappa, bcs_vel, time, &
       rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bcs_vel
    type(time_state_t), intent(in) :: time
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, k
    real(kind=rp) :: t, c
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
         temp_p, temp_u, temp_v, temp_w, temp_ruvw
    integer :: temp_indices(30)
    type(field_list_t) :: k_rho, k_m_x, k_m_y, k_m_z, k_E

    n = p%dof%size()
    s = rk_scheme%order
    call neko_scratch_registry%request_field(k_rho_1, temp_indices(1), .true.)
    call neko_scratch_registry%request_field(k_rho_2, temp_indices(2), .true.)
    call neko_scratch_registry%request_field(k_rho_3, temp_indices(3), .true.)
    call neko_scratch_registry%request_field(k_rho_4, temp_indices(4), .true.)
    call neko_scratch_registry%request_field(k_m_x_1, temp_indices(5), .true.)
    call neko_scratch_registry%request_field(k_m_x_2, temp_indices(6), .true.)
    call neko_scratch_registry%request_field(k_m_x_3, temp_indices(7), .true.)
    call neko_scratch_registry%request_field(k_m_x_4, temp_indices(8), .true.)
    call neko_scratch_registry%request_field(k_m_y_1, temp_indices(9), .true.)
    call neko_scratch_registry%request_field(k_m_y_2, temp_indices(10), .true.)
    call neko_scratch_registry%request_field(k_m_y_3, temp_indices(11), .true.)
    call neko_scratch_registry%request_field(k_m_y_4, temp_indices(12), .true.)
    call neko_scratch_registry%request_field(k_m_z_1, temp_indices(13), .true.)
    call neko_scratch_registry%request_field(k_m_z_2, temp_indices(14), .true.)
    call neko_scratch_registry%request_field(k_m_z_3, temp_indices(15), .true.)
    call neko_scratch_registry%request_field(k_m_z_4, temp_indices(16), .true.)
    call neko_scratch_registry%request_field(k_E_1, temp_indices(17), .true.)
    call neko_scratch_registry%request_field(k_E_2, temp_indices(18), .true.)
    call neko_scratch_registry%request_field(k_E_3, temp_indices(19), .true.)
    call neko_scratch_registry%request_field(k_E_4, temp_indices(20), .true.)
    call neko_scratch_registry%request_field(temp_rho, temp_indices(21), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_x, temp_indices(22), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_y, temp_indices(23), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_z, temp_indices(24), &
         .false.)
    call neko_scratch_registry%request_field(temp_E, temp_indices(25), .false.)
    call neko_scratch_registry%request_field(temp_p, temp_indices(26), .false.)
    call neko_scratch_registry%request_field(temp_u, temp_indices(27), .false.)
    call neko_scratch_registry%request_field(temp_v, temp_indices(28), .false.)
    call neko_scratch_registry%request_field(temp_w, temp_indices(29), .false.)
    call neko_scratch_registry%request_field(temp_ruvw, temp_indices(30), &
         .false.)

    call k_rho%init(4)
    call k_rho%assign(1, k_rho_1)
    call k_rho%assign(2, k_rho_2)
    call k_rho%assign(3, k_rho_3)
    call k_rho%assign(4, k_rho_4)
    call k_m_x%init(4)
    call k_m_x%assign(1, k_m_x_1)
    call k_m_x%assign(2, k_m_x_2)
    call k_m_x%assign(3, k_m_x_3)
    call k_m_x%assign(4, k_m_x_4)
    call k_m_y%init(4)
    call k_m_y%assign(1, k_m_y_1)
    call k_m_y%assign(2, k_m_y_2)
    call k_m_y%assign(3, k_m_y_3)
    call k_m_y%assign(4, k_m_y_4)
    call k_m_z%init(4)
    call k_m_z%assign(1, k_m_z_1)
    call k_m_z%assign(2, k_m_z_2)
    call k_m_z%assign(3, k_m_z_3)
    call k_m_z%assign(4, k_m_z_4)
    call k_E%init(4)
    call k_E%assign(1, k_E_1)
    call k_E%assign(2, k_E_2)
    call k_E%assign(3, k_E_3)
    call k_E%assign(4, k_E_4)

    compressible_res_device_add_physical_flux = &
         any(mu%x .ne. 0.0_rp) .or. any(kappa%x .ne. 0.0_rp)
    compressible_res_device_add_physical_stress = any(mu%x .ne. 0.0_rp)

    ! Runge-Kutta stages
    do i = 1, s
       call device_copy(temp_rho%x_d, rho_field%x_d, n)
       call device_copy(temp_m_x%x_d, m_x%x_d, n)
       call device_copy(temp_m_y%x_d, m_y%x_d, n)
       call device_copy(temp_m_z%x_d, m_z%x_d, n)
       call device_copy(temp_E%x_d, E%x_d, n)

       do j = 1, i-1
#ifdef HAVE_HIP
          call compressible_res_part_rk_sum_hip(temp_rho%x_d, temp_m_x%x_d, &
               temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
               k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_CUDA
          call compressible_res_part_rk_sum_cuda(temp_rho%x_d, temp_m_x%x_d, &
               temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
               k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_OPENCL
          call compressible_res_part_rk_sum_opencl(temp_rho%x_d, temp_m_x%x_d, &
               temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
               k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_METAL
          call compressible_res_part_rk_sum_metal(temp_rho%x_d, &
              temp_m_x%x_d, temp_m_y%x_d, temp_m_z%x_d, temp_E%x_d, &
              k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, &
              k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#endif
       end do

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
            Ax_stress, coef, gs, h, artificial_visc, mu, kappa)
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

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine advance_primitive_variables_device

  subroutine evaluate_rhs_device(rhs_rho_field, rhs_m_x, rhs_m_y, &
       rhs_m_z, rhs_E, rho_field, &
       m_x, m_y, m_z, E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, h, artificial_visc, mu, kappa)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, &
         rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n, i
    type(field_t), pointer :: temp, f_x, f_y, f_z, &
         visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: temp_indices(8)

    n = coef%dof%size()
    call neko_scratch_registry%request_field(f_x, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(f_y, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(f_z, temp_indices(3), .false.)

    !> rho = rho - dt * div(m)
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

    !> m = m - dt * div(rho * u * u^T + p*I)
    ! m_x
#ifdef HAVE_HIP
    call inviscid_res_part_mx_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call inviscid_res_part_mx_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call inviscid_res_part_mx_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_METAL
    call inviscid_res_part_mx_flux_metal(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#endif
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_y
#ifdef HAVE_HIP
    call inviscid_res_part_my_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call inviscid_res_part_my_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call inviscid_res_part_my_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_METAL
    call inviscid_res_part_my_flux_metal(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#endif
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_z
#ifdef HAVE_HIP
    call inviscid_res_part_mz_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call inviscid_res_part_mz_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call inviscid_res_part_mz_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#elif HAVE_METAL
    call inviscid_res_part_mz_flux_metal(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
#endif
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

    !> E = E - dt * div(u * (E + p))
    ! Inviscid energy flux for both NS and monolithic paths.
    ! Viscous energy diffusion is handled entirely by Ax(E) below.
#ifdef HAVE_HIP
    call inviscid_res_part_E_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, E%x_d, n)
#elif HAVE_CUDA
    call inviscid_res_part_E_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, E%x_d, n)
#elif HAVE_OPENCL
    call inviscid_res_part_E_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, E%x_d, n)
#elif HAVE_METAL
    call inviscid_res_part_E_flux_metal(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, E%x_d, n)
#endif
    call div(rhs_E%x_d, f_x%x_d, f_y%x_d, f_z%x_d, coef)

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, 1, coef)
    call gs%op(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, n, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD)

#ifdef HAVE_HIP
    call compressible_res_part_coef_mult_hip(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#elif HAVE_CUDA
    call compressible_res_part_coef_mult_cuda(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#elif HAVE_OPENCL
    call compressible_res_part_coef_mult_opencl(rhs_rho_field%x_d, &
         rhs_m_x%x_d, rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#elif HAVE_METAL
    call compressible_res_part_coef_mult_metal(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#endif

    call neko_scratch_registry%request_field(visc_rho, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_m_x, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(visc_m_y, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(visc_m_z, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(visc_E, temp_indices(8), .false.)

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

    ! Reset h1 coefficient back to 1.0 for other operations
    call device_rone(coef%h1_d, n)

    call gs%op(visc_rho, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x_d, visc_m_y%x_d, visc_m_z%x_d, 1, coef)
    call gs%op(visc_m_x%x, visc_m_y%x, visc_m_z%x, n, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x_d, visc_m_y%x_d, visc_m_z%x_d, 0, coef)
    call gs%op(visc_E, GS_OP_ADD)

    ! Apply artificial viscosity - the coefficient is already in the Laplacian
    ! rhs = -rhs - Binv * visc_lap
    call device_col2(visc_rho%x_d, coef%Binv_d, n)
    call device_col2(visc_m_x%x_d, coef%Binv_d, n)
    call device_col2(visc_m_y%x_d, coef%Binv_d, n)
    call device_col2(visc_m_z%x_d, coef%Binv_d, n)
    call device_col2(visc_E%x_d, coef%Binv_d, n)

    call device_cmult(rhs_rho_field%x_d, -1.0_rp, n)
    call device_sub2(rhs_rho_field%x_d, visc_rho%x_d, n)

    call device_cmult(rhs_m_x%x_d, -1.0_rp, n)
    call device_sub2(rhs_m_x%x_d, visc_m_x%x_d, n)

    call device_cmult(rhs_m_y%x_d, -1.0_rp, n)
    call device_sub2(rhs_m_y%x_d, visc_m_y%x_d, n)

    call device_cmult(rhs_m_z%x_d, -1.0_rp, n)
    call device_sub2(rhs_m_z%x_d, visc_m_z%x_d, n)

    call device_cmult(rhs_E%x_d, -1.0_rp, n)
    call device_sub2(rhs_E%x_d, visc_E%x_d, n)

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

    call grad(dudx%x, dudy%x, dudz%x, u%x, coef)
    call grad(dvdx%x, dvdy%x, dvdz%x, v%x, coef)
    call grad(dwdx%x, dwdy%x, dwdz%x, w%x, coef)

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
