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
  use operators, only : div, rotate_cyc
  use field_math, only : field_cmult
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  use device_math, only : device_copy, device_rone, &
       device_col2, device_cmult, device_sub2

  type, public, extends(euler_rhs_t) :: euler_res_device_t
   contains
     procedure, nopass :: step => advance_primitive_variables_device
     procedure, nopass :: evaluate_rhs_device
  end type euler_res_device_t

#ifdef HAVE_HIP
  interface
     subroutine euler_res_part_visc_hip(rhs_field_d, Binv_d, field_d, &
          effective_visc_d, n) &
          bind(c, name = 'euler_res_part_visc_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, effective_visc_d
       integer(c_int) :: n
     end subroutine euler_res_part_visc_hip
  end interface

  interface
     subroutine euler_res_part_mx_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mx_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mx_flux_hip
  end interface

  interface
     subroutine euler_res_part_my_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_my_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_my_flux_hip
  end interface

  interface
     subroutine euler_res_part_mz_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mz_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mz_flux_hip
  end interface

  interface
     subroutine euler_res_part_E_flux_hip(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'euler_res_part_E_flux_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine euler_res_part_E_flux_hip
  end interface

  interface
     subroutine euler_res_part_coef_mult_hip(rhs_rho_field_d, rhs_m_x_d, &
          rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'euler_res_part_coef_mult_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine euler_res_part_coef_mult_hip
  end interface

  interface
     subroutine euler_res_part_rk_sum_hip(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, b_i, n) &
          bind(c, name = 'euler_res_part_rk_sum_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, b_i
       integer(c_int) :: n
     end subroutine euler_res_part_rk_sum_hip
  end interface
#elif HAVE_CUDA
  interface
     subroutine euler_res_part_visc_cuda(rhs_field_d, Binv_d, field_d, &
          effective_visc_d, n) &
          bind(c, name = 'euler_res_part_visc_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, effective_visc_d
       integer(c_int) :: n
     end subroutine euler_res_part_visc_cuda
  end interface

  interface
     subroutine euler_res_part_mx_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mx_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mx_flux_cuda
  end interface

  interface
     subroutine euler_res_part_my_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_my_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_my_flux_cuda
  end interface

  interface
     subroutine euler_res_part_mz_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mz_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mz_flux_cuda
  end interface

  interface
     subroutine euler_res_part_E_flux_cuda(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'euler_res_part_E_flux_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine euler_res_part_E_flux_cuda
  end interface

  interface
     subroutine euler_res_part_coef_mult_cuda(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'euler_res_part_coef_mult_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine euler_res_part_coef_mult_cuda
  end interface

  interface
     subroutine euler_res_part_rk_sum_cuda(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, c, n) &
          bind(c, name = 'euler_res_part_rk_sum_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, c
       integer(c_int) :: n
     end subroutine euler_res_part_rk_sum_cuda
  end interface
#elif HAVE_OPENCL
  interface
     subroutine euler_res_part_visc_opencl(rhs_field_d, Binv_d, field_d, &
          effective_visc_d, n) &
          bind(c, name = 'euler_res_part_visc_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rhs_field_d, Binv_d, field_d, effective_visc_d
       integer(c_int) :: n
     end subroutine euler_res_part_visc_opencl
  end interface

  interface
     subroutine euler_res_part_mx_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mx_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mx_flux_opencl
  end interface

  interface
     subroutine euler_res_part_my_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_my_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_my_flux_opencl
  end interface

  interface
     subroutine euler_res_part_mz_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, n) &
          bind(c, name = 'euler_res_part_mz_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p
       integer(c_int) :: n
     end subroutine euler_res_part_mz_flux_opencl
  end interface

  interface
     subroutine euler_res_part_E_flux_opencl(f_x, f_y, f_z, &
          m_x, m_y, m_z, rho_field, p, E, n) &
          bind(c, name = 'euler_res_part_E_flux_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, E
       integer(c_int) :: n
     end subroutine euler_res_part_E_flux_opencl
  end interface

  interface
     subroutine euler_res_part_coef_mult_opencl(rhs_rho_field_d, &
          rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
          rhs_E_d, mult_d, n) &
          bind(c, name = 'euler_res_part_coef_mult_opencl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rhs_rho_field_d, rhs_m_x_d, rhs_m_y_d, rhs_m_z_d, &
            rhs_E_d, mult_d
       integer(c_int) :: n
     end subroutine euler_res_part_coef_mult_opencl
  end interface

  interface
     subroutine euler_res_part_rk_sum_opencl(rho, m_x, m_y, m_z, E, &
          k_rho_i, k_m_x_i, k_m_y_i, &
          k_m_z_i, k_E_i, &
          dt, c, n) &
          bind(c, name = 'euler_res_part_rk_sum_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: rho, m_x, m_y, m_z, E, &
            k_rho_i, k_m_x_i, k_m_y_i, &
            k_m_z_i, k_E_i
       real(c_rp) :: dt, c
       integer(c_int) :: n
     end subroutine euler_res_part_rk_sum_opencl
  end interface
#endif

contains
  subroutine advance_primitive_variables_device(rho_field, &
       m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, effective_visc, rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, effective_visc
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, k
    real(kind=rp) :: t, c
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E
    integer :: temp_indices(25)
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
    call neko_scratch_registry%request_field(temp_rho, temp_indices(21), .false.)
    call neko_scratch_registry%request_field(temp_m_x, temp_indices(22), .false.)
    call neko_scratch_registry%request_field(temp_m_y, temp_indices(23), .false.)
    call neko_scratch_registry%request_field(temp_m_z, temp_indices(24), .false.)
    call neko_scratch_registry%request_field(temp_E, temp_indices(25), .false.)

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

    ! Runge-Kutta stages
    do i = 1, s
       call device_copy(temp_rho%x_d, rho_field%x_d, n)
       call device_copy(temp_m_x%x_d, m_x%x_d, n)
       call device_copy(temp_m_y%x_d, m_y%x_d, n)
       call device_copy(temp_m_z%x_d, m_z%x_d, n)
       call device_copy(temp_E%x_d, E%x_d, n)

       do j = 1, i-1
#ifdef HAVE_HIP
          call euler_res_part_rk_sum_hip(temp_rho%x_d, temp_m_x%x_d, temp_m_y%x_d, &
               temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_CUDA
          call euler_res_part_rk_sum_cuda(temp_rho%x_d, temp_m_x%x_d, temp_m_y%x_d, &
               temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#elif HAVE_OPENCL
          call euler_res_part_rk_sum_opencl(temp_rho%x_d, temp_m_x%x_d, temp_m_y%x_d, &
               temp_m_z%x_d, temp_E%x_d, &
               k_rho%items(j)%ptr%x_d, k_m_x%items(j)%ptr%x_d, k_m_y%items(j)%ptr%x_d, &
               k_m_z%items(j)%ptr%x_d, k_E%items(j)%ptr%x_d, &
               dt, rk_scheme%coeffs_A(i, j), n)
#endif
       end do

       ! Compute f(U) = rhs(U) for the intermediate values
       call evaluate_rhs_device(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            p, u, v, w, Ax, &
            coef, gs, h, effective_visc)
    end do

    ! Update the solution
    do i = 1, s
#ifdef HAVE_HIP
       call euler_res_part_rk_sum_hip(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#elif HAVE_CUDA
       call euler_res_part_rk_sum_cuda(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#elif HAVE_OPENCL
       call euler_res_part_rk_sum_opencl(rho_field%x_d, &
            m_x%x_d, m_y%x_d, m_z%x_d, E%x_d, &
            k_rho%items(i)%ptr%x_d, k_m_x%items(i)%ptr%x_d, k_m_y%items(i)%ptr%x_d, &
            k_m_z%items(i)%ptr%x_d, k_E%items(i)%ptr%x_d, &
            dt, rk_scheme%coeffs_b(i), n)
#endif
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine advance_primitive_variables_device

  subroutine evaluate_rhs_device(rhs_rho_field, rhs_m_x, rhs_m_y, &
       rhs_m_z, rhs_E, rho_field, &
       m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, effective_visc)
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, effective_visc
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n, i
    real(kind=rp) :: visc_coeff
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
    call euler_res_part_mx_flux_hip(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_CUDA
    call euler_res_part_mx_flux_cuda(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
#elif HAVE_OPENCL
    call euler_res_part_mx_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, rho_field%x_d, p%x_d, n)
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
    call euler_res_part_my_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
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
    call euler_res_part_mz_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, n)
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
    call euler_res_part_E_flux_opencl(f_x%x_d, f_y%x_d, f_z%x_d, &
         m_x%x_d, m_y%x_d, m_z%x_d, &
         rho_field%x_d, p%x_d, E%x_d, n)
#endif
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
    call gs%op(rhs_m_x, GS_OP_ADD)
    call gs%op(rhs_m_y, GS_OP_ADD)
    call gs%op(rhs_m_z, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD)

#ifdef HAVE_HIP
    call euler_res_part_coef_mult_hip(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#elif HAVE_CUDA
    call euler_res_part_coef_mult_cuda(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#elif HAVE_OPENCL
    call euler_res_part_coef_mult_opencl(rhs_rho_field%x_d, rhs_m_x%x_d, &
         rhs_m_y%x_d, rhs_m_z%x_d, &
         rhs_E%x_d, coef%mult_d, n)
#endif

    call neko_scratch_registry%request_field(visc_rho, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_m_x, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(visc_m_y, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(visc_m_z, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(visc_E, temp_indices(8), .false.)

    ! Set h1 coefficient to the effective viscosity for the Laplacian operator
    call device_copy(coef%h1_d, effective_visc%x_d, n)

    ! Calculate artificial diffusion with variable viscosity
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    ! Reset h1 coefficient back to 1.0 for other operations
    call device_rone(coef%h1_d, n)

    call gs%op(visc_rho, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 1, coef)
    call gs%op(visc_m_x, GS_OP_ADD)
    call gs%op(visc_m_y, GS_OP_ADD)
    call gs%op(visc_m_z, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 0, coef)
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

end module euler_res_device
