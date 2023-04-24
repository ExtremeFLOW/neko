! Copyright (c) 2022, The Neko Authors
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
module rhs_maker_device
  use rhs_maker
  use device
  use utils
  use field_series, only : field_series_t
  use field, only : field_t
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, extends(rhs_maker_sumab_t) :: rhs_maker_sumab_device_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_sumab_device
  end type rhs_maker_sumab_device_t

  type, public, extends(rhs_maker_ext_t) ::  rhs_maker_ext_device_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_ext_device
     procedure, nopass :: compute_scalar => scalar_rhs_maker_ext_device
  end type rhs_maker_ext_device_t

  type, public, extends(rhs_maker_bdf_t) :: rhs_maker_bdf_device_t
   contains
     procedure, nopass :: compute_fluid => rhs_maker_bdf_device
     procedure, nopass :: compute_scalar => scalar_rhs_maker_bdf_device
  end type rhs_maker_bdf_device_t

#ifdef HAVE_HIP
    interface
     subroutine rhs_maker_sumab_hip(u_d, v_d, w_d, uu_d, vv_d, ww_d, &
          uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2, ab1, ab2, ab3, nab, n)&
          bind(c, name='rhs_maker_sumab_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, uu_d, vv_d, ww_d
       type(c_ptr), value :: uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2
       real(c_rp) :: ab1, ab2, ab3
       integer(c_int) :: nab, n
     end subroutine rhs_maker_sumab_hip
  end interface

  interface
     subroutine rhs_maker_ext_hip(abx1_d, aby1_d, abz1_d, &
                                  abx2_d, aby2_d, abz2_d, &
                                  bfx_d, bfy_d, bfz_d, &
                                  rho, ab1, ab2, ab3, n) &
                                  bind(c, name='rhs_maker_ext_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: abx1_d, aby1_d, abz1_d 
       type(c_ptr), value :: abx2_d, aby2_d, abz2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d
       real(c_rp) :: rho, ab1, ab2, ab3
       integer(c_int) :: n
     end subroutine rhs_maker_ext_hip
  end interface

  interface
     subroutine scalar_rhs_maker_ext_hip(fs_lag_d, fs_laglag_d, fs_d, rho, &
                                          ext1, ext2, ext3, n) &
                                          bind(c, name='scalar_rhs_maker_ext_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: fs_lag_d, fs_laglag_d, fs_d 
       real(c_rp) :: rho, ext1, ext2, ext3
       integer(c_int) :: n
     end subroutine scalar_rhs_maker_ext_hip
  end interface

  interface
     subroutine rhs_maker_bdf_hip(ulag1_d, ulag2_d, vlag1_d, vlag2_d, &
          wlag1_d, wlag2_d, bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) bind(c, name='rhs_maker_bdf_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: ulag1_d, ulag2_d, vlag1_d
       type(c_ptr), value :: vlag2_d, wlag1_d, wlag2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine rhs_maker_bdf_hip
  end interface 

  interface
     subroutine scalar_rhs_maker_bdf_hip(s_lag_d, s_laglag_d, fs_d, s_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) &
          bind(c, name='scalar_rhs_maker_bdf_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: s_lag_d, s_laglag_d
       type(c_ptr), value :: fs_d, s_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine scalar_rhs_maker_bdf_hip
  end interface  
#elif HAVE_CUDA
  interface
     subroutine rhs_maker_sumab_cuda(u_d, v_d, w_d, uu_d, vv_d, ww_d, &
          uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2, ab1, ab2, ab3, nab, n)&
          bind(c, name='rhs_maker_sumab_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, uu_d, vv_d, ww_d
       type(c_ptr), value :: uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2
       real(c_rp) :: ab1, ab2, ab3
       integer(c_int) :: nab, n
     end subroutine rhs_maker_sumab_cuda
  end interface

  interface
     subroutine rhs_maker_ext_cuda(abx1_d, aby1_d, abz1_d, &
                                   abx2_d, aby2_d, abz2_d, &
                                   bfx_d, bfy_d, bfz_d, &
                                   rho, ab1, ab2, ab3, n) &
                                   bind(c, name='rhs_maker_ext_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: abx1_d, aby1_d, abz1_d 
       type(c_ptr), value :: abx2_d, aby2_d, abz2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d
       real(c_rp) :: rho, ab1, ab2, ab3
       integer(c_int) :: n
     end subroutine rhs_maker_ext_cuda
  end interface

  interface
     subroutine scalar_rhs_maker_ext_cuda(fs_lag_d, fs_laglag_d, fs_d, rho, &
                                          ext1, ext2, ext3, n) &
                                          bind(c, name='scalar_rhs_maker_ext_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: fs_lag_d, fs_laglag_d, fs_d 
       real(c_rp) :: rho, ext1, ext2, ext3
       integer(c_int) :: n
     end subroutine scalar_rhs_maker_ext_cuda
  end interface

  interface
     subroutine rhs_maker_bdf_cuda(ulag1_d, ulag2_d, vlag1_d, vlag2_d, &
          wlag1_d, wlag2_d, bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) bind(c, name='rhs_maker_bdf_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: ulag1_d, ulag2_d, vlag1_d
       type(c_ptr), value :: vlag2_d, wlag1_d, wlag2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine rhs_maker_bdf_cuda
  end interface  
  
  interface
     subroutine scalar_rhs_maker_bdf_cuda(s_lag_d, s_laglag_d, fs_d, s_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) &
          bind(c, name='scalar_rhs_maker_bdf_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: s_lag_d, s_laglag_d
       type(c_ptr), value :: fs_d, s_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine scalar_rhs_maker_bdf_cuda
  end interface  
#elif HAVE_OPENCL
  interface
     subroutine rhs_maker_sumab_opencl(u_d, v_d, w_d, uu_d, vv_d, ww_d, &
          uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2, ab1, ab2, ab3, nab, n)&
          bind(c, name='rhs_maker_sumab_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, uu_d, vv_d, ww_d
       type(c_ptr), value :: uulag1, uulag2, vvlag1, vvlag2, wwlag1, wwlag2
       real(c_rp) :: ab1, ab2, ab3
       integer(c_int) :: nab, n
     end subroutine rhs_maker_sumab_opencl
  end interface

  interface
     subroutine rhs_maker_ext_opencl(abx1_d, aby1_d, abz1_d, &
                                     abx2_d, aby2_d, abz2_d, &
                                     bfx_d, bfy_d, bfz_d, &
                                     rho, ab1, ab2, ab3, n) &
                                     bind(c, name='rhs_maker_ext_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: abx1_d, aby1_d, abz1_d 
       type(c_ptr), value :: abx2_d, aby2_d, abz2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d
       real(c_rp) :: rho, ab1, ab2, ab3
       integer(c_int) :: n
     end subroutine rhs_maker_ext_opencl
  end interface

  interface
     subroutine scalar_rhs_maker_ext_opencl(fs_lag_d, fs_laglag_d, fs_d, rho, &
                                          ext1, ext2, ext3, n) &
                                          bind(c, name='scalar_rhs_maker_ext_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: fs_lag_d, fs_laglag_d, fs_d 
       real(c_rp) :: rho, ext1, ext2, ext3
       integer(c_int) :: n
     end subroutine scalar_rhs_maker_ext_opencl
  end interface

  interface
     subroutine rhs_maker_bdf_opencl(ulag1_d, ulag2_d, vlag1_d, vlag2_d, &
          wlag1_d, wlag2_d, bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) bind(c, name='rhs_maker_bdf_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: ulag1_d, ulag2_d, vlag1_d
       type(c_ptr), value :: vlag2_d, wlag1_d, wlag2_d
       type(c_ptr), value :: bfx_d, bfy_d, bfz_d, u_d, v_d, w_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine rhs_maker_bdf_opencl
  end interface 

  interface
     subroutine scalar_rhs_maker_bdf_opencl(s_lag_d, s_laglag_d, fs_d, s_d, B_d, &
          rho, dt, bd2, bd3, bd4, nbd, n) &
          bind(c, name='scalar_rhs_maker_bdf_opencl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: s_lag_d, s_laglag_d
       type(c_ptr), value :: fs_d, s_d, B_d
       reaL(c_rp) :: rho, dt, bd2, bd3, bd4
       integer(c_int) :: nbd, n
     end subroutine scalar_rhs_maker_bdf_opencl
  end interface  
#endif

contains

    subroutine rhs_maker_sumab_device(u, v, w, uu, vv, ww, uulag, vvlag, wwlag, ab, nab)
    type(field_t), intent(inout) :: u,v, w
    type(field_t), intent(inout) :: uu, vv, ww
    type(field_series_t), intent(inout) :: uulag, vvlag, wwlag
    real(kind=rp), dimension(3), intent(in) :: ab
    integer, intent(in) :: nab

#ifdef HAVE_HIP
    call rhs_maker_sumab_hip(u%x_d, v%x_d, w%x_d, uu%x_d, vv%x_d, ww%x_d, &
         uulag%lf(1)%x_d, uulag%lf(2)%x_d, vvlag%lf(1)%x_d, vvlag%lf(2)%x_d, &
         wwlag%lf(1)%x_d, wwlag%lf(2)%x_d, ab(1), ab(2), ab(3), nab, &
         uu%dof%size())
#elif HAVE_CUDA    
    call rhs_maker_sumab_cuda(u%x_d, v%x_d, w%x_d, uu%x_d, vv%x_d, ww%x_d, &
         uulag%lf(1)%x_d, uulag%lf(2)%x_d, vvlag%lf(1)%x_d, vvlag%lf(2)%x_d, &
         wwlag%lf(1)%x_d, wwlag%lf(2)%x_d, ab(1), ab(2), ab(3), nab, &
         uu%dof%size())
#elif HAVE_OPENCL
    call rhs_maker_sumab_opencl(u%x_d, v%x_d, w%x_d, uu%x_d, vv%x_d, ww%x_d, &
         uulag%lf(1)%x_d, uulag%lf(2)%x_d, vvlag%lf(1)%x_d, vvlag%lf(2)%x_d, &
         wwlag%lf(1)%x_d, wwlag%lf(2)%x_d, ab(1), ab(2), ab(3), nab, &
         uu%dof%size())
#endif
    
  end subroutine rhs_maker_sumab_device

  subroutine rhs_maker_ext_device(temp1, temp2, temp3, fx_lag, fy_lag, fz_lag, &
                           fx_laglag, fy_laglag, fz_laglag, fx, fy, fz, &
                           rho, ext_coeffs, n)
    type(field_t), intent(inout) :: temp1, temp2, temp3
    type(field_t), intent(inout) :: fx_lag, fy_lag, fz_lag
    type(field_t), intent(inout) :: fx_laglag, fy_laglag, fz_laglag
    real(kind=rp), intent(inout) :: rho, ext_coeffs(4)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: fx(n), fy(n), fz(n)
    type(c_ptr) :: fx_d, fy_d, fz_d

    fx_d = device_get_ptr(fx)
    fy_d = device_get_ptr(fy)
    fz_d = device_get_ptr(fz)

#ifdef HAVE_HIP
    call rhs_maker_ext_hip(fx_lag%x_d, fy_lag%x_d, fz_lag%x_d, &
                           fx_laglag%x_d, fy_laglag%x_d, fz_laglag%x_d, &
                           fx_d, fy_d, fz_d, rho, &
                           ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#elif HAVE_CUDA
    call rhs_maker_ext_cuda(fx_lag%x_d, fy_lag%x_d, fz_lag%x_d, &
                            fx_laglag%x_d, fy_laglag%x_d, fz_laglag%x_d, &
                            fx_d, fy_d, fz_d, rho, &
                            ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#elif HAVE_OPENCL
    call rhs_maker_ext_opencl(fx_lag%x_d, fy_lag%x_d, fz_lag%x_d, &
                              fx_laglag%x_d, fy_laglag%x_d, fz_laglag%x_d, &
                              fx_d, fy_d, fz_d, rho, &
                              ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#endif
    
  end subroutine rhs_maker_ext_device

  subroutine scalar_rhs_maker_ext_device(temp1, fs_lag, fs_laglag, fs, &
                           rho, ext_coeffs, n)
    type(field_t), intent(inout) :: temp1
    type(field_t), intent(inout) :: fs_lag
    type(field_t), intent(inout) :: fs_laglag
    real(kind=rp), intent(inout) :: rho, ext_coeffs(4)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: fs(n)
    type(c_ptr) :: fs_d

    fs_d = device_get_ptr(fs)

#ifdef HAVE_HIP
    call scalar_rhs_maker_ext_hip(fs_lag%x_d, fs_laglag%x_d, fs_d, rho, &
                              ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#elif HAVE_CUDA
    call scalar_rhs_maker_ext_cuda(fs_lag%x_d, fs_laglag%x_d, fs_d, rho, &
                              ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#elif HAVE_OPENCL
    call scalar_rhs_maker_ext_opencl(fs_lag%x_d, fs_laglag%x_d, fs_d, rho, &
                              ext_coeffs(1), ext_coeffs(2), ext_coeffs(3), n)
#endif
    
  end subroutine scalar_rhs_maker_ext_device

  subroutine rhs_maker_bdf_device(ta1, ta2, ta3, tb1, tb2, tb3, &
                               ulag, vlag, wlag, bfx, bfy, bfz, &
                               u, v, w, B, rho, dt, bd, nbd, n)    
    integer, intent(in) :: n, nbd
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: tb1, tb2, tb3
    type(field_series_t), intent(in) :: ulag, vlag, wlag        
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, rho, bd(4)
    type(c_ptr) :: bfx_d, bfy_d, bfz_d, B_d

    bfx_d = device_get_ptr(bfx)
    bfy_d = device_get_ptr(bfy)
    bfz_d = device_get_ptr(bfz)
    B_d = device_get_ptr(B)
    
#ifdef HAVE_HIP
    call rhs_maker_bdf_hip(ulag%lf(1)%x_d, ulag%lf(2)%x_d, &
                           vlag%lf(1)%x_d, vlag%lf(2)%x_d, &
                           wlag%lf(1)%x_d, wlag%lf(2)%x_d, &
                           bfx_d, bfy_d, bfz_d, u%x_d, v%x_d, w%x_d, &
                           B_d, rho, dt, bd(2), bd(3), bd(4), nbd, n)
#elif HAVE_CUDA
    call rhs_maker_bdf_cuda(ulag%lf(1)%x_d, ulag%lf(2)%x_d, &
                            vlag%lf(1)%x_d, vlag%lf(2)%x_d, &
                            wlag%lf(1)%x_d, wlag%lf(2)%x_d, &
                            bfx_d, bfy_d, bfz_d, u%x_d, v%x_d, w%x_d, &
                            B_d, rho, dt, bd(2), bd(3), bd(4), nbd, n)
#elif HAVE_OPENCL
    call rhs_maker_bdf_opencl(ulag%lf(1)%x_d, ulag%lf(2)%x_d, &
                              vlag%lf(1)%x_d, vlag%lf(2)%x_d, &
                              wlag%lf(1)%x_d, wlag%lf(2)%x_d, &
                              bfx_d, bfy_d, bfz_d, u%x_d, v%x_d, w%x_d, &
                              B_d, rho, dt, bd(2), bd(3), bd(4), nbd, n)
#endif

  end subroutine rhs_maker_bdf_device

  subroutine scalar_rhs_maker_bdf_device(temp1, temp2, s_lag, fs, s, B, rho, dt, &
                                         bd, nbd, n)    
    integer, intent(in) :: n, nbd
    type(field_t), intent(inout) :: temp1
    type(field_t), intent(in) :: s 
    type(field_t), intent(inout) :: temp2
    type(field_series_t), intent(in) :: s_lag
    real(kind=rp), intent(inout) :: fs(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, rho, bd(4)
    type(c_ptr) :: fs_d, B_d

    fs_d = device_get_ptr(fs)
    B_d = device_get_ptr(B)
    
#ifdef HAVE_HIP
    call scalar_rhs_maker_bdf_hip(s_lag%lf(1)%x_d, s_lag%lf(2)%x_d, &
                           fs_d, s%x_d, B_d, rho, dt, bd(2), bd(3), bd(4), &
                           nbd, n)
#elif HAVE_CUDA
    call scalar_rhs_maker_bdf_cuda(s_lag%lf(1)%x_d, s_lag%lf(2)%x_d, &
                           fs_d, s%x_d, B_d, rho, dt, bd(2), bd(3), bd(4), &
                           nbd, n)
#elif HAVE_OPENCL
    call scalar_rhs_maker_bdf_opencl(s_lag%lf(1)%x_d, s_lag%lf(2)%x_d, &
                           fs_d, s%x_d, B_d, rho, dt, bd(2), bd(3), bd(4), &
                           nbd, n)
#endif

  end subroutine scalar_rhs_maker_bdf_device
  
end module rhs_maker_device
