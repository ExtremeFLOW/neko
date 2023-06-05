! Copyright (c) 2018-2022, The Neko Authors
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
!> Routines to generate the right-hand sides for the convection-diffusion
!! equation. Employs the EXT/BDF time integration schemes to compute
!! the contributions coming from the explicitly extrapolated convective term
!! and the BDF scheme applied to the time derivative.
!! Inheritance is used to define implementation for different backends.
module rhs_maker
  use num_types
  use field_series, only : field_series_t
  use field, only : field_t
  implicit none 

  !> Abstract type to compute extrapolated velocity field for the pressure equation
  type, abstract :: rhs_maker_sumab_t
   contains
     procedure(rhs_maker_sumab), nopass, deferred :: compute_fluid

  end type rhs_maker_sumab_t

  !> Abstract type to sum up contributions to kth order extrapolation scheme
  type, abstract :: rhs_maker_ext_t
   contains
     procedure(rhs_maker_ext), nopass, deferred :: compute_fluid
     procedure(scalar_rhs_maker_ext), nopass, deferred :: compute_scalar
  end type rhs_maker_ext_t

  !> Abstract type to add contributions to F from lagged BD terms
  type, abstract :: rhs_maker_bdf_t
   contains
     procedure(rhs_maker_bdf), nopass, deferred :: compute_fluid
     procedure(scalar_rhs_maker_bdf), nopass, deferred :: compute_scalar
  end type rhs_maker_bdf_t

  abstract interface
     subroutine rhs_maker_sumab(u, v, w, uu, vv, ww, uulag, vvlag, wwlag, ab, nab)
       import field_t
       import field_series_t
       import rp
       type(field_t), intent(inout) :: u, v, w
       type(field_t), intent(inout) :: uu, vv, ww
       type(field_series_t), intent(inout) :: uulag, vvlag, wwlag
       real(kind=rp), dimension(3), intent(in) :: ab
       integer, intent(in) :: nab
     end subroutine rhs_maker_sumab
  end interface

  abstract interface
     subroutine rhs_maker_ext(temp1, temp2, temp3, fx_lag, fy_lag, fz_lag, &
          fx_laglag, fy_laglag, fz_laglag, fx, fy, fz, &
          rho, ext_coeffs, n)
       import field_t
       import rp
       type(field_t), intent(inout) :: temp1, temp2, temp3
       type(field_t), intent(inout) :: fx_lag, fy_lag, fz_lag
       type(field_t), intent(inout) :: fx_laglag, fy_laglag, fz_laglag
       real(kind=rp), intent(inout) :: rho, ext_coeffs(4)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: fx(n), fy(n), fz(n)
     end subroutine rhs_maker_ext
  end interface

  abstract interface
     subroutine scalar_rhs_maker_ext(temp1, fs_lag, fs_laglag, fs, rho, &
          ext_coeffs, n)
       import field_t
       import rp
       type(field_t), intent(inout) :: temp1
       type(field_t), intent(inout) :: fs_lag
       type(field_t), intent(inout) :: fs_laglag
       real(kind=rp), intent(inout) :: rho, ext_coeffs(4)
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: fs(n)
     end subroutine scalar_rhs_maker_ext
  end interface

  abstract interface
     subroutine rhs_maker_bdf(ta1, ta2, ta3, tb1, tb2, tb3, &
          ulag, vlag, wlag, bfx, bfy, bfz, &
          u, v, w, B, rho, dt, bd, nbd, n)
       import field_series_t
       import field_t
       import rp
       integer, intent(in) :: n, nbd
       type(field_t), intent(inout) :: ta1, ta2, ta3
       type(field_t), intent(in) :: u, v, w
       type(field_t), intent(inout) :: tb1, tb2, tb3
       type(field_series_t), intent(in) :: ulag, vlag, wlag        
       real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
       real(kind=rp), intent(in) :: B(n)
       real(kind=rp), intent(in) :: dt, rho, bd(4)
     end subroutine rhs_maker_bdf
  end interface

  abstract interface
     subroutine scalar_rhs_maker_bdf(temp1, temp2, s_lag, fs, s, B, rho, dt,&
          bd, nbd, n)
       import field_series_t
       import field_t
       import rp
       integer, intent(in) :: n, nbd
       type(field_t), intent(inout) :: temp1, temp2
       type(field_t), intent(in) :: s
       type(field_series_t), intent(in) :: s_lag
       real(kind=rp), intent(inout) :: fs(n)
       real(kind=rp), intent(in) :: B(n)
       real(kind=rp), intent(in) :: dt, rho, bd(4)
     end subroutine scalar_rhs_maker_bdf
  end interface

end module rhs_maker
