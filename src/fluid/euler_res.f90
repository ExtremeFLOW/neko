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
module euler_residual
  use gather_scatter, only : gs_t
  use ax_product, only : Ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use num_types, only : rp
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  implicit none
  private

  !> Abstract type to compute rhs
  type, public, abstract :: euler_rhs_t
   contains
     procedure(euler_rhs), nopass, deferred :: step
  end type euler_rhs_t

  !> Abstract interface to evaluate rhs
  abstract interface
     subroutine euler_rhs(rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
          coef, gs, h, c_avisc_low, rk_scheme, dt)
       import field_t
       import Ax_t
       import gs_t
       import coef_t
       import rp
       import runge_kutta_time_scheme_t
       type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
       type(field_t), intent(in) :: p, u, v, w, h
       class(Ax_t), intent(inout) :: Ax
       type(coef_t), intent(inout) :: coef
       type(gs_t), intent(inout) :: gs
       real(kind=rp) :: c_avisc_low
       class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
       real(kind=rp), intent(in) :: dt
     end subroutine euler_rhs
  end interface

  !> Abstract interface to choose bcknd for rhs evaluation
  interface
     module subroutine euler_rhs_factory(object)
       class(euler_rhs_t), allocatable, intent(inout) :: object
     end subroutine euler_rhs_factory
  end interface

  public :: euler_rhs_factory

end module euler_residual
