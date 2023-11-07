! Copyright (c) 2020-2021, The Neko Authors
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
!> Defines a Matrix-vector product
module ax_product
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use field, only : field_t
  use mesh, only : mesh_t
  implicit none
  private

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, public, abstract :: ax_t
   contains
     procedure(ax_compute), nopass, deferred :: compute
  end type ax_t

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w Vector of size @a (lx,ly,lz,nelv).
  !! @param u Vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  abstract interface
  subroutine ax_compute(w, u, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       type(space_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh       
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute
  end interface
  
end module ax_product
