! Copyright (c) 2024, The Neko Authors
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
module ax_helm
  use ax_product, only : ax_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  use field_list
  implicit none
  private

  !> Matrix-vector product for a Helmholtz problem.
  type, public, abstract, extends(ax_t) :: ax_helm_t
   contains
     !> Compute the product for 3 fields.
     procedure, pass(this) :: compute3 => ax_helm_compute3
  end type ax_helm_t

contains
  subroutine ax_helm_compute3(this, ax1, ax2, ax3, x1, x2, x3, coef, msh, Xh)
    implicit none
    class(ax_helm_t), intent(in) :: this
    type(space_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: ax1(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: ax2(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: ax3(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: x1(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: x2(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: x3(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call this%compute(ax1, x1, coef, msh, Xh)
    call this%compute(ax2, x2, coef, msh, Xh)
    call this%compute(ax3, x3, coef, msh, Xh)
  end subroutine ax_helm_compute3

end module ax_helm
