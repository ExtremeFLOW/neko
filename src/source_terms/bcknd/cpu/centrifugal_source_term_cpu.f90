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
!> Implements the cpu kernel for the `centrifugal_source_term_t` type.
!! Maintainer: Adam Peplinski.

module centrifugal_source_term_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use math, only : vcross
  use field, only : field_t
  use dofmap, only : dofmap_t
  implicit none
  private

  public :: centrifugal_source_term_compute_cpu

contains

  !> Computes the centrifugal source term on the cpu.
  !! @param omega The rotation vector.
  !! @param ref_point The reference point on the rotation axis.
  !! @param fields The right-hand side, which should be the velocity components.
  subroutine centrifugal_source_term_compute_cpu(omega, ref_point, fields)
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: ref_point(3)
    type(field_list_t), intent(inout) :: fields
    type(dofmap_t), pointer :: dof
    integer :: i, n
    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: rx, ry, rz, cx, cy, cz

    dof => fields%dof(1)
    n = fields%item_size(1)

    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)
    fw => fields%get_by_index(3)

    do concurrent (i = 1:n)
       ! displacement with respect to reference point
       rx = dof%x(i,1,1,1) - ref_point(1)
       ry = dof%y(i,1,1,1) - ref_point(2)
       rz = dof%z(i,1,1,1) - ref_point(3)

       ! Omega x r
       cx = (omega(2) * rz - omega(3) * ry)
       cy = (omega(3) * rx - omega(1) * rz)
       cz = (omega(1) * ry - omega(2) * rx)

       ! - Omega x ( Omega x r)
       fu%x(i,1,1,1) = fu%x(i,1,1,1) - (omega(2) * cz - omega(3) * cy)
       fv%x(i,1,1,1) = fv%x(i,1,1,1) - (omega(3) * cx - omega(1) * cz)
       fw%x(i,1,1,1) = fw%x(i,1,1,1) - (omega(1) * cy - omega(2) * cx)
    end do

  end subroutine centrifugal_source_term_compute_cpu

end module centrifugal_source_term_cpu
