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
!> Implements the cpu kernel for the `coriolis_source_term_t` type.
!! Maintainer: Timofey Mukha.

module coriolis_source_term_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use math, only : vcross
  use field, only : field_t
  implicit none
  private

  public :: coriolis_source_term_compute_cpu

contains

  !> Computes the generic Coriolis source term on the cpu.
  !! @param fields The right-hand side, which should be the velocity components.
  !! @param omega The rotation vector.
  !! @param omega The geostrophic wind.
  subroutine coriolis_source_term_compute_cpu(fields, omega, u_geo)
    type(field_list_t), intent(inout) :: fields
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: u_geo(3)
    integer :: i, n
    type(field_t), pointer :: u, v, w
    real(kind=rp) :: ui, vi, wi


    n = fields%item_size(1)

    u => fields%get_by_index(1)
    v => fields%get_by_index(2)
    w => fields%get_by_index(3)

    do concurrent (i = 1:n)
       ui = u%x(i,1,1,1) - u_geo(1)
       vi = v%x(i,1,1,1) - u_geo(2)
       wi = w%x(i,1,1,1) - u_geo(3)

       u%x(i,1,1,1) = u%x(i,1,1,1) - 2.0_rp * (omega(2) * wi - omega(3) * vi)
       v%x(i,1,1,1) = v%x(i,1,1,1) - 2.0_rp * (omega(3) * ui - omega(1) * wi)
       w%x(i,1,1,1) = w%x(i,1,1,1) - 2.0_rp * (omega(1) * vi - omega(2) * ui)
    end do

  end subroutine coriolis_source_term_compute_cpu

end module coriolis_source_term_cpu
