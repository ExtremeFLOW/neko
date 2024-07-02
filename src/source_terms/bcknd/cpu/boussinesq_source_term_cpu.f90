! Copyright (c) 2023, The Neko Authors
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
!> Implements the cpu kernel for the  `boussinesq_source_term_t` type.
module boussinesq_source_term_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use math, only : add2s2, cadd
  implicit none
  private

  public :: boussinesq_source_term_compute_cpu

contains

  !> Computes the Boussinesq source term on the cpu.
  !! @param fields The right-hand side.
  !! @param s The scalar field
  !! @param ref_value The reference value of the scalar field.
  !! @param g The gravity vector.
  !! @param beta The thermal expansion coefficient.
  subroutine boussinesq_source_term_compute_cpu(fields, s, ref_value, g, beta)
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: ref_value
    real(kind=rp), intent(in) :: g(3)
    real(kind=rp), intent(in) :: beta
    integer :: n_fields, i, n

    n_fields = fields%size()
    n = fields%item_size(1)

    do i=1, n_fields
       call add2s2(fields%items(i)%ptr%x, s%x, g(i)*beta, n)
       call cadd(fields%items(i)%ptr%x, -g(i)*beta*ref_value, n)
    end do
  end subroutine boussinesq_source_term_compute_cpu

end module boussinesq_source_term_cpu
