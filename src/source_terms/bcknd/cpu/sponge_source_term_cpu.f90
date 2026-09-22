! Copyright (c) 2026, The Neko Authors
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
!> Implements the cpu kernel for the `sponge_source_term_t` type.

module sponge_source_term_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  implicit none
  private

  public :: sponge_source_term_compute_cpu

contains

  !> Computes the generic sponge source term on the cpu.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param fringe The fringe field
  !! @param a_x The amplitude in the x-direction
  !! @param a_y The amplitude in the y-direction
  !! @param a_z The amplitude in the z-direction
  subroutine sponge_source_term_compute_cpu(fields, u, v, w, &
       u_bf, v_bf, w_bf, fringe, a_x, a_y, a_z)
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(in) :: u, v, w, fringe, u_bf, v_bf, w_bf
    real(kind=rp), intent(in) :: a_x, a_y, a_z
    integer :: i, n
    type(field_t), pointer :: fu, fv, fw

    n = fields%item_size(1)

    fu => fields%get(1)
    fv => fields%get(2)
    fw => fields%get(3)

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do private(i)
    do i = 1, n
       fu%x(i,1,1,1) = fu%x(i,1,1,1) + &
            a_x * fringe%x(i,1,1,1) * (u_bf%x(i,1,1,1) - u%x(i,1,1,1))
       fv%x(i,1,1,1) = fv%x(i,1,1,1) + &
            a_y * fringe%x(i,1,1,1) * (v_bf%x(i,1,1,1) - v%x(i,1,1,1))
       fw%x(i,1,1,1) = fw%x(i,1,1,1) + &
            a_z * fringe%x(i,1,1,1) * (w_bf%x(i,1,1,1) - w%x(i,1,1,1))
    end do
    !$omp end parallel do

  end subroutine sponge_source_term_compute_cpu

end module sponge_source_term_cpu
