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
!> Implements the device kernel for the `brinkman_source_term_t` type.
module brinkman_source_term_device
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use device_math, only : device_subcol3
  use field_registry, only : neko_field_registry
  implicit none
  private

  public :: brinkman_source_term_compute_device

contains

  !> Computes the Brinkman source term on the device.
  !! @param fields The right-hand side.
  !! @param values The values of the source components.
  subroutine brinkman_source_term_compute_device(fields, brinkman)
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(in) :: brinkman
    integer :: n
    type(field_t), pointer :: u, v, w

    n = fields%items(1)%ptr%dof%size()

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    call device_subcol3(fields%x_d(1), u%x_d, brinkman%x_d, n)
    call device_subcol3(fields%x_d(2), v%x_d, brinkman%x_d, n)
    call device_subcol3(fields%x_d(3), w%x_d, brinkman%x_d, n)

  end subroutine brinkman_source_term_compute_device

end module brinkman_source_term_device
