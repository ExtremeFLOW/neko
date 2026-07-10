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
!> Defines compressible residual factory.
submodule (compressible_residual) compressible_res_fctry
  use neko_config, only : NEKO_BCKND_DEVICE
  use compressible_res_cpu, only : compressible_res_cpu_t, &
       compressible_res_cpu_gamma
  use compressible_res_device, only : compressible_res_device_t, &
       compressible_res_device_gamma
  implicit none

contains

  !> Factory for the compressible residual computation routine.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  !! @param gamma Ratio of specific heats
  module subroutine compressible_rhs_factory(object, gamma)
    class(compressible_rhs_t), allocatable, intent(inout) :: object
    real(kind=rp), intent(in) :: gamma

    if (allocated(object)) then
       deallocate(object)
    end if

    !> TODO: Add support for SX
    if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(compressible_res_device_t::object)
    else
       allocate(compressible_res_cpu_t::object)
    end if

    ! Set gamma in the object
    object%gamma = gamma

    ! Sync module-level variables for CPU/GPU backends
    if (NEKO_BCKND_DEVICE .eq. 1) then
       compressible_res_device_gamma = gamma
    else
       compressible_res_cpu_gamma = gamma
    end if

  end subroutine compressible_rhs_factory

end submodule compressible_res_fctry
