! Copyright (c) 2022-2024, The Neko Authors
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
!> Defines Pressure residual factory for the Pn-Pn formulation
submodule (scalar_residual) scalar_residual_fctry
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_DEVICE
  use scalar_residual_device, only : scalar_residual_device_t
  use scalar_residual_cpu, only : scalar_residual_cpu_t
  use scalar_residual_sx, only : scalar_residual_sx_t
  implicit none

contains

  !> Factory for the scalar advection-diffusion residual.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine scalar_residual_factory(object)
    class(scalar_residual_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(scalar_residual_sx_t::object)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(scalar_residual_device_t::object)
    else
       allocate(scalar_residual_cpu_t::object)
    end if
    
  end subroutine scalar_residual_factory
  
end submodule scalar_residual_fctry
