! Copyright (c) 2022, The Neko Authors
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
module pnpn_res_fctry
  use neko_config
  use utils
  use pnpn_residual
  use pnpn_res_device, only : pnpn_prs_res_device_t, pnpn_vel_res_device_t
  use pnpn_res_cpu, only : pnpn_prs_res_cpu_t, pnpn_vel_res_cpu_t, &
                           pnpn_scalar_res_cpu_t
  use pnpn_res_sx, only : pnpn_prs_res_sx_t, pnpn_vel_res_sx_t
  implicit none

contains

  subroutine pnpn_prs_res_factory(prs_res)
    class(pnpn_prs_res_t), allocatable, intent(inout) :: prs_res

    if (allocated(prs_res)) then
       deallocate(prs_res)
    end if

    
    if (NEKO_BCKND_SX .eq. 1) then
       allocate(pnpn_prs_res_sx_t::prs_res)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(pnpn_prs_res_device_t::prs_res)
    else
       allocate(pnpn_prs_res_cpu_t::prs_res)
    end if
    
  end subroutine pnpn_prs_res_factory
  
  subroutine pnpn_vel_res_factory(vel_res)
    class(pnpn_vel_res_t), allocatable, intent(inout) :: vel_res

    if (allocated(vel_res)) then
       deallocate(vel_res)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(pnpn_vel_res_sx_t::vel_res)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(pnpn_vel_res_device_t::vel_res)
    else
       allocate(pnpn_vel_res_cpu_t::vel_res)
    end if
       
    
  end subroutine pnpn_vel_res_factory
  
  subroutine pnpn_scalar_res_factory(scalar_res)
    class(pnpn_scalar_res_t), allocatable, intent(inout) :: scalar_res

    if (allocated(scalar_res)) then
       deallocate(scalar_res)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       call neko_error("Not implemented")
!       allocate(pnpn_scalar_res_sx_t::scalar_res)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call neko_error("Not implemented")
!       allocate(pnpn_scalar_res_device_t::scalar_res)
    else
       allocate(pnpn_scalar_res_cpu_t::scalar_res)
    end if
       
    
  end subroutine pnpn_scalar_res_factory
end module pnpn_res_fctry
