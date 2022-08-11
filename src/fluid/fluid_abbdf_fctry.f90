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
!> Fluid abbdf factory for the Pn-Pn formulation
module fluid_abbdf_fctry
  use fluid_abbdf
  use fluid_abbdf_cpu
  use fluid_abbdf_sx
  use fluid_abbdf_device
  use neko_config
  implicit none

contains

  subroutine fluid_sumab_fctry(sumab)
    class(fluid_sumab_t), allocatable, intent(inout) :: sumab

    if (allocated(sumab)) then
       deallocate(sumab)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_sumab_sx_t::sumab)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &         
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_sumab_device_t::sumab)
    else
       allocate(fluid_sumab_cpu_t::sumab)
    end if
    
  end subroutine fluid_sumab_fctry

  subroutine fluid_makeabf_fctry(makeabf)
    class(fluid_makeabf_t), allocatable, intent(inout) :: makeabf

    if (allocated(makeabf)) then
       deallocate(makeabf)
    end if
    
    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_makeabf_sx_t::makeabf)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_makeabf_device_t::makeabf)
    else
       allocate(fluid_makeabf_cpu_t::makeabf)
    end if
    
  end subroutine fluid_makeabf_fctry

  subroutine fluid_makebdf_fctry(makebdf)
    class(fluid_makebdf_t), allocatable, intent(inout) :: makebdf

    if (allocated(makebdf)) then
       deallocate(makebdf)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_makebdf_sx_t::makebdf)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_makebdf_device_t::makebdf)
    else       
       allocate(fluid_makebdf_cpu_t::makebdf)
    end if
    
  end subroutine fluid_makebdf_fctry
  
end module fluid_abbdf_fctry
