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
module rhs_maker_fctry
  use rhs_maker
  use rhs_maker_cpu
  use rhs_maker_sx
  use rhs_maker_device
  use neko_config
  implicit none

contains

  subroutine rhs_maker_sumab_fctry(sumab)
    class(rhs_maker_sumab_t), allocatable, intent(inout) :: sumab

    if (allocated(sumab)) then
       deallocate(sumab)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_sumab_sx_t::sumab)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_sumab_device_t::sumab)
    else
       allocate(rhs_maker_sumab_cpu_t::sumab)
    end if

  end subroutine rhs_maker_sumab_fctry

  subroutine rhs_maker_ext_fctry(makeabf)
    class(rhs_maker_ext_t), allocatable, intent(inout) :: makeabf

    if (allocated(makeabf)) then
       deallocate(makeabf)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_ext_sx_t::makeabf)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_ext_device_t::makeabf)
    else
       allocate(rhs_maker_ext_cpu_t::makeabf)
    end if

  end subroutine rhs_maker_ext_fctry

  subroutine rhs_maker_bdf_fctry(makebdf)
    class(rhs_maker_bdf_t), allocatable, intent(inout) :: makebdf

    if (allocated(makebdf)) then
       deallocate(makebdf)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_bdf_sx_t::makebdf)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_bdf_device_t::makebdf)
    else       
       allocate(rhs_maker_bdf_cpu_t::makebdf)
    end if

  end subroutine rhs_maker_bdf_fctry

end module rhs_maker_fctry
