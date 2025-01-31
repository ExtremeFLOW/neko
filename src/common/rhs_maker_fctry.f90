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
!> Fluid right-hand-side factory for the Pn-Pn formulation
submodule (rhs_maker) rhs_maker_fctry
  use rhs_maker_cpu, only : rhs_maker_bdf_cpu_t, rhs_maker_ext_cpu_t, &
                            rhs_maker_sumab_cpu_t, rhs_maker_oifs_cpu_t
  use rhs_maker_sx, only : rhs_maker_bdf_sx_t, rhs_maker_ext_sx_t, &
                           rhs_maker_sumab_sx_t, rhs_maker_oifs_sx_t
  use rhs_maker_device, only : rhs_maker_bdf_device_t, &
       rhs_maker_ext_device_t, rhs_maker_sumab_device_t
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX

contains

  !> Factory routine for computing the extrapolated velocity values used in
  !! the pressure equation for the PnPn fluid scheme.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine rhs_maker_sumab_fctry(object)
    class(rhs_maker_sumab_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_sumab_sx_t::object)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_sumab_device_t::object)
    else
       allocate(rhs_maker_sumab_cpu_t::object)
    end if

  end subroutine rhs_maker_sumab_fctry

  !> Factory routine for computing the explicit-in-time contribution to the RHS.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine rhs_maker_ext_fctry(object)
    class(rhs_maker_ext_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_ext_sx_t::object)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_ext_device_t::object)
    else
       allocate(rhs_maker_ext_cpu_t::object)
    end if

  end subroutine rhs_maker_ext_fctry

  !> Factory routine for computing the RHS contributions from the BDF scheme.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine rhs_maker_bdf_fctry(object)
    class(rhs_maker_bdf_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_bdf_sx_t::object)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(rhs_maker_bdf_device_t::object)
    else
       allocate(rhs_maker_bdf_cpu_t::object)
    end if

  end subroutine rhs_maker_bdf_fctry

  !> Factory routine for computing the RHS contributions from the OIFS scheme.
  !! @details Only selects the compute backend.
  !! @param object The object to be allocated by the factory.
  module subroutine rhs_maker_oifs_fctry(object)
    class(rhs_maker_oifs_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       deallocate(object)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(rhs_maker_oifs_sx_t::object)
    else
       allocate(rhs_maker_oifs_cpu_t::object)
    end if

  end subroutine rhs_maker_oifs_fctry

end submodule rhs_maker_fctry
