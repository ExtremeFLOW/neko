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
!
!> Defines a factory subroutine for boundary conditions.
module bc_fctry
  use bc, only : bc_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use logger, only : neko_log
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use coefs, only : coef_t
  use facet_zone, only : facet_zone_t
  use mesh, only : NEKO_MSH_MAX_ZLBLS
  implicit none
  private

  public :: bc_factory

contains

  !> Boundary condition factory. Both constructs and initializes the object.
  !! Will mark a mesh zone for the bc and finalize.
  !! @param[in] coef SEM coefficients.
  !! @param[inout] json JSON object for initializing the bc.
  subroutine bc_factory(object, json, coef, zones)
    !class(bc_t), allocatable, intent(inout) :: object
    class(bc_t), pointer, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(facet_zone_t), intent(inout) :: zones(NEKO_MSH_MAX_ZLBLS)
    character(len=:), allocatable :: type
    integer :: zone_index

    call json_get(json, "zone_index", zone_index)
    call json_get(json, "type", type)

    if (trim(type) .eq. "dirichlet") then
       allocate(dirichlet_t::object)
    else if (trim(type) .eq. "neumann") then
       allocate(neumann_t::object)
    end if

    call object%init(coef, json)
    call object%mark_zone(zones(zone_index))
    call object%finalize()

  end subroutine bc_factory

end module bc_fctry
