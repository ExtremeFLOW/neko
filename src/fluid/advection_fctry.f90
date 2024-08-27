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
!> Contains the factory routine for `advection_t` children.
submodule (advection) advection_fctry
  use json_utils, only : json_get, json_get_or_default

  ! Advection and derivatives
  use adv_dealias, only : adv_dealias_t
  use adv_no_dealias, only : adv_no_dealias_t


contains

  !> A factory for \ref advection_t decendants. Both creates and initializes
  !! the object.
  !! @param object The object allocated by the factory.
  !! @param json The parameter file.
  !! @param coef The coefficients of the (space, mesh) pair.
  module subroutine advection_factory(object, json, coef)
    class(advection_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout), target :: coef
    logical :: dealias
    integer :: lxd, order

    ! Read the parameters from the json file
    call json_get(json, 'case.numerics.dealias', dealias)
    call json_get(json, 'case.numerics.polynomial_order', order)

    call json_get_or_default(json, 'case.numerics.dealiased_polynomial_order', &
                             lxd, ( 3 * (order + 1) ) / 2)

    ! Free allocatables if necessary
    if (allocated(object)) then
       call object%free
       deallocate(object)
    end if

    if (dealias) then
       allocate(adv_dealias_t::object)
    else
       allocate(adv_no_dealias_t::object)
    end if

    select type(adv => object)
      type is(adv_dealias_t)
       call adv%init(lxd, coef)
      type is(adv_no_dealias_t)
       call adv%init(coef)
    end select

  end subroutine advection_factory


end submodule advection_fctry
