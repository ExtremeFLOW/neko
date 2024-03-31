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
module advection_fctry
  use num_types, only : rp
  use advection, only : advection_t, adv_no_dealias_t, adv_dealias_t
  use coefs, only : coef_t
  use device, only : device_free
  use, intrinsic :: iso_c_binding, only : c_associated
  use json_utils, only : json_get
  use json_module, only : json_file
  implicit none
  private

  public :: advection_factory

contains

  !> A factory for \ref advection_t decendants.
  !! @param this Polymorphic object of class \ref advection_t.
  !! @param json The parameter file.
  !! @param coef The coefficients of the (space, mesh) pair.
  !! @note The factory both allocates and initializes `this`.
  subroutine advection_factory(this, json, coef)
    implicit none
    class(advection_t), allocatable, intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), target :: coef
    logical :: dealias, found
    integer :: lxd, order

    call json_get(json, 'case.numerics.dealias', dealias)
    call json%get('case.numerics.dealiased_polynomial_order', lxd, found)
    if (.not. found) then
       call json_get(json, 'case.numerics.polynomial_order', order)
       ! Note, assumes odd polynomial order
       lxd =  3.0_rp / 2.0_rp * (order + 1)
    end if

    ! Free allocatables if necessary
    if (allocated(this)) then
       call this%free
       deallocate(this)
    end if

    if (dealias) then
       allocate(adv_dealias_t::this)
    else
       allocate(adv_no_dealias_t::this)
    end if

    select type(adv => this)
    type is(adv_dealias_t)
       if (lxd .gt. 0) then
          call adv%init(lxd, coef)
       else
          call adv%init(coef%Xh%lx * 3/2,  coef)
       end if
    type is(adv_no_dealias_t)
       call adv%init(coef)
    end select

  end subroutine advection_factory


end module advection_fctry
