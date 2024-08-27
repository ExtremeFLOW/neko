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
  use coefs, only : coef_t
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use field_series, only: field_series_t
  use time_scheme_controller, only: time_scheme_controller_t
  use time_interpolator, only: time_interpolator_t
  use num_types, only: rp
  ! Advection and derivatives
  use advection, only : advection_t
  use adv_dealias, only : adv_dealias_t
  use adv_no_dealias, only : adv_no_dealias_t
  use adv_oifs, only : adv_oifs_t

  implicit none
  private

  public :: advection_factory

contains

  !> A factory for \ref advection_t decendants. Both creates and initializes the
  !! object.
  !! @param object The object allocated by the factory.
  !! @param json The parameter file.
  !! @param coef The coefficients of the (space, mesh) pair.
  !! @param ulag, vlag, wlag The lagged velocity fields.
  !! @param dtlag The lagged time steps.
  !! @param tlag The lagged times.
  !! @param time_scheme The bdf-ext time scheme used in the method.
  !! @param slag The lagged scalar field.
  !! @note The factory both allocates and initializes `object`.
  subroutine advection_factory(object, json, coef, ulag, vlag, wlag, &
                               dtlag, tlag, time_scheme, slag)
    implicit none
    class(advection_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(coef_t), target :: coef
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    real(kind=rp), target, intent(in) :: dtlag(10)
    real(kind=rp), target, intent(in) :: tlag(10)
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_series_t), target, optional :: slag

    logical :: dealias, oifs
    real(kind=rp) :: ctarget
    integer :: lxd, order

    ! Read the parameters from the json file
    call json_get(json, 'case.numerics.dealias', dealias)
    call json_get(json, 'case.numerics.polynomial_order', order)
    call json_get_or_default(json, 'case.numerics.oifs', oifs, .false.)

    call json_get_or_default(json, 'case.numerics.dealiased_polynomial_order', &
                             lxd, ( 3 * (order + 1) ) / 2)

    call json_get_or_default(json, 'case.numerics.target_cfl', ctarget, 1.9_rp)


    ! Free allocatables if necessary
    if (allocated(object)) then
       call object%free
       deallocate(object)
    end if

    if (oifs) then
      allocate(adv_oifs_t::object)
    else
      if (dealias) then
         allocate(adv_dealias_t::object)
      else
         allocate(adv_no_dealias_t::object)
      end if
    end if

    select type (adv => object)
      type is (adv_dealias_t)
       call adv%init(lxd, coef)
      type is (adv_no_dealias_t)
       call adv%init(coef)
      type is (adv_oifs_t)
       if (present(slag)) then
          call adv%init(lxd, coef, ctarget, ulag, vlag, wlag, &
                        dtlag, tlag, time_scheme, slag)
       else
          call adv%init(lxd, coef, ctarget, ulag, vlag, wlag, &
                        dtlag, tlag, time_scheme)
       end if
    end select

  end subroutine advection_factory


end module advection_fctry
