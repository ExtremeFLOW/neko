! Copyright (c) 2025-2026, The Neko Authors
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
!> Implements the viscous regularization factory.
submodule(viscous_regularization) viscous_regularization_fctry
  use artificial_viscosity, only : artificial_viscosity_t
  use utils, only : neko_error
  implicit none

contains

  !> Allocate and initialize a viscous regularization method.
  !! @param object The regularization object to allocate.
  !! @param type_name The regularization type name.
  !! @param json The regularization configuration.
  !! @param coef The SEM coefficients.
  !! @param dof The SEM map of degrees of freedom.
  module subroutine viscous_regularization_factory(object, type_name, json, &
       coef, dof)
    class(viscous_regularization_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ('artificial_viscosity')
       allocate(artificial_viscosity_t::object)
    case default
       call neko_error('Unknown viscous_regularization type: ' &
            // trim(type_name))
    end select

    call object%init(json, coef, dof)

  end subroutine viscous_regularization_factory

end submodule viscous_regularization_fctry
