! Copyright (c) 2025, The Neko Authors
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
submodule(regularization) regularization_fctry
  use entropy_viscosity, only: entropy_viscosity_t
  use utils, only: neko_error
  implicit none

contains

  module subroutine regularization_factory(object, type_name, json, &
                                           coef, dof, reg_coeff)
    class(regularization_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof
    type(field_t), intent(in), target :: reg_coeff

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ('entropy', 'entropy_viscosity')
       allocate(entropy_viscosity_t::object)
    case default
       call neko_error('Unknown regularization type: ' // trim(type_name))
    end select

    call object%init(json, coef, dof, reg_coeff)

  end subroutine regularization_factory

end submodule regularization_fctry

