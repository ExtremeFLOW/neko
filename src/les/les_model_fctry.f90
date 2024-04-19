
! Copyright (c) 2021-2022, The Neko Authors
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
module les_model_fctry
  use les_model, only : les_model_t
  use vreman, only : vreman_t
  use smagorinsky, only : smagorinsky_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use json_module, only : json_file
  implicit none
  private

  public :: les_model_factory

contains
  !> LES model factory. Both constructs and initializes the object.
  !! @param les_model The object to be allocated.
  !! @param name The name of the LES model.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine les_model_factory(les_model, name, dofmap, coef, json)
    class(les_model_t), allocatable, target, intent(inout) :: les_model
    character(len=*), intent(in) :: name
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json

    if (allocated(les_model)) then
       deallocate(les_model)
    end if

    if (trim(name) .eq. 'vreman') then
       allocate(vreman_t::les_model)
    end if

    if (trim(name) .eq. 'smagorinsky') then
      allocate(smagorinsky_t::les_model)
   end if

    call les_model%init(dofmap, coef, json)

  end subroutine les_model_factory

end module les_model_fctry