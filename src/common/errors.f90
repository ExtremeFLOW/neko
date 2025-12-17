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
!> Errors that can be thrown by Neko.
!! @details Various error handling routines that can be used in the code.
!! They use the throw_error and throw_warning mechanisms defined in utils.f90
!! in order to allow pFUnit to catch them during testing.
submodule (utils) errors
  implicit none

contains

  !> Reports a type selection error and stops execution.
  !! @param base_type The base type for which the selection was attempted.
  !! @param wrong_type The type that was attempted to be selected.
  !! @param known_types An array of known valid types.
  module subroutine neko_type_error(base_type, wrong_type, known_types)
    character(len=*), intent(in) :: base_type
    character(len=*), intent(in) :: wrong_type
    character(len=*), intent(in) :: known_types(:)
    integer :: i

    write(error_unit, *) '*** ERROR WHEN SELECTING TYPE ***'
    write(error_unit, *) 'Type ', wrong_type, ' does not exist for ', base_type
    write(error_unit, *) 'Valid types are:'
    do i = 1, size(known_types)
       write(error_unit, *) "    ", known_types(i)
    end do

    call throw_error('errors.f90', -1, message='')
  end subroutine neko_type_error

  !> Reports a type registration error and stops execution.
  !! @param base_type The base type for which the registration was attempted.
  !! @param wrong_type The type that was attempted to be registered.
  !! @param known Whether the conflicting type is known (standard) or custom.
  module subroutine neko_type_registration_error(base_type, wrong_type, known)
    character(len=*), intent(in) :: base_type
    character(len=*),intent(in) :: wrong_type
    logical, intent(in) :: known

    write(error_unit, *) '*** ERROR WHEN REGISTERING TYPE ***'
    write(error_unit, *) 'Type name ', wrong_type, &
         ' conflicts with and already existing ', base_type, " type"
    if (known) then
       write(error_unit, *) 'Please rename your custom type.'
    else
       write(error_unit, *) 'The already existing type is also custom.' // &
            ' Make all custom type names unique!'
    end if

    call throw_error('errors.f90', -1, message='')
  end subroutine neko_type_registration_error

  !> Reports a warning to standard output
  !! @param warning_msg The warning message to report.
  module subroutine neko_warning(warning_msg)
    character(len=*) :: warning_msg
    write(output_unit, *) '*** WARNING: ', warning_msg, ' ***'

    call throw_warning('errors.f90', -1, message='')
  end subroutine neko_warning

  !> Reports an error and stops execution.
  !! @param[optional] error_code The error code to report.
  module subroutine neko_error_plain(error_code)
    integer, optional, intent(in) :: error_code

    if (present(error_code)) then
       write(error_unit, *) '*** ERROR ***', error_code
    else
       write(error_unit, *) '*** ERROR ***'
    end if

    call throw_error('errors.f90', -1, message='')
  end subroutine neko_error_plain

  !> Reports an error and stops execution.
  !! @param error_msg The error message to report.
  module subroutine neko_error_msg(error_msg)
    character(len=*), intent(in) :: error_msg
    write(error_unit, *) '*** ERROR: ', error_msg, ' ***'

    call throw_error('errors.f90', -1, message='')
  end subroutine neko_error_msg


end submodule errors
