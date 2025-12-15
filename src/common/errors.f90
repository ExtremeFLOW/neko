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
!! Also contans additional hooks for pFUnit to catch errors during testing.
module errors
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  implicit none
  private

  abstract interface
     !> Interface for the throw procedure. Follows pFunit conventions.
     subroutine throw_intf(filename, line_number, message)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: line_number
        character(len=*), optional, intent(in) :: message
     end subroutine throw_intf
  end interface

  !> Pointer to the throw procedure.
  !! @details Ordinarily only raises an error stop. During testing,
  !! pFUnit will hijack this procedure to raise an excpeption that
  !! can be caught by the testing framework.
  procedure(throw_intf), public, pointer :: throw_error => &
       default_throw_error
  !> Same as above, but for warnings. Does nothing by default
  procedure(throw_intf), public, pointer :: throw_warning => &
       default_throw_warning

  public :: throw_intf, neko_error, neko_type_error, &
       neko_type_registration_error, neko_warning

contains
  !> Default throw method that stops execution.
  !! @details pFUnit will highjack this method during testing to catch
  !! exceptions.
  subroutine default_throw_error(filename, line_number, message)
     character(len=*), intent(in) :: filename
     integer, intent(in) :: line_number
     character(len=*), optional, intent(in) :: message

     error stop
  end subroutine default_throw_error

  !> Default throw method for warnings. Does nothing.
  subroutine default_throw_warning(filename, line_number, message)
     character(len=*), intent(in) :: filename
     integer, intent(in) :: line_number
     character(len=*), optional, intent(in) :: message
  end subroutine default_throw_warning

  !> Reports a type selection error and stops execution.
  !! @param base_type The base type for which the selection was attempted.
  !! @param wrong_type The type that was attempted to be selected.
  !! @param known_types An array of known valid types.
  subroutine neko_type_error(base_type, wrong_type, known_types)
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

  !> Reports an error and stops execution.
  !! @param error_msg Optional error message to report.
  !! @param error_code Optional error code to report.
  subroutine neko_error(error_msg, error_code)
    character(len=*), optional :: error_msg
    integer, optional :: error_code
    character(len=:), allocatable :: msg

    if (present(error_code)) then
       write(error_unit, *) '*** ERROR ***', error_code
       error stop
    else if (present(error_msg)) then
       write(error_unit, *) '*** ERROR: ', error_msg, ' ***'
    else
       write(error_unit, *) '*** ERROR ***'
       error stop
    end if

    call throw_error('errors.f90', -1, message='')
  end subroutine neko_error

  !> Reports a type registration error and stops execution.
  !! @param base_type The base type for which the registration was attempted.
  !! @param wrong_type The type that was attempted to be registered.
  !! @param known Whether the conflicting type is known (standard) or custom.
  subroutine neko_type_registration_error(base_type, wrong_type, known)
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
  subroutine neko_warning(warning_msg)
    character(len=*) :: warning_msg
    write(output_unit, *) '*** WARNING: ', warning_msg, ' ***'

    call throw_warning('errors.f90', -1, message='')
  end subroutine neko_warning

end module errors
