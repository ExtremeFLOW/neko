! Copyright (c) 2021-2024, The Neko Authors
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
!> Logging routines
module logger
  use neko_config, only : NEKO_VERSION
  use comm, only : pe_rank
  use utils, only : neko_error, neko_warning
  use, intrinsic :: iso_fortran_env, only : stdout => output_unit, &
       stderr => error_unit
  implicit none
  private

  ! > Size of the log message buffer
  !! @note This adjust for the leading space applied by `write`. 80 character
  !! output log leaves 79 characters for the message.
  integer, public, parameter :: LOG_SIZE = 79

  !> Length of the section header
  integer, public, parameter :: SEC_HEAD_SIZE = 30

  type, public :: log_t
     integer, private :: indent_
     integer, private :: section_id_
     integer, private :: tab_size_
     integer, private :: level_
     integer, private :: unit_

     character(len=LOG_SIZE), private :: section_header = ""

   contains
     procedure, pass(this) :: init => log_init
     procedure, pass(this) :: free => log_free
     procedure, pass(this) :: begin => log_begin
     procedure, pass(this) :: end => log_end
     procedure, pass(this) :: indent => log_indent
     procedure, pass(this) :: newline => log_newline
     procedure, pass(this) :: message => log_message
     procedure, pass(this) :: section => log_section
     procedure, pass(this) :: header => log_header
     procedure, pass(this) :: error => log_error
     procedure, pass(this) :: warning => log_warning
     procedure, pass(this) :: deprecated => log_deprecated
     procedure, pass(this) :: end_section => log_end_section

     procedure, private, pass(this) :: print_section_header => &
          log_print_section_header
  end type log_t

  !> Global log stream
  type(log_t), public :: neko_log
  !> Always logged
  integer, public, parameter :: NEKO_LOG_QUIET = 0
  !> Default log level
  integer, public, parameter :: NEKO_LOG_INFO = 1
  !> Verbose log level
  integer, public, parameter :: NEKO_LOG_VERBOSE = 2
  !> Deprecation log level
  integer, public, parameter :: NEKO_LOG_DEPRECATION = 5
  !> Debug log level
  integer, public, parameter :: NEKO_LOG_DEBUG = 10

  !> List of already logged deprecated features
  character(len=50), dimension(:), allocatable :: deprecated_list

contains

  !> Initialize a log
  subroutine log_init(this)
    class(log_t), intent(inout) :: this
    character(len=255) :: log_level
    character(len=255) :: log_tab_size
    character(len=255) :: log_file
    integer :: envvar_len

    this%indent_ = 0
    this%section_id_ = 0

    call get_environment_variable("NEKO_LOG_TAB_SIZE", log_tab_size, envvar_len)
    if (envvar_len .gt. 0) then
       read(log_tab_size(1:envvar_len), *) this%tab_size_
    else
       this%tab_size_ = 1
    end if

    call get_environment_variable("NEKO_LOG_LEVEL", log_level, envvar_len)
    if (envvar_len .gt. 0) then
       read(log_level(1:envvar_len), *) this%level_
    else
       this%level_ = NEKO_LOG_INFO
    end if

    call get_environment_variable("NEKO_LOG_FILE", log_file, envvar_len)
    if (envvar_len .gt. 0) then
       open(newunit = this%unit_, file = trim(log_file), status = 'replace', &
            action = 'write')
    else
       this%unit_ = stdout
    end if

  end subroutine log_init

  !> Free a log
  subroutine log_free(this)
    class(log_t), intent(inout) :: this
    integer :: i

    if (this%section_id_ .ne. 0) then
       call neko_error("Log is unbalanced")
    end if

    if (allocated(deprecated_list)) then
       call this%section("Deprecated features summary", NEKO_LOG_DEPRECATION)

       do i = 1, size(deprecated_list)
          call this%message(trim(deprecated_list(i)), NEKO_LOG_DEPRECATION)
       end do
       call this%end_section()
    end if

    if (this%unit_ .ne. stdout) then
       close(this%unit_)
    end if

    this%indent_ = 0
    this%level_ = NEKO_LOG_INFO
    this%unit_ = -1

    if (allocated(deprecated_list)) then
       deallocate(deprecated_list)
    end if

  end subroutine log_free

  !> Increase indention level
  subroutine log_begin(this)
    class(log_t), intent(inout) :: this

    if (pe_rank .eq. 0) then
       this%section_id_ = this%section_id_ + 1
       this%indent_ = this%indent_ + this%tab_size_
    end if

  end subroutine log_begin

  !> Decrease indention level
  subroutine log_end(this)
    class(log_t), intent(inout) :: this

    if (pe_rank .eq. 0) then
       if (this%section_id_ .eq. 0) then
          call neko_error("Log is unbalanced")
       end if
       this%section_id_ = this%section_id_ - 1
       this%indent_ = this%indent_ - this%tab_size_
    end if

    this%section_header = ""

  end subroutine log_end

  !> Indent a log
  subroutine log_indent(this)
    class(log_t), intent(in) :: this

    if (pe_rank .eq. 0) then
       write(this%unit_, '(A)', advance = 'no') repeat(' ', this%indent_)
    end if

  end subroutine log_indent

  !> Write a new line to a log
  subroutine log_newline(this, lvl)
    class(log_t), intent(in) :: this
    integer, optional :: lvl

    integer :: lvl_

    if (present(lvl)) then
       lvl_ = lvl
    else
       lvl_ = NEKO_LOG_INFO
    end if

    if (lvl_ .gt. this%level_) then
       return
    end if

    if (pe_rank .eq. 0) then
       write(this%unit_, '(A)') ''
    end if

  end subroutine log_newline

  !> Write a message to a log
  subroutine log_message(this, msg, lvl)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in) :: msg
    integer, optional :: lvl
    integer :: lvl_

    if (present(lvl)) then
       lvl_ = lvl
    else
       lvl_ = NEKO_LOG_INFO
    end if

    if (lvl_ .gt. this%level_) then
       return
    end if

    if (len_trim(this%section_header) .gt. 0) then
       call this%print_section_header(lvl)
    end if

    if (pe_rank .eq. 0) then
       call this%indent()
       write(this%unit_, '(A)') trim(msg)
    end if

  end subroutine log_message

  !> Write the Neko header to a log
  subroutine log_header(this, version, build_info)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: version
    character(len=*), intent(in) :: build_info

    if (pe_rank .eq. 0) then
       write(this%unit_, '(A)') ''
       write(this%unit_, '(1X,A)') '   _  __  ____  __ __  ____  '
       write(this%unit_, '(1X,A)') '  / |/ / / __/ / //_/ / __ \ '
       write(this%unit_, '(1X,A)') ' /    / / _/  / ,<   / /_/ / '
       write(this%unit_, '(1X,A)') '/_/|_/ /___/ /_/|_|  \____/  '
       write(this%unit_, '(A)') ''
       write(this%unit_, '(1X,A,A,A)') '(version: ', trim(version), ')'
       write(this%unit_, '(1X,A)') trim(build_info)
       write(this%unit_, '(A)') ''
    end if

  end subroutine log_header

  !> Write an error message to a log
  subroutine log_error(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (pe_rank .eq. 0) then
       call this%indent()
       write(stderr, '(A,A,A)') '*** ERROR: ', trim(msg), '  ***'
    end if

  end subroutine log_error

  !> Write a warning message to a log
  subroutine log_warning(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (pe_rank .eq. 0) then
       call this%indent()
       write(this%unit_, '(A,A,A)') '*** WARNING: ', trim(msg), '  ***'
    end if

  end subroutine log_warning

  !> Write a deprecation warning to a log
  !! @param feature Name of the deprecated feature
  !! @param removal_version Optional version when the feature will be removed
  !! @param extra_info Optional additional message to print
  subroutine log_deprecated(this, feature, removal_version, extra_info)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in) :: feature
    character(len=*), intent(in) :: removal_version
    character(len=*), intent(in), optional :: extra_info
    character(len=50), dimension(:), allocatable :: tmp_list
    character(len=255) :: deprecation_error
    character(len=LOG_SIZE) :: msg
    integer :: i

    if (pe_rank .ne. 0) return

    if (this%level_ .ge. NEKO_LOG_DEPRECATION .or. &
         is_deprecated(removal_version)) then

       ! Check that the feature have not already been logged
       if (.not. allocated(deprecated_list)) then
          allocate(character(len=50) :: deprecated_list(1))
          deprecated_list = trim(feature)
       else
          do i = 1, size(deprecated_list)
             if (trim(deprecated_list(i)) .eq. trim(feature)) return
          end do

          ! Save the feature to the list of deprecated features
          call move_alloc(deprecated_list, tmp_list)
          allocate(character(len=50)::deprecated_list(size(tmp_list)+1))
          deprecated_list(1:size(tmp_list)) = tmp_list
          deprecated_list(size(tmp_list) + 1) = trim(feature)
          deallocate(tmp_list)
       end if

       ! Construct deprecation message
       write(msg, '(A,A)') '*** DEPRECATION: ', trim(feature)
       call this%message(msg, NEKO_LOG_DEPRECATION)
       write(msg, '(A,A,A)') 'The feature "', trim(feature), '" is deprecated.'
       call this%message(msg, NEKO_LOG_DEPRECATION)
       write(msg, '(A,A,A)') 'It will be removed in version ', &
            trim(removal_version), '.'
       call this%message(msg, NEKO_LOG_DEPRECATION)

       if (present(extra_info)) then
          call this%message(extra_info, NEKO_LOG_DEPRECATION)
       end if

       call this%message('***', NEKO_LOG_DEPRECATION)

       deprecation_error = ""
       call get_environment_variable("NEKO_DEPRECATION_ERROR", &
            deprecation_error)

       if (trim(deprecation_error) .eq. "1") then
          call neko_error('Deprecated feature "' // trim(feature) // &
               '" is scheduled for removal in version: ' // &
               trim(removal_version) // ' (current version: ' // &
               trim(NEKO_VERSION) // ').')
       else if (is_deprecated(removal_version)) then
          call neko_warning('Deprecated feature "' // trim(feature) // &
               '" is scheduled for removal in version: ' // &
               trim(removal_version) // ' (current version: ' // &
               trim(NEKO_VERSION) // ').')
       end if
    end if

  end subroutine log_deprecated

  !> Begin a new log section
  subroutine log_section(this, msg, lvl)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in) :: msg
    integer, optional :: lvl

    integer :: pre, pos

    if (len_trim(this%section_header) .gt. 0) then
       call this%print_section_header(lvl)
    end if

    call this%begin()

    if (pe_rank .eq. 0) then
       pre = (30 - len_trim(msg)) / 2
       pos = 30 - (len_trim(msg) + pre)

       if (pre .lt. 0 .or. pos .lt. 0) then
          pre = 1
          pos = 1
          write(this%section_header, '(A,A,A)') &
               repeat('-', pre), trim(msg(1: SEC_HEAD_SIZE - 2)), &
               repeat('-', pos)
       else
          write(this%section_header, '(A,A,A)') &
               repeat('-', pre), trim(msg), repeat('-', pos)
       end if
    end if

  end subroutine log_section

  !> Print a section header
  subroutine log_print_section_header(this, lvl)
    class(log_t), intent(inout) :: this
    integer, optional :: lvl
    integer :: lvl_

    if (present(lvl)) then
       lvl_ = lvl
    else
       lvl_ = NEKO_LOG_INFO
    end if

    if (lvl_ .gt. this%level_) then
       return
    end if

    if (pe_rank .eq. 0) then
       call this%newline(lvl)
       call this%indent()
       write(this%unit_, '(A)') trim(this%section_header)
       this%section_header = ""
    end if

  end subroutine log_print_section_header

  !> End a log section
  subroutine log_end_section(this, msg, lvl)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in), optional :: msg
    integer, optional :: lvl

    if (present(msg)) then
       call this%message(msg, lvl)
    end if

    call this%end()

  end subroutine log_end_section

  !
  ! Rudimentary C interface
  !

  !> Write a message to a log (from C)
  !! @note This assumes the global log stream @a neko_log
  subroutine log_message_c(c_msg) bind(c, name = 'log_message')
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len

    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do

       call neko_log%message(trim(msg(1:len)))
    end if

  end subroutine log_message_c

  !> Write an error message to a log (from C)
  !! @note This assumes the global log stream @a neko_log
  subroutine log_error_c(c_msg) bind(c, name = "log_error")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len

    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do

       call neko_log%indent()
       write(stderr, '(A,A,A)') '*** ERROR: ', trim(msg(1:len)), '  ***'
    end if

  end subroutine log_error_c

  !> Write a warning message to a log (from C)
  !! @note This assumes the global log stream @a neko_log
  subroutine log_warning_c(c_msg) bind(c, name = "log_warning")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len

    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do

       call neko_log%indent()
       write(neko_log%unit_, '(A,A,A)') &
            '*** WARNING: ', trim(msg(1:len)), '  ***'
    end if

  end subroutine log_warning_c

  !> Begin a new log section (from C)
  !! @note This assumes the global log stream @a neko_log
  subroutine log_section_c(c_msg) bind(c, name = "log_section")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len

    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do

       call neko_log%section(trim(msg(1:len)))
    end if

  end subroutine log_section_c

  !> End a log section (from C)
  !! @note This assumes the global log stream @a neko_log
  subroutine log_end_section_c() bind(c, name = "log_end_section")

    call neko_log%end_section()

  end subroutine log_end_section_c

  !> Compare version strings
  !! @param version_removal Version string when the feature will be removed
  !! @return .true. if the current version is newer than or equal to the removal
  !! version, .false. otherwise
  logical function is_deprecated(version_removal)
    character(len=*), intent(in) :: version_removal
    character(len=50) :: current_str, removal_str
    integer :: current_number(3), removal_number(3)
    integer :: i, current_size, removal_size
    integer :: iostat_current, iostat_removal
    logical :: versions_are_valid, is_newer_than_removal

    current_str = trim(NEKO_VERSION)
    removal_str = trim(version_removal)

    current_size = 1
    do i = 1, len_trim(current_str)
       if (current_str(i:i) .eq. '.') then
          current_str(i:i) = ' '
          current_size = current_size + 1
       end if
    end do

    removal_size = 1
    do i = 1, len_trim(removal_str)
       if (removal_str(i:i) .eq. '.') then
          removal_str(i:i) = ' '
          removal_size = removal_size + 1
       end if
    end do

    read(current_str, *, iostat = iostat_current) &
         current_number(1:current_size)
    read(removal_str, *, iostat = iostat_removal) &
         removal_number(1:removal_size)
    versions_are_valid = (iostat_current .eq. 0 .and. iostat_removal .eq. 0)

    if (.not. versions_are_valid) then
       call neko_error('Invalid version string in deprecation check: ' // &
            'NEKO_VERSION=' // trim(NEKO_VERSION) // ', ' // &
            'removal_version=' // trim(version_removal))
    end if

    is_deprecated = .true.
    do i = 1, current_size
       if (current_number(i) .gt. removal_number(i)) then
          is_deprecated = .true.
          exit
       else if (current_number(i) .lt. removal_number(i)) then
          is_deprecated = .false.
          exit
       end if
    end do

  end function is_deprecated
end module logger
