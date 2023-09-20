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
!> Logging routines
module logger
  use comm, only : pe_rank
  use num_types, only : rp
  implicit none
  private

  integer, public, parameter :: LOG_SIZE = 80
  
  type, public :: log_t
     integer :: indent_
     integer :: section_id_
     integer :: level_
   contains
     procedure, pass(this) :: init => log_init
     procedure, pass(this) :: begin => log_begin
     procedure, pass(this) :: end => log_end
     procedure, pass(this) :: indent => log_indent
     procedure, nopass :: newline => log_newline          
     procedure, pass(this) :: message => log_message
     procedure, pass(this) :: section => log_section
     procedure, pass(this) :: status => log_status
     procedure, pass(this) :: error => log_error
     procedure, pass(this) :: warning => log_warning
     procedure, pass(this) :: end_section => log_end_section
  end type log_t
  
  !> Global log stream
  type(log_t), public :: neko_log
  
contains

  !> Initialize a log
  subroutine log_init(this)
    class(log_t), intent(inout) :: this
    character(len=255) :: log_level
    integer :: envvar_len

    this%indent_ = 1
    this%section_id_ = 0

    call get_environment_variable("NEKO_LOG_LEVEL", log_level, envvar_len)
    if (envvar_len .gt. 0) then
       read(log_level(1:envvar_len), *) this%level_
    else
       this%level_ = 1
    end if
    
  end subroutine log_init

  !> Increase indention level
  subroutine log_begin(this)
    class(log_t), intent(inout) :: this

    if (pe_rank .eq. 0) then
       this%indent_ = this%indent_ + 1
    end if
    
  end subroutine log_begin

  !> Decrease indention level
  subroutine log_end(this)
    class(log_t), intent(inout) :: this

    if (pe_rank .eq. 0) then
       this%indent_ = this%indent_ - 1
    end if
    
  end subroutine log_end
    
  !> Indent a log 
  subroutine log_indent(this)
    class(log_t), intent(in) :: this
    integer :: i

    if (pe_rank .eq. 0) then
       do i = 1, this%indent_
          write(*,'(A)', advance='no') ' '        
       end do
    end if
    
  end subroutine log_indent

  !> Write a new line to a log
  subroutine log_newline

    if (pe_rank .eq. 0) then
       write(*,*) ' '
    end if
    
  end subroutine log_newline

  !> Write a message to a log
  subroutine log_message(this, msg, lvl)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg
    integer, optional, intent(in) :: lvl

    if (present(lvl)) then
       if (lvl .gt. this%level_) then
          return
       end if
    end if

    if (pe_rank .eq. 0) then
       call this%indent()
       write(*, '(A)') trim(msg)
    end if
    
  end subroutine log_message

  !> Write an error message to a log
  subroutine log_error(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (pe_rank .eq. 0) then
       call this%indent()
       write(*, '(A,A,A)') '*** ERROR: ', trim(msg),'  ***'       
    end if

  end subroutine log_error

  !> Write a warning message to a log
  subroutine log_warning(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (pe_rank .eq. 0) then
       call this%indent()
       write(*, '(A,A,A)') '*** WARNING: ', trim(msg),'  ***'       
    end if

  end subroutine log_warning

  !> Begin a new log section
  subroutine log_section(this, msg)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in) :: msg
    integer :: i, pre

    if (pe_rank .eq. 0) then

       this%indent_ = this%indent_ + this%section_id_
       this%section_id_ = this%section_id_ + 1

       pre = (30 - len_trim(msg)) / 2
       
       write(*,*) ' '
       call this%indent()
       do i = 1, pre
          write(*,'(A)', advance='no') '-'
       end do
       
       write(*,'(A)', advance='no') trim(msg)
       do i = 1, 30 - (len_trim(msg) + pre)
          write(*,'(A)', advance='no') '-'
       end do
       write(*,*) ' '
    end if
    
  end subroutine log_section

  !> End a log section
  subroutine log_end_section(this, msg)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
       call this%message(msg)
    end if
    
    if (pe_rank .eq. 0) then
       this%section_id_ = this%section_id_ - 1
       this%indent_ = this%indent_ - this%section_id_
    end if
    
  end subroutine log_end_Section
  
  !> Write status banner
  !! @todo move to a future Time module
  subroutine log_status(this, t, T_end)
    class(log_t), intent(in) :: this
    real(kind=rp), intent(in) :: t
    real(kind=rp), intent(in) :: T_end
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_prog

     t_prog = 100d0 * t / T_end

    call this%message('----------------------------------------------------------------')
    write(log_buf, '(A,E15.7,A,F6.2,A)') 't = ', t,&
         '                                  [ ',t_prog,'% ]'

    call this%message(log_buf)
    call this%message('----------------------------------------------------------------')
  end subroutine log_status
  
  !
  ! Rudimentary C interface
  !

  !> Write a message to a log (from C)
  !! @note This assumes the global log stream @a neko_log  
  subroutine log_message_c(c_msg) bind(c, name='log_message') 
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len, i
  
    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do
       
       call neko_log%indent()
       write(*, '(A)') trim(msg(1:len))
    end if
    
  end subroutine log_message_c

  !> Write an error message to a log (from C)
  !! @note This assumes the global log stream @a neko_log  
  subroutine log_error_c(c_msg) bind(c, name="log_error")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len, i
    
    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do
       
       call neko_log%indent()
       write(*, '(A,A,A)') '*** ERROR: ',trim(msg(1:len)),'  ***'       
    end if
    
  end subroutine log_error_c

  !> Write a warning message to a log (from C)
  !! @note This assumes the global log stream @a neko_log  
  subroutine log_warning_c(c_msg) bind(c, name="log_warning")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len, i
    
    if (pe_rank .eq. 0) then
       len = 0
       do
          if (c_msg(len+1) .eq. C_NULL_CHAR) exit
          len = len + 1
          msg(len:len) = c_msg(len)
       end do
       
       call neko_log%indent()
       write(*, '(A,A,A)') '*** WARNING: ',trim(msg(1:len)),'  ***'       
    end if
    
  end subroutine log_warning_c

  !> Begin a new log section (from C)
  !! @note This assumes the global log stream @a neko_log  
  subroutine log_section_c(c_msg) bind(c, name="log_section")
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(in) :: c_msg
    character(len=LOG_SIZE) :: msg
    integer :: len, i
    
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
  subroutine log_end_section_c() bind(c, name="log_end_section")
  
    call neko_log%end_section()
    
  end subroutine log_end_section_c
  
end module logger
