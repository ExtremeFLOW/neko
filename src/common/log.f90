!> Logging routines
module log
  use comm
  use num_types
  implicit none

  integer, parameter :: LOG_SIZE = 80
  
  type log_t
     integer :: indent_
     integer :: section_id_
   contains
     procedure, pass(this) :: init => log_init
     procedure, pass(this) :: begin => log_begin
     procedure, pass(this) :: end => log_end
     procedure, pass(this) :: indent => log_indent
     procedure, pass(this) :: newline => log_newline          
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
    this%indent_ = 1
    this%section_id_ = 0
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
  subroutine log_newline(this)
    class(log_t), intent(in) :: this
    integer :: i

    if (pe_rank .eq. 0) then
       write(*,*) ' '
    end if
    
  end subroutine log_newline

  !> Write a message to a log
  subroutine log_message(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

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
  
end module log
