module log
  use comm
  implicit none

  integer, parameter :: LOG_SIZE = 80
  
  type log_t
     integer :: indent_
     integer :: section_id_
   contains
     procedure, pass(this) :: init => log_init
     procedure, pass(this) :: indent => log_indent
     procedure, pass(this) :: newline => log_newline          
     procedure, pass(this) :: message => log_message
     procedure, pass(this) :: section => log_section
     procedure, pass(this) :: end_section => log_end_section
  end type log_t
  
  type(log_t), public :: neko_log
  
contains

  subroutine log_init(this)
    class(log_t), intent(inout) :: this
    this%indent_ = 1
    this%section_id_ = 0
  end subroutine log_init

  subroutine log_indent(this)
    class(log_t), intent(in) :: this
    integer :: i

    do i = 1, this%indent_
       write(*,'(A)', advance='no') ' '        
    end do
    
  end subroutine log_indent
  
  subroutine log_newline(this)
    class(log_t), intent(in) :: this
    integer :: i

    if (pe_rank .eq. 0) then
       write(*,*) ' '
    end if
    
  end subroutine log_newline

  subroutine log_message(this, msg)
    class(log_t), intent(in) :: this
    character(len=*), intent(in) :: msg

    if (pe_rank .eq. 0) then
       call this%indent()
       write(*, '(A)') trim(msg)
    end if
    
  end subroutine log_message

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

  subroutine log_end_section(this, msg)
    class(log_t), intent(inout) :: this
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
       call this%message(msg)
    end if
    
    this%section_id_ = this%section_id_ - 1
    this%indent_ = this%indent_ - this%section_id_
    
  end subroutine log_end_Section
  
  
end module log
