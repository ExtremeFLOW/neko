module file
  use utils
  use generic_file
  use rea, only:rea_file_t
  use vtk_file
  implicit none
  
  type file_t
     class(*), private, pointer :: file_
     class(generic_file_t), allocatable :: file_type
   contains
     procedure :: write => file_write
     procedure :: read => file_read
     procedure :: free => file_free
  end type file_t

  interface file_t
     module procedure file_init
  end interface file_t

contains

  !> File reader/writer constructor
  !! @param fname Filename
  function file_init(fname) result(this)
    character(len=*), intent(inout) :: fname
    type(file_t), target :: this
    integer :: fname_len
    character(len=80) :: suffix
    integer suffix_pos
    class(generic_file_t), pointer :: q
    
    fname_len = len_trim(fname)
    suffix_pos = scan(trim(fname), '.', back=.true.)
    suffix = trim(fname(suffix_pos + 1:fname_len))
    
    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if
    
    if (suffix .eq. "rea") then
       allocate(rea_file_t::this%file_type)
    else if (suffix .eq.  "vtk") then
       allocate(vtk_file_t::this%file_type)
    else
       call neko_error('Unknown file format')
    end if

    q => this%file_type
    call q%init(fname)

  end function file_init

  !> File operation destructor
  subroutine file_free(this)
    class(file_t), intent(inout) :: this


  end subroutine file_free

  !> Write @a data to a file
  !! @param data Data to be written
  subroutine file_write(this, data)
    class(file_t), target :: this
    class(*), intent(inout) :: data
    class(generic_file_t), pointer :: q
    
    q => this%file_type
    call q%write(data)
    
   end subroutine file_write
   
  !> Read @a data from a file
  !! @param data Read data
  subroutine file_read(this, data)
    class(file_t), target :: this
    class(*), intent(inout) :: data
    class(generic_file_t), pointer :: q

    q => this%file_type
    call q%read(data)

  end subroutine file_read

end module file
