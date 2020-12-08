module file
  use utils
  use generic_file
  use nmsh_file
  use map_file
  use rea_file
  use re2_file
  use fld_file
  use vtk_file
  implicit none
  
  type file_t
     class(generic_file_t), allocatable :: file_type
   contains
     procedure :: write => file_write
     procedure :: read => file_read
     final :: file_free
  end type file_t

  interface file_t
     module procedure file_init
  end interface file_t

contains

  !> File reader/writer constructor
  !! @param fname Filename
  function file_init(fname) result(this)
    character(len=*) :: fname
    type(file_t), target :: this
    character(len=80) :: suffix
    class(generic_file_t), pointer :: q
    
    call filename_suffix(fname, suffix)
    
    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if
    
    if (suffix .eq. "rea") then
       allocate(rea_file_t::this%file_type)
    else if (suffix .eq. "re2") then
       allocate(re2_file_t::this%file_type)
    else if (suffix .eq. "map") then
       allocate(map_file_t::this%file_type)
    else if (suffix .eq. "vtk") then
       allocate(vtk_file_t::this%file_type)
    else if (suffix .eq. "nmsh") then
       allocate(nmsh_file_t::this%file_type)
    else if (suffix .eq. "fld") then
       allocate(fld_file_t::this%file_type)
    else
       call neko_error('Unknown file format')
    end if

    call this%file_type%init(fname)

  end function file_init

  !> File operation destructor
  subroutine file_free(this)
    type(file_t), intent(inout) :: this

    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if

  end subroutine file_free

  !> Write @a data to a file
  !! @param data Data to be written
  subroutine file_write(this, data)
    class(file_t), intent(inout) :: this
    class(*), intent(inout) :: data

    call this%file_type%write(data)
    
  end subroutine file_write
   
  !> Read @a data from a file
  !! @param data Read data
  subroutine file_read(this, data)
    class(file_t), intent(in) :: this
    class(*), intent(inout) :: data

    call this%file_type%read(data)
    
  end subroutine file_read

end module file
