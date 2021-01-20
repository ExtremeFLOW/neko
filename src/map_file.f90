!> NEKTON map file
!! @details This module is used to read/write NEKTON vertex mapping data
module map_file
  use generic_File
  use utils
  use comm
  use map
  implicit none
  private

  !> Interface for NEKTON map files
  type, public, extends(generic_file_t) :: map_file_t
   contains
     procedure :: read => map_file_read
     procedure :: write => map_file_write
  end type map_file_t

contains

  !> Load NEKTON map file
  subroutine map_file_read(this, data)
    class(map_file_t) :: this
    class(*), target, intent(inout) :: data
    type(map_t), pointer :: nm
    integer :: j, k, neli, nnzi, ierr
    
    select type(data)
    type is (map_t)
       nm => data
    class default
       call neko_error("Invalid output data")
    end select

    open(unit=10, file=trim(this%fname), status='old', iostat=ierr)
    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Reading NEKTON map file ", this%fname
    end if

    read(10, *) neli, nnzi

    !> @todo Check if neli matches map%nel
    
    do j = 1, nm%nel
       read(10, *) nm%imap(j),(nm%vertex(k, j), k=1,nm%nlv)
    end do
    
    close(unit=10)
    
  end subroutine map_file_read

  subroutine map_file_write(this, data, t)
    class(map_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=dp), intent(in), optional :: t
    call neko_error("Not implemented yet!")
  end subroutine map_file_write
  
end module map_file
