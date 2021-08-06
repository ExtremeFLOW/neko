!> Neko checkpoint file format
!! @details this module defines interface to read/write Neko's ceckpoint files
module chkp_file
  use generic_file
  use checkpoint    
  use num_types
  use utils
  implicit none
  private

  !> Interface for Neko checkpoint files
  type, public, extends(generic_file_t) :: chkp_file_t
   contains
     procedure :: read => chkp_file_read
     procedure :: write => chkp_file_write
  end type chkp_file_t

contains
  
  !> Write a Neko checkpoint
  subroutine chkp_file_write(this, data, t)
    class(chkp_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t

    select type(data)
    type is (chkp_t)
    class default
       call neko_error('Invalid data')
    end select

  end subroutine chkp_file_write
  
  !> Load a checkpoint from file
  subroutine chkp_file_read(this, data)
    class(chkp_file_t) :: this
    class(*), target, intent(inout) :: data
    
    call neko_error('Not implemented yet!')
    
  end subroutine chkp_file_read
  
  
end module chkp_file
