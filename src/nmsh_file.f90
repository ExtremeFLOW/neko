!> Neko binary mesh data
module nmsh_file
  use generic_file
  use utils
  use nmsh
  implicit none
  
  private

  !> Interface for Neko nmsh files
  type, public, extends(generic_file_t) :: nmsh_file_t
   contains
     procedure :: read => nmsh_file_read
     procedure :: write => nmsh_file_write
  end type nmsh_file_t

contains

  !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_read(this, data)
    class(nmsh_file_t) :: this
    class(*), target, intent(inout) :: data    
  end subroutine nmsh_file_read

    !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_write(this, data)
    class(nmsh_file_t), intent(in) :: this
    class(*), target, intent(in) :: data    
  end subroutine nmsh_file_write
  
end module nmsh_file
  
