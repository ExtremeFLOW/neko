!> Master module
!
module neko
  use num_types
  use comm
  use utils
  use math
  use speclib
  use space
  use htable
  use generic_file
  use entity
  use point
  use element
  use quad
  use hex
  use mesh
  use mesh_field
  use map
  use nmsh
  use re2
  use rea
  use mxm_wrapper
  use rea_file
  use re2_file
  use map_file
  use vtk_file
  use nmsh_file
  use file
  use field  
  use mpi_types
contains

  subroutine neko_init
    call comm_init
    call mpi_types_init
  end subroutine neko_init

  subroutine neko_finalize
    call mpi_types_free
    call comm_free
  end subroutine neko_finalize

end module neko
