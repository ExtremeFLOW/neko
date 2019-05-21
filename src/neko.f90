!> Master module
!
module neko
  use num_types
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
  use re2
  use rea
  use rea_file
  use re2_file
  use vtk_file
  use file
  use field
  use mpi
  use mpi_types
contains

  subroutine neko_init
    integer :: ierr

    call MPI_Init(ierr)

    call mpi_types_init

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine neko_init

  subroutine neko_finalize
    integer :: ierr

    call mpi_types_Free
    call MPI_Finalize(ierr)    
  end subroutine neko_finalize

end module neko
