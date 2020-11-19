!> Master module
!
module neko
  use num_types
  use comm
  use utils
  use math
  use speclib
  use dofmap
  use space
  use htable
  use generic_file
  use entity
  use point
  use element
  use quad
  use hex
  use uset
  use stack
  use tuple
  use mesh
  use mesh_field
  use map
  use nmsh
  use re2
  use rea
  use mxm_wrapper
  use mxm_std
  use rea_file
  use re2_file
  use map_file
  use vtk_file
  use fld_file
  use nmsh_file
  use file
  use field
  use mpi_types
  use gather_scatter
  use coefs
  use bc
  use wall
  use dirichlet
  use krylov
  use cg
  use precon
  use ax_product
  use gmres
  use neko_config

contains

  subroutine neko_init
    write(*,*) ''
    write(*,*) 'N E K O'
    write(*,*) '(version: ', trim(NEKO_VERSION),')'
    write(*,*) trim(NEKO_BUILD_INFO)
    write(*,*) ''
    
    call comm_init
    call mpi_types_init
  end subroutine neko_init

  subroutine neko_finalize
    call mpi_types_free
    call comm_free
  end subroutine neko_finalize

end module neko
