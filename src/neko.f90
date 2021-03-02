!> Master module
!
module neko
  use num_types
  use parameters    
  use comm
  use utils
  use log
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
  use jacobi
  use neko_config
  use case
  use sampler
  use output
  use fluid_output  
  use simulation
  use ax_helm
  use operators
  use mathops
  use projection
  use user_intf
  use parmetis
  use structs
  use curve
contains

  subroutine neko_init(C)
    type(case_t), intent(inout), optional :: C
    character(len=NEKO_FNAME_LEN) :: case_file
    character(len=10) :: suffix
    integer :: argc

    call comm_init
    call mpi_types_init

    call neko_log%init()

    if (pe_rank .eq. 0) then
       write(*,*) ''
       write(*,*) '   _  __  ____  __ __  ____ '
       write(*,*) '  / |/ / / __/ / //_/ / __ \'
       write(*,*) ' /    / / _/  / ,<   / /_/ /'
       write(*,*) '/_/|_/ /___/ /_/|_|  \____/ '
       write(*,*) ''
       write(*,*) '(version: ', trim(NEKO_VERSION),')'
       write(*,*) trim(NEKO_BUILD_INFO)
       write(*,*) ''
    end if

    if (present(C)) then

       argc = command_argument_count()

       if ((argc .lt. 1) .or. (argc .gt. 1)) then
          if (pe_rank .eq. 0) write(*,*) 'Usage: ./neko <case file>'
          stop
       end if

       call get_command_argument(1, case_file)

       call filename_suffix(case_file, suffix)

       if (trim(suffix) .ne. 'case') then
          call neko_error('Invalid case file')
       end if
       
       call case_init(C, case_file)
       
    end if
    
  end subroutine neko_init

  subroutine neko_finalize(C)
    type(case_t), intent(inout), optional :: C

    if (present(C)) then
       call case_free(C)
    end if

    call mpi_types_free
    call comm_free
  end subroutine neko_finalize

end module neko
