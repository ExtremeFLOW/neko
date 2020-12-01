!> Defines a simulation
module case
  use num_types
  use fluid_schemes
  use parameters
  use file
  use utils
  use mesh
  implicit none

  type :: case_t
     type(mesh_t) :: msh
     type(param_t) :: params
     class(fluid_scheme_t), allocatable :: fluid
  end type case_t

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init(C, case_file)
    type(case_t), intent(inout) :: C
    character(len=*), intent(in) :: case_file

    ! Namelist for case description
    character(len=NEKO_FNAME_LEN) :: mesh_file = ''
    character(len=NEKO_FNAME_LEN) :: fluid_scheme  = ''
    character(len=80) :: solver_velocity = ''
    character(len=80) :: solver_pressure = ''
    integer :: lx = 0
    type(param_t) :: params
    namelist /NEKO_CASE/ mesh_file, fluid_scheme, lx, params, &
         solver_velocity, solver_pressure


    type(file_t) :: msh_file
   
    open(10, file=trim(case_file))
    read(10, nml=NEKO_CASE)
    close(10)
        
    msh_file = file_t(mesh_file)
    call msh_file%read(C%msh)

    C%params = params

    if (trim(fluid_scheme) .eq. 'plan1') then
       allocate(fluid_plan1_t::C%fluid)
    else if (trim(fluid_scheme) .eq. 'plan4') then
       allocate(fluid_plan4_t::C%fluid)
    else
       call neko_error('Invalid fluid scheme')
    end if

    call C%fluid%init(C%msh, lx, solver_velocity, solver_pressure)
    
  end subroutine case_init

  !> Deallocate a case 
  subroutine case_free(C)
    type(case_t), intent(inout) :: C

    if (allocated(C%fluid)) then
       call C%fluid%free()
       deallocate(C%fluid)
    end if

    call mesh_free(C%msh)
    
  end subroutine case_free
  
end module case
