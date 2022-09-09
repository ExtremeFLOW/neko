module json_case
  use neko
  use json_module
  implicit none

contains
  
  subroutine json_case_create_neko_case(neko_case, json_case)
    type(case_t), intent(inout) :: neko_case
    type(json_file) :: json_case
    character(len=:), allocatable :: mesh_file
    character(len=:), allocatable :: fluid_scheme
    character(len=:), allocatable :: source_term
    character(len=:), allocatable :: initial_condition
    integer :: lx
    type(param_io_t) :: params
    logical found_key
    type(file_t) :: msh_file

    call neko_log%section('Case')
    call neko_log%message('Reading case from pyNeko json file')

    ! Read case description

    call json_case%get('case.mesh_file', mesh_file, found_key)
    if (.not. found_key) then
       call neko_error('Key mesh_file not found')
    else
       msh_file = file_t(trim(mesh_file))
       call msh_file%read(neko_case%msh)
    end if

    call json_case%get('case.lx', lx, found_key)
    if (.not. found_key) then
       call neko_error('lx not defined')
    end if

    call json_case%get('case.fluid_scheme', fluid_scheme, found_key)
    if (.not. found_key) then
       call neko_error('Fluid scheme not defined')
    else       
       call fluid_scheme_factory(neko_case%fluid, trim(fluid_scheme))
       call neko_case%fluid%init(neko_case%msh, lx, neko_case%params)
    end if

    call json_case%get('case.source_term', source_term, found_key)
    if (.not. found_key) then
       call neko_error('Source term not defined')
    else
       call neko_case%fluid%set_source(trim(source_term))
    end if

    call json_case%get('case.initial_condition', initial_condition, found_key)
    if (.not. found_key) then
       call neko_error('Intiail condition not defined')
    else       
       call set_flow_ic(neko_case%fluid%u, neko_case%fluid%v, neko_case%fluid%w,&
       neko_case%fluid%p, neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh,&
       initial_condition, neko_case%params)
    end if

    
    call neko_log%end_section()

  end subroutine json_case_create_neko_case
  
end module json_case
