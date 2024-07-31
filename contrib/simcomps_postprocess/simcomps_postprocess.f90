program simcomps_postprocess
  use neko
  use dofmap, only: dofmap_t, dofmap_free
  use json_utils
  use mpi_f08
  implicit none

  ! For case file
  character(len=NEKO_FNAME_LEN) :: case_fname
  type(json_file) :: json
  character(len=:), allocatable :: json_buffer
  integer :: len_json_buf
  type(case_t) :: empty_case

  ! For fld file
  character(len=NEKO_FNAME_LEN) :: fld_fname
  type(file_t) :: fld_file
  type(fld_file_data_t) :: fld_data
  integer :: n

  ! For the simcomp
  type(json_core) :: core
  type(json_file) :: comp_subdict
  type(json_value), pointer :: simcomp_object
  character(len=:), allocatable :: comp_type
  class(simulation_component_t), allocatable :: simcomp
  type(probes_t) :: p
  logical :: is_user, found

  ! Mesh data just in case
  type(mesh_t), target :: mesh
  type(dofmap_t) :: dof
  type(space_t), target :: Xh

  ! other
  integer :: ierr, argc, i

  ! For profilin
  real(kind=dp) :: ts, te
  character(len=LOG_SIZE) :: log_buf

  call neko_init

  !
  ! Read Args: field0.fld case.case
  !
  argc = command_argument_count()

  if ((argc .le. 1) .or. (argc .gt. 2)) then
     call usage()
     stop
  end if

  ! --- Read case file name
  call get_command_argument(1, case_fname)

  ! Load json case file
  if (pe_rank .eq. 0) then
     call json%load_file(filename=trim(case_fname))
     call json%print_to_string(json_buffer)
     len_json_buf = len(json_buffer)
  end if

  call MPI_Bcast(len_json_buf, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
  if (pe_rank .ne. 0) allocate(character(len=len_json_buf)::json_buffer)
  call MPI_Bcast(json_buffer, len_json_buf, MPI_CHARACTER, 0, NEKO_COMM, ierr)
  call json%load_from_string(json_buffer)

  deallocate(json_buffer)
  ! --------------------

  ! --- Read field name
  call get_command_argument(2, fld_fname)
  fld_file = file_t(trim(fld_fname), precision = dp)
  call fld_data%init
  ts = MPI_Wtime()
  call fld_file%read(fld_data)
  te = MPI_Wtime()
  write (log_buf, *) "Total time reading: ", te - ts
  call neko_log%message(log_buf)
  n = fld_data%u%n
  ! --------------------

  ! --- Create empty objects with just what is necessary
  empty_case%end_time = 999999.0_rp
  call Xh%init(GLL, fld_data%lx, fld_data%ly, lz=fld_data%lz)
  ! --------------

  ! ----- Serch for simulation components
  call json%get_core(core)
  call json%get('case.simulation_components', simcomp_object, found)
  if (.not. found) call neko_error("Did not find any simcomps")
  call json_extract_item(core, simcomp_object, 1, comp_subdict)
  ! ---------------

  ! Allocate simulation component
  call json_get(comp_subdict, "type", comp_type)

  select case (trim(comp_type))
  case ("probes")
     allocate(probes_t::simcomp)
  case default
     call neko_warning(trim(comp_type) // " is not postprocessable! Skipping.")
  end select

  select type (s => simcomp)
  type is (probes_t)
     ts = MPI_Wtime()
     call s%init_post(comp_subdict, empty_case, Xh, fld_data)
     te = MPI_Wtime()
     write (log_buf, *) "Total time initializing probes: ", te - ts
     call neko_log%message(log_buf)
  class default
     call neko_error("Problem")
  end select
  ! -------------------

  ! --- Do a first interpolation
  ts = MPI_Wtime()
  call simcomp%compute_(fld_data%time, 1)
  te = MPI_Wtime()
  write (log_buf, *) "Total time compute: ", te - ts
  call neko_log%message(log_buf)
  ! -------------------

  ! --- Loop through all the fld files and compute probes
  do i = 2, fld_data%meta_nsamples
     ts = MPI_Wtime()
     call fld_file%read(fld_data)
     te = MPI_Wtime()
     write (log_buf, *) "Total time reading: ", te - ts
     call neko_log%message(log_buf)
     ts = MPI_Wtime()
     call simcomp%compute_(fld_data%time, i)
     te = MPI_Wtime()
     write (log_buf, *) "Total time compute: ", te - ts
     call neko_log%message(log_buf)
  end do
  ! ------------------

  call neko_log%message("Done.")
  call simcomp%free()
  call Xh%free()
  call dofmap_free(dof)
  call file_free(fld_file)
  call fld_data%free
  call neko_finalize

contains

  subroutine usage()
  call neko_log%message("simcomps_postprocess <case file> <field series>")
  call neko_log%message("This script will read from the case file and execute")
  call neko_log%message("the simulation components using input from the")
  call neko_log%message("field files provided. Currently supported simcomps:")
  call neko_log%message(" - probes")
  call neko_log%message("")
  call neko_log%message("Example:")
  call neko_log%message("---------------------------------------------------")
  call neko_log%message("   simcomps_postprocess postprocess.case field0.fld")
  end subroutine usage

end program
