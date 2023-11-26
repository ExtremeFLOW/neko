program interpolate_to_unstructured_mesh
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t

  implicit none

  !> Declare objects
  type(json_file) :: params
  type(fld_io_controller_t) :: fld_ioc
  type(probes_t) :: pb
  type(field_t) :: work_field
  type(field_t) :: work_field2

  !> declare variables
  integer :: i,j,k,n


  !> --------------------
  !> Initialization phase 
  !> --------------------

  !> Initialize neko 
  call neko_init 

  !> Initialize the parameters
  call init_params(params)
  
  !> Initialize the file interpolator object and read first file
  call fld_ioc%init(params)

  !> Initialize work fields and other misc items
  call work_field%init(fld_ioc%dof)
  call work_field2%init(fld_ioc%dof)
  call neko_field_registry%add_field(fld_ioc%dof, "t_rms")

  !> Initialize the data processors not related to reading io
  call fld_ioc%init_data_processors(params)

  !> Initialize the standalone simulation components
  call init_simcomps_standalone(params, pb, fld_ioc)

  !> --------------------
  !> Compute phase
  !> --------------------

  !> Average the fields in the registry for the first data set
  call fld_ioc%average_registry(fld_ioc%t, sync=.true.)

  do i=1, fld_ioc%number_of_files -1

     !> Read the data in the next file
     call fld_ioc%step()
  
     !> Average the fields in the registry for the first data set
     call fld_ioc%average_registry(fld_ioc%t, sync=.true.)
        
  end do

  !> Move the data from mean_fields to the registry so it can be interpolated
  call fld_ioc%put_averages_in_registry()

  !> get quantities from the global averages that have been put into the registry
  call get_supporting_quantities(params, fld_ioc, work_field, work_field2) 
  
  ! Interpolate
  call pb%compute_(fld_ioc%t, fld_ioc%file_counter)
  
  !> --------------------
  !> Finalization phase
  !> --------------------
  
  call work_field%free()
  call work_field2%free()

  !> Finalize the probes
  call pb%free()

  !> Finalize neko
  call neko_finalize

end program interpolate_to_unstructured_mesh

!> ----------------------------------------------------------------------------

subroutine get_supporting_quantities(params, fld_ioc, work_field, work_field2) 
  use neko
  use fld_io_controller, only : fld_io_controller_t
  use json_utils
  implicit none
  type(json_file), intent(inout) :: params
  type(fld_io_controller_t), intent(inout) :: fld_ioc
  type(field_t), intent(inout) :: work_field
  type(field_t), intent(inout) :: work_field2
  type(field_t), pointer :: s, ss, s_rms, eps_t, eps_k
  real(kind=rp) :: ra=1e8_rp
  real(kind=rp) :: pr=1_rp
  real(kind=rp) :: lambda
  real(kind=rp) :: mu, bar_eps_t, bar_eps_k, nu_eps_t, nu_eps_k
  !> declare variables
  integer :: i,j,k,n

  !> Get all pointers you need from the registry to make it easier
  s => neko_field_registry%get_field("t")
  ss => neko_field_registry%get_field("tt")
  s_rms => neko_field_registry%get_field("t_rms")
  eps_t => neko_field_registry%get_field("eps_t")
  eps_k => neko_field_registry%get_field("eps_k")

  !> Compute t_rms
  !  First get mean(t)^2
  call col3(work_field%x, s%x,  s%x, s%dof%size())
  ! then substract mean(tt) - mean(t)^2
  call sub3(s_rms%x, ss%x,  work_field%x, s%dof%size())
  !Do it also in the gpu to verufy that all is good
  !  First get mean(t)^2
  call device_col3(work_field%x_d, s%x_d,  s%x_d, s%dof%size())
  ! then substract mean(tt) - mean(t)^2
  call device_sub3(s_rms%x_d, ss%x_d,  work_field%x_d, s%dof%size())

  ! Run the field diagnostic to see what is happening
  call report_status_of_field(s_rms, work_field)
  
  !> Integrate thermal dissipation
  mu = sqrt(pr/ra)
  lambda = 1_rp/sqrt(ra*pr)
  
  bar_eps_t = 0_rp
  call average_from_weights(bar_eps_t, eps_t, &
                     work_field, fld_ioc%coef%B, fld_ioc%coef%volume)
  !! Calculate nusselt from thermal dissipaton
  nu_eps_t = 1_rp/lambda*bar_eps_t
  
  !> Integrate kinetic
  bar_eps_k = 0_rp
  call average_from_weights(bar_eps_k, eps_k, &
                     work_field, fld_ioc%coef%B, fld_ioc%coef%volume)
  
  !! Calculate nusselt from kinetic dissipaton
  nu_eps_k = 1_rp + bar_eps_k/((mu**3_rp)*ra/(pr**2_rp))

  write(*,*) " field in registry nu_eps_t and nu_eps_k are:", nu_eps_t, nu_eps_k

end subroutine get_supporting_quantities


subroutine init_params(params)
  use neko
  use json_utils
  implicit none
  
  type(json_file), intent(inout) :: params
  character(20) :: case_file = "inputs.json"
  integer :: ierr, integer_val
  character(len=:), allocatable :: json_buffer
  
  !> Read input file
    if (pe_rank .eq. 0) then
      write(*,*)  trim(case_file)
      call params%load_file(filename=trim(case_file))
      call params%print_to_string(json_buffer)
      integer_val = len(json_buffer)
    end if

    call MPI_Bcast(integer_val, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    if (pe_rank .ne. 0) allocate(character(len=integer_val)::json_buffer)
    call MPI_Bcast(json_buffer, integer_val, MPI_CHARACTER, 0, NEKO_COMM, ierr)
    call params%load_from_string(json_buffer)

    deallocate(json_buffer)

end subroutine init_params


subroutine init_simcomps_standalone(params,pb, fld_ioc)
  use neko
  use fld_io_controller, only : fld_io_controller_t
  implicit none
  type(json_file), intent(inout) :: params
  type(probes_t), intent(inout) :: pb
  type(fld_io_controller_t), intent(inout) :: fld_ioc

  integer :: n_simcomps, i
  type(json_core) :: core
  type(json_value), pointer :: simcomp_object, comp_pointer 
  type(json_file) :: comp_subdict
  character(len=:), allocatable :: buffer
  logical :: found

  if (params%valid_path('case.simulation_components')) then
    
     if (pe_rank .eq. 0) then
        write(*,*)  "Initializing simulation components"
     end if

     call params%info('case.simulation_components', n_children=n_simcomps)

     call params%get_core(core)
     call params%get('case.simulation_components', simcomp_object, found)
     do i=1, n_simcomps
        ! Create a new json containing just the subdict for this simcomp
        call core%get_child(simcomp_object, i, comp_pointer, found)
        call core%print_to_string(comp_pointer, buffer)
        call comp_subdict%load_from_string(buffer)
        call pb%init_standalone(comp_subdict, fld_ioc%dof)
     end do
  end if

end subroutine init_simcomps_standalone
  
subroutine average_from_weights(avrg, field, work_field, weights, sum_weights)
    use neko
    implicit none
    real(kind=rp), intent(inout) :: avrg
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    real(kind=rp), dimension(field%Xh%lxyz,field%msh%nelv), intent(inout) :: weights
    real(kind=rp), intent(inout) :: sum_weights
    integer :: n
    type(c_ptr) :: weights_d
    
    n = field%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then 
       weights_d = device_get_ptr(weights)
       call device_col3(work_field%x_d,field%x_d,weights_d,n)      
       avrg = device_glsum(work_field%x_d,n)                  
       avrg = avrg / abs(sum_weights)            
    else
       call col3(work_field%x,field%x,weights,n)      
       avrg = glsum(work_field%x,n)                  
       avrg = avrg / abs(sum_weights)            
    end if

end subroutine average_from_weights

function assert_gpu_cpu_synch(field, work_field, tol) result(report)
    use neko
    implicit none
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(inout) :: tol
    logical :: report
    integer :: i, n

    n = field%dof%size()

    !> Make sure the work_fields are zeroed
    call rzero(work_field%x,n)
    call device_rzero(work_field%x_d,n)

    !> Copy the contents of the gpu to the work field
    call device_copy(work_field%x_d , field%x_d, n)

    !> Copy the information from the gpu to the cpu
    call device_memcpy(work_field%x, work_field%x_d, field%dof%size(), DEVICE_TO_HOST)

    !> Now assert
    report = .true.
    do i=1, n
       if (abs(field%x(i,1,1,1)-work_field%x(i,1,1,1)) .gt. tol ) then
          report = .false.
       end if
    end do
end function assert_gpu_cpu_synch

function assert_cpu_count_entries_less_than_tol(field, work_field, tol) result(report)
    use neko
    implicit none
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(inout) :: tol
    integer :: report
    integer :: i, n

    n = field%dof%size()
    
    !> Make sure the work_fields are zeroed
    call rzero(work_field%x,n)
    call device_rzero(work_field%x_d,n)

    !> Copy the contents of the cpu to the work field
    call copy(work_field%x , field%x, n)

    !> Now assert
    report = 0
    do i=1, n
       if ( work_field%x(i,1,1,1) .lt. tol ) then
          report = report + 1
       end if
    end do

end function assert_cpu_count_entries_less_than_tol

function assert_gpu_count_entries_less_than_tol(field, work_field, tol) result(report)
    use neko
    implicit none
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    real(kind=rp), intent(inout) :: tol
    integer :: report
    integer :: i, n

    n = field%dof%size()
    
    !> Make sure the work_fields are zeroed
    call rzero(work_field%x,n)
    call device_rzero(work_field%x_d,n)

    !> Copy the contents of the gpu to the work field
    call device_copy(work_field%x_d , field%x_d, n)
    
    !> Copy the information from the gpu to the cpu
    call device_memcpy(work_field%x, work_field%x_d, n, DEVICE_TO_HOST)

    !> Now assert
    report = 0
    do i=1, n
       if ( work_field%x(i,1,1,1) .lt. tol ) then
          report = report + 1
       end if
    end do
end function assert_gpu_count_entries_less_than_tol

subroutine report_status_of_field(field, work_field)
    use neko
    implicit none
    type(field_t), intent(inout) :: field
    type(field_t), intent(inout) :: work_field
    logical :: sync
    integer :: count_cpu, count_gpu
    real(kind=rp) tol

    logical :: assert_gpu_cpu_synch
    integer :: assert_cpu_count_entries_less_than_tol
    integer :: assert_gpu_count_entries_less_than_tol

    tol = 1e-8
    sync = assert_gpu_cpu_synch(field, work_field, tol)
    tol = -1e-7
    count_cpu = assert_cpu_count_entries_less_than_tol(field, work_field, tol)
    count_gpu = assert_cpu_count_entries_less_than_tol(field, work_field, tol)

    if (sync) then
       write(*,*) "The field is syncrhonized in cpu/gpu"
    else
       write(*,*) "The field is NOT syncrhonized in cpu/gpu"
    end if
       
    write(*,*) "The number of entries with tolerance less than ", tol, "is:"
    write(*,*) "cpu: ", count_cpu
    write(*,*) "gpu: ", count_gpu
    write(*,*) "The number of entries in the field is ", field%dof%size()
    !write(*,*) "---------------------------------------------------------------"

end subroutine report_status_of_field
