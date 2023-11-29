module data_processor
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t

  implicit none

  !> Declare type
  type, public :: data_processor_t     
     !> Types
     type(json_file) :: params
     type(fld_io_controller_t), allocatable :: fld_ioc(:)
     type(probes_t), allocatable :: pb(:)
  end type data_processor_t
    
contains
  
  !> Initialize reading the case object
  subroutine init_params(params)
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
     
     if (pe_rank.eq.0) then
        write(*,*) "contents of json file"     
        write(*,*) json_buffer
     end if

     deallocate(json_buffer)

  end subroutine init_params

  !> Initialize the reader arrays
  subroutine init_readers(this)
     class(data_processor_t), intent(inout) :: this
     !> Internal variables
     integer :: n_readers, i
     type(json_core) :: core
     type(json_value), pointer :: reader_object, reader_pointer 
     type(json_file) :: reader_subdict
     character(len=:), allocatable :: buffer
     logical :: found

     if (this%params%valid_path('case.field_reader')) then
    
        if (pe_rank .eq. 0) then
           write(*,*)  "Initializing data readers"
        end if

        call this%params%info('case.field_reader', n_children=n_readers)

        call this%params%get_core(core)
        call this%params%get('case.field_reader', reader_object, found)

        !> Allocating reader object array
        allocate(this%fld_ioc(n_readers))

        !> Initialize all the reader objects
        do i=1, n_readers

           ! Create a new json containing just the subdict for this simcomp
           call core%get_child(reader_object, i, reader_pointer, found)
           call core%print_to_string(reader_pointer, buffer)
           call reader_subdict%load_from_string(buffer)
           
           !> Initialize the object
           call this%fld_ioc(i)%init(reader_subdict)
           !> Initialize the data processors not related to reading io
           call this%fld_ioc(i)%init_data_processors(reader_subdict)

        end do
     end if
     
  end subroutine init_readers

  !> Initialize probes_standalone
  subroutine init_probes(this)
     class(data_processor_t), intent(inout) :: this
     !> Internal variables
     integer :: n_simcomps, i
     type(json_core) :: core
     type(json_value), pointer :: simcomp_object, comp_pointer 
     type(json_file) :: comp_subdict
     character(len=:), allocatable :: buffer
     logical :: found

     if (this%params%valid_path('case.simulation_components')) then
    
        if (pe_rank .eq. 0) then
           write(*,*)  "Initializing simulation components"
        end if

        call this%params%info('case.simulation_components', n_children=n_simcomps)

        call this%params%get_core(core)
        call this%params%get('case.simulation_components', simcomp_object, found)
     
        !> Allocating reader object array
        allocate(this%pb(n_simcomps))
        
        do i=1, n_simcomps

           ! Create a new json containing just the subdict for this simcomp
           call core%get_child(simcomp_object, i, comp_pointer, found)
           call core%print_to_string(comp_pointer, buffer)
           call comp_subdict%load_from_string(buffer)
           
           !> Initialize the object
           call this%pb(i)%init_standalone(comp_subdict, this%fld_ioc(1)%dof)

        end do
     end if

  end subroutine init_probes

  subroutine get_supporting_quantities(params, fld_ioc, work_field, work_field2) 

     type(json_file), intent(inout) :: params
     type(fld_io_controller_t), intent(inout) :: fld_ioc
     type(field_t), intent(inout) :: work_field
     type(field_t), intent(inout) :: work_field2
     type(field_t), pointer :: s, ss, s_rms, eps_t, eps_k
     real(kind=rp) :: ra
     real(kind=rp) :: pr
     real(kind=rp) :: lambda
     real(kind=rp) :: mu, bar_eps_t, bar_eps_k, nu_eps_t, nu_eps_k
     !> declare variables
     integer :: i,j,k,n

     !> Read from case file
     call json_get(params, 'case.non_dimensional_quantities.Ra', ra)
     call json_get(params, 'case.non_dimensional_quantities.Pr', pr)
     
     !> Get all pointers you need from the registry to make it easier
     s => neko_field_registry%get_field("t")
     ss => neko_field_registry%get_field("tt")
     s_rms => neko_field_registry%get_field("t_rms")
     eps_t => neko_field_registry%get_field("eps_t")
     eps_k => neko_field_registry%get_field("eps_k")

     !> Compute t_rms
     !Do it also in the gpu to verufy that all is good
     if (NEKO_BCKND_DEVICE .eq. 1) then 
        !  First get mean(t)^2
        call device_col3(work_field%x_d, s%x_d,  s%x_d, s%dof%size())
        ! then substract mean(tt) - mean(t)^2
        call device_sub3(s_rms%x_d, ss%x_d,  work_field%x_d, s%dof%size())
     else
        !  First get mean(t)^2
        call col3(work_field%x, s%x,  s%x, s%dof%size())
        ! then substract mean(tt) - mean(t)^2
        call sub3(s_rms%x, ss%x,  work_field%x, s%dof%size())
     end if

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

  subroutine average_from_weights(avrg, field, work_field, weights, sum_weights)
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
     type(field_t), intent(inout) :: field
     type(field_t), intent(inout) :: work_field
     logical :: sync
     integer :: count_cpu, count_gpu
     real(kind=rp) tol

     !logical :: assert_gpu_cpu_synch
     !integer :: assert_cpu_count_entries_less_than_tol
     !integer :: assert_gpu_count_entries_less_than_tol

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

  subroutine sync_field(field)
    type(field_t), intent(inout) :: field
    integer :: n

    n = field%dof%size()
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(field%x, &
                          field%x_d, &
                          n,DEVICE_TO_HOST,sync=.true.)
    end if

  end subroutine sync_field


  subroutine compare_mesh_kolmogorov(params, fld_ioc, work_field, work_field2) 
     type(json_file), intent(inout) :: params
     type(fld_io_controller_t), intent(inout) :: fld_ioc
     type(field_t), intent(inout) :: work_field
     type(field_t), intent(inout) :: work_field2
     type(field_t), pointer :: eps_t, eps_k, dx, dy, dz
     type(file_t) :: file_obj
     type(field_list_t) :: field_list
     !
     real(kind=rp) :: ra
     real(kind=rp) :: pr
     real(kind=rp) :: lambda
     real(kind=rp) :: mu, bar_eps_t, bar_eps_k, eta_t, eta_k
     integer :: i,j,k,n
     
     !> Read from case file
     call json_get(params, 'case.non_dimensional_quantities.Ra', ra)
     call json_get(params, 'case.non_dimensional_quantities.Pr', pr)
     
     !> Get all pointers you need from the registry to make it easier
     eps_t => neko_field_registry%get_field("eps_t")
     eps_k => neko_field_registry%get_field("eps_k")

     !> Get more non dimensional parameters
     mu = sqrt(pr/ra)
     lambda = 1_rp/sqrt(ra*pr)
     
     !> ----------------------------------------------------------------
     !> get kolmogorov scale
     
     !> Integrate thermal dissipation
     bar_eps_t = 0_rp
     call average_from_weights(bar_eps_t, eps_t, &
                     work_field, fld_ioc%coef%B, fld_ioc%coef%volume)
  
     !> Integrate kinetic dissipation
     bar_eps_k = 0_rp
     call average_from_weights(bar_eps_k, eps_k, &
                     work_field, fld_ioc%coef%B, fld_ioc%coef%volume)
   
     !> Calculate global kolmogorov scale
     eta_k = 0_rp
     eta_k = ((mu**3_rp)/(bar_eps_k))**(0.25_rp)
    
     !> Calculate global batchelor scale
     eta_t = 0_rp
     eta_t = ((lambda**3_rp)/(bar_eps_k))**(0.25_rp)

     !> ----------------------------------------------------------------
     !> Calculate the derivatives of the mesh

     call neko_field_registry%add_field(fld_ioc%dof, "dx")
     call neko_field_registry%add_field(fld_ioc%dof, "dy")
     call neko_field_registry%add_field(fld_ioc%dof, "dz")
     dx => neko_field_registry%get_field("dx")
     dy => neko_field_registry%get_field("dy")
     dz => neko_field_registry%get_field("dz")
     
     !> Get the derivatives
     call dudxyz(dx%x, fld_ioc%dof%x, fld_ioc%coef%drdx, fld_ioc%coef%dsdx, &
                       fld_ioc%coef%dtdx, fld_ioc%coef)
     call dudxyz(dy%x, fld_ioc%dof%y, fld_ioc%coef%drdy, fld_ioc%coef%dsdy, &
                       fld_ioc%coef%dtdy, fld_ioc%coef)
     call dudxyz(dz%x, fld_ioc%dof%z, fld_ioc%coef%drdz, fld_ioc%coef%dsdz, &
                       fld_ioc%coef%dtdz, fld_ioc%coef)

    
     !> ----------------------------------------------------------------
     !> scale it with eta
     if (NEKO_BCKND_DEVICE .eq. 1) then 
        call device_cmult(dx%x_d, 1/eta_k, dx%dof%size())               
        call device_cmult(dy%x_d, 1/eta_k, dy%dof%size())               
        call device_cmult(dz%x_d, 1/eta_k, dz%dof%size())               
     else
        call cmult(dx%x, 1/eta_k, dx%dof%size())               
        call cmult(dy%x, 1/eta_k, dy%dof%size())               
        call cmult(dz%x, 1/eta_k, dz%dof%size())               
     end if

     !> Syncrhonize the fields with the cpu
     call sync_field(dx)
     call sync_field(dy)
     call sync_field(dz)

     !> ----------------------------------------------------------------
     !> Write the fields
     file_obj = file_t("mesh_eta.fld")
     allocate(field_list%fields(3))
     field_list%fields(1)%f => dx
     field_list%fields(2)%f => dy
     field_list%fields(3)%f => dz 
     call file_obj%write(field_list,fld_ioc%t)

  end subroutine compare_mesh_kolmogorov

end module data_processor

program post_process
  use neko
  use data_processor

  implicit none

  !> --------------------
  !> Declaration phase 
  !> --------------------
  
  !> Declare objects
  type(data_processor_t) :: pd

  !> Declare fields
  type(field_t) :: work_field
  type(field_t) :: work_field2
  
  !> declare variables
  integer :: i,j,k,n, reader, probe

  !> --------------------
  !> Initialization phase 
  !> --------------------

  !> Initialize neko 
  call neko_init 

  !> Initialize the parameters
  call init_params(pd%params)
  
  !> Initialize the file interpolator object and read first file
  call init_readers(pd)
  
  !> Initialize work fields and other misc items
  call work_field%init(pd%fld_ioc(1)%dof)
  call work_field2%init(pd%fld_ioc(1)%dof)
  call neko_field_registry%add_field(pd%fld_ioc(1)%dof, "t_rms")
  
  !> Initialize the standalone simulation components 
  !! Do this after all variables you want in registry have been added
  call init_probes(pd)

  !> --------------------
  !> Compute phase
  !> --------------------
  
  !!> --------------------
  !!> Average the quantities of all files in all readers
  !!> --------------------
  do reader = 1, size(pd%fld_ioc)
     !> Average the fields in the registry for the first data set
     call pd%fld_ioc(reader)%average_registry(pd%fld_ioc(reader)%t, sync=.true.)

     do i=1, pd%fld_ioc(reader)%number_of_files -1

        !> Read the data in the next file
        call pd%fld_ioc(reader)%step()
  
        !> Average the fields in the registry for the first data set
        call pd%fld_ioc(reader)%average_registry(pd%fld_ioc(reader)%t, sync=.true.)
        
     end do
  end do

  !> For all the readers, mode the data from the averaging object to their corresponding
  !! fields in the registry
  do reader = 1, size(pd%fld_ioc)
     !> Move the data from mean_fields to the registry so it can be interpolated
     call pd%fld_ioc(reader)%put_averages_in_registry()
  end do

  !> get quantities from the global averages that have been put into the registry
  call get_supporting_quantities(pd%params, pd%fld_ioc(1), work_field, work_field2) 
  
  !> Interpolate
  do probe = 1, size(pd%pb)
     ! here assume one time for all the probes interpolated
     call pd%pb(probe)%compute_(pd%fld_ioc(1)%t, pd%fld_ioc(1)%file_counter)
  end do
  
  !!> --------------------
  !!> Find the ratio between the mesh spacing and kolmogorov scale
  !!> --------------------

  !> get quantities from the global averages that have been put into the registry
  call compare_mesh_kolmogorov(pd%params, pd%fld_ioc(1), work_field, work_field2) 

  !> --------------------
  !> Finalization phase
  !> --------------------
  
  call work_field%free()
  call work_field2%free()

  !> Finalize the probes
  do probe = 1, size(pd%pb)
     call pd%pb(probe)%free()
  end do

  !> Finalize neko
  call neko_finalize

end program post_process
