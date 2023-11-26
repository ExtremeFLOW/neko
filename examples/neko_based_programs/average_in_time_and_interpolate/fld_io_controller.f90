module fld_io_controller
  use neko
  use json_utils, only : json_get
  use tensor, only : trsp
  use mean_field, only : mean_field_t
  implicit none
  private

  !> Provides tools to read fields from files
  !! @details
  !! no details
  type, public :: fld_io_controller_t
     !> mesh file name and name with the .nek5000 field counter
     character(len=:), allocatable  :: mesh_fname, field_fname
     !> type to load the actual files
     type(file_t) :: field_file, mesh_file
     
     !> object to contain all the data in the files
     type(fld_file_data_t) :: field_file_data
     !> List to point to the fields inside the file
     type(vector_ptr_t), allocatable :: fields(:)

     !> Case parameters
     type(coef_t) :: coef
     type(dofmap_t) :: dof
     type(space_t) :: Xh
     type(mesh_t) :: msh
     type(gs_t) :: gs_h

     !> Time in the currently read field
     real(kind=rp) :: t

     !> counter of the currently read file
     integer :: file_counter 

     !> number of files
     integer :: number_of_files

     !> number of fields
     integer :: number_of_fields

     !> Fields that should be added to the registry
     logical :: add_fields_to_registry = .false.
     integer :: n_fields_in_registry = 0
     character(len=20), allocatable  :: field_in_registry_name(:)
     integer, allocatable :: field_in_registry_position_in_file(:)
     type(field_list_t) :: field_in_registry_pointer

     !> Average fields in registry
     logical :: average_fields_in_registry = .false.
     real(kind=rp) :: average_starting_t
     real(kind=rp) :: average_last_t
     type(mean_field_t), allocatable :: mean_fields(:)

   contains
     
     !> Initialize object.
     procedure, pass(this) :: init => initialize
     
     !> Initialize object.
     procedure, pass(this) :: init_data_processors => initialize_data_processors
     
     !> Execute object.
     procedure, pass(this) :: step => step
    
     !> Syncronize the data read from a file into a field int he registry 
     procedure, pass(this) :: synchronize_registry => synchronize_registry
     
     !> Destructor
     procedure, pass(this) :: free => finalize

     !> Average fields that are in the registry and put the values in mean_field
     procedure, pass(this) :: average_registry => average_registry
  
     !> Move the data from mean_fields to the registry
     procedure, pass(this) :: put_averages_in_registry => put_averages_in_registry

  end type fld_io_controller_t

contains

  !> Constructor
  !! @param param parameter file
  subroutine initialize(this, params)
    class(fld_io_controller_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: i,j,k,lx, ierr

    !> Read from case file
    call json_get(params, 'case.field_reader.mesh_name', &
         this%mesh_fname)
    call json_get(params, 'case.field_reader.field_file_name', &
         this%field_fname)

    !!> Define the file names. The counters are defined in the .nek5000 file 
    this%mesh_file = file_t(trim(this%mesh_fname))
    this%field_file = file_t(trim(this%field_fname))

    !> Read the mesh
    call this%mesh_file%read(this%msh) 
    call this%field_file_data%init(this%msh%nelv,this%msh%offset_el)

    !> syncrhonize 
    !call MPI_Barrier(NEKO_COMM, ierr)

    !> Read the first field in sequence
    this%file_counter = 1
    if (pe_rank .eq. 0) write(*,*) '-------------------------------------------'
    if (pe_rank .eq. 0) write(*,*) 'Reading file: ', this%file_counter
    call this%field_file%read(this%field_file_data)
    if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', this%field_file_data%time 
    this%t = this%field_file_data%time
    
    !> determine the number of files
    this%number_of_files = this%field_file_data%meta_nsamples

    !> determine the number of fields
    this%number_of_fields = 5 + this%field_file_data%n_scalars

    !> Create dof map and adjust the mesh if deformed
    call adjust_mesh_deformation_and_create_space(this)
 
    !> Allocate list to contain the fields that are read
    allocate(this%fields(this%field_file_data%size()))
    
    !> Point the field list to the data
    call this%field_file_data%get_list(this%fields,this%field_file_data%size())

  end subroutine initialize
  
  !> Constructor 2
  !! @param param parameter file
  subroutine initialize_data_processors(this, params)
    class(fld_io_controller_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: i,j,k,lx, ierr

    !> See if some fields should be added to the registry
    if (params%valid_path('case.field_reader.add_to_registry')) then
       this%add_fields_to_registry = .true.
    
       call params%info('case.field_reader.add_to_registry.field_position_in_file', &
             n_children=this%n_fields_in_registry)
       call json_get(params, &
            'case.field_reader.add_to_registry.field_position_in_file', &
             this%field_in_registry_position_in_file)
       call json_get(params, &
            'case.field_reader.add_to_registry.field_name_in_registry', &
             this%field_in_registry_name)

       !> Add the fields to the registry
       do i = 1, this%n_fields_in_registry
          call neko_field_registry%add_field(this%dof, &
                trim(this%field_in_registry_name(i)))
       end do

       !> Have a list to pointers to the fields in the registry for easy access
       allocate(this%field_in_registry_pointer%fields(this%n_fields_in_registry))
       do i = 1, this%n_fields_in_registry
          this%field_in_registry_pointer%fields(i)%f => neko_field_registry%get_field(&
                                          trim(this%field_in_registry_name(i)))
       end do

       !> syncrhonize the file and the registry
       call this%synchronize_registry()
   
    end if    


    if (params%valid_path('case.field_reader.average_registry')) then
       this%average_fields_in_registry = .true.
     
       call json_get(params, &
            'case.field_reader.average_registry.start_time', &
             this%average_starting_t)
     
       this%average_last_t = this%average_starting_t

       !> Allocate and initialize the mean fields
       allocate(this%mean_fields(this%n_fields_in_registry))
       do i = 1, this%n_fields_in_registry
          call this%mean_fields(i)%init(neko_field_registry%get_field(&
                                   trim(this%field_in_registry_name(i))))
       end do

    end if    


  end subroutine initialize_data_processors

  subroutine synchronize_registry(this)
    class(fld_io_controller_t), intent(inout) :: this
    integer :: i,j,k,n

    n =  this%dof%size()   
    
    !> Copy the data from field_file_data object to the field object
    do i = 1, this%n_fields_in_registry

       !> Only syncrhonize the registry with the fields that are in the actual file
       if (this%field_in_registry_position_in_file(i) .gt. 0) then

          ! Copy the data from the file to our pointer
          call copy(this%field_in_registry_pointer%fields(i)%f%x, &
                 this%fields(this%field_in_registry_position_in_file(i))%v%x,n)

          ! Put it also on device
          if (NEKO_BCKND_DEVICE .eq. 1) then 
            call device_memcpy(this%field_in_registry_pointer%fields(i)%f%x, &
                            this%field_in_registry_pointer%fields(i)%f%x_d, &
                            n,HOST_TO_DEVICE)
          end if

       end if

    end do
     

  end subroutine synchronize_registry

  subroutine finalize(this)
    class(fld_io_controller_t), intent(inout) :: this
     
    if (allocated(this%field_in_registry_pointer%fields)) deallocate(this%field_in_registry_pointer%fields)
    
    call this%Xh%free()
    call this%msh%free()
    call this%gs_h%free()

  end subroutine finalize
  
  subroutine step(this)
    class(fld_io_controller_t), intent(inout) :: this
    integer :: i,j,k,e,lx 

    this%file_counter = this%file_counter + 1

    if (pe_rank .eq. 0) write(*,*) '-------------------------------------------'
    if (pe_rank .eq. 0) write(*,*) 'Reading file: ', this%file_counter
    call this%field_file%read(this%field_file_data)
    if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', this%field_file_data%time 
    this%t = this%field_file_data%time

    !> Point the field list to the data
    call this%field_file_data%get_list(this%fields,this%field_file_data%size())

    if (this%add_fields_to_registry .eqv. .true.) then
       !> syncrhonize the file and the registry
       call this%synchronize_registry()
    end if

  end subroutine step


  subroutine average_registry(this, t, sync)
    class(fld_io_controller_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    logical, intent(in) :: sync
    integer :: i,j,k,e,lx,n

    n =  this%dof%size()   
 
    !> Update the running average of the fields
    do i = 1, this%n_fields_in_registry
        call this%mean_fields(i)%update(t-this%average_last_t) 
        if (sync .eqv. .true. .and. NEKO_BCKND_DEVICE .eq. 1) then 
            call device_memcpy(this%mean_fields(i)%mf%x, &
                               this%mean_fields(i)%mf%x_d, &
                               n,DEVICE_TO_HOST)
        end if

    end do
    this%average_last_t = t

  end subroutine average_registry
  
  
  subroutine put_averages_in_registry(this)
    class(fld_io_controller_t), intent(inout) :: this
    integer :: i,j,k,n

    n =  this%dof%size()   
    
    !> Copy the data from field_file_data object to the field object
    do i = 1, this%n_fields_in_registry
       ! Copy the data from the file to our pointer
       call copy(this%field_in_registry_pointer%fields(i)%f%x, &
                 this%mean_fields(i)%mf%x,n)
       ! Put it also on device
       if (NEKO_BCKND_DEVICE .eq. 1) then 
         call device_copy(this%field_in_registry_pointer%fields(i)%f%x_d, &
                            this%mean_fields(i)%mf%x_d, &
                            n)
       end if
    end do
     

  end subroutine put_averages_in_registry


  subroutine adjust_mesh_deformation_and_create_space(this)
    class(fld_io_controller_t), intent(inout) :: this
    integer :: i,j,k,e,lx 

    !> Copy the mesh from the read field to ensure mesh deformation
    lx = this%field_file_data%lx
    !To make sure any deformation made in the user file is passed onto here as well
    do i = 1,this%msh%nelv
       this%msh%elements(i)%e%pts(1)%p%x(1) = this%field_file_data%x%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(1) = this%field_file_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(1) = this%field_file_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(1) = this%field_file_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(1) = this%field_file_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(1) = this%field_file_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(1) = this%field_file_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(1) = this%field_file_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

       this%msh%elements(i)%e%pts(1)%p%x(2) = this%field_file_data%y%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(2) = this%field_file_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(2) = this%field_file_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(2) = this%field_file_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(2) = this%field_file_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(2) = this%field_file_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(2) = this%field_file_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(2) = this%field_file_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

       this%msh%elements(i)%e%pts(1)%p%x(3) = this%field_file_data%z%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(3) = this%field_file_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(3) = this%field_file_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(3) = this%field_file_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(3) = this%field_file_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(3) = this%field_file_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(3) = this%field_file_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(3) = this%field_file_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
    end do
   
    !> Based on data read, initialize dofmap, gs, coef
    call this%Xh%init(GLL, this%field_file_data%lx, this%field_file_data%ly, this%field_file_data%lz)
    this%dof = dofmap_t(this%msh, this%Xh)

    !> Update dof map to account for non linear defformations of the elements
    !! %apply_deform is the last thing done in dofmat_init before copying to the gpu
    !! so now we will just override that with direct copy from the file
    do e = 1,this%msh%nelv
       do k = 1,lx
          do j = 1,lx
             do i = 1,lx
                this%dof%x(i,j,k,e) = this%field_file_data%x%x(linear_index(i,j,k,e,lx,lx,lx))  
                this%dof%y(i,j,k,e) = this%field_file_data%y%x(linear_index(i,j,k,e,lx,lx,lx))  
                this%dof%z(i,j,k,e) = this%field_file_data%z%x(linear_index(i,j,k,e,lx,lx,lx))   
             end do 
          end do
       end do
    end do
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%dof%x, this%dof%x_d, this%dof%size(), HOST_TO_DEVICE)
       call device_memcpy(this%dof%y, this%dof%y_d, this%dof%size(), HOST_TO_DEVICE)
       call device_memcpy(this%dof%z, this%dof%z_d, this%dof%size(), HOST_TO_DEVICE)
    end if

    write(*,*) "dof size is= ", this%dof%size()

    call this%gs_h%init(this%dof)
    call this%coef%init(this%gs_h)

  end subroutine adjust_mesh_deformation_and_create_space

end module fld_io_controller
