module field_file_interpolator
  use neko
  use json_utils, only : json_get
  use tensor, only : trsp
  implicit none
  private

  !> Provides tools to interpolate from data written to file
  !! @details
  !! no details
  type, public :: field_file_interpolator_t
     !> File names
     character(len=:), allocatable  :: mesh_fname, field_fname
     !> Files 
     type(file_t) :: field_file, mesh_file
     !> Data read from files
     type(fld_file_data_t) :: field_data
     !> Case parameters
     type(coef_t) :: coef
     type(dofmap_t) :: dof
     type(space_t) :: Xh
     type(mesh_t) :: msh
     type(gs_t) :: gs_h
     !> List to point to the fields inside the file
     type(vector_ptr_t), allocatable :: fields_in_file(:)
     real(kind=rp) :: t

     !> Fields to be probed
     type(field_list_t) :: sampled_fields
     character(len=20), allocatable  :: which_fields(:)
     !> Position in the file where relevant field is
     integer, allocatable :: field_pos(:)

     !> Probe type
     type(probes_t) :: pb
     !> Output variables
     type(file_t) :: fout
     type(matrix_t) :: mat_out
     !> Case IO parameters  
     integer            :: n_fields
     character(len=:), allocatable  :: output_file
     !> Output control
     logical :: write_output = .false.

   contains
     
     !> Initialize object.
     procedure, pass(this) :: init => initialize
     
     !> Execute object.
     procedure, pass(this) :: interpolate => interpolate_field
     
     !> Destructor
     procedure, pass(this) :: free => finalize

  end type field_file_interpolator_t

contains

  !> Constructor
  !! @param param parameter file
  subroutine initialize(this, params)
    class(field_file_interpolator_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: i,j,k,lx 

    !> Read from case file
    call params%info('case.probes.fields', n_children=this%n_fields)
    call json_get(params, 'case.probes.fields', this%which_fields) 
    call json_get(params, 'case.probes.field_pos_file', this%field_pos) 


    call json_get(params, 'case.probes.mesh_name', &
         this%mesh_fname)
    call json_get(params, 'case.probes.field_name', &
         this%field_fname)

    !!> Define the file names. The counters are defined in the .nek5000 file 
    this%mesh_file = file_t(trim(this%mesh_fname))
    this%field_file = file_t(trim(this%field_fname))

    !> Read the mesh
    call this%mesh_file%read(this%msh) 
    call this%field_data%init(this%msh%nelv,this%msh%offset_el)

    !> Read the first field in sequence
    if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
    call this%field_file%read(this%field_data)
    if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', this%field_data%time 
    this%t = this%field_data%time

    !> Create dof map and adjust the mesh if deformed
    call adjust_mesh_deformation_and_create_space(this)

    !> IMPORTANT! 
    !! Because the interpolator works with fields in the registry, create the ones with appropiate names
    !! based on the case file. Later when using the interpolator, make sure to copy the data to these fields
    do i = 1, this%n_fields
       call neko_field_registry%add_field(this%dof, trim(this%which_fields(i)))
    end do
    !! Then get a list with pointers to the fields
    allocate(this%sampled_fields%fields(this%n_fields))
    do i = 1, this%n_fields
       this%sampled_fields%fields(i)%f => neko_field_registry%get_field(&
                                          trim(this%which_fields(i)))
    end do
   
    !> Initialize the interpolator
    call initialize_interpolator(this, params)

    !> Allocate list to contain the fields that are read
    allocate(this%fields_in_file(this%field_data%size()))
    

  end subroutine initialize


  subroutine interpolate_field(this)
    class(field_file_interpolator_t), intent(inout) :: this
    integer :: i,j,k,n

    !> For the first field that is already in memory
    !> Copy the data from the field_data to the correct field and use the probe interpolator
    call use_interpolator(this)
    
    !> Read the rest of the fields in the sequence
    do i = 1, this%field_data%meta_nsamples-1
       
       if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
       call this%field_file%read(this%field_data)
       if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', this%field_data%time 
       this%t = this%field_data%time
    
       !> Copy the data from the field_data to the correct field and use the probe interpolator
       call use_interpolator(this)
       
    end do
     

  end subroutine interpolate_field


  subroutine finalize(this)
    class(field_file_interpolator_t), intent(inout) :: this
    
    if (allocated(this%sampled_fields%fields)) deallocate(this%sampled_fields%fields)
    call this%pb%free
    call this%mat_out%free
    call file_free(this%fout)
  
    call this%Xh%free()
    call this%msh%free()
    call this%gs_h%free()

  end subroutine finalize


  subroutine use_interpolator(this)
    class(field_file_interpolator_t), intent(inout) :: this
    integer :: i,j,k,n

    call this%field_data%get_list(this%fields_in_file,this%field_data%size())
    n =  this%dof%size()   

    do i = 1, this%n_fields
       !> Copy the data from the file to our pointer
       call copy(this%sampled_fields%fields(i)%f%x,this%fields_in_file(this%field_pos(i))%v%x,n)
       ! Put it also on device
       if (NEKO_BCKND_DEVICE .eq. 1) then 
         call device_memcpy(this%sampled_fields%fields(i)%f%x, &
                            this%sampled_fields%fields(i)%f%x_d, &
                            n,HOST_TO_DEVICE)
       end if
    end do

    !> Interpolate the desired fields
    call this%pb%interpolate(this%t,1, this%write_output)
    !! Write if the interpolate function returs write_output=.true.
    if (this%write_output) then
       this%mat_out%x = this%pb%out_fields
       call this%fout%write(this%mat_out, this%t)
       this%write_output = .false.
    end if

  end subroutine use_interpolator

  subroutine initialize_interpolator(this,params)
    class(field_file_interpolator_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(matrix_t) :: mat_coords

    !> Initialize the probes asumming t=0
    !> Read the output information
    call json_get(params, 'case.probes.output_file', this%output_file) 
    !! Read probe info and initialize the controller, arrays, etc.
    call this%pb%init(0.0_rp, params, this%coef%Xh)
    !! Perform the set up of gslib_findpts
    call this%pb%setup(this%coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call this%pb%map(this%coef)
    !> Write a summary of the probe info
    call this%pb%show()
    !> Initialize the output
    this%fout = file_t(trim(this%output_file))
    call this%mat_out%init(this%pb%n_probes, this%pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(this%pb%n_probes,3)
    !! Array them as rows
    call trsp(mat_coords%x, this%pb%n_probes, this%pb%xyz, 3)
    !! Write the data to the file
    call this%fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

  end subroutine initialize_interpolator
    
  subroutine adjust_mesh_deformation_and_create_space(this)
    class(field_file_interpolator_t), intent(inout) :: this
    integer :: i,j,k,e,lx 

    !> Copy the mesh from the read field to ensure mesh deformation
    lx = this%field_data%lx
    !To make sure any deformation made in the user file is passed onto here as well
    do i = 1,this%msh%nelv
       this%msh%elements(i)%e%pts(1)%p%x(1) = this%field_data%x%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(1) = this%field_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(1) = this%field_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(1) = this%field_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(1) = this%field_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(1) = this%field_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(1) = this%field_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(1) = this%field_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

       this%msh%elements(i)%e%pts(1)%p%x(2) = this%field_data%y%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(2) = this%field_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(2) = this%field_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(2) = this%field_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(2) = this%field_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(2) = this%field_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(2) = this%field_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(2) = this%field_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

       this%msh%elements(i)%e%pts(1)%p%x(3) = this%field_data%z%x(linear_index(1,1,1,i,lx,lx,lx))  
       this%msh%elements(i)%e%pts(2)%p%x(3) = this%field_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(3)%p%x(3) = this%field_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(4)%p%x(3) = this%field_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(5)%p%x(3) = this%field_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(6)%p%x(3) = this%field_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(7)%p%x(3) = this%field_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
       this%msh%elements(i)%e%pts(8)%p%x(3) = this%field_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
    end do
   
    !> Based on data read, initialize dofmap, gs, coef
    call this%Xh%init(GLL, this%field_data%lx, this%field_data%ly, this%field_data%lz)
    this%dof = dofmap_t(this%msh, this%Xh)

    !> Update dof map to account for non linear defformations of the elements
    !! %apply_deform is the last thing done in dofmat_init before copying to the gpu
    !! so now we will just override that with direct copy from the file
    do e = 1,this%msh%nelv
       do k = 1,lx
          do j = 1,lx
             do i = 1,lx
                this%dof%x(i,j,k,e) = this%field_data%x%x(linear_index(i,j,k,e,lx,lx,lx))  
                this%dof%y(i,j,k,e) = this%field_data%y%x(linear_index(i,j,k,e,lx,lx,lx))  
                this%dof%z(i,j,k,e) = this%field_data%z%x(linear_index(i,j,k,e,lx,lx,lx))   
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

end module field_file_interpolator
