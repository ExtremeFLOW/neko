module field_file_interpolator
  use neko
  implicit none
  private

  !> Provides tools to interpolate from data written to file
  !! @details
  !! no details
  type, public :: field_file_interpolator_t
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
     !> Pointer for fields to add to the registry
     type(field_t), pointer :: field_pointer
     !> List to point to the fields inside the file
     type(vector_ptr_t), allocatable :: fields(:)
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
  !! @param u u velocity field
  !! @param v v velocity field
  !! @param w w velocity field
  !! @param coef type with all geometrical variables
  subroutine initialize(this, params)
    class(field_file_interpolator_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
     
    !> Support variables for probes 
     type(matrix_t) :: mat_coords

     !> Define the file names. The counters are defined in the .nek5000 file 
     mesh_fname  = 'cylinder.nmsh'
     field_fname = 'mean_field0.fld'
     file_point = 1  
     mesh_file = file_t(trim(mesh_fname))
     field_file = file_t(trim(field_fname))

     !> Read the mesh
     call mesh_file%read(msh) 
     call field_data%init(msh%nelv,msh%offset_el)

     !> Read the first field in sequence
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
     call field_file%read(field_data)
     if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', field_data%time 
     t = field_data%time

     !> Copy the mesh from the read field to ensure mesh deformation
     lx = field_data%lx
     !To make sure any deformation made in the user file is passed onto here as well
     do i = 1,msh%nelv
        msh%elements(i)%e%pts(1)%p%x(1) = field_data%x%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(1) = field_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(1) = field_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(1) = field_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(1) = field_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(1) = field_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(1) = field_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(1) = field_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

        msh%elements(i)%e%pts(1)%p%x(2) = field_data%y%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(2) = field_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(2) = field_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(2) = field_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(2) = field_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(2) = field_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(2) = field_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(2) = field_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

        msh%elements(i)%e%pts(1)%p%x(3) = field_data%z%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(3) = field_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(3) = field_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(3) = field_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(3) = field_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(3) = field_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(3) = field_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(3) = field_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
     end do
   
     !> Based on data read, initialize dofmap, gs, coef
     call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)
     dof = dofmap_t(msh, Xh)
     call gs_h%init(dof)
     call coef%init(gs_h)
 
     !> Now create the field in the field registry such that probes can find it
     call neko_field_registry%add_field(dof, 'v1')
     v1 => neko_field_registry%get_field('v1')
     n = v1%dof%size()

    !> ========== Initialize Probes =================
    
    !> Read the output information
    call json_get(params, 'case.probes.output_file', output_file) 

    !> Probe set up
    !! Read probe info and initialize the controller, arrays, etc.
    call pb%init(0.0_rp, params, coef%Xh)
    !! Perform the set up of gslib_findpts
    call pb%setup(coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call pb%map(coef)
    !> Write a summary of the probe info
    call pb%show()

    !> Initialize the output
    fout = file_t(trim(output_file))
    call mat_out%init(pb%n_probes, pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(pb%n_probes,3)
    !! Array them as rows
    call trsp(mat_coords%x, pb%n_probes, pb%xyz, 3)
    !! Write the data to the file
    call fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

    !> ==============================================

     
     !> Allocate list to contain the fields that are read
     allocate(fields(field_data%size()))
     !> Get a list that contains pointers to the fields in the file
     call field_data%get_list(fields,field_data%size())




  end subroutine initialize


end module field_file_interpolator



module user
  use neko
  use time_based_controller, only : time_based_controller_t
  use tensor, only : trsp
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => user_check
  end subroutine user_setup

  
  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !> Variables for bias
    real(kind=rp) :: arg, bias
    
    arg  = 0
    bias = x*0.2*exp(arg)
    bias = 0

    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z + bias

  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_u_and_s(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    s => neko_field_registry%get_field('s')

    !> Initialize with zero velocity
    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())

    !> Initialize with rand perturbations on temperature
    call rzero(s%x,w%dof%size())
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) 
    end do
    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e+s%msh%offset_el,rp)*real(i*j*k,rp))
                r = sqrt(s%dof%x(i,j,k,e)**2+s%dof%y(i,j,k,e)**2)
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001*rand*s%dof%x(i,j,k,e)*&
                                                    sin(3*pi*r/0.05_rp)*sin(10*pi*z)
             end do
          end do
       end do
    end do
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine set_initial_conditions_for_u_and_s

  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    if (pe_rank.eq.0) write(*,*) 'Rayleigh Number is Ra=', Ra
    
  end subroutine user_initialize



  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param
 
  end subroutine user_finalize

  subroutine set_bousinesq_forcing_term(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine user_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

  
    call interpolate_from_file(params)

 
  end subroutine user_check


  subroutine interpolate_from_file(params)
     type(json_file), intent(inout) :: params
     character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, field_fname, hom_dir, output_fname
     type(file_t) :: field_file, output_file_fld, mesh_file
     real(kind=rp) :: start_time, el_h, el_dim(3,3), domain_height
     real(kind=rp), allocatable :: temp_el(:,:,:)
     type(fld_file_data_t) :: field_data
     type(coef_t) :: coef
     type(dofmap_t) :: dof
     type(space_t) :: Xh
     type(mesh_t) :: msh
     type(gs_t) :: gs_h
     type(field_t), pointer :: v1, avg_u, old_u, el_heights
     type(vector_ptr_t), allocatable :: fields(:)
     integer, allocatable :: hom_dir_el(:)
     integer :: argc, i, n, lx, j, e, n_levels
     integer :: file_point
     real(kind=rp) :: t
     !> ========== Needed for Probes =================
  
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

     !> =============================================
    
     !> Support variables for probes 
     type(matrix_t) :: mat_coords

     !> Define the file names. The counters are defined in the .nek5000 file 
     mesh_fname  = 'cylinder.nmsh'
     field_fname = 'mean_field0.fld'
     file_point = 1  
     mesh_file = file_t(trim(mesh_fname))
     field_file = file_t(trim(field_fname))

     !> Read the mesh
     call mesh_file%read(msh) 
     call field_data%init(msh%nelv,msh%offset_el)

     !> Read the first field in sequence
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
     call field_file%read(field_data)
     if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', field_data%time 
     t = field_data%time

     !> Copy the mesh from the read field to ensure mesh deformation
     lx = field_data%lx
     !To make sure any deformation made in the user file is passed onto here as well
     do i = 1,msh%nelv
        msh%elements(i)%e%pts(1)%p%x(1) = field_data%x%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(1) = field_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(1) = field_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(1) = field_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(1) = field_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(1) = field_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(1) = field_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(1) = field_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

        msh%elements(i)%e%pts(1)%p%x(2) = field_data%y%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(2) = field_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(2) = field_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(2) = field_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(2) = field_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(2) = field_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(2) = field_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(2) = field_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

        msh%elements(i)%e%pts(1)%p%x(3) = field_data%z%x(linear_index(1,1,1,i,lx,lx,lx))  
        msh%elements(i)%e%pts(2)%p%x(3) = field_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(3)%p%x(3) = field_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(4)%p%x(3) = field_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
        msh%elements(i)%e%pts(5)%p%x(3) = field_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(6)%p%x(3) = field_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(7)%p%x(3) = field_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
        msh%elements(i)%e%pts(8)%p%x(3) = field_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
     end do
   
     !> Based on data read, initialize dofmap, gs, coef
     call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)
     dof = dofmap_t(msh, Xh)
     call gs_h%init(dof)
     call coef%init(gs_h)
 
     !> Now create the field in the field registry such that probes can find it
     call neko_field_registry%add_field(dof, 'v1')
     v1 => neko_field_registry%get_field('v1')
     n = v1%dof%size()

    !> ========== Initialize Probes =================
    
    !> Read the output information
    call json_get(params, 'case.probes.output_file', output_file) 

    !> Probe set up
    !! Read probe info and initialize the controller, arrays, etc.
    call pb%init(0.0_rp, params, coef%Xh)
    !! Perform the set up of gslib_findpts
    call pb%setup(coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call pb%map(coef)
    !> Write a summary of the probe info
    call pb%show()

    !> Initialize the output
    fout = file_t(trim(output_file))
    call mat_out%init(pb%n_probes, pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(pb%n_probes,3)
    !! Array them as rows
    call trsp(mat_coords%x, pb%n_probes, pb%xyz, 3)
    !! Write the data to the file
    call fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

    !> ==============================================

     
     !> Allocate list to contain the fields that are read
     allocate(fields(field_data%size()))
     !> Get a list that contains pointers to the fields in the file
     call field_data%get_list(fields,field_data%size())
      
          
     !> Copy the data from the file to our pointer
     call copy(v1%x,fields(file_point)%v%x,v1%dof%size())
     ! Put it also on device
     if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(v1%x,v1%x_d, &
                      v1%dof%size(),HOST_TO_DEVICE)
     end if

     !> ========== Needed for Probes =================

     !> Interpolate the desired fields
     call pb%interpolate(t,1, write_output)
     !! Write if the interpolate function returs write_output=.true.
     if (write_output) then
        mat_out%x = pb%out_fields
        call fout%write(mat_out, t)
        write_output = .false.
     end if

     !> ==============================================
     

     !> Read the rest of the fields in the sequence
     do i = 1, field_data%meta_nsamples-1
        if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
        call field_file%read(field_data)
        if (pe_rank .eq. 0) write(*,*) 'time_in_file= ', field_data%time 
        t = field_data%time
        !> Get a list that contains pointers to the fields in the file
        call field_data%get_list(fields,field_data%size())
       
        !> Copy the data from the file to our pointer
        call copy(v1%x,fields(file_point)%v%x,v1%dof%size())
        ! Put it also on device
        if (NEKO_BCKND_DEVICE .eq. 1) then 
           call device_memcpy(v1%x,v1%x_d, &
                          v1%dof%size(),HOST_TO_DEVICE)
        end if
    
        !> ========== Needed for Probes =================

        !> Interpolate the desired fields
        call pb%interpolate(t,i+1, write_output)
        !! Write if the interpolate function returs write_output=.true.
        if (write_output) then
           mat_out%x = pb%out_fields
           call fout%write(mat_out, t)
           write_output = .false.
        end if

    !> ==============================================

     end do
    
     call pb%free
     call mat_out%free
     call file_free(fout)

  end subroutine interpolate_from_file



end module user
