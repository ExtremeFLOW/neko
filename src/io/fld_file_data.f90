!> Simple module to handle fld file series.
!! Provides an interface to the different fields sotred in a fld file
!! Also provides simple functions to scale and add different fld files.
!! An example of using this module is shown in contrib/average_fields.f90
!! The fld_file_data_t should dynamically update each time one reads a new fld
!! file.
!! Martin Karp 1/2-2023
module fld_file_data
  use num_types, only : rp
  use math, only: cmult, add2
  use vector, only : vector_t, vector_ptr_t
  use field, only: field_t
  use dofmap, only: dofmap_t
  use space, only: space_t, GLL
  use global_interpolation, only: global_interpolation_t
  use utils, only: neko_error
  use mesh, only: mesh_t
  implicit none
  private

  type, public :: fld_file_data_t
     type(vector_t) :: x !< x-coords
     type(vector_t) :: y !< y-coords
     type(vector_t) :: z !< z-coords
     type(vector_t) :: u !< x-velocity field
     type(vector_t) :: v !< y-velocity field
     type(vector_t) :: w !< z-velocity field
     type(vector_t) :: p !< pressure field
     type(vector_t) :: t !< temperature
     integer, allocatable :: idx(:) !< element idxs
     type(vector_t), allocatable :: s(:) !< Numbered scalar fields
     integer :: gdim !< spatial dimensions
     integer :: n_scalars = 0 !< number of numbered scalar fields
     real(kind=rp) :: time = 0.0 !< time of sample
     integer :: glb_nelv = 0 !< global number of elements
     integer :: nelv = 0  !< n elements on the pe
     integer :: offset_el = 0 !< element offset for this pe
     integer :: lx = 0 !< N GLL points in x
     integer :: ly = 0
     integer :: lz = 0
     integer :: t_counter = 0 !< counter of samples
     ! meta file information (if any)
     !> number of samples specified in .nek5000 file
     integer :: meta_nsamples = 0
     integer :: meta_start_counter = 0 !< number of first field
     !> name of fld series as specified in .nek5000 (meta) file
     character(len=1024) :: fld_series_fname
   contains
     procedure, pass(this) :: init => fld_file_data_init
     procedure, pass(this) :: free => fld_file_data_free
     procedure, pass(this) :: scale => fld_file_data_scale
     procedure, pass(this) :: add => fld_file_data_add
     procedure, pass(this) :: size => fld_file_data_size
     procedure, pass(this) :: get_list => fld_file_data_get_list
     procedure, pass(this) :: init_same => fld_file_data_init_same
     procedure, pass(this) :: init_n_fields => fld_file_data_init_n_fields
     !> Generates a global_interpolation object to interpolate the fld data.
     procedure, pass(this) :: generate_interpolator => &
          fld_file_data_generate_interpolator
  end type fld_file_data_t

contains

  !> Initialise a fld_file_data object with nelv elements with a offset_nel
  subroutine fld_file_data_init(this, nelv, offset_el)
    class(fld_file_data_t), intent(inout) :: this
    integer, intent(in), optional :: nelv, offset_el

    call this%free()
    if (present(nelv)) this%nelv = nelv
    if (present(offset_el)) this%offset_el = offset_el

  end subroutine fld_file_data_init
  !> Get number of fields in this fld file
  function fld_file_data_size(this) result(i)
    class(fld_file_data_t) :: this
    integer :: i
    i = 0
    if (this%u%n .gt. 0) i = i + 1
    if (this%v%n .gt. 0) i = i + 1
    if (this%w%n .gt. 0) i = i + 1
    if (this%p%n .gt. 0) i = i + 1
    if (this%t%n .gt. 0) i = i + 1
    i = i + this%n_scalars

  end function fld_file_data_size

  !> Genereate same fields as in another fld_file
  subroutine fld_file_data_init_same(this, fld_file, n)
    class(fld_file_data_t), target, intent(inout) :: this
    class(fld_file_data_t), target, intent(in) :: fld_file
    integer, intent(in) :: n
    integer :: i, j
    
    if(fld_file%u%n .gt. 0) then
       call this%u%init(n)
    end if
    if(fld_file%v%n .gt. 0) then
       call this%v%init(n)
    end if
    if(fld_file%w%n .gt. 0) then
       call this%w%init(n)
    end if
    if(fld_file%p%n .gt. 0) then
       call this%p%init(n)
    end if
    if(fld_file%t%n .gt. 0) then
       call this%t%init(n)
    end if
    this%n_scalars = fld_file%n_scalars
    allocate(this%s(fld_file%n_scalars))
    do j = 1, fld_file%n_scalars
       call this%s(j)%init(n)
    end do

  end subroutine fld_file_data_init_same

  !> Genereate same fields as in another fld_file
  subroutine fld_file_data_init_n_fields(this, n_fields, n)
    class(fld_file_data_t), target, intent(inout) :: this
    integer, intent(in) :: n, n_fields
    integer :: i, j


    if(n_fields .gt. 0) then
       call this%u%init(n)
    end if
    if(n_fields .gt. 1) then
       call this%v%init(n)
    end if
    if(n_fields .gt. 2) then
       call this%w%init(n)
    end if
    if(n_fields .gt. 3) then
       call this%p%init(n)
    end if
    if(n_fields .gt. 4) then
       call this%t%init(n)
    end if 
    if (n_fields .gt. 5) then
       this%n_scalars = n_fields-5
       allocate(this%s(this%n_scalars))
       do j = 1, this%n_scalars
          call this%s(j)%init(n)
       end do
    end if

  end subroutine fld_file_data_init_n_fields

  !> Get a list with pointers to the fields in the fld file
  subroutine fld_file_data_get_list(this, ptr_list, n)
    class(fld_file_data_t), target, intent(in) :: this
    integer, intent(in) :: n
    integer :: i, j
    type(vector_ptr_t), intent(inout) :: ptr_list(n)
    i = 1
    if (this%u%n .gt. 0) then
       ptr_list(i)%ptr => this%u
       i = i + 1
    end if
    if (this%v%n .gt. 0) then
       ptr_list(i)%ptr => this%v
       i = i + 1
    end if
    if (this%w%n .gt. 0) then
       ptr_list(i)%ptr => this%w
       i = i + 1
    end if
    if (this%p%n .gt. 0) then
       ptr_list(i)%ptr => this%p
       i = i + 1
    end if
    if (this%t%n .gt. 0) then
       ptr_list(i)%ptr => this%t
       i = i + 1
    end if
    do j = 1, this%n_scalars
       ptr_list(i)%ptr => this%s(j)
       i = i +1
    end do

  end subroutine fld_file_data_get_list



  !> Scale the values stored in this fld_file_data
  subroutine fld_file_data_scale(this, c)
    class(fld_file_data_t), intent(inout) :: this
    real(kind=rp), intent(in) :: c
    integer :: i

    if (this%u%n .gt. 0) call cmult(this%u%x, c, this%u%n)
    if (this%v%n .gt. 0) call cmult(this%v%x, c, this%v%n)
    if (this%w%n .gt. 0) call cmult(this%w%x, c, this%w%n)
    if (this%p%n .gt. 0) call cmult(this%p%x, c, this%p%n)
    if (this%t%n .gt. 0) call cmult(this%t%x, c, this%t%n)

    do i = 1, this%n_scalars
       if (this%s(i)%n .gt. 0) call cmult(this%s(i)%x, c, this%s(i)%n)
    end do

  end subroutine fld_file_data_scale

  !> Add the values in another fld file to this
  subroutine fld_file_data_add(this, fld_data_add)
    class(fld_file_data_t), intent(inout) :: this
    class(fld_file_data_t), intent(in) :: fld_data_add
    integer :: i

    if (this%u%n .gt. 0) call add2(this%u%x, fld_data_add%u%x, this%u%n)
    if (this%v%n .gt. 0) call add2(this%v%x, fld_data_add%v%x, this%v%n)
    if (this%w%n .gt. 0) call add2(this%w%x, fld_data_add%w%x, this%w%n)
    if (this%p%n .gt. 0) call add2(this%p%x, fld_data_add%p%x, this%p%n)
    if (this%t%n .gt. 0) call add2(this%t%x, fld_data_add%t%x, this%t%n)

    do i = 1, this%n_scalars
       if (this%s(i)%n .gt. 0) call add2(this%s(i)%x, fld_data_add%s(i)%x, &
            this%s(i)%n)
    end do
  end subroutine fld_file_data_add

  !> Deallocate fld file data type
  subroutine fld_file_data_free(this)
    class(fld_file_data_t), intent(inout) :: this
    integer :: i
    call this%x%free()
    call this%y%free()
    call this%z%free()
    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%p%free()
    call this%t%free()
    if (allocated(this%s)) then
       do i = 1, this%n_scalars
          call this%s(i)%free()
       end do
       deallocate(this%s)
    end if
    this%n_scalars = 0
    this%time = 0.0
    this%glb_nelv = 0
    this%nelv = 0
    this%offset_el = 0
    this%lx = 0
    this%ly = 0
    this%lz = 0
    this%t_counter = 0
    this%meta_nsamples = 0
    this%meta_start_counter = 0
    if(allocated(this%idx)) deallocate(this%idx)
  end subroutine fld_file_data_free

  !> Generates a global_interpolation object to interpolate the fld data.
  !! @param to_dof Dofmap on which to interpolate.
  !! @param to_msh Mesh on which to interpolate.
  !! @param tolerance Tolerance for the newton iterations.
  function fld_file_data_generate_interpolator(this, to_dof, &
       to_msh, tolerance) result(global_interp)
    class(fld_file_data_t), intent(in) :: this
    type(dofmap_t), intent(in), target :: to_dof
    type(mesh_t), intent(in), target :: to_msh
    real(kind=rp), intent(in) :: tolerance

    type(global_interpolation_t) :: global_interp

    ! --- variables for interpolation
    type(space_t) :: fld_Xh
    real(kind=rp), allocatable :: x_coords(:,:,:,:), y_coords(:,:,:,:), &
         z_coords(:,:,:,:)
    real(kind=rp) :: center_x,  center_y, center_z
    integer :: e, i
    ! ---

    type(space_t), pointer :: to_Xh
    to_Xh => to_dof%Xh

    ! Safeguard in case we didn't read mesh information
    if (.not. allocated(this%x%x) .or. &
         .not. allocated(this%y%x) .or. &
         .not. allocated(this%z%x)) call neko_error("Unable to retrieve &
&mesh information from fld data.")

    ! Create a space based on the fld data
    call fld_Xh%init(GLL, this%lx,  this%ly, this%lz)

    ! These are the coordinates of our current dofmap
    ! that we use for the interpolation
    allocate(x_coords(to_Xh%lx,  to_Xh%ly, to_Xh%lz, to_msh%nelv))
    allocate(y_coords(to_Xh%lx,  to_Xh%ly, to_Xh%lz, to_msh%nelv))
    allocate(z_coords(to_Xh%lx,  to_Xh%ly, to_Xh%lz, to_msh%nelv))

    !> To ensure that each point is within an element
    !! Remedies issue with points on the boundary
    !! Technically gives each point a slightly different value
    !! but still within the specified tolerance
    do e = 1, to_msh%nelv
       center_x = 0d0
       center_y = 0d0
       center_z = 0d0
       do i = 1, to_Xh%lxyz
          center_x = center_x + to_dof%x(i, 1, 1, e)
          center_y = center_y + to_dof%y(i, 1, 1, e)
          center_z = center_z + to_dof%z(i, 1, 1, e)
       end do
       center_x = center_x / to_Xh%lxyz
       center_y = center_y / to_Xh%lxyz
       center_z = center_z / to_Xh%lxyz
       do i = 1, to_Xh%lxyz
          x_coords(i, 1, 1, e) = to_dof%x(i, 1, 1, e) - &
               tolerance * (to_dof%x(i, 1, 1, e) - center_x)
          y_coords(i, 1, 1, e) = to_dof%y(i, 1, 1, e) - &
               tolerance * (to_dof%y(i, 1, 1, e) - center_y)
          z_coords(i, 1, 1, e) = to_dof%z(i, 1, 1, e) - &
               tolerance * (to_dof%z(i, 1, 1, e) - center_z)
       end do
    end do

    ! The initialization is done based on the variables created from
    ! fld data
    call global_interp%init(this%x%x, this%y%x, this%z%x, this%gdim, &
         this%nelv, fld_Xh, tol = tolerance)
    call global_interp%find_points(x_coords, y_coords, z_coords, &
         to_dof%size())

    deallocate(x_coords)
    deallocate(y_coords)
    deallocate(z_coords)

  end function fld_file_data_generate_interpolator

end module fld_file_data
