!> Simple module to handle fld file series.
!! Provides an interface to the different fields sotred in a fld file
!! Also provides simple functions to scale and add different fld files.
!! An example of using this module is shown in contrib/average_fields.f90
!! The fld_file_data_t should dynamically update each time one reads a new fld file
!! Martin Karp 1/2-2023
module fld_file_data
  use field
  use vector
  use math
  implicit none
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
     integer :: meta_nsamples = 0 !< number of samples specified in .nek5000 file
     integer :: meta_start_counter = 0 !< number of first field
     character(len=1024) :: fld_series_fname !< name of fld series as specified in .nek5000 (meta) file

   contains
     procedure, pass(this) :: init => fld_file_data_init
     procedure, pass(this) :: free => fld_file_data_free
     procedure, pass(this) :: scale => fld_file_data_scale
     procedure, pass(this) :: add => fld_file_data_add
     procedure, pass(this) :: size => fld_file_data_size
     procedure, pass(this) :: get_list => fld_file_get_list
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
    if(this%u%n .gt. 0) i = i + 1 
    if(this%v%n .gt. 0) i = i + 1 
    if(this%w%n .gt. 0) i = i + 1 
    if(this%p%n .gt. 0) i = i + 1 
    if(this%t%n .gt. 0) i = i + 1 
    i = i + this%n_scalars

  end function fld_file_data_size

  !> Get a list with pointers to the fields in the fld file
  subroutine fld_file_get_list(this, ptr_list, n) 
    class(fld_file_data_t), target, intent(in) :: this
    integer, intent(in) :: n
    integer :: i, j
    type(vector_ptr_t), intent(inout) :: ptr_list(n)
    i = 1
    if(this%u%n .gt. 0) then
       ptr_list(i)%v => this%u
       i = i + 1
    end if     
    if(this%v%n .gt. 0) then
       ptr_list(i)%v => this%v
       i = i + 1
    end if     
    if(this%w%n .gt. 0) then
       ptr_list(i)%v => this%w
       i = i + 1
    end if     
    if(this%p%n .gt. 0) then
       ptr_list(i)%v => this%p
       i = i + 1
    end if     
    if(this%t%n .gt. 0) then
       ptr_list(i)%v => this%t
       i = i + 1
    end if     
    do j = 1, this%n_scalars
       ptr_list(i)%v => this%s(j)
       i = i +1
    end do

  end subroutine fld_file_get_list



  !> Scale the values stored in this fld_file_data
  subroutine fld_file_data_scale(this, c)
    class(fld_file_data_t), intent(inout) :: this
    real(kind=rp), intent(in) :: c
    integer :: i

    if(this%u%n .gt. 0) call cmult(this%u%x,c,this%u%n)
    if(this%v%n .gt. 0) call cmult(this%v%x,c,this%v%n)
    if(this%w%n .gt. 0) call cmult(this%w%x,c,this%w%n)
    if(this%p%n .gt. 0) call cmult(this%p%x,c,this%p%n)
    if(this%t%n .gt. 0) call cmult(this%t%x,c,this%t%n)

    do i = 1, this%n_scalars
       if(this%s(i)%n .gt. 0) call cmult(this%s(i)%x,c,this%s(i)%n)
    end do

  end subroutine fld_file_data_scale

  !> Add the values in another fld file to this
  subroutine fld_file_data_add(this, fld_data_add)
    class(fld_file_data_t), intent(inout) :: this
    class(fld_file_data_t), intent(in) :: fld_data_add
    integer :: i

    if(this%u%n .gt. 0) call add2(this%u%x,fld_data_add%u%x,this%u%n)
    if(this%v%n .gt. 0) call add2(this%v%x,fld_data_add%v%x,this%v%n)
    if(this%w%n .gt. 0) call add2(this%w%x,fld_data_add%w%x,this%w%n)
    if(this%p%n .gt. 0) call add2(this%p%x,fld_data_add%p%x,this%p%n)
    if(this%t%n .gt. 0) call add2(this%t%x,fld_data_add%t%x,this%t%n)

    do i = 1, this%n_scalars
       if(this%s(i)%n .gt. 0) call add2(this%s(i)%x,fld_data_add%s(i)%x,this%s(i)%n)
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
  end subroutine fld_file_data_free

end module fld_file_data
