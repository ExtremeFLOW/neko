module fld_file_data
  use field
  use vector
  implicit none

  type, public :: fld_file_data_t
     type(vector_t) :: x
     type(vector_t) :: y
     type(vector_t) :: z
     type(vector_t) :: u
     type(vector_t) :: v
     type(vector_t) :: w
     type(vector_t) :: p
     type(vector_t) :: t
     integer, allocatable :: idx(:)
     type(vector_t), allocatable :: s(:)
     integer :: gdim
     integer :: n_scalars = 0
     real(kind=rp) :: time = 0.0
     integer :: glb_nelv = 0 
     integer :: nelv = 0 
     integer :: offset_nelv = 0
     integer :: lx = 0
     integer :: ly = 0
     integer :: lz = 0
     integer :: t_counter = 0
     ! meta file information (if any)
     integer :: meta_nsamples = 0
     integer :: meta_start_counter = 0
     character(len=1024) :: fld_series_fname

   contains
     procedure, pass(this) :: init => fld_file_data_init
     procedure, pass(this) :: free => fld_file_data_free
  end type fld_file_data_t

contains
  !> Initialise a fld_file_data object with nelv elements with a offset_nel
  subroutine fld_file_data_init(this, nelv, offset_nelv)
    class(fld_file_data_t), intent(inout) :: this
    integer, intent(in) :: nelv, offset_nelv
    call this%free()
    this%nelv = nelv
    this%offset_nelv = offset_nelv
    
  end subroutine fld_file_data_init

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
    this%offset_nelv = 0
    this%lx = 0
    this%ly = 0
    this%lz = 0
    this%t_counter = 0
    this%meta_nsamples = 0
    this%meta_start_counter = 0
  end subroutine fld_file_data_free

end module fld_file_data
