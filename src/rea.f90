!
!> NEKTON session data
!! @details This module is used to represent NEKTON session data
!
module rea
  use generic_file
  use num_types
  use utils
  use mesh
  implicit none
  private

  !> Interface for NEKTON ascii files
  type, extends(generic_file_t) :: rea_file_t
   contains
     procedure :: read => rea_file_read
     procedure :: write => rea_file_write
  end type rea_file_t

  !> NEKTON session data struct.
  !! @todo add missing data fields
  type rea_t
     type(mesh_t) :: msh                     !< Mesh (rep. as a Neko mesh)
     real(kind=dp), allocatable :: params(:) !< Parameters
  end type rea_t
  
  public :: rea_file_t, rea_t
  
contains
  
  !> Free a NEKTON session data
  subroutine rea_free(r)
    type(rea_t), intent(inout) :: r

    call mesh_free(r%msh)
    
    if (allocated(r%params)) then
       deallocate(r%params)
    end if

  end subroutine rea_free

  !> Load NEKTON session data from an ascii file
  subroutine rea_file_read(this, data)
    class(rea_file_t) :: this
    class(*), target, intent(inout) :: data
    type(mesh_t), pointer :: msh
    real(kind=dp), pointer :: params(:)
    integer :: ndim, nparam, nskip, nlogic
    integer :: nelgs, nelgv, i, j, ierr
    logical :: read_param

    select type(data)
    type is (rea_t)
       call rea_free(data)       
       msh => data%msh
       params => data%params
       read_param = .true.
    type is (mesh_t)    
       msh => data
       read_param = .false.
    class default
       call neko_error('Invalid output data')
    end select

    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    write(*, '(A,A)') ' Reading ', this%fname
    
    read(9, *)
    read(9, *)
    read(9, *) ndim
    read(9, *) nparam
    
    if (.not. read_param) then
       ! Skip parameters
       do i = 1, nparam
          read(9, *)
       end do
    else       
       allocate(params(nparam))
       do i = 1, nparam
          read(9, *) params(i)
       end do
    end if
    
    ! Skip passive scalars
    read(9, *) nskip
    do i = 1, nskip
       read(9, *)
    end do
    
    ! Skip logic switches
    read(9, *) nlogic
    do i = 1, nlogic
       read(9, *)
    end do
    
    ! Read mesh info
    read(9, *)
    read(9, *)
    read(9, *) nelgs,ndim, nelgv
    if (nelgs .lt. 0) then
       !> @todo Add support to load binary NEKTON meshes
       call neko_error('Binary NEKTON meshes are not supported')
    end if

    write(*,*) nelgs, ndim, nelgv
    
    
    call mesh_init_coordinates(msh, ndim, nelgv)       
    do i = 1, nelgv
       read(9, *)
       if (ndim .eq. 2) then
          read(9, *) (msh%xc(j, i),j=1,4)
          read(9, *) (msh%yc(j, i),j=1,4)
       else if (ndim .eq. 3) then
          read(9, *) (msh%xc(j, i),j=1,4)
          read(9, *) (msh%yc(j, i),j=1,4)
          read(9, *) (msh%zc(j, i),j=1,4)
          read(9, *) (msh%xc(j, i),j=5,8)
          read(9, *) (msh%yc(j, i),j=5,8)
          read(9, *) (msh%zc(j, i),j=5,8)
       end if
    end do
    !> @todo Add support for curved side data

    close(9)
    
  end subroutine rea_file_read


  subroutine rea_file_write(this, data)
    class(rea_file_t) :: this
    class(*), target, intent(in) :: data
  end subroutine rea_file_write
end module rea
