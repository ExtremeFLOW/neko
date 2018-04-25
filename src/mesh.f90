!> Defines a mesh
module mesh
  use num_types
  implicit none

  type mesh_t

     integer :: lelv
     integer :: dim 

     integer :: lx1
     integer :: ly1
     integer :: lz1

     real(kind=dp), allocatable :: xc(:,:)   !< X-coordinates
     real(kind=dp), allocatable :: yc(:,:)   !< Y-coordinates
     real(kind=dp), allocatable :: zc(:,:)   !< Z-coordinates

  end type mesh_t

contains 

  subroutine mesh_init_coordinates(m, ndim, nelv)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: ndim
    integer, intent(in) :: nelv
    
    if (.not. allocated(m%xc)) then
       allocate(m%xc(4, nelv))
    end if

    if (.not. allocated(m%yc)) then
       allocate(m%yc(4, nelv))
    end if

    if (ndim .gt. 2 .and. (.not. allocated(m%zc))) then
       allocate(m%zc(4, nelv))
    end if
    
  end subroutine mesh_init_coordinates
  
  subroutine mesh_free(m)
    type(mesh_t), intent(inout) :: m

    if (allocated(m%xc)) then
       deallocate(m%xc)
    end if

    if (allocated(m%yc)) then
       deallocate(m%yc)
    end if

    if (allocated(m%zc)) then
       deallocate(m%zc)
    end if

  end subroutine mesh_free


end module mesh
