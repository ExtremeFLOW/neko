!
!> NEKTON session data
!! @details This module is used to represent NEKTON session data
module rea
  use num_types
  use mesh
  implicit none
  private


  !> NEKTON session data struct.
  !! @todo add missing data fields
  type, public :: rea_t
     type(mesh_t) :: msh                     !< Mesh (rep. as a Neko mesh)
     real(kind=dp), allocatable :: params(:) !< Parameters
     character(len=3), allocatable :: cbc(:,:)
  end type rea_t

  
  public :: rea_free

contains
  
  !> Free a NEKTON session data
  subroutine rea_free(r)
    type(rea_t), intent(inout) :: r

    call mesh_free(r%msh)
    
    if (allocated(r%params)) then
       deallocate(r%params)
    end if

  end subroutine rea_free



end module rea
