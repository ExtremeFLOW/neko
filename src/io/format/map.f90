!> NEKTON map
!! @todo Figure out a better name for this module
module map
  use mesh
  implicit none

  !> NEKTON vertex mapping
  type :: map_t
     integer :: nel, nlv
     integer, allocatable :: imap(:)
     integer, allocatable :: vertex(:,:)
  end type map_t

  interface map_init
     module procedure map_init_nel_nelv, map_init_mesh
  end interface map_init
contains

  subroutine map_init_mesh(m, msh)
    type(map_t), intent(inout) :: m
    type(mesh_t), intent(in) :: msh

    call map_free(m)

    m%nel = msh%nelv
    m%nlv = msh%npts

    call map_init_common(m)

  end subroutine map_init_mesh

  subroutine map_init_nel_nelv(m, nel, nlv)
    type(map_t), intent(inout) :: m
    integer, intent(in) :: nel
    integer, intent(in) :: nlv

    call map_free(m)

    m%nel = nel
    m%nlv = nlv

    call map_init_common(m)

  end subroutine map_init_nel_nelv

  subroutine map_init_common(m)
    type(map_t), intent(inout) :: m

    allocate(m%imap(m%nel))

    allocate(m%vertex(m%nlv, m%nel))

  end subroutine map_init_common

  subroutine map_free(m)
    type(map_t), intent(inout) :: m

    if (allocated(m%imap)) then
       deallocate(m%imap)
    end if

    if (allocated(m%vertex)) then
       deallocate(m%vertex)
    end if
  end subroutine map_free

end module map
