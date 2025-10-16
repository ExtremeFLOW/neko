!> NEKTON map
!! @todo Figure out a better name for this module
module map
  use mesh, only : mesh_t
  implicit none

  !> NEKTON vertex mapping
  type :: map_t
     integer :: nel, nlv
     integer, allocatable :: imap(:)
     integer, allocatable :: vertex(:,:)
   contains
     !> Contrutctor
     generic :: init => init_nel_nelv, init_mesh
     procedure, private, pass(this) :: init_common => map_init_common
     procedure, private, pass(this) :: init_nel_nelv => map_init_nel_nelv
     procedure, private, pass(this) :: init_mesh => map_init_mesh
     procedure, pass(this) :: init_nel_ => map_free
     !> Destructor
     procedure, pass(this) :: free => map_free
  end type map_t

contains

  subroutine map_init_mesh(this, msh)
    class(map_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh

    call this%free()

    this%nel = msh%nelv
    this%nlv = msh%npts

    call this%init_common()

  end subroutine map_init_mesh

  subroutine map_init_nel_nelv(this, nel, nlv)
    class(map_t), intent(inout) :: this
    integer, intent(in) :: nel
    integer, intent(in) :: nlv

    call this%free()

    this%nel = nel
    this%nlv = nlv

    call this%init_common()

  end subroutine map_init_nel_nelv

  subroutine map_init_common(this)
    class(map_t), intent(inout) :: this

    allocate(this%imap(this%nel))

    allocate(this%vertex(this%nlv, this%nel))

  end subroutine map_init_common

  subroutine map_free(this)
    class(map_t), intent(inout) :: this

    if (allocated(this%imap)) then
       deallocate(this%imap)
    end if

    if (allocated(this%vertex)) then
       deallocate(this%vertex)
    end if
  end subroutine map_free

end module map
