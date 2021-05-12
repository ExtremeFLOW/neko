!> Defines a tetrahedral mesh
!! @details Mesh dervied from an existing hexahedral mesh via bisection
module tet_mesh
  use mesh
  use tet
  use point
  implicit none
  
  type :: tet_mesh_t
     type(tet_t), allocatable :: el(:)
     type(mesh_t), pointer :: msh
   contains
     procedure, pass(this) :: init => tet_mesh_init
     procedure, pass(this) :: free => tet_mesh_free
  end type tet_mesh_t

contains

  !> Initialise a tetrahedral mesh based on a hexahedral mesh @a msh
  subroutine tet_mesh_init(this, msh)
    class(tet_mesh_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    integer :: i, j
    type(point_t), pointer :: p1, p2, p3, p4

    call this%free()
    
    this%msh => msh

    allocate(this%el(msh%nelv * 8))


    ! Bisect hexahedral mesh into a tetrahedral mesh using
    ! one tet per vertex as described in,
    ! P. D. Bello-Maldonado and P. F. Fischer
    ! SIAM J. Sci. Comput., 41(5), S2–S18, 2019.
    j = 0
    do i = 1, msh%nelv

       j = j + 1
       p1 => msh%elements(i)%e%pts(1)%p
       p2 => msh%elements(i)%e%pts(2)%p
       p3 => msh%elements(i)%e%pts(5)%p
       p4 => msh%elements(i)%e%pts(4)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(1)%p
       p2 => msh%elements(i)%e%pts(2)%p
       p3 => msh%elements(i)%e%pts(3)%p
       p4 => msh%elements(i)%e%pts(6)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(4)%p
       p2 => msh%elements(i)%e%pts(3)%p
       p3 => msh%elements(i)%e%pts(2)%p
       p4 => msh%elements(i)%e%pts(7)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(4)%p
       p2 => msh%elements(i)%e%pts(3)%p
       p3 => msh%elements(i)%e%pts(8)%p
       p4 => msh%elements(i)%e%pts(1)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(1)%p
       p2 => msh%elements(i)%e%pts(5)%p
       p3 => msh%elements(i)%e%pts(6)%p
       p4 => msh%elements(i)%e%pts(8)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(2)%p
       p2 => msh%elements(i)%e%pts(5)%p
       p3 => msh%elements(i)%e%pts(6)%p
       p4 => msh%elements(i)%e%pts(7)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(3)%p
       p2 => msh%elements(i)%e%pts(8)%p
       p3 => msh%elements(i)%e%pts(7)%p
       p4 => msh%elements(i)%e%pts(6)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => msh%elements(i)%e%pts(4)%p
       p2 => msh%elements(i)%e%pts(8)%p
       p3 => msh%elements(i)%e%pts(7)%p
       p4 => msh%elements(i)%e%pts(5)%p
       call this%el(j)%init(j, p1, p2, p3, p4)

    end do
    
  end subroutine tet_mesh_init

  !> Deallocate a tetrahedral mesh
  subroutine tet_mesh_free(this)
    class(tet_mesh_t), intent(inout) :: this

    if (allocated(this%el)) then
       deallocate(this%el)
    end if

    nullify(this%msh)
    
  end subroutine tet_mesh_free
  
end module tet_mesh
