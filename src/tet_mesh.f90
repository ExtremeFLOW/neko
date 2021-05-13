!> Defines a tetrahedral mesh
!! @details Mesh dervied from an existing hexahedral mesh via bisection
module tet_mesh
  use mesh
  use tet
  use point
  use utils
  implicit none
  
  integer, parameter :: TET_MSH_OTPV = 1

  type :: tet_mesh_t
     type(tet_t), allocatable :: el(:) !< Tetrahedron elements
     type(mesh_t), pointer :: msh      !< Hexahedron mesh
     integer :: nelv                   !< Number of Tetrahedrons
   contains
     procedure, pass(this) :: init => tet_mesh_init
     procedure, pass(this) :: free => tet_mesh_free
  end type tet_mesh_t

  private :: tet_mesh_bisect_otpv

contains

  !> Initialise a tetrahedral mesh based on a hexahedral mesh @a msh
  subroutine tet_mesh_init(this, msh, mthd)
    class(tet_mesh_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    integer, intent(in) :: mthd

    call this%free()
    
    this%msh => msh

    if (mthd .eq. TET_MSH_OTPV) then
       this%nelv = msh%nelv * 8
       allocate(this%el(this%nelv))
       call tet_mesh_bisect_otpv(this)
    else
       call neko_error('Invalid bisection strategy')
    end if

  end subroutine tet_mesh_init

  !> Deallocate a tetrahedral mesh
  subroutine tet_mesh_free(this)
    class(tet_mesh_t), intent(inout) :: this

    if (allocated(this%el)) then
       deallocate(this%el)
    end if

    nullify(this%msh)
    
  end subroutine tet_mesh_free

  ! Bisect hexahedral mesh into a tetrahedral mesh using
  ! one tet per vertex as described in,
  ! P. D. Bello-Maldonado and P. F. Fischer
  ! SIAM J. Sci. Comput., 41(5), S2–S18, 2019.
  subroutine tet_mesh_bisect_otpv(tet_msh)
    type(tet_mesh_t), intent(inout) :: tet_msh
    integer :: i, j
    type(point_t), pointer :: p1, p2, p3, p4
    
    j = 0
    do i = 1, tet_msh%msh%nelv

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(5)%p
       p4 => tet_msh%msh%elements(i)%e%pts(4)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(2)%p
       p3 => tet_msh%msh%elements(i)%e%pts(3)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(3)%p
       p3 => tet_msh%msh%elements(i)%e%pts(2)%p
       p4 => tet_msh%msh%elements(i)%e%pts(7)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(3)%p
       p3 => tet_msh%msh%elements(i)%e%pts(8)%p
       p4 => tet_msh%msh%elements(i)%e%pts(1)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(1)%p
       p2 => tet_msh%msh%elements(i)%e%pts(5)%p
       p3 => tet_msh%msh%elements(i)%e%pts(6)%p
       p4 => tet_msh%msh%elements(i)%e%pts(8)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(2)%p
       p2 => tet_msh%msh%elements(i)%e%pts(5)%p
       p3 => tet_msh%msh%elements(i)%e%pts(6)%p
       p4 => tet_msh%msh%elements(i)%e%pts(7)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(3)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(6)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

       j = j + 1
       p1 => tet_msh%msh%elements(i)%e%pts(4)%p
       p2 => tet_msh%msh%elements(i)%e%pts(8)%p
       p3 => tet_msh%msh%elements(i)%e%pts(7)%p
       p4 => tet_msh%msh%elements(i)%e%pts(5)%p
       call tet_msh%el(j)%init(j, p1, p2, p3, p4)

    end do
    
  end subroutine tet_mesh_bisect_otpv
  
end module tet_mesh
