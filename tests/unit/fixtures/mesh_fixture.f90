module mesh_fixture
  use mesh, only : mesh_t
  use point, only : point_t
  implicit none
  private

  public :: single_unit_hex_mesh
  public :: single_reference_element_mesh
  public :: two_adjacent_unit_hex_mesh
  public :: three_hex_right_triangle_tip_mesh

contains

  subroutine single_unit_hex_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(8)

    ! Single unit hex extruded in z from the unit square in the xy plane:
    !
    !   p3 ----- p4
    !    |       |
    !    |  e1   |
    !    |       |
    !   p1 ----- p2
    !
    !   z = 0: p1, p2, p3, p4
    !   z = 1: p5, p6, p8, p7

    call p(1)%init(0d0, 0d0, 0d0)
    call p(1)%set_id(1)
    call p(2)%init(1d0, 0d0, 0d0)
    call p(2)%set_id(2)
    call p(3)%init(0d0, 1d0, 0d0)
    call p(3)%set_id(4)
    call p(4)%init(1d0, 1d0, 0d0)
    call p(4)%set_id(3)
    call p(5)%init(0d0, 0d0, 1d0)
    call p(5)%set_id(5)
    call p(6)%init(1d0, 0d0, 1d0)
    call p(6)%set_id(6)
    call p(7)%init(1d0, 1d0, 1d0)
    call p(7)%set_id(7)
    call p(8)%init(0d0, 1d0, 1d0)
    call p(8)%set_id(8)

    call msh%init(3, 1)
    call msh%add_element(1, 1, p(1), p(2), p(4), p(3), &
         p(5), p(6), p(7), p(8))
    call msh%generate_conn()
  end subroutine single_unit_hex_mesh

  subroutine single_reference_element_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(8)

    ! Canonical reference hexahedron on [-1, 1]^3 with local node ordering
    ! matching src/mesh/hex.f90.

    call p(1)%init(-1d0, -1d0, -1d0)
    call p(1)%set_id(1)
    call p(2)%init(1d0, -1d0, -1d0)
    call p(2)%set_id(2)
    call p(3)%init(-1d0, 1d0, -1d0)
    call p(3)%set_id(3)
    call p(4)%init(1d0, 1d0, -1d0)
    call p(4)%set_id(4)
    call p(5)%init(-1d0, -1d0, 1d0)
    call p(5)%set_id(5)
    call p(6)%init(1d0, -1d0, 1d0)
    call p(6)%set_id(6)
    call p(7)%init(-1d0, 1d0, 1d0)
    call p(7)%set_id(7)
    call p(8)%init(1d0, 1d0, 1d0)
    call p(8)%set_id(8)

    call msh%init(3, 1)
    call msh%add_element(1, 1, p(1), p(2), p(3), p(4), &
         p(5), p(6), p(7), p(8))
    call msh%generate_conn()
  end subroutine single_reference_element_mesh

  subroutine two_adjacent_unit_hex_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(12)

    ! Two adjacent unit hexes extruded in z from two unit squares in the xy
    ! plane sharing the vertical edge at x = 1:
    !
    !   p3 ----- p4 ----- p6
    !    |   e1  |   e2  |
    !    |       |       |
    !   p1 ----- p2 ----- p5
    !
    ! The shared face is:
    !   element 1 face x = 1 <-> element 2 face x = 1

    call p(1)%init(0d0, 0d0, 0d0)
    call p(1)%set_id(1)
    call p(2)%init(1d0, 0d0, 0d0)
    call p(2)%set_id(2)
    call p(3)%init(0d0, 1d0, 0d0)
    call p(3)%set_id(4)
    call p(4)%init(1d0, 1d0, 0d0)
    call p(4)%set_id(3)
    call p(5)%init(2d0, 0d0, 0d0)
    call p(5)%set_id(5)
    call p(6)%init(2d0, 1d0, 0d0)
    call p(6)%set_id(6)
    call p(7)%init(0d0, 0d0, 1d0)
    call p(7)%set_id(7)
    call p(8)%init(1d0, 0d0, 1d0)
    call p(8)%set_id(8)
    call p(9)%init(1d0, 1d0, 1d0)
    call p(9)%set_id(9)
    call p(10)%init(0d0, 1d0, 1d0)
    call p(10)%set_id(10)
    call p(11)%init(2d0, 0d0, 1d0)
    call p(11)%set_id(11)
    call p(12)%init(2d0, 1d0, 1d0)
    call p(12)%set_id(12)

    call msh%init(3, 2)
    call msh%add_element(1, 1, p(1), p(2), p(4), p(3), &
         p(7), p(8), p(9), p(10))
    call msh%add_element(2, 2, p(2), p(5), p(6), p(4), &
         p(8), p(11), p(12), p(9))
    call msh%generate_conn()
  end subroutine two_adjacent_unit_hex_mesh

  subroutine three_hex_right_triangle_tip_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(18)

    ! Three unit hexes extruded in z from the following 2D layout:
    !
    !   p7 ----- p8 ----- p9      y = 1
    !    |   e2   |
    !   p4 ----- p5 ----- p6      y = 0
    !    |   e3   |   e1
    !   p1 ----- p2 ----- p3      y = -1
    !
    !   x = -1    x = 0    x = 1

    call p(1)%init(-1d0, -1d0, 0d0)
    call p(1)%set_id(1)
    call p(2)%init(0d0, -1d0, 0d0)
    call p(2)%set_id(2)
    call p(3)%init(1d0, -1d0, 0d0)
    call p(3)%set_id(3)
    call p(4)%init(-1d0, 0d0, 0d0)
    call p(4)%set_id(4)
    call p(5)%init(0d0, 0d0, 0d0)
    call p(5)%set_id(5)
    call p(6)%init(1d0, 0d0, 0d0)
    call p(6)%set_id(6)
    call p(7)%init(-1d0, 1d0, 0d0)
    call p(7)%set_id(7)
    call p(8)%init(0d0, 1d0, 0d0)
    call p(8)%set_id(8)
    call p(9)%init(1d0, 1d0, 0d0)
    call p(9)%set_id(9)

    call p(10)%init(-1d0, -1d0, 1d0)
    call p(10)%set_id(10)
    call p(11)%init(0d0, -1d0, 1d0)
    call p(11)%set_id(11)
    call p(12)%init(1d0, -1d0, 1d0)
    call p(12)%set_id(12)
    call p(13)%init(-1d0, 0d0, 1d0)
    call p(13)%set_id(13)
    call p(14)%init(0d0, 0d0, 1d0)
    call p(14)%set_id(14)
    call p(15)%init(1d0, 0d0, 1d0)
    call p(15)%set_id(15)
    call p(16)%init(-1d0, 1d0, 1d0)
    call p(16)%set_id(16)
    call p(17)%init(0d0, 1d0, 1d0)
    call p(17)%set_id(17)
    call p(18)%init(1d0, 1d0, 1d0)
    call p(18)%set_id(18)

    call msh%init(3, 3)

    ! Element adjacent to the x-axis cathetus.
    call msh%add_element(1, 1, p(2), p(3), p(6), p(5), &
         p(11), p(12), p(15), p(14))

    ! Element adjacent to the y-axis cathetus.
    call msh%add_element(2, 2, p(4), p(5), p(8), p(7), &
         p(13), p(14), p(17), p(16))

    ! Connector element touching the triangle only at the tip.
    call msh%add_element(3, 3, p(1), p(2), p(5), p(4), &
         p(10), p(11), p(14), p(13))

    call msh%generate_conn()
  end subroutine three_hex_right_triangle_tip_mesh

end module mesh_fixture
