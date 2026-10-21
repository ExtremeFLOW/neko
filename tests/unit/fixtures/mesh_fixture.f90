module mesh_fixture
  use num_types, only : dp
  use mesh, only : mesh_t
  use point, only : point_t
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use datadist, only : linear_dist_t
  implicit none
  private

  public :: single_unit_hex_mesh
  public :: single_reference_element_mesh
  public :: single_skewed_element_mesh
  public :: two_adjacent_unit_hex_mesh
  public :: three_hex_right_triangle_tip_mesh
  public :: stacked_unit_hex_column_mesh

contains

  subroutine single_unit_hex_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(8)

    ! Single unit hex on [0, 1]^3. The points are listed in the NEKTON
    ! symmetric ordering of src/mesh/hex.f90, i.e. x fastest, then y, then z.
    ! Passing them counterclockwise instead twists the element and makes the
    ! metric lose positive definiteness.
    !
    !   z = 0:  p3 ----- p4      z = 1:  p7 ----- p8
    !            |       |                |       |
    !            |  e1   |                |  e1   |
    !            |       |                |       |
    !           p1 ----- p2              p5 ----- p6

    call p(1)%init(0d0, 0d0, 0d0)
    call p(1)%set_id(1)
    call p(2)%init(1d0, 0d0, 0d0)
    call p(2)%set_id(2)
    call p(3)%init(0d0, 1d0, 0d0)
    call p(3)%set_id(3)
    call p(4)%init(1d0, 1d0, 0d0)
    call p(4)%set_id(4)
    call p(5)%init(0d0, 0d0, 1d0)
    call p(5)%set_id(5)
    call p(6)%init(1d0, 0d0, 1d0)
    call p(6)%set_id(6)
    call p(7)%init(0d0, 1d0, 1d0)
    call p(7)%set_id(7)
    call p(8)%init(1d0, 1d0, 1d0)
    call p(8)%set_id(8)

    call msh%init(3, 1)
    call msh%add_element(1, 1, p(1), p(2), p(3), p(4), &
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

  subroutine single_skewed_element_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(8)
    real(kind=kind(1.0d0)), parameter :: tilt = 2.0d-2

    ! Affine image of the reference hexahedron with z = t + tilt * r.
    ! The t-normal therefore has a small x component while remaining close
    ! enough to the z-axis to exercise the near-z tangent construction.

    call p(1)%init(-1d0, -1d0, -1d0 - tilt)
    call p(1)%set_id(1)
    call p(2)%init(1d0, -1d0, -1d0 + tilt)
    call p(2)%set_id(2)
    call p(3)%init(-1d0, 1d0, -1d0 - tilt)
    call p(3)%set_id(3)
    call p(4)%init(1d0, 1d0, -1d0 + tilt)
    call p(4)%set_id(4)
    call p(5)%init(-1d0, -1d0, 1d0 - tilt)
    call p(5)%set_id(5)
    call p(6)%init(1d0, -1d0, 1d0 + tilt)
    call p(6)%set_id(6)
    call p(7)%init(-1d0, 1d0, 1d0 - tilt)
    call p(7)%set_id(7)
    call p(8)%init(1d0, 1d0, 1d0 + tilt)
    call p(8)%set_id(8)

    call msh%init(3, 1)
    call msh%add_element(1, 1, p(1), p(2), p(3), p(4), &
         p(5), p(6), p(7), p(8))
    call msh%generate_conn()
  end subroutine single_skewed_element_mesh

  subroutine two_adjacent_unit_hex_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(12)

    ! Two adjacent unit hexes extruded in z from two unit squares in the xy
    ! plane sharing the vertical face at x = 1. The points are listed in the
    ! NEKTON symmetric ordering of src/mesh/hex.f90, i.e. x fastest, then y,
    ! then z. Passing them counterclockwise instead twists both elements and
    ! makes the metric lose positive definiteness.
    !
    !   z = 0:   p4 ---- p5 ---- p6     z = 1:  p10 --- p11 --- p12
    !             |  e1  |  e2  |                |  e1  |  e2  |
    !             |      |      |                |      |      |
    !            p1 ---- p2 ---- p3               p7 ---- p8 ---- p9
    !
    ! The shared face is:
    !   element 1 face x = 1 <-> element 2 face x = 1

    call p(1)%init(0d0, 0d0, 0d0)
    call p(1)%set_id(1)
    call p(2)%init(1d0, 0d0, 0d0)
    call p(2)%set_id(2)
    call p(3)%init(2d0, 0d0, 0d0)
    call p(3)%set_id(3)
    call p(4)%init(0d0, 1d0, 0d0)
    call p(4)%set_id(4)
    call p(5)%init(1d0, 1d0, 0d0)
    call p(5)%set_id(5)
    call p(6)%init(2d0, 1d0, 0d0)
    call p(6)%set_id(6)
    call p(7)%init(0d0, 0d0, 1d0)
    call p(7)%set_id(7)
    call p(8)%init(1d0, 0d0, 1d0)
    call p(8)%set_id(8)
    call p(9)%init(2d0, 0d0, 1d0)
    call p(9)%set_id(9)
    call p(10)%init(0d0, 1d0, 1d0)
    call p(10)%set_id(10)
    call p(11)%init(1d0, 1d0, 1d0)
    call p(11)%set_id(11)
    call p(12)%init(2d0, 1d0, 1d0)
    call p(12)%set_id(12)

    call msh%init(3, 2)
    call msh%add_element(1, 1, p(1), p(2), p(4), p(5), &
         p(7), p(8), p(10), p(11))
    call msh%add_element(2, 2, p(2), p(3), p(5), p(6), &
         p(8), p(9), p(11), p(12))
    call msh%generate_conn()
  end subroutine two_adjacent_unit_hex_mesh

  subroutine three_hex_right_triangle_tip_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    type(point_t) :: p(18)

    ! Three unit hexes extruded in z from the following 2D layout, with the
    ! points listed in the NEKTON symmetric ordering of src/mesh/hex.f90,
    ! i.e. x fastest, then y, then z:
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
    call msh%add_element(1, 1, p(2), p(3), p(5), p(6), &
         p(11), p(12), p(14), p(15))

    ! Element adjacent to the y-axis cathetus.
    call msh%add_element(2, 2, p(4), p(5), p(7), p(8), &
         p(13), p(14), p(16), p(17))

    ! Connector element touching the triangle only at the tip.
    call msh%add_element(3, 3, p(1), p(2), p(4), p(5), &
         p(10), p(11), p(13), p(14))

    call msh%generate_conn()
  end subroutine three_hex_right_triangle_tip_mesh

  !> Column of `n` unit hexes stacked in z on [0, 1]^2 x [z0, z0 + n], where
  !! `z0` is the optional `z_offset` and defaults to 0.
  !!
  !! Element `g` in global numbering occupies z in [z0 + g - 1, z0 + g]. The
  !! elements
  !! are distributed linearly over the ranks of NEKO_COMM, so the fixture can
  !! be used by tests that run on several ranks. The point ids are global,
  !! which is what generate_conn() uses to find the faces shared between
  !! elements on different ranks. The points are listed in the NEKTON
  !! symmetric ordering of src/mesh/hex.f90, i.e. x fastest, then y, then z.
  !!
  !!   z = z0 + g:  p7 ----- p8      z = z0 + g - 1:  p3 ----- p4
  !!                 |       |                         |       |
  !!                 |  e_g  |                         |  e_g  |
  !!                 |       |                         |       |
  !!                p5 ----- p6                       p1 ----- p2
  subroutine stacked_unit_hex_column_mesh(msh, n, z_offset)
    type(mesh_t), intent(inout) :: msh
    integer, intent(in) :: n
    real(kind=dp), intent(in), optional :: z_offset
    type(linear_dist_t) :: dist
    type(point_t) :: p(8)
    integer :: e, g
    real(kind=dp) :: z0

    dist = linear_dist_t(n, pe_rank, pe_size, NEKO_COMM)
    call msh%init(3, dist)

    do e = 1, dist%num_local()
       g = dist%start_idx() + e
       z0 = real(g - 1, dp)
       if (present(z_offset)) z0 = z0 + z_offset

       ! The four points of layer g - 1 have ids 4 * (g - 1) + 1 .. 4 * g
       ! and the four points of layer g have ids 4 * g + 1 .. 4 * (g + 1),
       ! so consecutive elements share the ids of the layer between them.
       call p(1)%init(0d0, 0d0, z0)
       call p(1)%set_id(4 * (g - 1) + 1)
       call p(2)%init(1d0, 0d0, z0)
       call p(2)%set_id(4 * (g - 1) + 2)
       call p(3)%init(0d0, 1d0, z0)
       call p(3)%set_id(4 * (g - 1) + 3)
       call p(4)%init(1d0, 1d0, z0)
       call p(4)%set_id(4 * (g - 1) + 4)
       call p(5)%init(0d0, 0d0, z0 + 1d0)
       call p(5)%set_id(4 * g + 1)
       call p(6)%init(1d0, 0d0, z0 + 1d0)
       call p(6)%set_id(4 * g + 2)
       call p(7)%init(0d0, 1d0, z0 + 1d0)
       call p(7)%set_id(4 * g + 3)
       call p(8)%init(1d0, 1d0, z0 + 1d0)
       call p(8)%set_id(4 * g + 4)

       call msh%add_element(e, g, p(1), p(2), p(3), p(4), &
            p(5), p(6), p(7), p(8))
    end do
    call msh%generate_conn()
  end subroutine stacked_unit_hex_column_mesh

end module mesh_fixture
