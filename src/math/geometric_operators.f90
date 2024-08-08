! Lets have a library for the geometric operators, such as signed distance and
! boolean containment tests. This will be useful across a wide part of the code.

! Further more lets try to set it up in a proper way with submodules and module
! operators.

!> \module geometric_operators
!! This module contains geometric operators that are useful for the mesh
!! generation and manipulation. The operators are implemented in a way that
!! they can be used in a modular way, and can be extended with new operators
!! easily.
!!
!! Please note that this module only contains the interfaces, the actual
!! implementation is located in approproate submodules. The submodules are
!! located in the backend directories.
module geometric_operators
  use num_types, only: dp
  use point, only: point_t
  use tri, only: tri_t
  use tet, only: tet_t
  use quad, only: quad_t
  use hex, only: hex_t
  implicit none

  ! ========================================================================== !
  ! Distance for elements from the mesh structure.

  !> \interface distance_point
  !! This interface defines the distance operator for a point.
  !! The distance operator calculates the distance between a query point and a
  !! given point.
  !!
  !! @param p Query point
  !! @param point Point
  !! @return Distance value
  interface distance_point
     module procedure distance_point_real
     module procedure distance_point_t
  end interface distance_point

  ! ========================================================================== !
  ! Distance to geometric primitives.

  interface

     !> Distance to a point
     module function distance_point_real(p, point)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: point
       real(kind=dp) :: distance_point_real
     end function distance_point_real

     !> Distance to the infinite line defined by a point and a direction.
     module function distance_line(p, point, direction)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: point
       real(kind=dp), dimension(3), intent(in) :: direction
       real(kind=dp) :: distance_line
     end function distance_line

     !> Distance to the infinite ray starting at a point and going in a
     !! direction.
     module function distance_line_ray(p, point, direction)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: point
       real(kind=dp), dimension(3), intent(in) :: direction
       real(kind=dp) :: distance_line_ray
     end function distance_line_ray

     !> Distance to a line segment defined by two points.
     module function distance_line_segment(p, point_0, point_1)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: point_0
       real(kind=dp), dimension(3), intent(in) :: point_1
       real(kind=dp) :: distance_line_segment
     end function distance_line_segment

     !> Distance to a plane defined by a point and a normal.
     module function distance_plane(p, point, normal)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: point
       real(kind=dp), dimension(3), intent(in) :: normal
       real(kind=dp) :: distance_plane
     end function distance_plane

     !> Distance to a sphere defined by a center and a radius.
     module function distance_sphere(p, sphere_center, sphere_radius)
       real(kind=dp), dimension(3), intent(in) :: p
       real(kind=dp), dimension(3), intent(in) :: sphere_center
       real(kind=dp), intent(in) :: sphere_radius
       real(kind=dp) :: distance_sphere
     end function distance_sphere

  end interface

  ! ========================================================================== !
  ! Distance to elements.

  !> \interface distance_elements
  !! This interface defines the distance operator. The distance
  !! operator is used to calculate the distance between a point and a geometric
  !! object.
  !! The distance is signed for volumetric objects, meaning that the distance is
  !! positive if the point is outside the object, and negative if the point is
  !! inside the object.
  !! For non-volumetric objects, the distance is always positive.
  interface distance_element

     !> Distance to a point
     module function distance_point_t(p, point)
       real(kind=dp), dimension(3), intent(in) :: p
       type(point_t), intent(in) :: point
       real(kind=dp) :: distance_point_t
     end function distance_point_t

     !> Distance to a triangle
     module function distance_triangle(p, triangle)
       real(kind=dp), dimension(3), intent(in) :: p
       type(tri_t), intent(in) :: triangle
       real(kind=dp) :: distance_triangle
     end function distance_triangle

     !> Distance to tetrahedron
     module function distance_tetrahedron(p, tetrahedron)
       real(kind=dp), dimension(3), intent(in) :: p
       type(tet_t), intent(in) :: tetrahedron
       real(kind=dp) :: distance_tetrahedron
     end function distance_tetrahedron

  end interface distance_element

  ! ========================================================================== !
  ! Local coordinate operators.

  !> \interface linear coordinates
  !! This interface defines the linear coordinates operator. The linear
  !! coordinates operator is used to calculate the linear coordinates of a point
  !! with respect to a geometric object. The linear coordinates are bilinear for
  !! a quad, trilinear for a hexahedron and baricentric for a triangle and
  !! tetrahedron.
  interface linear_coordinate

     !> Barycentric coordinates for a point and a triangle
     module function barycentric_coordinate_triangle(p, triangle)
       real(kind=dp), dimension(3), intent(in) :: p
       type(tri_t), intent(in) :: triangle
       real(kind=dp), dimension(3) :: barycentric_coordinate_triangle
     end function barycentric_coordinate_triangle

     !> Barycentric coordinates for a point and a tetrahedron
     module function barycentric_coordinate_tetrahedron(p, tetrahedron)
       real(kind=dp), dimension(3), intent(in) :: p
       type(tet_t), intent(in) :: tetrahedron
       real(kind=dp), dimension(4) :: barycentric_coordinate_tetrahedron
     end function barycentric_coordinate_tetrahedron

     !> Linear coordinates for a point and a quadrilateral.
     module function bilinear_coordinate_quadrilateral(p, quadrilateral)
       real(kind=dp), dimension(3), intent(in) :: p
       type(quad_t), intent(in) :: quadrilateral
       real(kind=dp), dimension(2) :: bilinear_coordinate_quadrilateral
     end function bilinear_coordinate_quadrilateral

     !> Linear coordinates for a point and a hexahedron.
     module function trilinear_coordinate_hexahedron(p, hexahedron)
       real(kind=dp), dimension(3), intent(in) :: p
       type(hex_t), intent(in) :: hexahedron
       real(kind=dp), dimension(3) :: trilinear_coordinate_hexahedron
     end function trilinear_coordinate_hexahedron

  end interface linear_coordinate

  !> \interface euclidean coordinates
  !! This interface defines the euclidean coordinates operator. The euclidean
  !! coordinates operator is used to calculate the euclidean coordinates of a
  !! point with respect to a geometric object based on a given set of linear
  !! coordinates.
  interface euclidean_coordinate

     !> Map barycentric coordinates to euclidean coordinates based on a given
     !! triangle.
     module function euclidean_coordinate_triangle(barycoord, triangle)
       real(kind=dp), dimension(3), intent(in) :: barycoord
       type(tri_t), intent(in) :: triangle
       real(kind=dp), dimension(3) :: euclidean_coordinate_triangle
     end function euclidean_coordinate_triangle

     !> Map barycentric coordinates to euclidean coordinates based on a given
     !! tetrahedron.
     module function euclidean_coordinate_tetrahedron(barycoord, tetrahedron)
       real(kind=dp), dimension(4), intent(in) :: barycoord
       type(tet_t), intent(in) :: tetrahedron
       real(kind=dp), dimension(3) :: euclidean_coordinate_tetrahedron
     end function euclidean_coordinate_tetrahedron

     !> Map linear coordinates to euclidean coordinates based on a given
     !! quadrilateral.
     module function euclidean_coordinate_quadrilateral(linearcoord, quadrilateral)
       real(kind=dp), dimension(2), intent(in) :: linearcoord
       type(quad_t), intent(in) :: quadrilateral
       real(kind=dp), dimension(3) :: euclidean_coordinate_quadrilateral
     end function euclidean_coordinate_quadrilateral

     !> Map linear coordinates to euclidean coordinates based on a given
     !! hexahedron.
     module function euclidean_coordinate_hexahedron(linearcoord, hexahedron)
       real(kind=dp), dimension(3), intent(in) :: linearcoord
       type(hex_t), intent(in) :: hexahedron
       real(kind=dp), dimension(3) :: euclidean_coordinate_hexahedron
     end function euclidean_coordinate_hexahedron

  end interface euclidean_coordinate

end module geometric_operators
