! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> @brief Module containing Signed Distance Functions.
module signed_distance
  use num_types, only: rp
  use field, only: field_t
  use tri, only: tri_t
  use tri_mesh, only: tri_mesh_t
  use aabb_tree, only: aabb_tree_t

  implicit none
  private
  public :: signed_distance_field

contains

  !> @brief Signed distance field
  !! @details This routine computes the signed distance field for a given
  !! object.
  !!
  !! Currently supported objects are:
  !! - Triangular mesh (tri_mesh_t)
  !!
  !! @param[in,out] field_data Field data
  !! @param[in] object Object
  !! @param[in,optional] max_distance Maximum distance outside the mesh
  subroutine signed_distance_field(field_data, object, max_distance)
    USE utils, ONLY: neko_error
    implicit none

    type(field_t), intent(inout) :: field_data
    class(*), intent(in) :: object
    real(kind=rp), intent(in), optional :: max_distance

    real(kind=rp) :: max_dist

    if (present(max_distance)) then
       max_dist = max_distance
    else
       max_dist = huge(0.0_rp)
    end if

    select type(object)
      type is (tri_mesh_t)
       call signed_distance_field_tri_mesh(field_data, object, max_dist)

      class default
       call neko_error("signed_distance_field: Object type not supported.")
    end select

  end subroutine signed_distance_field

  !> @brief Signed distance field for a triangular mesh
  !! @details This routine computes the signed distance field for a given
  !! triangular mesh. The algorithm utilizes an AABB tree to accelerate the
  !! earch for potential elements. The signed distance is computed using the
  !! brute force approach, where we compute the signed distance to each element
  !! found through the AABB tree, and return the minimum distance.
  !!
  !! @param[in,out] field_data Field data
  !! @param[in] mesh Triangular mesh
  !! @param[in] max_distance Maximum distance outside the mesh
  subroutine signed_distance_field_tri_mesh(field_data, mesh, max_distance)
    type(field_t), intent(inout) :: field_data
    type(tri_mesh_t), intent(in) :: mesh
    real(kind=rp), intent(in) :: max_distance

    integer :: total_size
    integer :: id, el_id
    type(aabb_tree_t) :: search_tree
    real(kind=rp), dimension(3) :: p

    ! Zero the field
    field_data%x = 0.0_rp
    total_size = field_data%dof%size()

    print *, "Building the search tree."

    ! Build the search tree
    call search_tree%init(mesh%nelv)
    do el_id = 1, mesh%nelv
       call search_tree%insert_object(mesh%el(el_id), el_id)
    end do

    print *, "Computing the signed distance field."

    do id = 1, total_size
       p(1) = field_data%dof%x(id, 1, 1, 1)
       p(2) = field_data%dof%y(id, 1, 1, 1)
       p(3) = field_data%dof%z(id, 1, 1, 1)

       field_data%x(id, 1, 1, 1) = tri_mesh_aabb_tree(search_tree, mesh%el, p, max_distance)
       !  field_data%x(id, 1, 1, 1) = tri_mesh_brute_force(mesh, p)
    end do

  end subroutine signed_distance_field_tri_mesh

  !> @brief Signed distance function
  !! @deprecated This routine is deprecated and will be removed in the future.
  !! @details This routine computes the signed distance function for the
  !! boundary mesh, to a given point (x, y, z). The algorithm is a
  !! brute force approach, where we compute the signed distance to each
  !! element in the mesh, and return the minimum distance.
  !!
  !! @param p Point
  !! @param mesh Boundary mesh
  !! @return Signed distance value
  function tri_mesh_brute_force(mesh, p, max_distance) result(distance)
    use tri, only: tri_t
    use point, only: point_t
    use num_types, only: rp

    implicit none

    type(tri_mesh_t), intent(in) :: mesh
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(in) :: max_distance

    integer :: id
    real(kind=rp) :: distance, weighted_sign
    real(kind=rp) :: cd, cs
    real(kind=rp) :: tol = 1e-6_rp

    distance = 1e10_rp
    weighted_sign = 0.0_rp

    do id = 1, mesh%nelv
       call element_distance(mesh%el(id), p, cd, cs)

       ! Update the weighted sign, if the relative difference is small
       if (abs(cd - distance) / distance .lt. tol) then
          weighted_sign = weighted_sign + cs
       else if (cd .lt. distance) then
          weighted_sign = cs
       end if

       distance = min(cd, distance)
    end do

    distance = sign(distance, weighted_sign)

  end function tri_mesh_brute_force

  !> @brief Signed distance function using an AABB tree
  !! @details This routine computes the signed distance function for the
  !! boundary mesh, to a given point (x, y, z). The algorithm utilizes an
  !! AABB tree to accelerate the search for potential elements. The signed
  !! distance is computed using the brute force approach, where we compute the
  !! signed distance to each element found through the AABB tree, and return
  !! the minimum distance.
  !!
  !! @param tree AABB tree
  !! @param object_list List of objects
  !! @param p Point
  !! @param max_distance Maximum distance outside the mesh
  !! @return Signed distance value
  function tri_mesh_aabb_tree(tree, object_list, p, max_distance) result(distance)
    use aabb, only: aabb_t
    use aabb_tree, only: aabb_node_t
    use stack, only: stack_i4_t
    implicit none

    class(aabb_tree_t), intent(in) :: tree
    class(tri_t), dimension(:), intent(in) :: object_list
    real(kind=rp), dimension(3), intent(in) :: p
    real(kind=rp), intent(in) :: max_distance

    real(kind=rp) :: distance
    real(kind=rp) :: weighted_sign

    real(kind=rp), parameter :: tol = 1.0e-6_rp

    type(stack_i4_t) :: simple_stack
    integer :: current_index

    type(aabb_node_t) :: current_node
    type(aabb_t) :: current_aabb
    integer :: current_object_index
    real(kind=rp) :: current_distance
    real(kind=rp) :: current_sign

    type(aabb_node_t) :: left_node
    type(aabb_node_t) :: right_node

    type(aabb_t) :: search_box

    integer :: root_index, left_index, right_index

    root_index = tree%get_root_index()

    call simple_stack%init(size(object_list) * 2)
    call simple_stack%push(root_index)

    ! Assign the initial distance as the distance to the center a random object
    ! in the list, with the diagonal of the object added to the distance.
    current_object_index = floor(rand() * size(object_list) + 1)
    call element_distance(object_list(current_object_index), p, distance, weighted_sign)
    distance = distance + object_list(current_object_index)%diameter()

    call search_box%init(p - min(distance, max_distance), &
                         p + min(distance, max_distance))

    do while (.not. simple_stack%is_empty())
       current_index = simple_stack%pop()

       current_node = tree%get_node(current_index)
       current_aabb = current_node%get_aabb()

       if (current_node%is_leaf()) then
          if (distance .lt. current_aabb%min_distance(p)) then
             cycle
          end if

          current_object_index = current_node%get_object_index()
          call element_distance(object_list(current_object_index), p, current_distance, current_sign)

          ! Update the weighted sign, if the relative difference is small
          if (abs(current_distance - distance) / distance .lt. tol) then
             weighted_sign = weighted_sign + current_sign
          else if (current_distance .lt. distance) then
             weighted_sign = current_sign
          end if

          distance = min(distance, current_distance)

          ! Update the search box to the new distance
          if (distance .gt. current_aabb%get_diagonal()) then
             call search_box%init(p - distance, p + distance)
          end if

       else if (current_aabb%overlaps(search_box)) then

          left_node = tree%get_left_node(current_index)
          right_node = tree%get_right_node(current_index)

          if (left_node%aabb%overlaps(search_box)) then
             left_index = tree%get_left_index(current_index)
             call simple_stack%push(left_index)
          end if

          if (right_node%aabb%overlaps(search_box)) then
             right_index = tree%get_right_index(current_index)
             call simple_stack%push(right_index)
          end if

       end if

    end do

    if (distance .gt. max_distance) then
       distance = max_distance
    else
       distance = sign(distance, weighted_sign)
    end if

  end function tri_mesh_aabb_tree

  ! ========================================================================== !
  ! Element distance functions
  ! ========================================================================== !

  !> @brief Main interface for the signed distance function for an element.
  subroutine element_distance(element, p, distance, weighted_sign)
    class(*), intent(in) :: element
    real(kind=rp), dimension(3), intent(in) :: p
    real(kind=rp), intent(out) :: distance
    real(kind=rp), intent(out), optional :: weighted_sign

    select type(element)
      type is (tri_t)
       call element_distance_triangle(element, p, distance, weighted_sign)
      class default
       print *, "Error: Element type not supported."
       stop
    end select
  end subroutine element_distance

  ! -------------------------------------------------------------------------- !
  ! Specific element distance functions

  !> @brief Signed distance function for a triangle
  !! @details This routine computes the signed distance function for the current
  !! triangle. We compute the barycentric coordinate of the point projected onto
  !! the triangle. If the projection is inside the triangle, the distance is the
  !! distance to the plane. If the projection is outside the triangle, the
  !! distance is the distance to the nearest edge or vertex.
  !!
  !! In order to improve precision of the sign estimation, we also compute the
  !! weighted sign, which is the perpendicular distance to the plane divided by
  !! the distance to the nearest point.
  !!
  !! @Note The returned distance is signed if the weighted_sign is not present.
  !!
  !! @param this Triangle
  !! @param p Point
  !! @return Distance value
  !! @return[optional] Weighted sign
  subroutine element_distance_triangle(triangle, p, distance, weighted_sign)

    type(tri_t), intent(in) :: triangle
    real(kind=rp), dimension(3), intent(in) :: p

    real(kind=rp), intent(out) :: distance
    real(kind=rp), intent(out), optional :: weighted_sign

    real(kind=rp), dimension(3) :: v1, v2, v3
    real(kind=rp), dimension(3) :: normal
    real(kind=rp) :: normal_length
    real(kind=rp) :: b1, b2, b3

    real(kind=rp), dimension(3) :: projection
    real(kind=rp), dimension(3) :: edge
    real(kind=rp) :: tol = 1e-10_rp

    real(kind=rp) :: face_distance
    real(kind=rp) :: t

    ! Get vertices and the normal vector
    v1 = triangle%pts(1)%p%x
    v2 = triangle%pts(2)%p%x
    v3 = triangle%pts(3)%p%x

    normal = cross(v2 - v1, v3 - v1)
    normal_length = norm2(normal)
    if (normal_length .gt. tol)then
       normal = normal / normal_length
    end if

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! triangular prism, along an edge or by a face.
    face_distance = dot_product(p - v1, normal)

    projection = p - normal * face_distance
    b1 = dot_product(normal, cross(v2 - v1, projection - v1)) / normal_length
    b2 = dot_product(normal, cross(v3 - v2, projection - v2)) / normal_length
    b3 = 1.0_rp - b1 - b2

    if (b1 .le. tol) then
       edge = v2 - v1
       t = dot_product(p - v1, edge) / norm2(edge)
       t = max(0.0_rp, min(1.0_rp, t))

       projection = v1 + t * edge
    else if (b2 .le. tol) then
       edge = v3 - v2
       t = dot_product(p - v2, edge) / norm2(edge)
       t = max(0.0_rp, min(1.0_rp, t))

       projection = v2 + t * edge
    else if (b3 .le. tol) then
       edge = v1 - v3
       t = dot_product(p - v3, edge) / norm2(edge)
       t = max(0.0_rp, min(1.0_rp, t))

       projection = v3 + t * edge
    end if

    if (present(weighted_sign)) then
       distance = norm2(projection - p)
       weighted_sign = face_distance / distance
    else
       distance = sign(norm2(projection - p), face_distance)
    end if

  end subroutine element_distance_triangle

  !> Compute cross product of two vectors
  !> @param[in] a First vector
  !> @param[in] b Second vector
  !> @return Cross product \f$ a \times b \f$
  pure function cross(a, b) result(c)
    use num_types, only: rp
    implicit none

    real(kind=rp), dimension(3), intent(in) :: a
    real(kind=rp), dimension(3), intent(in) :: b
    real(kind=rp), dimension(3) :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

end module signed_distance
