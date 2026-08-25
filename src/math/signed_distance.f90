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
  use num_types, only : dp, rp
  use field, only : field_t
  use tri, only : tri_t
  use tri_mesh, only : tri_mesh_t
  use aabb, only : aabb_t
  use aabb_tree, only : aabb_tree_t, aabb_node_t, AABB_NULL_NODE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error, neko_warning
  use math, only : rzero
  !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num

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
    type(field_t), intent(inout) :: field_data
    class(*), intent(in) :: object
    real(kind=dp), intent(in), optional :: max_distance

    real(kind=dp) :: max_dist

    if (present(max_distance)) then
       max_dist = max_distance
    else
       max_dist = huge(0.0_dp)
    end if

    select type (object)
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
    real(kind=dp), intent(in) :: max_distance

    integer :: total_size
    integer :: id
    type(aabb_tree_t) :: search_tree
    real(kind=dp), dimension(3) :: p
    real(kind=dp) :: distance
    integer :: thread_id, n_threads, stack_size
    integer, allocatable :: stack(:,:)

    ! Zero the field
    call rzero(field_data%x, field_data%size())
    total_size = field_data%dof%size()

    call search_tree%init(mesh%nelv)
    call search_tree%build(mesh%el)

    if (search_tree%get_size() .ne. mesh%nelv) then
       call neko_error("signed_distance_field_tri_mesh: " // &
            "Error building the search tree.")
    end if

    n_threads = 1
    !$ n_threads = omp_get_max_threads()
    stack_size = 2 * mesh%nelv
    allocate(stack(stack_size, n_threads))

    !$omp parallel private(id, p, distance, thread_id)
    thread_id = 1
    !$ thread_id = omp_get_thread_num() + 1
    !$omp do
    do id = 1, total_size
       p(1) = field_data%dof%x(id, 1, 1, 1)
       p(2) = field_data%dof%y(id, 1, 1, 1)
       p(3) = field_data%dof%z(id, 1, 1, 1)

       distance = tri_mesh_aabb_tree(search_tree, mesh%el, p, &
            max_distance, stack(:, thread_id))

       field_data%x(id, 1, 1, 1) = real(distance, kind=rp)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(stack)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_warning("signed_distance_field_tri_mesh:&
       & Device version not implemented.")
       call device_memcpy(field_data%x, field_data%x_d, field_data%size(), &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine signed_distance_field_tri_mesh

  !> @brief Signed distance function
  !! @details This routine computes the signed distance function for the
  !! boundary mesh, to a given point (x, y, z). The algorithm is a
  !! brute force approach, where we compute the signed distance to each
  !! element in the mesh, and return the minimum distance.
  !!
  !! @note This routine is here for completeness, and is not optimized for
  !! performance. It is recommended to use the AABB tree version for large meshes.
  !!
  !! @param p Point
  !! @param mesh Boundary mesh
  !! @return Signed distance value
  function tri_mesh_brute_force(mesh, p, max_distance) result(distance)
    type(tri_mesh_t), intent(in) :: mesh
    real(kind=dp), intent(in) :: p(3)
    real(kind=dp), intent(in) :: max_distance

    integer :: id
    real(kind=dp) :: distance, weighted_sign
    real(kind=dp) :: cd, cs
    real(kind=dp) :: tol = 1e-6_dp

    distance = huge(0.0_dp)
    weighted_sign = 0.0_dp

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

    distance = sign(min(distance, max_distance), weighted_sign)

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
  function tri_mesh_aabb_tree(tree, object_list, p, max_distance, stack) &
       result(distance)
    class(aabb_tree_t), intent(in) :: tree
    class(tri_t), contiguous, dimension(:), intent(in) :: object_list
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), intent(in) :: max_distance
    integer, dimension(:), intent(inout) :: stack

    real(kind=dp) :: distance
    real(kind=dp) :: weighted_sign

    real(kind=dp), parameter :: tol = 1.0e-6_dp

    integer :: current_index, stack_top

    type(aabb_node_t) :: current_node
    integer :: current_object_index
    real(kind=dp) :: current_distance
    real(kind=dp) :: current_sign

    type(aabb_node_t) :: left_node
    type(aabb_node_t) :: right_node

    type(aabb_t) :: root_box
    type(aabb_t) :: search_box

    integer :: root_index, left_index, right_index

    ! Initialize the stack and the search box.
    stack_top = 0
    call search_box%init(p - max_distance, p + max_distance)

    ! Check if the root node overlaps the search box.
    root_index = tree%get_root_index()
    root_box = tree%get_aabb(root_index)

    if (.not. root_box%overlaps(search_box)) then
       distance = max_distance
       return
    end if

    distance = max_distance
    weighted_sign = 1.0_dp
    call push_node(stack, stack_top, root_index)

    ! Traverse the tree and compute the signed distance to the elements
    do while (stack_top .gt. 0)
       current_index = stack(stack_top)
       stack_top = stack_top - 1
       if (current_index .eq. AABB_NULL_NODE) cycle

       current_node = tree%get_node(current_index)

       if (aabb_distance(current_node%get_aabb(), p) .gt. distance) then
          cycle
       end if

       if (current_node%is_leaf()) then

          current_object_index = current_node%get_object_index()
          call element_distance(object_list(current_object_index), p, &
               current_distance, current_sign)

          ! Update the weighted sign, if the relative difference is small
          if (current_distance .lt. distance) then
             weighted_sign = current_sign
          else if (distance .le. tol) then
             if (current_distance .le. tol) then
                weighted_sign = weighted_sign + current_sign
             end if
          else if (abs(current_distance - distance) / distance .lt. tol) then
             weighted_sign = weighted_sign + current_sign
          end if

          distance = min(distance, current_distance)

          ! Update the search box to the new distance
          call search_box%init(p - distance, p + distance)
       else

          left_index = tree%get_left_index(current_index)
          if (left_index .ne. AABB_NULL_NODE) then
             left_node = tree%get_left_node(current_index)
             if (left_node%aabb%overlaps(search_box)) then
                call push_node(stack, stack_top, left_index)
             end if
          end if

          right_index = tree%get_right_index(current_index)
          if (right_index .ne. AABB_NULL_NODE) then
             right_node = tree%get_right_node(current_index)
             if (right_node%aabb%overlaps(search_box)) then
                call push_node(stack, stack_top, right_index)
             end if
          end if
       end if
    end do

    if (distance .gt. max_distance) then
       distance = max_distance
    end if
    distance = sign(distance, weighted_sign)

  end function tri_mesh_aabb_tree

  !> Push a node index onto a caller-provided stack.
  subroutine push_node(stack, top, node_index)
    integer, dimension(:), intent(inout) :: stack
    integer, intent(inout) :: top
    integer, intent(in) :: node_index

    if (top .ge. size(stack)) then
       call neko_error("signed_distance_field_tri_mesh: AABB stack overflow.")
    end if
    top = top + 1
    stack(top) = node_index

  end subroutine push_node

  !> Euclidean distance from a point to an axis-aligned bounding box.
  function aabb_distance(box, p) result(distance)
    type(aabb_t), intent(in) :: box
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp) :: distance
    real(kind=dp), dimension(3) :: box_min, box_max, delta

    box_min = box%get_min()
    box_max = box%get_max()
    delta = max(max(box_min - p, 0.0_dp), p - box_max)
    distance = norm2(delta)

  end function aabb_distance

  ! ========================================================================== !
  ! Element distance functions
  ! ========================================================================== !

  !> @brief Main interface for the signed distance function for an element.
  subroutine element_distance(element, p, distance, weighted_sign)
    class(*), intent(in) :: element
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp), intent(out) :: distance
    real(kind=dp), intent(out), optional :: weighted_sign

    select type (element)
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
    real(kind=dp), dimension(3), intent(in) :: p

    real(kind=dp), intent(out) :: distance
    real(kind=dp), intent(out), optional :: weighted_sign

    real(kind=dp), dimension(3) :: v1, v2, v3
    real(kind=dp), dimension(3) :: normal
    real(kind=dp) :: normal_length
    real(kind=dp) :: b1, b2, b3

    real(kind=dp), dimension(3) :: projection
    real(kind=dp), dimension(3) :: edge
    real(kind=dp) :: tol = 1e-10_dp

    real(kind=dp) :: face_distance
    real(kind=dp) :: t

    ! Get vertices and the normal vector
    v1 = triangle%pts(1)%p%x
    v2 = triangle%pts(2)%p%x
    v3 = triangle%pts(3)%p%x

    normal = cross(v2 - v1, v3 - v1)
    normal_length = norm2(normal)

    if (normal_length .lt. tol) then
       distance = huge(1.0_dp)
       weighted_sign = 0.0_dp
       return
    end if
    normal = normal / normal_length

    ! Compute Barycentric coordinates to determine if the point is inside the
    ! triangular prism, of along an edge or by a face.
    face_distance = dot_product(p - v1, normal)

    projection = p - normal * face_distance
    b1 = dot_product(normal, cross(v2 - v1, projection - v1)) / normal_length
    b2 = dot_product(normal, cross(v3 - v2, projection - v2)) / normal_length
    b3 = dot_product(normal, cross(v1 - v3, projection - v3)) / normal_length

    if (b1 .le. tol) then
       edge = v2 - v1
       t = dot_product(p - v1, edge) / norm2(edge)
       t = max(0.0_dp, min(1.0_dp, t))

       projection = v1 + t * edge
    else if (b2 .le. tol) then
       edge = v3 - v2
       t = dot_product(p - v2, edge) / norm2(edge)
       t = max(0.0_dp, min(1.0_dp, t))

       projection = v2 + t * edge
    else if (b3 .le. tol) then
       edge = v1 - v3
       t = dot_product(p - v3, edge) / norm2(edge)
       t = max(0.0_dp, min(1.0_dp, t))

       projection = v3 + t * edge
    end if

    distance = norm2(projection - p)
    if (present(weighted_sign)) then
       weighted_sign = face_distance / distance
    else
       distance = sign(distance, face_distance)
    end if

  end subroutine element_distance_triangle

  !> Compute cross product of two vectors
  !> @param[in] a First vector
  !> @param[in] b Second vector
  !> @return Cross product \f$ a \times b \f$
  pure function cross(a, b) result(c)
    real(kind=dp), dimension(3), intent(in) :: a
    real(kind=dp), dimension(3), intent(in) :: b
    real(kind=dp), dimension(3) :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

end module signed_distance
