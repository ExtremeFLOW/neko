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
! ============================================================================ !
! Original C++ Implementation from:
! https://github.com/JamesRandall/SimpleVoxelEngine/blob/master/voxelEngine/src/AABBTree.h
!
! Translated to Fortran by:
! @author Tim Felle Olsen
! @date 9 Feb 2024
!
! C++ Code License:
! The MIT License (MIT)
!
! Copyright (c) 2017 James Randall
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of
! this software and associated documentation files (the "Software"), to deal in
! the Software without restriction, including without limitation the rights to
! use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
! the Software, and to permit persons to whom the Software is furnished to do so,
! subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
! FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
! COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! ============================================================================ !

!> @brief Axis Aligned Bounding Box (aabb) Tree data structure.
!! @details
!! This is a Fortran implementation of an Axis Aligned Bounding Box Tree
!! data structure.
!! The purpose of this is to accelerate a Signed Distance Function and other
!! spatial computations.
module aabb_tree
  use aabb
  use tri, only: tri_t
  use num_types, only: rp, dp

  implicit none
  private

  integer, parameter, public :: AABB_NULL_NODE = -1

  ! ========================================================================== !
  ! Type definitions
  ! ========================================================================== !

  !> @brief Node type for the Axis Aligned Bounding Box (aabb) Tree
  type, public :: aabb_node_t
     private
     type(aabb_t), public :: aabb
     integer :: object_index = -1

     ! tree links
     integer :: parent_node_index = AABB_NULL_NODE
     integer :: left_node_index = AABB_NULL_NODE
     integer :: right_node_index = AABB_NULL_NODE

     ! node linked list link
     integer :: next_node_index = AABB_NULL_NODE

   contains
     procedure, pass(this), public :: init => aabb_node_init

     ! Getters
     procedure, pass(this), public :: get_aabb => aabb_node_get_aabb
     procedure, pass(this), public :: get_object_index => &
          aabb_node_get_object_index
     procedure, pass(this), public :: get_parent_index => &
          aabb_node_get_parent_index
     procedure, pass(this), public :: get_left_index => &
          aabb_node_get_left_index
     procedure, pass(this), public :: get_right_index => &
          aabb_node_get_right_index

     ! Unary operations
     procedure, pass(this), public :: min_distance => aabb_node_min_distance

     ! Boolean operators
     procedure, pass(this), public :: is_leaf => aabb_node_is_leaf
     procedure, pass(this), public :: is_valid => aabb_node_is_valid

     ! Comparison operators
     generic :: operator(.lt.) => less
     generic :: operator(.gt.) => greater

     procedure, pass(this) :: less => aabb_node_less
     procedure, pass(this) :: greater => aabb_node_greater

  end type aabb_node_t

  !> @brief Axis Aligned Bounding Box (aabb) Tree
  type, public :: aabb_tree_t
     private
     type(aabb_node_t), allocatable :: nodes(:)
     integer :: root_node_index = AABB_NULL_NODE
     integer :: allocated_node_count = 0
     integer :: next_free_node_index = AABB_NULL_NODE
     integer :: node_capacity = 0
     integer :: growth_size = 1

   contains

     ! Initializers
     procedure, pass(this), public :: init => aabb_tree_init
     procedure, pass(this), public :: build => aabb_tree_build_tree
     procedure, pass(this), public :: insert_object => &
          aabb_tree_insert_object

     ! Getters
     procedure, pass(this), public :: get_size => aabb_tree_get_size

     procedure, pass(this), public :: get_root_index => &
          aabb_tree_get_root_index
     procedure, pass(this), public :: get_parent_index => &
          aabb_tree_get_parent_index
     procedure, pass(this), public :: get_left_index => &
          aabb_tree_get_left_index
     procedure, pass(this), public :: get_right_index => &
          aabb_tree_get_right_index

     procedure, pass(this), public :: get_node => aabb_tree_get_node
     procedure, pass(this), public :: get_root_node => &
          aabb_tree_get_root_node
     procedure, pass(this), public :: get_parent_node => &
          aabb_tree_get_parent_node
     procedure, pass(this), public :: get_left_node => &
          aabb_tree_get_left_node
     procedure, pass(this), public :: get_right_node => &
          aabb_tree_get_right_node

     procedure, pass(this), public :: get_aabb => aabb_tree_get_aabb

     procedure, pass(this), public :: query_overlaps => &
          aabb_tree_query_overlaps

     procedure, pass(this), public :: print => aabb_tree_print

     ! ----------------------------------------------------------------------- !
     ! Internal methods

     procedure, pass(this) :: allocate_node => aabb_tree_allocate_node
     procedure, pass(this) :: deallocate_node => aabb_tree_deallocate_node
     procedure, pass(this) :: resize_node_pool => aabb_tree_resize_node_pool
     procedure, pass(this) :: insert_leaf => aabb_tree_insert_leaf

     procedure, pass(this) :: fix_upwards_tree => aabb_tree_fix_upwards_tree

     procedure, pass(this) :: valid_tree => aabb_tree_valid_tree

  end type aabb_tree_t

contains

  ! ========================================================================== !
  ! Definitions of node methods
  ! ========================================================================== !

  !> @brief Initializes the AABB node.
  subroutine aabb_node_init(this)
    class(aabb_node_t), intent(inout) :: this

    this%object_index = -1
    this%parent_node_index = AABB_NULL_NODE
    this%left_node_index = AABB_NULL_NODE
    this%right_node_index = AABB_NULL_NODE
    this%next_node_index = AABB_NULL_NODE
  end subroutine aabb_node_init

  ! -------------------------------------------------------------------------- !
  ! Getters

  !> @brief Returns the Axis Aligned Bounding Box (aabb) of the node.
  pure function aabb_node_get_aabb(this) result(res)
    class(aabb_node_t), intent(in) :: this
    type(aabb_t) :: res

    res = this%aabb
  end function aabb_node_get_aabb

  !> @brief Returns the object index of the node.
  pure function aabb_node_get_object_index(this) result(object_index)
    class(aabb_node_t), intent(in) :: this
    integer :: object_index

    object_index = this%object_index
  end function aabb_node_get_object_index

  !> @brief Returns the parent index of the node.
  pure function aabb_node_get_parent_index(this) result(parent_index)
    class(aabb_node_t), intent(in) :: this
    integer :: parent_index

    parent_index = this%parent_node_index
  end function aabb_node_get_parent_index

  !> @brief Returns the left index of the node.
  pure function aabb_node_get_left_index(this) result(left_index)
    class(aabb_node_t), intent(in) :: this
    integer :: left_index

    left_index = this%left_node_index
  end function aabb_node_get_left_index

  !> @brief Returns the right index of the node.
  pure function aabb_node_get_right_index(this) result(right_index)
    class(aabb_node_t), intent(in) :: this
    integer :: right_index

    right_index = this%right_node_index
  end function aabb_node_get_right_index

  !> @brief Get the minimum possible distance from the aabb to a point.
  function aabb_node_min_distance(this, p) result(distance)
    class(aabb_node_t), intent(in) :: this
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp) :: distance

    distance = 0.5_rp * this%aabb%get_diameter() &
      - norm2(this%aabb%get_center() - p)
  end function aabb_node_min_distance

  ! -------------------------------------------------------------------------- !
  ! Boolean operators

  !> @brief Returns true if the node is a leaf node.
  pure function aabb_node_is_leaf(this) result(res)
    class(aabb_node_t), intent(in) :: this
    logical :: res

    res = this%left_node_index == AABB_NULL_NODE .and. &
      this%right_node_index == AABB_NULL_NODE
  end function aabb_node_is_leaf

  !> @brief Returns true if the node is a valid node.
  pure function aabb_node_is_valid(this) result(valid)
    class(aabb_node_t), intent(in) :: this
    logical :: valid

    if (this%is_leaf()) then
       valid = &
         & this%left_node_index .eq. AABB_NULL_NODE .and. &
         & this%right_node_index .eq. AABB_NULL_NODE .and. &
         & this%object_index .gt. 0
    else
       valid = &
         & this%left_node_index .ne. AABB_NULL_NODE .and. &
         & this%right_node_index .ne. AABB_NULL_NODE .and. &
         & this%object_index .eq. -1
    end if

  end function aabb_node_is_valid

  ! -------------------------------------------------------------------------- !
  ! Comparison operators

  !> @brief Returns true if the node is less than the other node.
  pure function aabb_node_less(this, other) result(res)
    class(aabb_node_t), intent(in) :: this
    class(aabb_node_t), intent(in) :: other
    logical :: res

    res = this%aabb .lt. other%aabb

  end function aabb_node_less

  !> @brief Returns true if the node is greater than the other node.
  pure function aabb_node_greater(this, other) result(res)
    class(aabb_node_t), intent(in) :: this
    class(aabb_node_t), intent(in) :: other
    logical :: res

    res = this%aabb .gt. other%aabb

  end function aabb_node_greater

  ! ========================================================================== !
  ! Definitions of tree methods
  ! ========================================================================== !

  !> @brief Initializes the AABB tree.
  subroutine aabb_tree_init(this, initial_capacity)
    class(aabb_tree_t), intent(inout) :: this
    integer, intent(in) :: initial_capacity

    integer :: i

    this%root_node_index = AABB_NULL_NODE
    this%allocated_node_count = 0
    this%next_free_node_index = 1
    this%node_capacity = initial_capacity
    this%growth_size = initial_capacity

    if (allocated(this%nodes)) deallocate(this%nodes)
    allocate(this%nodes(initial_capacity))

    do i = 1, initial_capacity
       this%nodes(i)%next_node_index = i + 1
    end do
    this%nodes(initial_capacity)%next_node_index = AABB_NULL_NODE
  end subroutine aabb_tree_init

  !> @brief Builds the tree.
  subroutine aabb_tree_build_tree(this, objects, padding)
    use utils, only: neko_error
    implicit none

    class(aabb_tree_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: objects
    real(kind=dp), optional, intent(in) :: padding

    integer :: i_obj, i_node, i
    logical :: done

    integer :: start_layer, end_layer

    type(aabb_t), dimension(:), allocatable :: box_list
    integer, dimension(:), allocatable :: sorted_indices

    real(kind=dp) :: aabb_padding

    call this%init(size(objects) * 2)

    ! ------------------------------------------------------------------------ !
    ! Start by sorting the list of objects, then build a balanced binary tree
    ! from the sorted list

    allocate(box_list(size(objects)))

    if (present(padding)) then
       aabb_padding = padding
    else
       aabb_padding = 0.0_dp
    end if
       
    do i_obj = 1, size(objects)
       box_list(i_obj) = get_aabb(objects(i_obj), aabb_padding)
    end do
    sorted_indices = sort(box_list)

    do i = 1, size(sorted_indices)
       i_obj = sorted_indices(i)
       i_node = this%allocate_node()
       this%nodes(i_node)%aabb = get_aabb(objects(i_obj))
       this%nodes(i_node)%object_index = i_obj
    end do


    start_layer = 1
    end_layer = size(objects)
    done = .false.
    do while (.not. done)

       ! build the next layer
       do i = start_layer, end_layer - 1, 2
          i_node = this%allocate_node()

          this%nodes(i_node)%aabb = merge(this%nodes(i)%aabb, &
               this%nodes(i + 1)%aabb)

          this%nodes(i_node)%left_node_index = i
          this%nodes(i_node)%right_node_index = i + 1

          this%nodes(i)%parent_node_index = i_node
          this%nodes(i + 1)%parent_node_index = i_node
       end do

       ! if the number of nodes is odd, we need to create a new node to hold the
       ! last node
       if (mod(end_layer - start_layer, 2) .eq. 0) then
          i_node = this%allocate_node()
          this%nodes(i_node)%aabb = this%nodes(end_layer)%aabb
          this%nodes(i_node)%left_node_index = end_layer
          this%nodes(i_node)%right_node_index = AABB_NULL_NODE

          this%nodes(end_layer)%parent_node_index = i_node
       end if

       ! move to the next layer
       start_layer = end_layer + 1
       end_layer = this%allocated_node_count

       ! If there is only one node left, we are done
       done = start_layer .eq. end_layer
    end do

    ! The last node allocated is the root node
    this%root_node_index = this%allocated_node_count

    if (this%get_size() .ne. size(objects)) then
       print *, "this%get_size() = ", this%get_size()
       print *, "size(objects) = ", size(objects)
       call neko_error("Invalid tree size")
    end if

  end subroutine aabb_tree_build_tree

  function sort(array) result(indices)
    type(aabb_t), dimension(:), intent(in) :: array
    integer, dimension(:), allocatable :: indices
    logical, dimension(:), allocatable :: visited

    integer :: i, imin
    integer :: minidx

    allocate(indices(size(array)))
    allocate(visited(size(array)))

    visited = .false.
    indices = 0
    do i = 1, size(array)
       minidx = -1
       do imin = 1, size(array)
          if (.not. visited(imin) .and. minidx .eq. -1) minidx = imin

          if (visited(imin) .and. array(imin) .lt. array(minidx)) minidx = imin
       end do

       indices(i) = minidx
       visited(minidx) = .true.
    end do

  end function sort

  ! -------------------------------------------------------------------------- !
  ! Getters

  !> @brief Returns the size of the tree, in number of leaves.
  function aabb_tree_get_size(this) result(size)
    use stack, only: stack_i4_t
    use utils, only: neko_error
    class(aabb_tree_t), intent(in) :: this
    integer :: size

    type(stack_i4_t) :: simple_stack
    integer :: idx, tmp

    call simple_stack%init(this%allocated_node_count)
    size = 0
    tmp = this%get_root_index()
    if (tmp .ne. AABB_NULL_NODE) then
       call simple_stack%push(tmp)
    end if

    do while (.not. simple_stack%is_empty())
       idx = simple_stack%pop()
       if (idx .eq. AABB_NULL_NODE) cycle

       if (this%nodes(idx)%is_leaf()) then
          size = size + 1
       else
          tmp = this%get_left_index(idx)
          call simple_stack%push(tmp)
          tmp = this%get_right_index(idx)
          call simple_stack%push(tmp)
       end if
    end do

  end function aabb_tree_get_size

  ! -------------------------------------------------------------------------- !
  ! Get index of nodes

  !> @brief Returns the index of the root node.
  pure function aabb_tree_get_root_index(this) result(root_index)
    class(aabb_tree_t), intent(in) :: this
    integer :: root_index

    root_index = this%root_node_index
  end function aabb_tree_get_root_index

  !> @brief Returns the index of the parent node of the node at the given index.
  pure function aabb_tree_get_parent_index(this, node_index) &
       result(parent_index)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    integer :: parent_index

    parent_index = this%nodes(node_index)%parent_node_index
  end function aabb_tree_get_parent_index

  !> @brief Returns the index of the left node of the node at the given index.
  pure function aabb_tree_get_left_index(this, node_index) &
       result(left_index)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    integer :: left_index

    left_index = this%nodes(node_index)%left_node_index
  end function aabb_tree_get_left_index

  !> @brief Returns the index of the right node of the node at the given index.
  pure function aabb_tree_get_right_index(this, node_index) &
       result(right_index)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    integer :: right_index

    right_index = this%nodes(node_index)%right_node_index
  end function aabb_tree_get_right_index

  ! -------------------------------------------------------------------------- !
  ! Get nodes

  !> @brief Returns the node at the given index.
  pure function aabb_tree_get_node(this, node_index) result(node)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    type(aabb_node_t) :: node

    node = this%nodes(node_index)
  end function aabb_tree_get_node

  !> @brief Returns the root node of the tree.
  pure function aabb_tree_get_root_node(this) result(root_node)
    class(aabb_tree_t), intent(in) :: this
    type(aabb_node_t) :: root_node

    root_node = this%nodes(this%root_node_index)
  end function aabb_tree_get_root_node

  !> @brief Returns the parent node of the node at the given index.
  pure function aabb_tree_get_parent_node(this, node_index) &
       result(parent_node)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    type(aabb_node_t) :: parent_node

    parent_node = this%nodes(this%nodes(node_index)%parent_node_index)
  end function aabb_tree_get_parent_node

  !> @brief Returns the left node of the node at the given index.
  pure function aabb_tree_get_left_node(this, node_index) result(left_node)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    type(aabb_node_t) :: left_node

    left_node = this%nodes(this%nodes(node_index)%left_node_index)
  end function aabb_tree_get_left_node

  !> @brief Returns the right node of the node at the given index.
  pure function aabb_tree_get_right_node(this, node_index) &
       result(right_node)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    type(aabb_node_t) :: right_node

    right_node = this%nodes(this%nodes(node_index)%right_node_index)
  end function aabb_tree_get_right_node

  pure function aabb_tree_get_aabb(this, node_index) result(out_box)
    class(aabb_tree_t), intent(in) :: this
    integer, intent(in) :: node_index
    type(aabb_t) :: out_box

    out_box = this%nodes(node_index)%aabb
  end function aabb_tree_get_aabb

  ! -------------------------------------------------------------------------- !

  !> @brief Inserts an object into the tree.
  subroutine aabb_tree_insert_object(this, object, object_index)
    class(aabb_tree_t), intent(inout) :: this
    class(*), intent(in) :: object
    integer, intent(in) :: object_index

    integer :: node_index

    node_index = this%allocate_node()
    this%nodes(node_index)%aabb = get_aabb(object)
    this%nodes(node_index)%object_index = object_index

    call this%insert_leaf(node_index)
  end subroutine aabb_tree_insert_object

  !> @brief Queries the tree for overlapping objects.
  subroutine aabb_tree_query_overlaps(this, object, object_index, overlaps)
    use stack, only: stack_i4_t
    implicit none

    class(aabb_tree_t), intent(in) :: this
    class(*), intent(in) :: object
    integer, intent(in) :: object_index
    type(stack_i4_t), intent(inout) :: overlaps

    type(stack_i4_t) :: simple_stack
    type(aabb_t) :: object_box

    integer :: root_index, left_index, right_index

    integer :: node_index, tmp_index

    object_box = get_aabb(object)
    root_index = this%get_root_index()

    call simple_stack%init()
    call simple_stack%push(root_index)

    do while (.not. simple_stack%is_empty())
       node_index = simple_stack%pop()

       if (node_index == AABB_NULL_NODE) cycle

       if (this%nodes(node_index)%aabb%overlaps(object_box)) then
          if (this%nodes(node_index)%is_leaf()) then
             if (this%nodes(node_index)%object_index .ne. object_index) then
                tmp_index = this%nodes(node_index)%object_index
                call overlaps%push(tmp_index)
             end if
          else
             left_index = this%get_left_index(node_index)
             call simple_stack%push(left_index)
             right_index = this%get_right_index(node_index)
             call simple_stack%push(right_index)
          end if
       end if
    end do
  end subroutine aabb_tree_query_overlaps

  ! -------------------------------------------------------------------------- !
  ! Internal methods

  !> @brief Allocates a new node in the tree.
  function aabb_tree_allocate_node(this) result(node_index)
    class(aabb_tree_t), intent(inout) :: this
    integer :: node_index

    if (this%next_free_node_index == AABB_NULL_NODE) then
       call this%resize_node_pool(this%node_capacity + this%growth_size)
    end if

    node_index = this%next_free_node_index

    associate(new_node => this%nodes(node_index))
      this%next_free_node_index = new_node%next_node_index

      new_node%parent_node_index = AABB_NULL_NODE
      new_node%left_node_index = AABB_NULL_NODE
      new_node%right_node_index = AABB_NULL_NODE

      this%next_free_node_index = new_node%next_node_index
      this%allocated_node_count = this%allocated_node_count + 1

    end associate
  end function aabb_tree_allocate_node

  !> @brief Deallocates a node in the tree.
  subroutine aabb_tree_deallocate_node(this, node_index)
    class(aabb_tree_t), intent(inout) :: this
    integer, intent(in) :: node_index

    this%nodes(node_index)%next_node_index = this%next_free_node_index
    this%next_free_node_index = node_index
    this%allocated_node_count = this%allocated_node_count - 1
  end subroutine aabb_tree_deallocate_node

  !> @brief Inserts a leaf into the tree.
  subroutine aabb_tree_insert_leaf(this, leaf_node_index)
    class(aabb_tree_t), intent(inout) :: this
    integer, intent(in) :: leaf_node_index

    integer :: tree_node_index

    real(kind=rp) :: cost_left
    real(kind=rp) :: cost_right

    type(aabb_node_t) :: leaf_node
    type(aabb_node_t) :: tree_node
    type(aabb_node_t) :: left_node
    type(aabb_node_t) :: right_node

    type(aabb_t) :: combined_aabb
    real(kind=rp) :: new_parent_node_cost
    real(kind=rp) :: minimum_push_down_cost
    type(aabb_t) :: new_left_aabb
    type(aabb_t) :: new_right_aabb
    integer :: leaf_sibling_index
    type(aabb_node_t) :: leaf_sibling
    integer :: old_parent_index
    integer :: new_parent_index
    type(aabb_node_t) :: new_parent
    type(aabb_node_t) :: old_parent

    ! make sure were inserting a new leaf
    leaf_node = this%nodes(leaf_node_index)

    ! if the tree is empty then we make the root the leaf
    if (this%root_node_index .eq. AABB_NULL_NODE) then
       this%root_node_index = leaf_node_index
       leaf_node%parent_node_index = AABB_NULL_NODE
       leaf_node%left_node_index = AABB_NULL_NODE
       leaf_node%right_node_index = AABB_NULL_NODE

       return
    end if

    ! search for the best place to put the new leaf in the tree
    ! we use surface area and depth as search heuristics
    tree_node_index = this%root_node_index
    tree_node = this%get_node(tree_node_index)
    do while (.not. tree_node%is_leaf())

       ! because of the test in the while loop above we know we are never a
       ! leaf inside it
       left_node = this%get_left_node(tree_node_index)
       right_node = this%get_right_node(tree_node_index)

       ! ------------------------------------------------------------------- !

       combined_aabb = merge(tree_node%aabb, leaf_node%get_aabb())

       new_parent_node_cost = 2.0_rp * combined_aabb%get_surface_area()
       minimum_push_down_cost = 2.0_rp * ( &
         & combined_aabb%get_surface_area() &
         & - tree_node%aabb%get_surface_area()&
         & )

       ! use the costs to figure out whether to create a new parent here or
       ! descend
       if (left_node%is_leaf()) then
          new_left_aabb = merge(leaf_node%aabb, left_node%get_aabb())
          cost_left = new_left_aabb%get_surface_area() + minimum_push_down_cost
       else
          new_left_aabb = merge(leaf_node%aabb, left_node%get_aabb())
          cost_left = ( &
            & new_left_aabb%get_surface_area() &
            & - left_node%aabb%get_surface_area()&
            & ) + minimum_push_down_cost
       end if

       if (right_node%is_leaf()) then
          new_right_aabb = merge(leaf_node%aabb, right_node%aabb)
          cost_right = new_right_aabb%get_surface_area() + &
               minimum_push_down_cost
       else
          new_right_aabb = merge(leaf_node%aabb, right_node%aabb)
          cost_right = ( &
            & new_right_aabb%get_surface_area() &
            & - right_node%aabb%get_surface_area() &
            & ) + minimum_push_down_cost
       end if

       ! if the cost of creating a new parent node here is less than descending
       ! in either direction then we know we need to create a new parent node,
       ! errrr, here and attach the leaf to that
       if (new_parent_node_cost < cost_left .and. &
            new_parent_node_cost < cost_right) then
          exit
       end if

       ! otherwise descend in the cheapest direction
       if (cost_left .lt. cost_right) then
          tree_node_index = tree_node%get_left_index()
       else
          tree_node_index = tree_node%get_right_index()
       end if

       ! ------------------------------------------------------------------- !
       ! Update the node and continue the loop
       tree_node = this%get_node(tree_node_index)
    end do

    ! the leafs sibling is going to be the node we found above and we are
    ! going to create a new parent node and attach the leaf and this item
    leaf_sibling_index = tree_node_index
    leaf_sibling = this%nodes(leaf_sibling_index)
    old_parent_index = this%get_parent_index(leaf_sibling_index)
    new_parent_index = this%allocate_node()
    new_parent = this%nodes(new_parent_index)
    new_parent%parent_node_index = old_parent_index
    new_parent%aabb = merge(leaf_node%aabb, leaf_sibling%aabb)

    if (leaf_node .lt. leaf_sibling) then
       new_parent%left_node_index = leaf_node_index
       new_parent%right_node_index = leaf_sibling_index
    else
       new_parent%left_node_index = leaf_sibling_index
       new_parent%right_node_index = leaf_node_index
    end if

    leaf_node%parent_node_index = new_parent_index
    leaf_sibling%parent_node_index = new_parent_index

    if (old_parent_index .eq. AABB_NULL_NODE) then
       ! the old parent was the root and so this is now the root
       this%root_node_index = new_parent_index
    else
       ! the old parent was not the root and so we need to patch the left or
       ! right index to point to the new node
       old_parent = this%nodes(old_parent_index)
       if (old_parent%left_node_index .eq. leaf_sibling_index) then
          old_parent%left_node_index = new_parent_index
       else
          old_parent%right_node_index = new_parent_index
       end if
       this%nodes(old_parent_index) = old_parent
    end if

    this%nodes(leaf_node_index) = leaf_node
    this%nodes(leaf_sibling_index) = leaf_sibling
    this%nodes(new_parent_index) = new_parent

    ! finally we need to walk back up the tree fixing heights and areas
    tree_node_index = leaf_node%parent_node_index

    call this%fix_upwards_tree(tree_node_index)

  end subroutine aabb_tree_insert_leaf

  !> @brief Validates the tree.
  function aabb_tree_valid_tree(this) result(valid)
    use stack, only: stack_i4_t
    implicit none

    class(aabb_tree_t), intent(in) :: this
    logical :: valid

    type(stack_i4_t) :: simple_stack
    integer :: current_index
    integer :: root_index, left_index, right_index

    valid = .true.
    if (this%root_node_index .eq. AABB_NULL_NODE) then
       valid = .false.
    end if

    root_index = this%get_root_index()

    call simple_stack%init(this%node_capacity)
    call simple_stack%push(root_index)

    do while (.not. simple_stack%is_empty())
       current_index = simple_stack%pop()
       if (current_index == AABB_NULL_NODE) cycle

       valid = valid .and. this%nodes(current_index)%is_valid()

       if (.not. this%nodes(current_index)%is_leaf()) then
          left_index = this%get_left_index(current_index)
          right_index = this%get_right_index(current_index)

          call simple_stack%push(left_index)
          call simple_stack%push(right_index)
       end if
    end do
  end function aabb_tree_valid_tree

  !> @brief Fixes the tree upwards.
  !! @details This method is used to fix the tree upwards after an insertion.
  !! It is used to expand the nodes of the tree to fit the new leaf node.
  subroutine aabb_tree_fix_upwards_tree(this, tree_start_index)
    class(aabb_tree_t), intent(inout) :: this
    integer, intent(in) :: tree_start_index

    type(aabb_node_t) :: left_node
    type(aabb_node_t) :: right_node
    integer :: tree_node_index

    tree_node_index = tree_start_index
    do while (tree_node_index .ne. AABB_NULL_NODE)
       left_node = this%get_left_node(tree_node_index)
       right_node = this%get_right_node(tree_node_index)

       this%nodes(tree_node_index)%aabb = merge(left_node%aabb, right_node%aabb)

       tree_node_index = this%get_parent_index(tree_node_index)
    end do
  end subroutine aabb_tree_fix_upwards_tree

  !> @brief Prints the tree.
  subroutine aabb_tree_print(this)
    use stack, only: stack_i4_t
    class(aabb_tree_t), intent(inout) :: this
    type(stack_i4_t) :: simple_stack

    integer :: current_index
    integer :: root_index, left_index, right_index

    root_index = this%get_root_index()
    call simple_stack%init(this%node_capacity)
    call simple_stack%push(root_index)

    do while (.not. simple_stack%is_empty())
       current_index = simple_stack%pop()
       if (current_index .eq. AABB_NULL_NODE) cycle

       left_index = this%get_left_index(current_index)
       right_index = this%get_right_index(current_index)

       call simple_stack%push(this%nodes(current_index)%left_node_index)
       call simple_stack%push(this%nodes(current_index)%right_node_index)

       write(*, *) "i = ", current_index
       write(*, *) "  Parent  : ", this%get_parent_index(current_index)
       write(*, *) "  Children: ", this%get_left_index(current_index), &
            this%get_right_index(current_index)

       write(*, *) "  object_index = ", this%nodes(current_index)%object_index
    end do

  end subroutine aabb_tree_print

  !> @brief Resizes the node pool.
  subroutine aabb_tree_resize_node_pool(this, new_capacity)
    class(aabb_tree_t), intent(inout) :: this
    integer, intent(in) :: new_capacity

    type(aabb_node_t), dimension(:), allocatable :: temp
    integer :: i

    allocate(temp(new_capacity))
    temp(:this%node_capacity) = this%nodes(:this%node_capacity)

    do i = this%allocated_node_count, new_capacity
       temp(i)%next_node_index = i + 1
    end do
    temp(new_capacity)%next_node_index = AABB_NULL_NODE

    call move_alloc(temp, this%nodes)

    this%node_capacity = new_capacity
    this%next_free_node_index = this%allocated_node_count + 1

  end subroutine aabb_tree_resize_node_pool

end module aabb_tree
