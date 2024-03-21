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
! https://github.com/JamesRandall/SimpleVoxelEngine/blob/master/voxelEngine/include/AABB.h
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

!> @brief Axis Aligned Bounding Box (aabb) implementation in Fortran
!! @details
!! This is a Fortran implementation of an Axis Aligned Bounding Box (aabb) data
!! structure. The aabb is a box that is aligned to the x, y and z axes. It is
!! defined by two points, the lower left front corner and the upper right back
!! corner. This is the base data structure for the aabb_Tree, which is used to
!! accelerate a Signed Distance Function.
module aabb
  use num_types, only: dp
  use element, only: element_t
  use point, only: point_t
  use tri, only: tri_t
  use quad, only: quad_t
  use tet, only: tet_t
  use hex, only: hex_t
  use utils, only: neko_error

  implicit none
  private
  public :: aabb_t, get_aabb, merge, intersection

  ! ========================================================================== !
  ! Public interface for free functions
  ! ========================================================================== !

  !> @brief Merge two aabbs.
  interface merge
     module procedure merge_aabb
  end interface merge

  !> @brief Intersect two aabbs.
  interface intersection
     module procedure intersection_aabb
  end interface intersection

  !> @brief Axis Aligned Bounding Box (aabb) data structure.
  !! @details The aabb is a box that is aligned to the x, y and z axes. It is
  !! defined by two points, the lower left front corner and the upper right back
  !! corner. The purpose of this is to accelerate a Signed Distance Function,
  !! through an aabb_Tree.
  type :: aabb_t
     private

     logical :: initialized = .false.
     real(kind=dp) :: box_min(3)
     real(kind=dp) :: box_max(3)
     real(kind=dp) :: center(3)
     real(kind=dp) :: diameter
     real(kind=dp) :: surface_area

   contains

     ! Initializers
     procedure, pass(this), public :: init => aabb_init

     ! Getters
     procedure, pass(this), public :: get_width => aabb_get_width
     procedure, pass(this), public :: get_height => aabb_get_height
     procedure, pass(this), public :: get_depth => aabb_get_depth
     procedure, pass(this), public :: get_diameter => aabb_get_diameter
     procedure, pass(this), public :: get_surface_area => aabb_get_surface_area
     procedure, pass(this), public :: get_center => aabb_get_center
     procedure, pass(this), public :: get_diagonal => aabb_get_diagonal

     ! Unary operations
     procedure, pass(this), public :: min_distance => aabb_min_distance

     ! Comparison operators
     generic :: operator(.lt.) => less
     generic :: operator(.gt.) => greater

     !> @brief Check if two aabbs are overlapping.
     procedure, pass(this), public :: overlaps => aabb_overlaps
     !> @brief Check if this aabb fully contains another aabb.
     procedure, pass(this), public :: contains => aabb_contains_other
     !> @brief Check if this aabb contains a point.
     procedure, pass(this), public :: contains_point => aabb_contains_point

     ! Private comparison operators
     procedure, pass(this) :: less => aabb_less
     procedure, pass(this) :: greater => aabb_greater

     ! Private operations
     procedure, pass(this), private :: calculate_surface_area

  end type aabb_t

contains

  ! ========================================================================== !
  ! Constructors
  ! ========================================================================== !

  !> @brief Construct the aabb of a predefined object.
  !!
  !! @details This function is used to get the aabb of a predefined object.
  !! Optionally, the user can define the padding of the aabb, which is a
  !! multiple of the diameter of the aabb. This is used to avoid numerical
  !! issues when the object itself it axis aligned.
  !!
  !! Current support:
  !! - Triangle (tri_t)
  !! - Quadrilateral (quad_t)
  !! - Tetrahedron (tet_t)
  !! - Hexahedron (hex_t)
  !!
  !! @param[in] object The object to get the aabb of.
  !! @param[in] padding The padding of the aabb.
  !! @return The aabb of the object.
  function get_aabb(object, padding) result(box)
    use utils, only: neko_error
    implicit none

    class(*), intent(in) :: object
    real(kind=dp), intent(in), optional :: padding
    type(aabb_t) :: box

    select type(object)
      type is (tri_t)
       box = get_aabb_element(object, padding)
      type is (hex_t)
       box = get_aabb_element(object, padding)
      type is (tet_t)
       box = get_aabb_element(object, padding)
      type is (quad_t)
       box = get_aabb_element(object, padding)

      class default
       call neko_error("get_aabb: Unsupported object type")
    end select

  end function get_aabb

  !> @brief Get the aabb of an arbitrary element.
  !!
  !! @details This function calculates the aabb of an element. The aabb is
  !! defined by the lower left front corner and the upper right back corner.
  !! The aabb is calculated by finding the minimum and maximum x, y and z
  !! coordinate for all points in the arbitrary element type.
  !!
  !! @param element The arbitrary element to get the aabb of.
  !! @return The aabb of the element.
  function get_aabb_element(object, padding) result(box)
    class(element_t), intent(in) :: object
    real(kind=dp), intent(in), optional :: padding
    type(aabb_t) :: box

    real(kind=dp), dimension(3) :: box_min, box_max
    integer :: i
    type(point_t), pointer :: pi

    do i = 1, object%n_points()
       pi => object%p(i)
       box_min = min(box_min, pi%x)
       box_max = max(box_max, pi%x)
    end do

    if (present(padding)) then
       box_min = box_min - padding * (box_max - box_min)
       box_max = box_max + padding * (box_max - box_min)
    end if

    call box%init(box_min, box_max)
  end function get_aabb_element


  ! ========================================================================== !
  ! Initializers
  ! ========================================================================== !

  !> @brief Initialize the aabb.
  !! @param lower_left_front The lower left front corner of the aabb.
  !! @param upper_right_back The upper right back corner of the aabb.
  subroutine aabb_init(this, lower_left_front, upper_right_back)
    class(aabb_t), intent(inout) :: this
    real(kind=dp), dimension(3), intent(in) :: lower_left_front
    real(kind=dp), dimension(3), intent(in) :: upper_right_back

    this%box_min = lower_left_front
    this%box_max = upper_right_back
    this%center = (this%box_min + this%box_max) / 2.0_dp
    this%diameter = norm2(this%box_max - this%box_min)

    this%initialized = .true.
    this%surface_area = this%calculate_surface_area()

  end subroutine aabb_init

  ! ========================================================================== !
  ! Getters
  ! ========================================================================== !

  !> @brief Get the width of the aabb. Also known as the x-axis length.
  pure function aabb_get_width(this) result(width)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: width

    width = this%box_max(1) - this%box_min(1)
  end function aabb_get_width

  !> @brief Get the depth of the aabb. Also known as the y-axis length.
  pure function aabb_get_depth(this) result(depth)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: depth

    depth = this%box_max(2) - this%box_min(2)
  end function aabb_get_depth

  !> @brief Get the height of the aabb. Also known as the z-axis length.
  pure function aabb_get_height(this) result(height)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: height

    height = this%box_max(3) - this%box_min(3)
  end function aabb_get_height

  !> @brief Get the diameter length of the aabb.
  pure function aabb_get_diameter(this) result(diameter)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: diameter

    diameter = this%diameter
  end function aabb_get_diameter

  !> @brief Get the surface area of the aabb.
  pure function aabb_get_surface_area(this) result(surface_area)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: surface_area

    surface_area = this%surface_area
  end function aabb_get_surface_area

  !> @brief Get the center of the aabb.
  pure function aabb_get_center(this) result(center)
    class(aabb_t), intent(in) :: this
    real(kind=dp), dimension(3) :: center

    center = this%center
  end function aabb_get_center

  !> @brief Get the diagonal of the aabb.
  pure function aabb_get_diagonal(this) result(diagonal)
    class(aabb_t), intent(in) :: this
    real(kind=dp), dimension(3) :: diagonal

    diagonal = this%box_max - this%box_min
  end function aabb_get_diagonal

  ! ========================================================================== !
  ! Operations
  ! ========================================================================== !

  !> @brief Check if two aabbs are overlapping.
  function aabb_overlaps(this, other) result(is_overlapping)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: is_overlapping

    if (.not. this%initialized .or. .not. other%initialized) then
       !  call neko_error("aabb_overlaps: One or both aabbs are not initialized")
       is_overlapping = .false.
    else

       is_overlapping = all(this%box_min .le. other%box_max) .and. &
         all(this%box_max .ge. other%box_min)
    end if

  end function aabb_overlaps

  !> @brief Check if this aabb contains another aabb.
  function aabb_contains_other(this, other) result(is_contained)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: is_contained

    ! if (.not. this%initialized .or. .not. other%initialized) then
    !  call neko_error("aabb_contains: One or both aabbs are not initialized")
    ! end if

    is_contained = all(this%box_min .le. other%box_min) .and. &
      all(this%box_max .ge. other%box_max)

  end function aabb_contains_other

  !> @brief Check if this aabb contains a point.
  function aabb_contains_point(this, p) result(is_contained)
    class(aabb_t), intent(in) :: this
    real(kind=dp), dimension(3), intent(in) :: p
    logical :: is_contained

    ! if (.not. this%initialized) then
    !  call neko_error("aabb_contains_point: One or both aabbs are not initialized")
    ! end if

    is_contained = all(p .ge. this%box_min) .and. all(p .le. this%box_max)
  end function aabb_contains_point

  !> @brief Get the minimum possible distance from the aabb to a point.
  function aabb_min_distance(this, p) result(distance)
    class(aabb_t), intent(in) :: this
    real(kind=dp), dimension(3), intent(in) :: p
    real(kind=dp) :: distance

    if (.not. this%initialized) then
       distance = huge(0.0_dp)
    end if

    distance = this%get_diameter() / 2.0_dp - norm2(this%get_center() - p)
  end function aabb_min_distance

  ! ========================================================================== !
  ! Binary operations
  ! ========================================================================== !

  !> @brief Merge two aabbs.
  function merge_aabb(box1, box2) result(merged)
    class(aabb_t), intent(in) :: box1
    class(aabb_t), intent(in) :: box2
    type(aabb_t) :: merged

    real(kind=dp), dimension(3) :: box_min, box_max

    box_min = min(box1%box_min, box2%box_min)
    box_max = max(box1%box_max, box2%box_max)

    call merged%init(box_min, box_max)
  end function merge_aabb

  !> @brief Get the intersection of two aabbs.
  function intersection_aabb(box1, box2) result(intersected)
    class(aabb_t), intent(in) :: box1
    class(aabb_t), intent(in) :: box2
    type(aabb_t) :: intersected

    real(kind=dp), dimension(3) :: box_min, box_max

    box_min = max(box1%box_min, box2%box_min)
    box_max = min(box1%box_max, box2%box_max)

    call intersected%init(box_min, box_max)
  end function intersection_aabb

  ! ========================================================================== !
  ! Private operations
  ! ========================================================================== !

  !> @brief Calculate the surface area of the aabb.
  pure function calculate_surface_area(this) result(surface_area)
    class(aabb_t), intent(in) :: this
    real(kind=dp) :: surface_area

    surface_area = 2.0 * (&
      & this%get_width() * this%get_height() &
      & + this%get_width() * this%get_depth() &
      & + this%get_height() * this%get_depth() &
      &)
  end function calculate_surface_area

  ! ========================================================================== !
  ! Comparison operators
  ! ========================================================================== !

  !> @brief Less than comparison operator.
  pure function aabb_less(this, other)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: aabb_less
    logical :: equal

    if (.not. this%initialized .or. .not. other%initialized) then
       aabb_less = .false.
       return
    end if

    aabb_less = this%box_min(1) .lt. other%box_min(1)
    equal = this%box_min(1) .le. other%box_min(1)

    if (.not. aabb_less .and. equal) then
       aabb_less = this%box_min(2) .lt. other%box_min(2)
       equal = this%box_min(2) .le. other%box_min(2)
    end if
    if (.not. aabb_less .and. equal) then
       aabb_less = this%box_min(3) .lt. other%box_min(3)
       equal = this%box_min(3) .le. other%box_min(3)
    end if
    if (.not. aabb_less .and. equal) then
       aabb_less = this%box_max(1) .lt. other%box_max(1)
       equal = this%box_max(1) .le. other%box_max(1)
    end if
    if (.not. aabb_less .and. equal) then
       aabb_less = this%box_max(2) .lt. other%box_max(2)
       equal = this%box_max(2) .le. other%box_max(2)
    end if
    if (.not. aabb_less .and. equal) then
       aabb_less = this%box_max(3) .lt. other%box_max(3)
    end if

  end function aabb_less

  !> @brief Greater than comparison operator.
  pure function aabb_greater(this, other)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: aabb_greater
    logical :: equal

    if (.not. this%initialized .or. .not. other%initialized) then
       aabb_greater = .false.
       return
    end if

    aabb_greater = this%box_min(1) .gt. other%box_min(1)
    equal = this%box_min(1) .ge. other%box_min(1)

    if (.not. aabb_greater .and. equal) then
       aabb_greater = this%box_min(2) .gt. other%box_min(2)
       equal = this%box_min(2) .ge. other%box_min(2)
    end if
    if (.not. aabb_greater .and. equal) then
       aabb_greater = this%box_min(3) .gt. other%box_min(3)
       equal = this%box_min(3) .ge. other%box_min(3)
    end if
    if (.not. aabb_greater .and. equal) then
       aabb_greater = this%box_max(1) .gt. other%box_max(1)
       equal = this%box_max(1) .ge. other%box_max(1)
    end if
    if (.not. aabb_greater .and. equal) then
       aabb_greater = this%box_max(2) .gt. other%box_max(2)
       equal = this%box_max(2) .ge. other%box_max(2)
    end if
    if (.not. aabb_greater .and. equal) then
       aabb_greater = this%box_max(3) .gt. other%box_max(3)
    end if

  end function aabb_greater

end module aabb
