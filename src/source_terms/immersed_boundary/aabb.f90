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
  use num_types, only: rp
  use tri, only: tri_t

  implicit none
  private
  public :: aabb_t, get_aabb

  !> @brief Axis Aligned Bounding Box (aabb) data structure.
  !! @details The aabb is a box that is aligned to the x, y and z axes. It is
  !! defined by two points, the lower left front corner and the upper right back
  !! corner. The purpose of this is to accelerate a Signed Distance Function,
  !! through an aabb_Tree.
  type :: aabb_t
     private

     real(kind=rp) :: box_min(3) = [huge(0.0_rp), huge(0.0_rp), huge(0.0_rp)]
     real(kind=rp) :: box_max(3) = [-huge(0.0_rp), -huge(0.0_rp), -huge(0.0_rp)]
     real(kind=rp) :: center(3)
     real(kind=rp) :: diagonal
     real(kind=rp) :: surface_area

   contains

     ! Initializers
     procedure, pass(this), public :: init => aabb_init

     ! Getters
     procedure, pass(this), public :: get_width => aabb_get_width
     procedure, pass(this), public :: get_height => aabb_get_height
     procedure, pass(this), public :: get_depth => aabb_get_depth
     procedure, pass(this), public :: get_diagonal => aabb_get_diagonal
     procedure, pass(this), public :: get_surface_area => aabb_get_surface_area
     procedure, pass(this), public :: get_center => aabb_get_center

     ! Binary operations
     procedure, pass(this), public :: merge => aabb_merge
     procedure, pass(this), public :: intersection => aabb_intersection

     ! Unary operations
     procedure, pass(this), public :: min_distance => aabb_min_distance

     ! Comparison operators
     generic :: operator(.lt.) => less
     generic :: operator(.gt.) => greater

     procedure, pass(this), public :: overlaps => aabb_overlaps
     procedure, pass(this), public :: contains => aabb_contains_other
     procedure, pass(this), public :: contains_point => aabb_contains_point

     ! Private comparison operators
     procedure, pass(this) :: less => aabb_less
     procedure, pass(this) :: greater => aabb_greater

  end type aabb_t

contains

  ! ========================================================================== !
  ! Constructors
  ! ========================================================================== !

  !> @brief Construct the aabb of a predefined object.
  !! @details This function is used to get the aabb of a predefined object.
  !! Optionally, the user can define the padding of the aabb, which is a
  !! multiple of the diagonal of the aabb. This is used to avoid numerical
  !! issues when the object itself it axis aligned.
  !!
  !! Current support:
  !! - Triangle (tri_t)
  !!
  !! @param[in] object The object to get the aabb of.
  !! @param[in] padding The padding of the aabb.
  !! @return The aabb of the object.
  function get_aabb(object, padding) result(box)
    use utils, only: neko_error
    implicit none

    class(*), intent(in) :: object
    real(kind=rp), intent(in), optional :: padding
    type(aabb_t) :: box

    select type(object)
      type is (tri_t)
       box = get_aabb_triangle(object)

      class default
       print *, "Error: get_aabb not implemented for this type"
       stop
    end select

    if (present(padding)) then
       box%box_min = box%box_min - padding * box%diagonal
       box%box_max = box%box_max + padding * box%diagonal
    end if

  end function get_aabb

  !> @brief Get the aabb of a triangle.
  !! @details This function calculates the aabb of a triangle. The padding is a
  !! multiple of the diagonal of the aabb, and is used to avoid numerical issues
  !! when the triangle itself it axis aligned.
  !! @param triangle The triangle to get the aabb of.
  !! @return The aabb of the triangle.
  function get_aabb_triangle(triangle) result(aabb)
    type(tri_t), intent(in) :: triangle
    type(aabb_t) :: aabb

    real(kind=rp), dimension(3) :: box_min, box_max

    associate(pts => triangle%pts)
      box_min = min(pts(1)%p%x, pts(2)%p%x, pts(3)%p%x)
      box_max = max(pts(1)%p%x, pts(2)%p%x, pts(3)%p%x)
    end associate

    call aabb%init(box_min, box_max)
  end function get_aabb_triangle

  ! ========================================================================== !
  ! Initializers
  ! ========================================================================== !

  !> @brief Initialize the aabb.
  !! @param lower_left_front The lower left front corner of the aabb.
  !! @param upper_right_back The upper right back corner of the aabb.
  subroutine aabb_init(this, lower_left_front, upper_right_back)
    class(aabb_t), intent(inout) :: this
    real(kind=rp), dimension(3), intent(in) :: lower_left_front
    real(kind=rp), dimension(3), intent(in) :: upper_right_back

    if (.not. all(upper_right_back >= lower_left_front)) then
       print *, "Error: upper_right_back must be greater than lower_left_front"
       print *, "lower_left_front: ", lower_left_front
       print *, "upper_right_back: ", upper_right_back

       stop
    end if

    this%box_min = lower_left_front
    this%box_max = upper_right_back
    this%center = (this%box_min + this%box_max) / 2.0_rp
    this%diagonal = norm2(this%box_max - this%box_min)
    this%surface_area = calculate_surface_area(this)
  end subroutine aabb_init

  ! ========================================================================== !
  ! Getters
  ! ========================================================================== !

  !> @brief Get the width of the aabb. Also known as the x-axis length.
  pure function aabb_get_width(this) result(width)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: width

    width = this%box_max(1) - this%box_min(1)
  end function aabb_get_width

  !> @brief Get the depth of the aabb. Also known as the y-axis length.
  pure function aabb_get_depth(this) result(depth)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: depth

    depth = this%box_max(2) - this%box_min(2)
  end function aabb_get_depth

  !> @brief Get the height of the aabb. Also known as the z-axis length.
  pure function aabb_get_height(this) result(height)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: height

    height = this%box_max(3) - this%box_min(3)
  end function aabb_get_height

  !> @brief Get the diagonal length of the aabb.
  pure function aabb_get_diagonal(this) result(diagonal)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: diagonal

    diagonal = this%diagonal
  end function aabb_get_diagonal

  !> @brief Get the surface area of the aabb.
  pure function aabb_get_surface_area(this) result(surface_area)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: surface_area

    surface_area = this%surface_area
  end function aabb_get_surface_area

  !> @brief Get the center of the aabb.
  pure function aabb_get_center(this) result(center)
    class(aabb_t), intent(in) :: this
    real(kind=rp), dimension(3) :: center

    center = this%center
  end function aabb_get_center

  ! ========================================================================== !
  ! Operations
  ! ========================================================================== !

  !> @brief Check if two aabbs are overlapping.
  pure function aabb_overlaps(this, other) result(is_overlapping)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: is_overlapping

    is_overlapping = this%box_max(1) >= other%box_min(1) .and. &
      this%box_min(1) <= other%box_max(1) .and. &
      this%box_max(2) >= other%box_min(2) .and. &
      this%box_min(2) <= other%box_max(2) .and. &
      this%box_max(3) >= other%box_min(3) .and. &
      this%box_min(3) <= other%box_max(3)
  end function aabb_overlaps

  !> @brief Check if this aabb contains another aabb.
  pure function aabb_contains_other(this, other) result(is_contained)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    logical :: is_contained

    is_contained = other%box_min(1) >= this%box_min(1) .and. &
      other%box_max(1) <= this%box_max(1) .and. &
      other%box_min(2) >= this%box_min(2) .and. &
      other%box_max(2) <= this%box_max(2) .and. &
      other%box_min(3) >= this%box_min(3) .and. &
      other%box_max(3) <= this%box_max(3)
  end function aabb_contains_other

  !> @brief Check if this aabb contains a point.
  pure function aabb_contains_point(this, p) result(is_contained)
    class(aabb_t), intent(in) :: this
    real(kind=rp), dimension(3), intent(in) :: p
    logical :: is_contained

    is_contained = &
      p(1) .ge. this%box_min(1) .and. &
      p(1) .le. this%box_max(1) .and. &
      p(2) .ge. this%box_min(2) .and. &
      p(2) .le. this%box_max(2) .and. &
      p(3) .ge. this%box_min(3) .and. &
      p(3) .le. this%box_max(3)
  end function aabb_contains_point

  !> @brief Merge two aabbs.
  function aabb_merge(this, other) result(merged)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    type(aabb_t) :: merged

    call merged%init( &
      & min(this%box_min, other%box_min), &
      & max(this%box_max, other%box_max) &
      & )
  end function aabb_merge

  !> @brief Get the intersection of two aabbs.
  function aabb_intersection(this, other) result(intersection)
    class(aabb_t), intent(in) :: this
    class(aabb_t), intent(in) :: other
    type(aabb_t) :: intersection

    call intersection%init( &
      & min(this%box_min, other%box_min), &
      & max(this%box_max, other%box_max) &
      & )
  end function aabb_intersection

  !> @brief Get the minimum possible distance from the aabb to a point.
  pure function aabb_min_distance(this, p) result(distance)
    class(aabb_t), intent(in) :: this
    real(kind=rp), dimension(3), intent(in) :: p
    real(kind=rp) :: distance

    distance = this%get_diagonal() / 2.0_rp - norm2(this%get_center() - p)
  end function aabb_min_distance

  ! ========================================================================== !
  ! Private operations
  ! ========================================================================== !

  !> @brief Calculate the surface area of the aabb.
  pure function calculate_surface_area(this) result(surface_area)
    class(aabb_t), intent(in) :: this
    real(kind=rp) :: surface_area

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
