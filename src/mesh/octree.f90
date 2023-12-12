! Copyright (c) 2022, The Neko Authors
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
!> Implements an Octree
!! @details Fast search/lookup of points
module octree
  use num_types
  use point
  use utils
  implicit none
  private

  type oct_ptr_t
     type(oct_t), pointer :: ptr => null()
  end type oct_ptr_t

  !> Defines an octree octant
  type oct_t
     type(oct_ptr_t) :: oct(8)
     integer :: level
     type(point_t) :: origin
     real(kind=dp) :: width
     type(point_t) :: point
     logical valid
  end type oct_t

  !> Defines an octree
  type, public ::  octree_t
     type(oct_t), pointer :: root => null()
   contains
     procedure, pass(t) :: init => octree_init
     procedure, pass(t) :: free => octree_free
     procedure, pass(t) :: insert => octree_insert
     procedure, pass(t) :: find => octree_find
  end type octree_t

contains

  !> Initialize an octree
  subroutine octree_init(t, width)
    class(octree_t), intent(inout) :: t
    real(kind=dp), intent(in) :: width
    type(point_t) :: origin
    integer, parameter :: top_level = 0

    call octree_free(t)

    origin = (/ 0d0, 0d0, 0d0 /)
    call octree_oct_init(t%root, origin, width, top_level)

  end subroutine octree_init

  !> Destroy an octree
  subroutine octree_free(t)
    class(octree_t), intent(inout) :: t

    call octree_free_oct(t%root)

  end subroutine octree_free

  !> Insert a point @a p into the octree
  subroutine octree_insert(t, p)
    class(octree_t), intent(inout) :: t
    type(point_t), intent(in) :: p

    call octree_oct_insert(t%root, p)

  end subroutine octree_insert

  !> Find a point @a p in an octree
  function octree_find(t, p) result(rcode)
    class(octree_t), intent(in) :: t
    type(point_t), intent(in) :: p
    logical rcode

    rcode = (octree_oct_find(t%root, p) .eq. 0)

  end function octree_find


  !> Insert a point @a p into the octree rooted at @a o
  recursive subroutine octree_oct_insert(o, p)
    type(oct_t), intent(inout), pointer :: o
    type(point_t), intent(in) :: p
    type(point_t) :: tmp_pt, offset, new_origin
    integer :: i

    if (.not. associated(o%oct(1)%ptr)) then
       if (.not. o%valid) then
          o%point = p
          o%valid = .true.
          return
       else
          if (o%point .eq. p) then
             return
          else
             tmp_pt = o%point
             o%valid = .false.

             do i = 1, 8
                offset = (/ -0.5d0, -0.5d0, -0.5d0 /)
                if (iand((i - 1), 4) .gt. 0) offset%x(1) = 0.5d0
                if (iand((i - 1), 2) .gt. 0) offset%x(2) = 0.5d0
                if (iand((i - 1), 1) .gt. 0) offset%x(3) = 0.5d0

                new_origin = o%origin%x + (o%width * offset%x)
                call octree_oct_init(o%oct(i)%ptr, new_origin, &
                     o%width * 0.5d0, o%level + 1)
             end do

             call octree_oct_insert(o%oct(octree_oct(o, tmp_pt))%ptr, tmp_pt)
             call octree_oct_insert(o%oct(octree_oct(o, p))%ptr, p)

          end if
       end if
    else
       call octree_oct_insert(o%oct(octree_oct(o, p))%ptr, p)
    end if
  end subroutine octree_oct_insert

  !> Find the octant containing a point @a p
  recursive function octree_oct_find(o, p) result(rcode)
    type(oct_t), pointer, intent(in) :: o
    type(point_t), intent(in) :: p
    integer :: rcode
    integer :: oct_idx

    rcode = 1

    if (.not. associated(o%oct(1)%ptr)) then
       if (o%valid .and. octree_oct_inside(o, p)) then
          rcode = 0
          return
       end if
    else
       oct_idx = octree_oct(o, p)
       rcode = octree_oct_find(o%oct(oct_idx)%ptr, p)
    end if

  end function octree_oct_find

  !> Initialize an octant width a given width, origin and level
  subroutine octree_oct_init(o, origin, width, level)
    type(oct_t), pointer, intent(inout) :: o
    type(point_t), intent(in) :: origin
    real(kind=dp), intent(in) :: width
    integer, intent(in) :: level
    integer :: i

    if (associated(o)) then
       call neko_error('Octree octant already initialized')
    else
       allocate(o)
       o%origin = origin
       o%width = width
       o%level = level
       o%valid = .false.

       do i = 1, 8
          nullify(o%oct(i)%ptr)
       end do
    end if

  end subroutine octree_oct_init

  !> Deallocate an oct in an octree
  recursive subroutine octree_free_oct(o)
    type(oct_t), pointer, intent(inout) :: o
    integer :: i

    if (.not. associated(o)) then
       return
    else if (.not. associated(o%oct(1)%ptr)) then
       deallocate(o)
       nullify(o)
       return
    else
       do i = 1, 8
          call octree_free_oct(o%oct(i)%ptr)
       end do
       deallocate(o)
       nullify(o)
    end if


  end subroutine octree_free_oct

  !> Return the octant for a given point
  pure function octree_oct(oct, point) result(oct_idx)
    type(oct_t), pointer, intent(in) :: oct
    type(point_t), intent(in) :: point
    integer :: oct_idx

    oct_idx = 0
    if (point%x(1) .ge. oct%origin%x(1)) oct_idx = ior(oct_idx, 4)
    if (point%x(2) .ge. oct%origin%x(2)) oct_idx = ior(oct_idx, 2)
    if (point%x(3) .ge. oct%origin%x(3)) oct_idx = ior(oct_idx, 1)
    oct_idx = oct_idx + 1

  end function octree_oct

  !> Return if a point is inside an octant
  pure function octree_oct_inside(oct, point) result(inside)
    type(oct_t), pointer, intent(in) :: oct
    type(point_t), intent(in) :: point
    logical :: inside

    inside = .false.
    if ((point%x(1) .ge. (oct%origin%x(1) - oct%width) .and. &
         point%x(1) .le. (oct%origin%x(1) + oct%width)) .and. &
         (point%x(2) .ge. (oct%origin%x(2) - oct%width) .and. &
         point%x(2) .le. (oct%origin%x(2) + oct%width)) .and. &
         (point%x(3) .ge. (oct%origin%x(3) - oct%width) .and. &
         point%x(3) .le. (oct%origin%x(3) + oct%width))) then
       inside = .true.
    end if
  end function octree_oct_inside

end module octree
