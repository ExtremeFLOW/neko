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
! Implements a geometry subset that combines different zones.
module combine_point_zone
  use point_zone, only: point_zone_t, point_zone_wrapper_t
  use box_point_zone, only: box_point_zone_t
  use sphere_point_zone, only: sphere_point_zone_t
  use cylinder_point_zone, only: cylinder_point_zone_t
  use num_types, only: rp
  use json_utils, only: json_get, json_get_or_default
  use json_module, only: json_file, json_core, json_value
  use utils, only: neko_error, concat_string_array
  use logger, only: neko_log
  implicit none
  private

  ! List of all possible types that one can combine
  character(len=20) :: KNOWN_TYPES(3) = [character(len=20) :: &
       "box", &
       "sphere", &
       "cylinder"]

  !> A point zone that combines different point zones.
  type, public, extends(point_zone_t) :: combine_point_zone_t
     !> List of sub-zones to construct.
     type(point_zone_wrapper_t), allocatable :: internal_zones(:)
     !> Number of zones to construct.
     integer :: n_zones = 0
     !> Operator with which to combine the point zones (AND, OR, XOR)
     character(len=:), allocatable :: operator
   contains
     !> Constructor from json object file.
     procedure, pass(this) :: init => combine_point_zone_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => combine_point_zone_free
     !> Defines the criterion of selection of a GLL point in the combine point 
     !! zone.
     procedure, pass(this) :: criterion => combine_point_zone_criterion
  end type combine_point_zone_t

contains

  !> Constructor from json object file. Reads.
  !! @param json Json object file.
  !! @param size Size with which to initialize the stack
  subroutine combine_point_zone_init_from_json(this, json, size)
    class(combine_point_zone_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(in) :: size

    ! Json low-level manipulator.
    type(json_core) :: core
    ! Pointer to the source_terms JSON object and the individual sources.
    type(json_value), pointer :: source_object, source_pointer
    ! Buffer for serializing the json.
    character(len=:), allocatable :: buffer
    ! A single source term as its own json_file.
    type(json_file) :: source_subdict
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    character(len=:), allocatable :: str_read
    integer :: i, n_zones
    logical :: found

    call json_get(json, "name", str_read)
    call this%init_base(size, trim(str_read))

    call json%get_core(core)
    call json%get('subsets', source_object, found)

    if (.not. found) call neko_error("No subsets found")

    n_zones = core%count(source_object)
    this%n_zones = n_zones

    ! Allocate arrays if we found things
    if (n_zones .gt. 0) allocate(this%internal_zones(n_zones))

    do i = 1, n_zones

       ! Create a new json containing just the subdict for this source.
       call core%get_child(source_object, i, source_pointer, found)
       call core%print_to_string(source_pointer, buffer)
       call source_subdict%load_from_string(buffer)

       call json_get(source_subdict, "geometry", type_name)

       if (trim(type_name) .eq. "box") then
          allocate(box_point_zone_t::this%internal_zones(i)%pz)
       else if (trim(type_name) .eq. "sphere") then
          allocate(sphere_point_zone_t::this%internal_zones(i)%pz)
       else if (trim(type_name) .eq. "cylinder") then
          allocate(cylinder_point_zone_t::this%internal_zones(i)%pz)
       else
          type_string = concat_string_array(KNOWN_TYPES, &
               NEW_LINE('A') // "-  ", .true.)
          call neko_error("Unknown point zone type: " &
               // trim(type_name) // ".  Known types are: " &
               // type_string)
       end if

       ! Initialize the sub-point zone but do not map anything, here we
       ! initialize the stack with just size 1.
       call this%internal_zones(i)%pz%init(source_subdict, 1)

    end do

    ! Chcek that we got the proper operator
    call json_get_or_default(json, "operator", this%operator, "OR")
    select case (trim(this%operator))
    case ("OR")
    case ("AND")
    case ("XOR")
    case default
       call neko_error("Unknown operator " // trim(this%operator))
    end select

    call json_get_or_default(json, "invert", this%inverse, .false.)

  end subroutine combine_point_zone_init_from_json

  !> Destructor.
  subroutine combine_point_zone_free(this)
    class(combine_point_zone_t), intent(inout) :: this

    integer :: i

    if (allocated(this%internal_zones)) then
       do i = 1, this%n_zones
          call this%internal_zones(i)%pz%free
       end do
       deallocate(this%internal_zones)
    end if

    this%n_zones = 0
    call this%free_base()

  end subroutine combine_point_zone_free

  !> Defines the criterion of selection of a GLL point in the combined point
  !! zone.
  !! @param x x-coordinate of the GLL point.
  !! @param y y-coordinate of the GLL point.
  !! @param z z-coordinate of the GLL point.
  !! @param j 1st nonlinear index of the GLL point.
  !! @param k 2nd nonlinear index of the GLL point.
  !! @param l 3rd nonlinear index of the GLL point.
  !! @param e element index of the GLL point.
  pure function combine_point_zone_criterion(this, x, y, z, j, k, l, e) &
       result(is_inside)
    class(combine_point_zone_t), intent(in) :: this
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    logical :: is_inside

    integer :: i

    is_inside = this%internal_zones(1)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%internal_zones(1)%pz%inverse

    do i = 2, this%n_zones
       select case (trim(this%operator))
       case ("OR")
          is_inside = is_inside .or. (this%internal_zones(i)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%internal_zones(i)%pz%inverse)

       case ("AND")
          is_inside = is_inside .and. (this%internal_zones(i)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%internal_zones(i)%pz%inverse)

       case ("XOR")
          is_inside = is_inside .neqv. (this%internal_zones(i)%pz%criterion(x, &
               y, z, j, k, l, e).neqv. this%internal_zones(i)%pz%inverse)

       case default
       end select
    end do

  end function combine_point_zone_criterion

end module combine_point_zone
