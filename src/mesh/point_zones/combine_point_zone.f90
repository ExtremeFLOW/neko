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
  use point_zone, only : point_zone_t, point_zone_pointer_t, &
       point_zone_wrapper_t, point_zone_factory
  use box_point_zone, only : box_point_zone_t
  use sphere_point_zone, only : sphere_point_zone_t
  use cylinder_point_zone, only : cylinder_point_zone_t
  use num_types, only : rp
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file, json_core, json_value
  use utils, only : neko_error, concat_string_array
  use logger, only : neko_log
  implicit none
  private

  !> A point zone that combines different point zones.
  type, public, extends(point_zone_t) :: combine_point_zone_t
     !> List of all the sub zones.
     type(point_zone_pointer_t), allocatable :: zones(:)
     !> List of the sub-zones to be created internally
     type(point_zone_wrapper_t), allocatable :: internal_zones(:)
     !> List of the names of the sub-zones to construct.
     character(len=80), allocatable :: names(:)
     !> Number of total zones.
     integer :: n_zones = 0
     !> Number of external zones to be filled by the registry.
     integer :: n_external_zones = 0
     !> Number of internal zone, to be created inside init.
     integer :: n_internal_zones = 0

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
    ! Pointer to the point_zones JSON object and the individual sources.
    type(json_value), pointer :: source_object, source_pointer
    ! Buffer for serializing the json.
    character(len=:), allocatable :: buffer
    ! A single source term as its own json_file.
    type(json_file) :: source_subdict
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    character(len=:), allocatable :: str_read
    integer :: i, n_zones, i_internal, i_external
    logical :: found, invert

    call json_get(json, "name", str_read)
    call json_get_or_default(json, "invert", invert, .false.)
    call this%init_base(size, trim(str_read), invert)

    call json%get_core(core)
    call json%get('subsets', source_object, found)

    if (.not. found) call neko_error("No subsets found")

    this%n_zones = core%count(source_object)

    ! Allocate arrays if we found things
    if (this%n_zones .gt. 0) then
       allocate(this%zones(this%n_zones))
    end if

    ! First, count how many external zones we have (external = only "name",
    ! to be retrieved by the register later).
    do i = 1, this%n_zones
       ! Create a new json containing just the subdict for this source.
       call core%get_child(source_object, i, source_pointer, found)
       call core%print_to_string(source_pointer, buffer)
       call source_subdict%load_from_string(buffer)

       if (.not. source_subdict%valid_path("geometry")) then
          this%n_external_zones = this%n_external_zones + 1
       end if
    end do

    this%n_internal_zones = this%n_zones - this%n_external_zones
    if (this%n_external_zones .gt. 0) &
         allocate(this%names(this%n_external_zones))
    if (this%n_internal_zones .gt. 0) &
         allocate(this%internal_zones(this%n_internal_zones))

    i_internal = 1
    i_external = 1

    ! now go through everything again and either construct a point zone or
    ! save its name for the registry to fill it in later
    do i = 1, this%n_zones

       ! Create a new json containing just the subdict for this source.
       call core%get_child(source_object, i, source_pointer, found)
       call core%print_to_string(source_pointer, buffer)
       call source_subdict%load_from_string(buffer)

       if (source_subdict%valid_path("geometry")) then
          call point_zone_factory(this%internal_zones(i_internal)%pz, &
               source_subdict)
          call assign_point_zone(this%zones(i_internal)%pz, &
               this%internal_zones(i_internal)%pz)
          i_internal = i_internal + 1
       else
          call json_get(source_subdict, "name", type_name)
          this%names(i_external) = trim(type_name)
          i_external = i_external + 1
       end if

    end do

    ! Chcek that we got the proper operator
    call json_get(json, "operator", this%operator)
    select case (trim(this%operator))
    case ("OR")
    case ("AND")
    case ("XOR")
    case default
       call neko_error("Unknown operator " // trim(this%operator))
    end select

  end subroutine combine_point_zone_init_from_json

  subroutine assign_point_zone(pt, tgt)
    class(point_zone_t), intent(inout), pointer :: pt
    class(point_zone_t), intent(inout), target :: tgt

    pt => tgt

  end subroutine assign_point_zone

  !> Destructor.
  subroutine combine_point_zone_free(this)
    class(combine_point_zone_t), intent(inout) :: this

    integer :: i

    if (allocated(this%zones)) then
       do i = 1, this%n_zones
          nullify(this%zones(i)%pz)
       end do
       deallocate(this%zones)
    end if

    if (allocated(this%internal_zones)) then
       do i = 1, this%n_internal_zones
          call this%internal_zones(i)%pz%free
       end do
       deallocate(this%internal_zones)
    end if

    if (allocated(this%names)) deallocate(this%names)

    this%n_zones = 0
    this%n_internal_zones = 0
    this%n_external_zones = 0
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

    is_inside = this%zones(1)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%zones(1)%pz%invert

    do i = 2, this%n_zones
       select case (trim(this%operator))
       case ("OR")
          is_inside = is_inside .or. (this%zones(i)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%zones(i)%pz%invert)

       case ("AND")
          is_inside = is_inside .and. (this%zones(i)%pz%criterion(x, &
               y, z, j, k, l, e) .neqv. this%zones(i)%pz%invert)

       case ("XOR")
          is_inside = is_inside .neqv. (this%zones(i)%pz%criterion(x, &
               y, z, j, k, l, e).neqv. this%zones(i)%pz%invert)

       case default
       end select
    end do

  end function combine_point_zone_criterion

end module combine_point_zone
