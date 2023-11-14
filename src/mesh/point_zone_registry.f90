! Copyright (c) 2019-2021, The Neko Authors
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
! Implements a point zone registry for storing point zones.
module point_zone_registry
  use point_zone, only : point_zone_t, point_zone_wrapper_t
  use point_zone_fctry, only: point_zone_factory
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use utils, only: neko_error
  use json_utils, only: json_get
  use json_module, only: json_file, json_core, json_value
  implicit none
  private

  type :: point_zone_registry_t
     !> List of point_zones stored.
     type(point_zone_wrapper_t), allocatable :: point_zones(:)
     !> Number of registered point_zones.
     integer, private :: n = 0
     !> The size the point_zones array is increased by upon reallocation.
     integer, private :: expansion_size
   contains
     !> Expand the point_zones array so as to accomodate more point_zones.
     procedure, private, pass(this) :: expand
     !> Constructor, reading from json point zones.
     procedure, pass(this) :: init => point_zone_registry_init
     !> Destructor.
     procedure, pass(this) :: free => point_zone_registry_free
     !> Adds a point zone object to the registry from a json object.
     procedure, pass(this) :: add_point_zone_from_json
     !> Returns the number of point zones in the registry.
     procedure, pass(this) :: n_point_zones
     !> Retrieves a point zone in the registry by its index in the 
     !! `point_zones` array.
     procedure, pass(this) :: get_point_zone_by_index
     !> Retrieves a point zone in the registry by its name.
     procedure, pass(this) :: get_point_zone_by_name
     !> Returns the expansion size with which the `point_zone_registry_t`
     !! was initialized.
     procedure, pass(this) :: get_expansion_size
     !> Returns the total size of the `point_zones` array (not the number of 
     !! point zones in the registry!).
     procedure, pass(this) :: get_size
     !> Checks if a point zone exists in the registry.
     procedure, pass(this) :: point_zone_exists
     generic :: get_point_zone => get_point_zone_by_index, get_point_zone_by_name
     generic :: add_point_zone => add_point_zone_from_json
  end type point_zone_registry_t

  !> Global point_zone registry
  type(point_zone_registry_t), public, target :: neko_point_zone_registry

contains
  !> Constructor, reading from json point zones.
  !! @param json Json file object.
  !! @param dof Dofmap to map the point zone from GLL points.
  !! @param size Size of the point zone registry.
  !! @param expansion_size Expansion size for the point zone registry.
  !! @note At this stage, the point_zone registry is only allocated
  !! if we find anything in the `case.point_zones` json path. Any
  !! point_zones that are not defined in that way will need to be added
  !! using the `add_point_zone` subroutine.
  subroutine point_zone_registry_init(this, json, dof, expansion_size)
    class(point_zone_registry_t), intent(inout):: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(inout) :: dof
    integer, optional, intent(in) :: expansion_size

    ! Json low-level manipulator.
    type(json_core) :: core
    ! Pointer to the source_terms JSON object and the individual sources.
    type(json_value), pointer :: source_object, source_pointer
    ! Buffer for serializing the json.
    character(len=:), allocatable :: buffer
    ! A single source term as its own json_file.
    type(json_file) :: source_subdict
    character(len=:), allocatable :: type
    logical :: found
    integer :: n_zones, i

    call this%free()

     if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 10
    end if

    this%n = 0

    !
    ! Count if there are any point zones defined in the json
    !
    if(json%valid_path('case.point_zones')) then

       call json%get_core(core)
       call json%get('case.point_zones', source_object, found)

       n_zones = core%count(source_object)
       this%n = n_zones

       allocate(this%point_zones(n_zones))

       ! Initialize every point zone
       do i = 1, n_zones
          ! Create a new json containing just the subdict for this source.
          call core%get_child(source_object, i, source_pointer, found)
          call core%print_to_string(source_pointer, buffer)
          call source_subdict%load_from_string(buffer)

          call point_zone_factory(this%point_zones(i)%pz, source_subdict, dof)
       end do
    end if

  end subroutine point_zone_registry_init

  !> Destructor.
  subroutine point_zone_registry_free(this)
    class(point_zone_registry_t), intent(inout):: this
    integer :: i

    if (allocated(this%point_zones)) then

       do i=1, this%n_point_zones()
          call this%point_zones(i)%pz%free()
       end do

       deallocate(this%point_zones)

    end if

    this%n = 0
    this%expansion_size = 0
  end subroutine point_zone_registry_free

  !> Expand the point_zones array so as to accomodate more point_zones.
  subroutine expand(this)
    class(point_zone_registry_t), intent(inout) :: this
    type(point_zone_wrapper_t), allocatable :: temp(:)
    integer :: i

    allocate(temp(this%n + this%expansion_size))
    temp(1:this%n) = this%point_zones(1:this%n)
    call move_alloc(temp, this%point_zones)

  end subroutine expand

  !> Adds a point zone object to the registry from a json object.
  !! @param json Json object from which to initialize the point zone.
  !! @param dof Dofmap from which to map the point zone.
  subroutine add_point_zone_from_json(this, json, dof)
    class(point_zone_registry_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), target, intent(inout) :: dof
!    type(h_cptr_t) :: key
    character(len=:), allocatable :: str_read
    integer :: i

    !
    ! Allocate the point zones array as it was not necessarily done
    ! in init
    !
    if (.not. allocated(this%point_zones)) then
       allocate(this%point_zones(this%expansion_size))
    end if

    call json_get(json, "name", str_read)

    ! Check if point zone exists with the input name
    if (this%point_zone_exists(trim(str_read))) then
       call neko_error("Field with name " // trim(str_read) // &
            " is already registered")
    end if

    !
    ! This will always be true if point_zones was allocated in
    ! init.
    !
    if (this%n_point_zones() .eq. this%get_size()) then
      call this%expand()
    end if

    this%n = this%n + 1

    ! initialize the point_zone at the appropriate index
    call point_zone_factory(this%point_zones(this%n)%pz, json, dof)

    ! generate a key for the name lookup map and assign it to the index
    !    key%ptr = c_loc(fld_name)
    !    call this%name_index_map%set(key, this%n)

    !    write(*,*) "HTABLE DATA, ", this%name_index_map%get(key, i)
  end subroutine add_point_zone_from_json

  !> Returns the number of point zones in the registry.
  pure function n_point_zones(this) result(n)
    class(point_zone_registry_t), intent(in) :: this
    integer :: n

    n = this%n
  end function n_point_zones

  !> Returns the total size of the `point_zones` array (not the number of 
  !! point zones in the registry!).
  !! @note Use `n_point_zones()` to retrieve the actual number of point
  !! zones in the registry.
  pure function get_size(this) result(n)
    class(point_zone_registry_t), intent(in) :: this
    integer :: n

    n = size(this%point_zones)
  end function get_size

  !> Returns the expansion size with which the `point_zone_registry_t`
  !! was initialized.
  pure function get_expansion_size(this) result(n)
    class(point_zone_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Retrieves a point zone in the registry by its index in the 
  !! `point_zones` array.
  !! @param i Index in the `point_zones` array.
  function get_point_zone_by_index(this, i) result(pz)
    class(point_zone_registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    class(point_zone_t), pointer :: pz

    if (i < 1) then
       call neko_error("Field index must be > 1")
    else if (i > this%n_point_zones()) then
       call neko_error("Field index exceeds number of stored point_zones")
    endif

    pz => this%point_zones(i)%pz
  end function get_point_zone_by_index

  !> Retrieves a point zone in the registry by its name.
  !! @param name Name of the point zone.
  function get_point_zone_by_name(this, name) result(pz)
    class(point_zone_registry_t), target, intent(in) :: this
    character(len=*), intent(in) :: name
    class(point_zone_t), pointer :: pz
    logical :: found
    integer :: i

    found = .false.
    do i=1, this%n_point_zones()
       if (trim(this%point_zones(i)%pz%name) .eq. trim(name)) then
          pz => this%point_zones(i)%pz
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       call neko_error("Point zone " // trim(name) // &
            " could not be found in the registry")
    end if
  end function get_point_zone_by_name

  !> Checks if a point zone exists in the registry.
  !! @param name Name of the point zone.
  function point_zone_exists(this, name) result(found)
    class(point_zone_registry_t), target, intent(in) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i=1, this%n_point_zones()
       if (trim(this%point_zones(i)%pz%name) .eq. trim(name)) then
          found = .true.
          exit
       end if
    end do
  end function point_zone_exists

end module point_zone_registry
