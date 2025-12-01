! Copyright (c) 2018-2025, The Neko Authors
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
!> Defines a registry for storing solution fields
module field_registry
  use, intrinsic :: iso_fortran_env, only: error_unit
  use field, only : field_t
  use registry_entry, only : registry_entry_t
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use htable, only : h_cptr_t
  use utils, only: neko_error
  use comm, only : pe_rank
  use json_module, only : json_file
  use json_utils, only : json_get
  implicit none
  private

  type, public :: field_registry_t
     !> List of fields stored.
     type(registry_entry_t), private, allocatable :: entries(:)
     !> List of aliases to fields stored.
     type(json_file), private :: aliases
     !> Number of registered fields.
     integer, private :: n_fields_ = 0
     !> Number of registered entries
     integer, private :: n_entries_ = 0
     !> Number of aliases.
     integer, private :: n_aliases_ = 0
     !> The size the fields array is increased by upon reallocation.
     integer, private :: expansion_size_
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_registry_init
     !> Destructor.
     procedure, pass(this) :: free => field_registry_free
     !> Expand the fields array so as to accommodate more fields.
     procedure, private, pass(this) :: expand => field_registry_expand


     !> Add a field to the registry.
     procedure, pass(this) :: add_field => field_registry_add_field
     !> Add an alias to a field in the registry.
     procedure, pass(this) :: add_alias => field_registry_add_alias

     !> Get pointer to a stored field by index.
     procedure, pass(this) :: get_field_by_index => registry_get_field_by_index
     !> Get pointer to a stored field by name.
     procedure, pass(this) :: get_field_by_name => registry_get_field_by_name

     !> Check if a field with a given name is already in the registry.
     procedure, pass(this) :: field_exists => field_registry_field_exists

     !> Generic field getter
     generic :: get_field => get_field_by_index, get_field_by_name

     !> Get total allocated size of `fields`.
     procedure, pass(this) :: get_size => field_registry_get_size
     !> Get number of registered entries.
     procedure, pass(this) :: n_entries => field_registry_n_entries
     !> Get the number of fields in the registry.
     procedure, pass(this) :: n_fields => field_registry_n_fields
     !> Get the number of aliases in the registry.
     procedure, pass(this) :: n_aliases => field_registry_n_aliases
     !> Get the `expansion_size`
     procedure, pass(this) :: get_expansion_size => registry_get_expansion_size
  end type field_registry_t

  !> Global field registry
  type(field_registry_t), public, target :: neko_field_registry

contains
  ! ========================================================================== !
  ! Constructors/Destructors

  !> Constructor
  !! @param size The allocation size of `fields` on init.
  !! @param expansion_size The number of entries added to `fields` on expansion.
  subroutine field_registry_init(this, size, expansion_size)
    class(field_registry_t), intent(inout):: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size

    call this%free()

    if (present(size)) then
       allocate(this%entries(size))
    else
       allocate(this%entries(25))
    end if

    call this%aliases%initialize()

    if (present(expansion_size)) then
       this%expansion_size_ = expansion_size
    else
       this%expansion_size_ = 5
    end if

  end subroutine field_registry_init

  !> Destructor
  subroutine field_registry_free(this)
    class(field_registry_t), intent(inout):: this
    integer :: i
    if (allocated(this%entries)) then
       do i = 1, this%n_entries()
          call this%entries(i)%free()
       end do
       deallocate(this%entries)
    end if

    call this%aliases%destroy()

    this%n_fields_ = 0
    this%n_entries_ = 0
    this%n_aliases_ = 0
    this%expansion_size_ = 0
  end subroutine field_registry_free

  !> Expand the fields array so as to accommodate more fields.
  subroutine field_registry_expand(this)
    class(field_registry_t), intent(inout) :: this
    type(registry_entry_t), allocatable :: temp(:)

    allocate(temp(this%n_fields_ + this%expansion_size_))
    temp(1:this%n_fields_) = this%entries(1:this%n_fields_)
    call move_alloc(temp, this%entries)
  end subroutine field_registry_expand

  ! ========================================================================== !
  ! Methods for adding objects to the registry

  !> Add a field to the registry.
  !! @param dof The map of degrees of freedom.
  !! @param fld_name The name of the field.
  !! @param ignore_existing If true, will do nothing if the field is already in
  !! the registry. If false, will throw an error. Optional, defaults to false.
  subroutine field_registry_add_field(this, dof, fld_name, ignore_existing)
    class(field_registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), target, intent(in) :: fld_name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%field_exists(fld_name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Field with name " // fld_name // &
               " is already registered")
       end if
    end if

    if (this%n_fields() .eq. size(this%entries)) then
       call this%expand()
    end if

    this%n_fields_ = this%n_fields_ + 1

    ! initialize the field at the appropriate index
    call this%entries(this%n_fields_)%init_field(dof, fld_name)

  end subroutine field_registry_add_field

  !> Add an alias for an existing field in the registry.
  !! @param alias The alias.
  !! @param fld_name The name of the field.
  subroutine field_registry_add_alias(this, alias, fld_name)
    class(field_registry_t), intent(inout) :: this
    character(len=*), intent(in) :: alias
    character(len=*), intent(in) :: fld_name

    if (this%field_exists(alias)) then
       call neko_error("Cannot create alias. Field " // alias // &
            " already exists in the registry")
    end if

    if (this%field_exists(fld_name)) then
       this%n_aliases_ = this%n_aliases_ + 1
       call this%aliases%add(trim(alias), trim(fld_name))
    else
       call neko_error("Cannot create alias. Field " // fld_name // &
            " could not be found in the registry")
    end if
  end subroutine field_registry_add_alias

  ! ========================================================================== !
  ! Methods for retrieving objects from the registry

  !> Get pointer to a stored field by index.
  function registry_get_field_by_index(this, i) result(f)
    class(field_registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(field_t), pointer :: f

    if (i < 1) then
       call neko_error("Field index must be > 1")
    else if (i > this%n_fields()) then
       call neko_error("Field index exceeds number of stored fields")
    endif

    f => this%entries(i)%field_ptr
  end function registry_get_field_by_index

  !> Get pointer to a stored field by field name.
  recursive function registry_get_field_by_name(this, name) result(f)
    class(field_registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(field_t), pointer :: f
    logical :: found, is_alias
    integer :: i

    do i = 1, this%n_fields()
       if (this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%field_ptr
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_field_by_name(alias_target)
       return
    end if

    if (pe_rank .eq. 0) then
       write(error_unit, *) "Current field_registry contents:"

       do i = 1, this%n_fields()
          write(error_unit, *) "- ", this%entries(i)%get_name()
       end do
    end if
    call neko_error("Field " // name // " could not be found in the registry")

  end function registry_get_field_by_name

  !> Check if a field with a given name is already in the registry.
  function field_registry_field_exists(this, name) result(found)
    class(field_registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_fields()
       if (this%entries(i)%get_name() .eq. name) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function field_registry_field_exists

  ! ========================================================================== !
  ! Generic component accessor methods

  !> Get number of registered entries.
  pure function field_registry_n_entries(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%n_entries_
  end function field_registry_n_entries

  !> Get the number of fields stored in the registry
  pure function field_registry_n_fields(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%n_fields_
  end function field_registry_n_fields

  !> Get the number of aliases stored in the registry
  pure function field_registry_n_aliases(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%n_aliases_
  end function field_registry_n_aliases

  !> Get the size of the fields array.
  pure function field_registry_get_size(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = size(this%entries)
  end function field_registry_get_size

  !> Get the expansion size.
  pure function registry_get_expansion_size(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size_
  end function registry_get_expansion_size

end module field_registry
