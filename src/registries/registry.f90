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
module registry
  use, intrinsic :: iso_fortran_env, only: error_unit
  use num_types, only : rp
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t
  use registry_entry, only : registry_entry_t
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use htable, only : h_cptr_t
  use utils, only: neko_error
  use comm, only : pe_rank
  use json_module, only : json_file
  use json_utils, only : json_get
  use logger, only : neko_log, LOG_SIZE
  implicit none
  private

  type, public :: registry_t
     !> List of entries stored.
     type(registry_entry_t), private, allocatable :: entries(:)
     !> List of aliases to entries stored.
     type(json_file), private :: aliases
     !> Number of registered entries
     integer, private :: n_entries_ = 0
     !> Number of aliases.
     integer, private :: n_aliases_ = 0
     !> The size the entries array is increased by upon reallocation.
     integer, private :: expansion_size_ = 5
   contains
     !> Constructor.
     procedure, pass(this) :: init => registry_init
     !> Destructor.
     procedure, pass(this) :: free => registry_free
     !> Expand the array of entries so as to accommodate more entries.
     procedure, private, pass(this) :: expand => registry_expand

     !> Add a field to the registry.
     procedure, pass(this) :: add_field => registry_add_field
     !> Add a vector to the registry.
     procedure, pass(this) :: add_vector => registry_add_vector
     !> Add a matrix to the registry.
     procedure, pass(this) :: add_matrix => registry_add_matrix
     !> Add a scalar to the registry.
     procedure, pass(this) :: add_scalar => registry_add_scalar
     !> Add an alias to a field in the registry.
     procedure, pass(this) :: add_alias => registry_add_alias

     !> Get pointer to a stored field by index.
     procedure, pass(this) :: get_field_by_index => registry_get_field_by_index
     !> Get pointer to a stored vector by index.
     procedure, pass(this) :: get_vector_by_index => &
          registry_get_vector_by_index
     !> Get pointer to a stored matrix by index.
     procedure, pass(this) :: get_matrix_by_index => &
          registry_get_matrix_by_index
     !> Get pointer to a stored scalar by index.
     procedure, pass(this) :: get_scalar_by_index => &
          registry_get_scalar_by_index

     !> Get pointer to a stored field by name.
     procedure, pass(this) :: get_field_by_name => registry_get_field_by_name
     !> Get pointer to a stored vector by name.
     procedure, pass(this) :: get_vector_by_name => registry_get_vector_by_name
     !> Get pointer to a stored matrix by name.
     procedure, pass(this) :: get_matrix_by_name => registry_get_matrix_by_name
     !> Get pointer to a stored scalar by name.
     procedure, pass(this) :: get_scalar_by_name => registry_get_scalar_by_name

     !> Generic field getter
     generic :: get_field => get_field_by_index, get_field_by_name
     !> Generic vector getter
     generic :: get_vector => get_vector_by_index, get_vector_by_name
     !> Generic matrix getter
     generic :: get_matrix => get_matrix_by_index, get_matrix_by_name
     !> Generic scalar getter
     generic :: get_scalar => get_scalar_by_index, get_scalar_by_name

     !> Check if an entry with a given name is already in the registry.
     procedure, pass(this) :: entry_exists => registry_entry_exists
     !> Check if a field with a given name is already in the registry.
     procedure, pass(this) :: field_exists => registry_field_exists
     !> Check if a vector with a given name is already in the registry.
     procedure, pass(this) :: vector_exists => registry_vector_exists
     !> Check if a matrix with a given name is already in the registry.
     procedure, pass(this) :: matrix_exists => registry_matrix_exists
     !> Check if a scalar with a given name is already in the registry.
     procedure, pass(this) :: scalar_exists => registry_scalar_exists

     !> Get total allocated size of `fields`.
     procedure, pass(this) :: get_size => registry_get_size
     !> Get number of registered entries.
     procedure, pass(this) :: n_entries => registry_n_entries
     !> Get the number of fields in the registry.
     procedure, pass(this) :: n_fields => registry_n_fields
     !> Get the number of vectors in the registry.
     procedure, pass(this) :: n_vectors => registry_n_vectors
     !> Get the number of matrices in the registry.
     procedure, pass(this) :: n_matrices => registry_n_matrices
     !> Get the number of scalars in the registry.
     procedure, pass(this) :: n_scalars => registry_n_scalars
     !> Get the number of aliases in the registry.
     procedure, pass(this) :: n_aliases => registry_n_aliases
     !> Get the `expansion_size`
     procedure, pass(this) :: get_expansion_size => registry_get_expansion_size
     !> Print registry contents optionally filtered by type.
     procedure, pass(this) :: print_contents => registry_print_contents
  end type registry_t

  !> Global field registry
  type(registry_t), public, target :: neko_registry

contains
  ! ========================================================================== !
  ! Constructors/Destructors

  !> Constructor
  !! @param size The allocation size of `entries` on init.
  !! @param expansion_size The number of entries added to `entries` on expansion.
  subroutine registry_init(this, size, expansion_size)
    class(registry_t), intent(inout):: this
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
    end if

  end subroutine registry_init

  !> Destructor
  subroutine registry_free(this)
    class(registry_t), intent(inout):: this
    integer :: i

    if (allocated(this%entries)) then
       do i = 1, this%n_entries()
          call this%entries(i)%free()
       end do
       deallocate(this%entries)
    end if

    call this%aliases%destroy()

    this%n_entries_ = 0
    this%n_aliases_ = 0
    this%expansion_size_ = 5
  end subroutine registry_free

  !> Expand the fields array so as to accommodate more fields.
  subroutine registry_expand(this)
    class(registry_t), intent(inout) :: this
    type(registry_entry_t), allocatable :: temp(:)

    allocate(temp(this%n_entries_ + this%expansion_size_))
    temp(1:this%n_entries_) = this%entries(1:this%n_entries_)
    call move_alloc(temp, this%entries)
  end subroutine registry_expand

  ! ========================================================================== !
  ! Methods for adding objects to the registry

  !> Add a field to the registry.
  !! @param dof The map of degrees of freedom.
  !! @param name The name of the field.
  !! @param ignore_existing If true, will do nothing if the field is already in
  !! the registry. If false, will throw an error. Optional, defaults to false.
  subroutine registry_add_field(this, dof, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%field_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Field with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! initialize the field at the appropriate index
    call this%entries(this%n_entries_)%init_field(dof, name)

  end subroutine registry_add_field

  !> Add a vector to the registry.
  !! @param n The size of the vector.
  !! @param name The name of the vector.
  !! @param ignore_existing If true, will do nothing if the vector is already in
  !! the registry. If false, will throw an error. Optional, defaults to false.
  subroutine registry_add_vector(this, n, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    integer, intent(in) :: n
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%vector_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Vector with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named vector at the appropriate index
    call this%entries(this%n_entries_)%init_vector(n, name)

  end subroutine registry_add_vector

  !> Add a matrix to the registry.
  !! @param n The size of the matrix.
  !! @param name The name of the matrix.
  !! @param ignore_existing If true, will do nothing if the matrix is already in
  !! the registry. If false, will throw an error. Optional, defaults to false.
  subroutine registry_add_matrix(this, nrows, ncols, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    integer, intent(in) :: nrows, ncols
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%matrix_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Vector with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named matrix at the appropriate index
    call this%entries(this%n_entries_)%init_matrix(nrows, ncols, name)

  end subroutine registry_add_matrix

  !> Add a scalar to the registry.
  !! @param value The scalar value.
  !! @param name The name of the scalar.
  !! @param ignore_existing If true, skip if scalar already registered.
  subroutine registry_add_scalar(this, value, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    real(kind=rp), intent(in) :: value
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%scalar_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Scalar with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named scalar at the appropriate index
    call this%entries(this%n_entries_)%init_scalar(value, name)

  end subroutine registry_add_scalar

  !> Add an alias for an existing entry in the registry.
  !! @param alias The alias.
  !! @param name The name of the entry.
  subroutine registry_add_alias(this, alias, name)
    class(registry_t), intent(inout) :: this
    character(len=*), intent(in) :: alias
    character(len=*), intent(in) :: name

    if (this%entry_exists(alias)) then
       call neko_error("Cannot create alias. Entry " // alias // &
            " already exists in the registry")
    end if

    if (this%entry_exists(name)) then
       this%n_aliases_ = this%n_aliases_ + 1
       call this%aliases%add(trim(alias), trim(name))
    else
       call neko_error("Cannot create alias. Entry " // name // &
            " could not be found in the registry")
    end if
  end subroutine registry_add_alias

  ! ========================================================================== !
  ! Methods for retrieving objects from the registry by index

  !> Get pointer to a stored field by index.
  function registry_get_field_by_index(this, i) result(f)
    class(registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(field_t), pointer :: f
    character(len=:), allocatable :: buffer

    if (i < 1) then
       call neko_error("Field index must be > 1")
    else if (i > this%n_entries()) then
       call neko_error("Field index exceeds number of stored fields")
    endif

    if (this%entries(i)%get_type() .ne. 'field') then
       write(buffer, *) "Requested index ", i, " is not a field, but a ", &
            this%entries(i)%get_type()
       call neko_error(buffer)
    end if

    f => this%entries(i)%get_field()
  end function registry_get_field_by_index

  !> Get pointer to a stored vector by index.
  function registry_get_vector_by_index(this, i) result(f)
    class(registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(vector_t), pointer :: f
    character(len=:), allocatable :: buffer

    if (i < 1) then
       call neko_error("Vector index must be > 1")
    else if (i > this%n_entries()) then
       call neko_error("Vector index exceeds number of stored vectors")
    endif

    if (this%entries(i)%get_type() .ne. 'vector') then
       write(buffer, *) "Requested index ", i, " is not a vector, but a ", &
            this%entries(i)%get_type()
       call neko_error(buffer)
    end if

    f => this%entries(i)%get_vector()
  end function registry_get_vector_by_index

  !> Get pointer to a stored matrix by index.
  function registry_get_matrix_by_index(this, i) result(f)
    class(registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(matrix_t), pointer :: f
    character(len=:), allocatable :: buffer

    if (i < 1) then
       call neko_error("Matrix index must be > 1")
    else if (i > this%n_entries()) then
       call neko_error("Matrix index exceeds number of stored matrices")
    endif

    if (this%entries(i)%get_type() .ne. 'matrix') then
       write(buffer, *) "Requested index ", i, " is not a matrix, but a ", &
            this%entries(i)%get_type()
       call neko_error(buffer)
    end if

    f => this%entries(i)%get_matrix()
  end function registry_get_matrix_by_index

  !> Get pointer to a stored scalar by index.
  function registry_get_scalar_by_index(this, i) result(s)
    class(registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp), pointer :: s
    character(len=:), allocatable :: buffer

    if (i < 1) then
       call neko_error("Scalar index must be > 1")
    else if (i > this%n_entries()) then
       call neko_error("Scalar index exceeds number of stored scalars")
    endif

    if (this%entries(i)%get_type() .ne. 'scalar') then
       write(buffer, *) "Requested index ", i, " is not a scalar, but a ", &
            this%entries(i)%get_type()
       call neko_error(buffer)
    end if

    s => this%entries(i)%get_scalar()
  end function registry_get_scalar_by_index

  ! ========================================================================== !
  ! Methods for retrieving objects from the registry by name

  !> Get pointer to a stored field by field name.
  recursive function registry_get_field_by_name(this, name) result(f)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(field_t), pointer :: f
    logical :: found
    integer :: i

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'field' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%get_field()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_field_by_name(alias_target)
       return
    end if

    if (pe_rank .eq. 0) then
       write(error_unit, *) "Current registry contents:"

       do i = 1, this%n_entries()
          write(error_unit, *) "- ", this%entries(i)%get_name()
       end do
    end if
    call neko_error("Field " // name // " could not be found in the registry")

  end function registry_get_field_by_name


  !> Get pointer to a stored vector by name.
  recursive function registry_get_vector_by_name(this, name) result(f)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(vector_t), pointer :: f
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'vector' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%get_vector()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_vector_by_name(alias_target)
       return
    end if

    if (pe_rank .eq. 0) then
       write(error_unit, *) "Current registry contents:"

       do i = 1, this%n_entries()
          write(error_unit, *) "- ", this%entries(i)%get_name()
       end do
    end if
    call neko_error("Vector " // name // " could not be found in the registry")

  end function registry_get_vector_by_name

  !> Get pointer to a stored matrix by name.
  recursive function registry_get_matrix_by_name(this, name) result(f)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(matrix_t), pointer :: f
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'matrix' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%get_matrix()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_matrix_by_name(alias_target)
       return
    end if

    if (pe_rank .eq. 0) then
       write(error_unit, *) "Current registry contents:"

       do i = 1, this%n_entries()
          write(error_unit, *) "- ", this%entries(i)%get_name()
       end do
    end if
    call neko_error("Matrix " // name // " could not be found in the registry")

  end function registry_get_matrix_by_name

  !> Get pointer to a stored scalar by name.
  recursive function registry_get_scalar_by_name(this, name) result(s)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    real(kind=rp), pointer :: s
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'scalar' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          s => this%entries(i)%get_scalar()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       s => this%get_scalar_by_name(alias_target)
       return
    end if

    if (pe_rank .eq. 0) then
       write(error_unit, *) "Current registry contents:"

       do i = 1, this%n_entries()
          write(error_unit, *) "- ", this%entries(i)%get_name()
       end do
    end if
    call neko_error("Scalar " // name // " could not be found in the registry")

  end function registry_get_scalar_by_name

  ! ========================================================================== !
  ! Methods for checking existence of objects in the registry

  !> Check if a field with a given name is already in the registry.
  function registry_entry_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (trim(this%entries(i)%get_name()) .eq. trim(name)) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_entry_exists

  !> Check if a field with a given name is already in the registry.
  function registry_field_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'field' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_field_exists

  !> Check if a vector with a given name is already in the registry.
  function registry_vector_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'vector' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_vector_exists

  !> Check if a matrix with a given name is already in the registry.
  function registry_matrix_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'matrix' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_matrix_exists

  !> Check if a scalar with a given name is already in the registry.
  function registry_scalar_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'scalar' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          found = .true.
          return
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_scalar_exists

  ! ========================================================================== !
  ! Generic component accessor methods

  !> Get number of registered entries.
  pure function registry_n_entries(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n

    n = this%n_entries_
  end function registry_n_entries

  !> Get the number of fields stored in the registry
  pure function registry_n_fields(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'field') then
          n = n + 1
       end if
    end do
  end function registry_n_fields

  !> Get the number of vector stored in the registry
  pure function registry_n_vectors(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'vector') then
          n = n + 1
       end if
    end do
  end function registry_n_vectors

  !> Get the number of matrix stored in the registry
  pure function registry_n_matrices(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'matrix') then
          n = n + 1
       end if
    end do
  end function registry_n_matrices

  !> Get the number of scalars stored in the registry
  pure function registry_n_scalars(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'scalar') then
          n = n + 1
       end if
    end do
  end function registry_n_scalars

  !> Get the number of aliases stored in the registry
  pure function registry_n_aliases(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n

    n = this%n_aliases_
  end function registry_n_aliases

  !> Get the size of the fields array.
  pure function registry_get_size(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n

    if (allocated(this%entries)) then
       n = size(this%entries)
    else
       n = 0
    end if
  end function registry_get_size

  !> Get the expansion size.
  pure function registry_get_expansion_size(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size_
  end function registry_get_expansion_size

  !> Print the contents of the registry to standard output.
  subroutine registry_print(this)
    class(registry_t), intent(in) :: this
    character(len=LOG_SIZE), allocatable :: buffer
    integer :: i

    call neko_log%section("Field Registry Contents")
    do i = 1, this%n_entries()
       write(buffer, '(A,I4,A,A)') "- [", i, "] ", &
            this%entries(i)%get_type(), ": ", this%entries(i)%get_name()
       call neko_log%message(trim(buffer))
    end do

    call neko_log%end_section()
  end subroutine registry_print

  !> Print the registry contents grouped by entity type.
  subroutine registry_print_contents(this, type)
    class(registry_t), intent(in) :: this
    character(len=*), optional, intent(in) :: type
    character(len=:), allocatable :: filter_type
    character(len=6), parameter :: types(4) = (/ 'field ', 'vector', 'matrix', &
         'scalar' /)
    logical :: filter_active
    integer :: i
    logical :: known_type

    filter_active = .false.
    if (present(type)) then
       filter_type = trim(type)
       filter_active = .true.
       known_type = .false.
       do i = 1, size(types)
          if (filter_type == types(i)) then
             known_type = .true.
             exit
          end if
       end do
       if (.not. known_type) then
          call neko_error("registry::print_contents: Unsupported type " &
               // trim(filter_type))
       end if
    end if

    call neko_log%section("Registry Contents")
    do i = 1, size(types)
       if (filter_active .and. (filter_type .ne. types(i))) cycle
       call registry_print_section(this, types(i))
    end do
    call neko_log%end_section()
  end subroutine registry_print_contents

  !> Print a single section of the registry for the given type.
  subroutine registry_print_section(this, entity_type)
    class(registry_t), intent(in) :: this
    character(len=*), intent(in) :: entity_type
    integer :: i
    logical :: found
    character(len=LOG_SIZE) :: buffer

    call neko_log%message("  "//trim(entity_type)//" entries:")
    found = .false.
    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. entity_type) then
          found = .true.
          write(buffer, '(A,I4,A,A)') "    [", i, "] ", &
               trim(this%entries(i)%get_name())
          call neko_log%message(trim(buffer))
       end if
    end do
    if (.not. found) then
       call neko_log%message("    <none>")
    end if
  end subroutine registry_print_section

end module registry
