! Copyright (c) 2018-2026, The Neko Authors
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
  use num_types, only : rp
  use field, only : field_t
  use vector, only : vector_t
  use matrix, only : matrix_t
  use tensor3, only : tensor3_t
  use tensor4, only : tensor4_t
  use registry_entry, only : registry_entry_t
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use json_module, only : json_file
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
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
     !> Add a tensor3 to the registry.
     procedure, pass(this) :: add_tensor3 => registry_add_tensor3
     !> Add a tensor4 to the registry.
     procedure, pass(this) :: add_tensor4 => registry_add_tensor4
     !> Add a real scalar to the registry.
     procedure, pass(this) :: add_real_scalar => registry_add_real_scalar
     !> Add an integer scalar to the registry.
     procedure, pass(this) :: add_integer_scalar => registry_add_integer_scalar
     !> Add an alias to a field in the registry.
     procedure, pass(this) :: add_alias => registry_add_alias

     !> Get pointer to a stored field by name.
     procedure, pass(this) :: get_field => registry_get_field
     !> Get pointer to a stored vector by name.
     procedure, pass(this) :: get_vector => registry_get_vector
     !> Get pointer to a stored matrix by name.
     procedure, pass(this) :: get_matrix => registry_get_matrix
     !> Get pointer to a stored tensor3 by name.
     procedure, pass(this) :: get_tensor3 => registry_get_tensor3
     !> Get pointer to a stored tensor4 by name.
     procedure, pass(this) :: get_tensor4 => registry_get_tensor4
     !> Get pointer to a stored real scalar by name.
     procedure, pass(this) :: get_real_scalar => registry_get_real_scalar
     !> Get pointer to a stored integer scalar by name.
     procedure, pass(this) :: get_integer_scalar => registry_get_integer_scalar

     ! Just to retain the old API for backwards compatibility.
     generic :: get_field_by_name => get_field
     generic :: get_vector_by_name => get_vector
     generic :: get_matrix_by_name => get_matrix
     generic :: get_tensor3_by_name => get_tensor3
     generic :: get_tensor4_by_name => get_tensor4
     generic :: get_real_scalar_by_name => get_real_scalar
     generic :: get_integer_scalar_by_name => get_integer_scalar

     !> Check if an entry with a given name is already in the registry.
     procedure, pass(this) :: entry_exists => registry_entry_exists
     !> Check if a field with a given name is already in the registry.
     procedure, pass(this) :: field_exists => registry_field_exists
     !> Check if a vector with a given name is already in the registry.
     procedure, pass(this) :: vector_exists => registry_vector_exists
     !> Check if a matrix with a given name is already in the registry.
     procedure, pass(this) :: matrix_exists => registry_matrix_exists
     !> Check if a tensor3 with a given name is already in the registry.
     procedure, pass(this) :: tensor3_exists => registry_tensor3_exists
     !> Check if a tensor4 with a given name is already in the registry.
     procedure, pass(this) :: tensor4_exists => registry_tensor4_exists
     !> Check if a real scalar with a given name is already in the registry.
     procedure, pass(this) :: real_scalar_exists => registry_real_scalar_exists
     !> Check if an integer scalar with a given name is already in the registry.
     procedure, pass(this) :: integer_scalar_exists => &
          registry_integer_scalar_exists
     !> Backwards compatible scalar existence check (real).
     procedure, pass(this) :: scalar_exists => registry_real_scalar_exists

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
     !> Get the number of tensor3s in the registry.
     procedure, pass(this) :: n_tensor3s => registry_n_tensor3s
     !> Get the number of tensor4s in the registry.
     procedure, pass(this) :: n_tensor4s => registry_n_tensor4s
     !> Get the number of real scalars in the registry.
     procedure, pass(this) :: n_real_scalars => registry_n_real_scalars
     !> Get the number of integer scalars in the registry.
     procedure, pass(this) :: n_integer_scalars => registry_n_integer_scalars
     !> Backwards compatible scalar count (real).
     procedure, pass(this) :: n_scalars => registry_n_real_scalars
     !> Get the number of aliases in the registry.
     procedure, pass(this) :: n_aliases => registry_n_aliases
     !> Get the `expansion_size`
     procedure, pass(this) :: get_expansion_size => registry_get_expansion_size
     !> Print registry contents optionally filtered by type.
     procedure, pass(this) :: print_contents => registry_print_contents
  end type registry_t

  !> Global field registry
  type(registry_t), public, target :: neko_registry

  !> This registry is used to store user-defined scalars and vectors, provided
  !! under the `constants` section of the case file. These are separated
  !! from the global registry to prevent name clashes with registered objects
  !! used by Neko itself.
  type(registry_t), public, target :: neko_const_registry

contains
  ! ========================================================================== !
  ! Constructors/Destructors

  !> Constructor
  !! @param size The allocation size of `entries` on init.
  !! @param expansion_size The number of entries added to `entries` on
  !! expansion.
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
    integer :: n, i

    n = this%get_size()

    if (n .gt. 0) then
       call move_alloc(this%entries, temp)
    end if

    allocate(this%entries(n + this%expansion_size_))

    if (n .gt. 0) then
       do i = 1, n
          call this%entries(i)%move_from(temp(i))
          call temp(i)%free()
       end do
    end if

    if (allocated(temp)) deallocate(temp)
  end subroutine registry_expand

  ! ========================================================================== !
  ! Methods for adding objects to the registry

  !> Add a field to the registry.
  !! @param dof The map of degrees of freedom.
  !! @param name The name of the field.
  !! @param ignore_existing If true, will do nothing if the field is already
  !! in the registry. If false, will throw an error. Optional, default is false.
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

    call neko_log%message("Field " // trim(name) // " added to the registry", &
         lvl=NEKO_LOG_DEBUG)

  end subroutine registry_add_field

  !> Add a vector to the registry.
  !! @param n The size of the vector.
  !! @param name The name of the vector.
  !! @param ignore_existing If true, will do nothing if the vector is already
  !! in the registry. If false, will throw an error. Optional, default is false.
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

    call neko_log%message("Vector " // trim(name) // " added to the registry", &
         lvl=NEKO_LOG_DEBUG)

  end subroutine registry_add_vector

  !> Add a matrix to the registry.
  !! @param n The size of the matrix.
  !! @param name The name of the matrix.
  !! @param ignore_existing If true, will do nothing if the matrix is already
  !! in the registry. If false, will throw an error. Optional, default is false.
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
          call neko_error("Matrix with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named matrix at the appropriate index
    call this%entries(this%n_entries_)%init_matrix(nrows, ncols, name)

    call neko_log%message("Matrix " // trim(name) // " added to the registry", &
         lvl=NEKO_LOG_DEBUG)

  end subroutine registry_add_matrix

  !> Add a tensor3 to the registry.
  !! @param n The size of the tensor3.
  !! @param name The name of the tensor3.
  !! @param ignore_existing If true, will do nothing if the tensor3 is already
  !! in the registry. If false, will throw an error. Optional, default is false.
  subroutine registry_add_tensor3(this, n, m, k, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    integer, intent(in) :: n, m, k
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%tensor3_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Tensor3 with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named tensor3 at the appropriate index
    call this%entries(this%n_entries_)%init_tensor3(n, m, k, name)

    call neko_log%message("Tensor3 " // trim(name) // &
         " added to the registry", lvl=NEKO_LOG_DEBUG)

  end subroutine registry_add_tensor3

  !> Add a tensor4 to the registry.
  !! @param n The size of the tensor4.
  !! @param name The name of the tensor4.
  !! @param ignore_existing If true, will do nothing if the tensor4 is already
  !! in the registry. If false, will throw an error. Optional, default is false.
  subroutine registry_add_tensor4(this, n, m, k, l, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    integer, intent(in) :: n, m, k, l
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%tensor4_exists(name)) then
       if (ignore_existing_) then
          return
       else
          call neko_error("Tensor4 with name " // name // &
               " is already registered")
       end if
    end if

    if (this%n_entries() .eq. this%get_size()) then
       call this%expand()
    end if

    this%n_entries_ = this%n_entries_ + 1

    ! Initialize the named tensor4 at the appropriate index
    call this%entries(this%n_entries_)%init_tensor4(n, m, k, l, name)

    call neko_log%message("Tensor4 " // trim(name) // &
         " added to the registry", lvl=NEKO_LOG_DEBUG)

  end subroutine registry_add_tensor4

  !> Add a real scalar to the registry.
  !! @param value The scalar value.
  !! @param name The name of the scalar.
  !! @param ignore_existing If true, skip if scalar already registered.
  subroutine registry_add_real_scalar(this, value, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    real(kind=rp), intent(in) :: value
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%real_scalar_exists(name)) then
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
    call this%entries(this%n_entries_)%init_real_scalar(value, name)

  end subroutine registry_add_real_scalar

  !> Add an integer scalar to the registry.
  !! @param value The scalar value.
  !! @param name The name of the scalar.
  !! @param ignore_existing If true, skip if scalar already registered.
  subroutine registry_add_integer_scalar(this, value, name, ignore_existing)
    class(registry_t), intent(inout) :: this
    integer, intent(in) :: value
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing)) then
       ignore_existing_ = ignore_existing
    end if

    if (this%integer_scalar_exists(name)) then
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
    call this%entries(this%n_entries_)%init_integer_scalar(value, name)

  end subroutine registry_add_integer_scalar

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
  ! Methods for retrieving objects from the registry by name

  !> Get pointer to a stored field by field name.
  recursive function registry_get_field(this, name) result(f)
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
       f => this%get_field(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Field " // name // " could not be found in the registry")

  end function registry_get_field


  !> Get pointer to a stored vector by name.
  recursive function registry_get_vector(this, name) result(f)
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
       f => this%get_vector(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Vector " // name // " could not be found in the registry")

  end function registry_get_vector

  !> Get pointer to a stored matrix by name.
  recursive function registry_get_matrix(this, name) result(f)
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
       f => this%get_matrix(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Matrix " // name // " could not be found in the registry")

  end function registry_get_matrix

  !> Get pointer to a stored tensor3 by name.
  recursive function registry_get_tensor3(this, name) result(f)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(tensor3_t), pointer :: f
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'tensor3' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%get_tensor3()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_tensor3(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Tensor3 " // name // " could not be found in the registry")

  end function registry_get_tensor3

  !> Get pointer to a stored tensor4 by name.
  recursive function registry_get_tensor4(this, name) result(f)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    type(tensor4_t), pointer :: f
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'tensor4' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          f => this%entries(i)%get_tensor4()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       f => this%get_tensor4(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Tensor4 " // name // " could not be found in the registry")

  end function registry_get_tensor4

  !> Get pointer to a stored real scalar by name.
  recursive function registry_get_real_scalar(this, name) result(s)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    real(kind=rp), pointer :: s
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'real_scalar' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          s => this%entries(i)%get_real_scalar()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       s => this%get_real_scalar(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Real scalar " // name // " could not be found in the registry")

  end function registry_get_real_scalar

  !> Get pointer to a stored integer scalar by name.
  recursive function registry_get_integer_scalar(this, name) result(s)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias_target
    integer, pointer :: s
    logical :: found
    integer :: i

    found = .false.

    do i = 1, this%n_entries()
       if (this%entries(i)%get_type() .eq. 'integer_scalar' .and. &
            this%entries(i)%get_name() .eq. trim(name)) then
          s => this%entries(i)%get_integer_scalar()
          return
       end if
    end do

    call this%aliases%get(name, alias_target, found)
    if (found) then
       s => this%get_integer_scalar(alias_target)
       return
    end if

    call this%print_contents()
    call neko_error("Integer scalar " // name // &
         " could not be found in the registry")

  end function registry_get_integer_scalar

  ! ========================================================================== !
  ! Methods for checking existence of objects in the registry

  !> Check if an entry with a given name is already in the registry.
  !! @param name The name of the entry.
  !! @param type The type of the entry. Optional, if not provided, will check
  !!        all types.
  function registry_entry_exists(this, name, type) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: type
    logical :: found
    integer :: i

    found = .false.
    do i = 1, this%n_entries()
       if (trim(this%entries(i)%get_name()) .eq. trim(name)) then
          if (present(type)) then
             if (trim(this%entries(i)%get_type()) .eq. trim(type)) then
                found = .true.
                return
             end if
          else
             found = .true.
             return
          end if
       end if
    end do

    found = this%aliases%valid_path(name)
  end function registry_entry_exists

  !> Check if a field with a given name is already in the registry.
  function registry_field_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'field')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_field_exists

  !> Check if a vector with a given name is already in the registry.
  function registry_vector_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'vector')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_vector_exists

  !> Check if a matrix with a given name is already in the registry.
  function registry_matrix_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'matrix')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_matrix_exists

  !> Check if a tensor3 with a given name is already in the registry.
  function registry_tensor3_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'tensor3')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_tensor3_exists

  !> Check if a tensor4 with a given name is already in the registry.
  function registry_tensor4_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'tensor4')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_tensor4_exists

  !> Check if a real scalar with a given name is already in the registry.
  function registry_real_scalar_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'real_scalar')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_real_scalar_exists

  !> Check if an integer scalar with a given name is already in the registry.
  function registry_integer_scalar_exists(this, name) result(found)
    class(registry_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    logical :: found

    found = this%entry_exists(name, 'integer_scalar')
    if (.not. found) found = this%aliases%valid_path(name)

  end function registry_integer_scalar_exists

  ! ========================================================================== !
  ! Generic component accessor methods

  !> Get number of registered entries.
  pure function registry_n_entries(this, type) result(n)
    class(registry_t), intent(in) :: this
    character(len=*), intent(in), optional :: type
    integer :: n, i

    if (present(type)) then
       n = 0
       do i = 1, this%n_entries_
          if (this%entries(i)%get_type() .eq. trim(type)) then
             n = n + 1
          end if
       end do
    else
       n = this%n_entries_
    end if

  end function registry_n_entries

  !> Get the number of fields stored in the registry
  pure function registry_n_fields(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('field')
  end function registry_n_fields

  !> Get the number of vector stored in the registry
  pure function registry_n_vectors(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('vector')
  end function registry_n_vectors

  !> Get the number of matrix stored in the registry
  pure function registry_n_matrices(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('matrix')
  end function registry_n_matrices

  !> Get the number of tensor3 stored in the registry
  pure function registry_n_tensor3s(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('tensor3')
  end function registry_n_tensor3s

  !> Get the number of tensor4 stored in the registry
  pure function registry_n_tensor4s(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('tensor4')
  end function registry_n_tensor4s

  !> Get the number of real scalars stored in the registry
  pure function registry_n_real_scalars(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('real_scalar')
  end function registry_n_real_scalars

  !> Get the number of integer scalars stored in the registry
  pure function registry_n_integer_scalars(this) result(n)
    class(registry_t), intent(in) :: this
    integer :: n, i
    n = this%n_entries('integer_scalar')
  end function registry_n_integer_scalars

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
    character(len=14), parameter :: types(7) = [ &
         'field         ', &
         'vector        ', &
         'matrix        ', &
         'tensor3       ', &
         'tensor4       ', &
         'real_scalar   ', &
         'integer_scalar' ]
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
