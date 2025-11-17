! Copyright (c) 2025, The Neko Authors
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
!> Defines the `vector_registry_t` type.

module vector_registry
  use, intrinsic :: iso_fortran_env, only: error_unit
  use vector, only : vector_t
  use utils, only: neko_error
  use comm, only : pe_rank
  use json_module, only : json_file
  use json_utils, only : json_get
  implicit none
  private

  !> A named `vector_t`.
  type :: named_vector_t
     character(len=:), allocatable :: name
     type(vector_t) :: vector
   contains
     !> Constructor
     procedure, pass(this) :: init => named_vector_init
     !> Destructor
     procedure, pass(this) :: free => named_vector_free
  end type named_vector_t


  !! A registry of named `vector_t` objects.
  !! @details Same as field_registry_t but for vectors. Stores named
  !! vectors, but raw vectors are returned by the `get_` methods.
  !! Supports adding aliases to vectors as well, so that the same vector can be
  !! accessed by different names.
  type, public :: vector_registry_t
     !> List of vectors stored.
     type(named_vector_t), private, allocatable :: vectors(:)
     !> List of aliases to vectors stored.
     type(json_file), private, allocatable :: aliases(:)
     !> Number of registered vectors.
     integer, private :: n_vectors_
     !> Number of aliases.
     integer, private :: n_aliases_
     !> The size the vectors array is increased by upon reallocation.
     integer, private :: expansion_size
   contains
     !> Expand the vectors array so as to accomodate more vectors.
     procedure, private, pass(this) :: expand => vector_registry_expand
     !> Expand the aliases array so as to accomodate more aliases.
     procedure, private, pass(this) :: expand_aliases => &
          vector_registry_expand_aliases
     !> Constructor.
     procedure, pass(this) :: init => vector_registry_init
     !> Destructor.
     procedure, pass(this) :: free => vector_registry_free
     !> Add a vector to the registry.
     procedure, pass(this) :: add_vector => vector_registry_add_vector
     !> Add an alias to a vector in the registry.
     procedure, pass(this) :: add_alias => vector_registry_add_alias
     !> Get the number of vectors in the registry.
     procedure, pass(this) :: n_vectors => vector_registry_n_vectors
     !> Get the number of aliases in the registry.
     procedure, pass(this) :: n_aliases => vector_registry_n_aliases
     !> Get pointer to a stored vector by index.
     procedure, pass(this) :: get_vector_by_index => &
          vector_registry_get_vector_by_index
     !> Get pointer to a stored vector by name.
     procedure, pass(this) :: get_vector_by_name => &
          vector_registry_get_vector_by_name
     !> Get the `expansion_size`
     procedure, pass(this) :: get_expansion_size => &
          vector_registry_get_expansion_size
     !> Get total allocated size of `vectors`.
     procedure, pass(this) :: get_size => vector_registry_get_size
     !> Check if a vector with a given name is already in the registry.
     procedure, pass(this) :: vector_exists => vector_registry_vector_exists
     generic :: get_vector => get_vector_by_index, get_vector_by_name
  end type vector_registry_t

  !> Global vector registry
  type(vector_registry_t), public, target :: neko_vector_registry

contains
  !> Constructor
  subroutine named_vector_init(this, name, n)
    class(named_vector_t), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    integer, intent(in) :: n

    this%name = name
    call this%vector%init(n)
  end subroutine named_vector_init

  !> Destructor
  subroutine named_vector_free(this)
    class(named_vector_t), intent(inout) :: this

    call this%vector%free()
    if (allocated(this%name)) then
       deallocate(this%name)
    end if
  end subroutine named_vector_free

  !> Constructor
  !! @param size The allocation size of `vectors` on init.
  !! @param expansion_size The number of entries added to `vectors` on
  !! expansion.
  subroutine vector_registry_init(this, size, expansion_size)
    class(vector_registry_t), intent(inout):: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size

    call this%free()

    if (present(size)) then
       allocate (this%vectors(size))
       allocate (this%aliases(size))
    else
       allocate (this%vectors(50))
       allocate (this%aliases(50))
    end if

    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 50
    end if

    this%n_vectors_ = 0
    this%n_aliases_ = 0
  end subroutine vector_registry_init

  !> Destructor
  subroutine vector_registry_free(this)
    class(vector_registry_t), intent(inout):: this
    integer :: i
    if (allocated(this%vectors)) then
       do i = 1, this%n_vectors()
          call this%vectors(i)%free()
       end do
       deallocate(this%vectors)
    end if

    if (allocated(this%aliases)) then
       deallocate(this%aliases)
    end if

    this%n_vectors_ = 0
    this%n_aliases_ = 0
    this%expansion_size = 0
  end subroutine vector_registry_free

  !> Expand the vectors array so as to accomodate more vectors.
  subroutine vector_registry_expand(this)
    class(vector_registry_t), intent(inout) :: this
    type(named_vector_t), allocatable :: temp(:)

    allocate(temp(this%n_vectors_ + this%expansion_size))
    temp(1:this%n_vectors_) = this%vectors(1:this%n_vectors_)
    call move_alloc(temp, this%vectors)
  end subroutine vector_registry_expand

  !> Expand the aliases array so as to accomodate more aliases.
  subroutine vector_registry_expand_aliases(this)
    class(vector_registry_t), intent(inout) :: this
    type(json_file), allocatable :: temp(:)

    allocate(temp(this%n_aliases() + this%expansion_size))
    temp(1:this%n_aliases()) = this%aliases(1:this%n_vectors_)
    call move_alloc(temp, this%aliases)
  end subroutine vector_registry_expand_aliases

  !> Add a vector to the registry.
  !! @param n The size of the vector.
  !! @param name The name of the vector.
  !! @param ignore_existing If true, will do nothing if the vector is already in
  !! the registry. If false, will throw an error. Optional, defaults to false.
  subroutine vector_registry_add_vector(this, n, name, ignore_existing)
    class(vector_registry_t), intent(inout) :: this
    integer, intent(in) :: n
    character(len=*), target, intent(in) :: name
    logical, optional, intent(in) :: ignore_existing
    logical :: ignore_existing_

    ignore_existing_ = .false.
    if (present(ignore_existing_)) then
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

    if (this%n_vectors() == size(this%vectors)) then
       call this%expand()
    end if

    this%n_vectors_ = this%n_vectors_ + 1

    ! Initialize the named vector at the appropriate index
    call this%vectors(this%n_vectors_)%init(name, n)

  end subroutine vector_registry_add_vector

  !> Add an alias for an existing vector in the registry.
  !! @param alias The alias.
  !! @param name The name of the vector.
  subroutine vector_registry_add_alias(this, alias, name)
    class(vector_registry_t), intent(inout) :: this
    character(len=*), target, intent(in) :: alias
    character(len=*), target, intent(in) :: name

    if (this%vector_exists(alias)) then
       call neko_error("Cannot create alias. Vector " // alias // &
            " already exists in the registry")
    end if

    if (this%vector_exists(name)) then
       if (this%n_aliases_ == size(this%aliases)) then
          call this%expand_aliases()
       end if

       this%n_aliases_ = this%n_aliases_ + 1

       call this%aliases(this%n_aliases_)%initialize()
       call this%aliases(this%n_aliases_)%add("alias", alias)
       call this%aliases(this%n_aliases_)%add("target", name)
    else
       call neko_error("Cannot create alias. Vector " // name // &
            " could not be found in the registry")
    end if
  end subroutine vector_registry_add_alias

  !> Get the number of vectors stored in the registry
  pure function vector_registry_n_vectors(this) result(n)
    class(vector_registry_t), intent(in) :: this
    integer :: n

    n = this%n_vectors_
  end function vector_registry_n_vectors

  !> Get the number of aliases stored in the registry
  pure function vector_registry_n_aliases(this) result(n)
    class(vector_registry_t), intent(in) :: this
    integer :: n

    n = this%n_aliases_
  end function vector_registry_n_aliases

  !> Get the size of the vectors array.
  pure function vector_registry_get_size(this) result(n)
    class(vector_registry_t), intent(in) :: this
    integer :: n

    n = size(this%vectors)
  end function vector_registry_get_size

  !> Get the expansion size.
  pure function vector_registry_get_expansion_size(this) result(n)
    class(vector_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function vector_registry_get_expansion_size

  !> Get pointer to a stored vector by index.
  function vector_registry_get_vector_by_index(this, i) result(f)
    class(vector_registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(vector_t), pointer :: f

    if (i < 1) then
       call neko_error("Vector index must be > 1")
    else if (i > this%n_vectors()) then
       call neko_error("Vector index exceeds number of stored fields")
    endif

    f => this%vectors(i)%vector
  end function vector_registry_get_vector_by_index

  !> Get pointer to a stored vector by name.
  recursive function vector_registry_get_vector_by_name(this, name) result(f)
    class(vector_registry_t), target, intent(in) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias
    character(len=:), allocatable :: alias_target
    type(vector_t), pointer :: f
    logical :: found
    integer :: i
    type(json_file), pointer :: alias_json ! need this for some reason

    found = .false.

    do i = 1, this%n_vectors()
       if (this%vectors(i)%name == trim(name)) then
          f => this%vectors(i)%vector
          found = .true.
          exit
       end if
    end do

    do i = 1, this%n_aliases()
       alias_json => this%aliases(i)
       call json_get(alias_json, "alias", alias)
       if (alias == trim(name)) then
          call json_get(alias_json, "target", alias_target)
          f => this%get_vector_by_name(alias_target)
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       if (pe_rank .eq. 0) then
          write(error_unit,*) "Current vector_registry contents:"

          do i=1, this%n_vectors()
             write(error_unit,*) "- ", this%vectors(i)%name
          end do
       end if
       call neko_error("Vector " // name // &
            " could not be found in the registry")
    end if
  end function vector_registry_get_vector_by_name

  !> Check if a vector with a given name is already in the registry.
  function vector_registry_vector_exists(this, name) result(found)
    class(vector_registry_t), target, intent(in) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: alias
    logical :: found
    integer :: i
    type(json_file), pointer :: alias_json

    found = .false.
    do i=1, this%n_vectors()
       if (this%vectors(i)%name == name) then
          found = .true.
          exit
       end if
    end do

    do i=1, this%n_aliases()
       alias_json => this%aliases(i)
       call json_get(alias_json, "alias", alias)
       if (alias == name) then
          found = .true.
          exit
       end if
    end do
  end function vector_registry_vector_exists

end module vector_registry
