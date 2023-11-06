! Copyright (c) 2018-2023, The Neko Authors
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
!
module field_registry
  use field, only : field_t
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  use htable, only : h_cptr_t
  use utils, only: neko_error
  implicit none
  private
  
  type :: field_registry_t
     !> list of fields stored
     type(field_t), private, allocatable :: fields(:) 
     !> number of registered fields
     integer, private :: n                            
     !> the size the fields array is increased by upon reallocation
     integer, private :: expansion_size               
   contains
     procedure, private, pass(this) :: expand
     procedure, pass(this) :: init => field_registry_init
     procedure, pass(this) :: free => field_registry_free
     procedure, pass(this) :: add_field
     procedure, pass(this) :: n_fields
     procedure, pass(this) :: get_field_by_index
     procedure, pass(this) :: get_field_by_name
     procedure, pass(this) :: get_expansion_size
     procedure, pass(this) :: get_size
     procedure, pass(this) :: field_exists
     generic :: get_field => get_field_by_index, get_field_by_name
  end type field_registry_t

  !> Global field registry
  type(field_registry_t), public, target :: neko_field_registry

contains
  !> Constructor, optionally taking initial registry and expansion
  !> size as argument
  subroutine field_registry_init(this, size, expansion_size)
    class(field_registry_t), intent(inout):: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size

    call this%free()

    if (present(size)) then
       allocate (this%fields(size))
    else
       allocate (this%fields(50))
    end if

    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 50
    end if

    this%n = 0
  end subroutine field_registry_init

  !> Destructor
  subroutine field_registry_free(this)
    class(field_registry_t), intent(inout):: this
    integer :: i
    if (allocated(this%fields)) then
       do i=1, this%n_fields()
          call this%fields(i)%free()
       end do
       deallocate(this%fields)
    end if
    this%n = 0
    this%expansion_size = 0
  end subroutine field_registry_free

  !> expand the fields array so as to accomodate more fields
  subroutine expand(this)
    class(field_registry_t), intent(inout) :: this
    type(field_t), allocatable :: temp(:)  
    integer :: i

    allocate(temp(this%n + this%expansion_size))
    temp(1:this%n) = this%fields(1:this%n)
    call move_alloc(temp, this%fields)


  end subroutine expand

  subroutine add_field(this, dof, fld_name)
    class(field_registry_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=*), target, intent(in) :: fld_name 
!    type(h_cptr_t) :: key
    integer :: i

    if (this%field_exists(fld_name)) then
       call neko_error("Field with name " // fld_name // &
            " is already registered")
    end if

    if (this%n_fields() == size(this%fields)) then
       call this%expand()
    end if

    this%n = this%n + 1

    ! initialize the field at the appropraite index
    call this%fields(this%n)%init( dof, fld_name)

    ! generate a key for the name lookup map and assign it to the index
    !    key%ptr = c_loc(fld_name)
    !    call this%name_index_map%set(key, this%n)

    !    write(*,*) "HTABLE DATA, ", this%name_index_map%get(key, i)
  end subroutine add_field

  !> Get the number of fields stored in the registry
  pure function n_fields(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%n  
  end function n_fields

  !> Get the size of the fields array
  pure function get_size(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = size(this%fields)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(field_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Get pointer to a stored field by index
  function get_field_by_index(this, i) result(f)
    class(field_registry_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(field_t), pointer :: f

    if (i < 1) then
       call neko_error("Field index must be > 1")
    else if (i > this%n_fields()) then
       call neko_error("Field index exceeds number of stored fields")
    endif

    f => this%fields(i)
  end function get_field_by_index

  !> Get pointer to a stored field by field name
  function get_field_by_name(this, name) result(f)
    class(field_registry_t), target, intent(in) :: this
    character(len=*), intent(in) ::name 
    type(field_t), pointer :: f
    logical :: found
    integer :: i

    found = .false.
    do i=1, this%n_fields()
       if (this%fields(i)%name == name) then
          f => this%fields(i)
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       call neko_error("Field " // name // &
            " could not be found in the registry")
    end if
  end function get_field_by_name

  !> Check if a field with a given name is already in the registry
  function field_exists(this, name) result(found)
    class(field_registry_t), target, intent(in) :: this
    character(len=*), intent(in) ::name 
    logical :: found
    integer :: i

    found = .false.
    do i=1, this%n_fields()
       if (this%fields(i)%name == name) then
          found = .true.
          exit
       end if
    end do
  end function field_exists



end module field_registry
