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
!> Defines a registry for storing and requesting temporary fields
!! This can be used when you have a function that will be called
!! often and you don't want to create temporary fields (work arrays) inside
!! it on each call.
module scratch_registry
  use num_types
  use field
  use utils
  implicit none
  private
  
  type :: scratch_registry_t
     type(field_t), private, allocatable :: fields(:) !< list of scratch fields 
     logical, private, allocatable :: inuse(:)        !< Tracks which fields are used
     integer, private :: n                            !< number of registered fields
     integer, private :: expansion_size               !< the size the fields array is increased by upon reallocation
     type(dofmap_t), pointer :: dof                   !< Dofmap
   contains
     procedure, private, pass(this) :: expand
     procedure, pass(this) :: init => scratch_registry_init
     procedure, pass(this) :: free => scratch_registry_free
     procedure, pass(this) :: add_field
     procedure, pass(this) :: n_fields
     procedure, pass(this) :: get_expansion_size
     procedure, pass(this) :: get_size
     procedure, pass(this) :: get_field
  end type scratch_registry_t

  !> Global scratch registry
  type(scratch_registry_t), public, target :: neko_scratch_registry

contains
  !> Constructor, optionally taking initial registry and expansion
  !> size as argument
  subroutine scratch_registry_init(this, dof, size, expansion_size)
    class(scratch_registry_t), intent(inout):: this
    type(dofmap_t), target, intent(in) :: dof
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size
    
    this%dof => dof

    if (present(size)) then
       allocate (this%fields(size))
       allocate (this%inuse(size))
    else
       allocate (this%fields(10))
       allocate (this%inuse(10))
    end if

    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 10
    end if

    this%n = 0
  end subroutine scratch_registry_init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(scratch_registry_t), intent(inout):: this
    integer :: i

    do i=1, this%n_fields()
       call field_free(this%fields(i))
    end do
    deallocate(this%fields)
    deallocate(this%inuse)
  end subroutine scratch_registry_free

  !> expand the fields array so as to accomodate more fields
  subroutine expand(this)
    class(scratch_registry_t), intent(inout) :: this
    type(field_t), allocatable :: temp(:)  
    logical, allocatable :: temp2(:)  

    allocate(temp(this%n + this%expansion_size))
    allocate(temp2(this%n + this%expansion_size))
    temp(1:this%n) = this%fields(1:this%n)

    temp2(1:this%n) = this%inuse(1:this%n)
    temp2(this%n+1:) = .false.
    call move_alloc(temp, this%fields)
    this%inuse = temp2
  end subroutine expand

  subroutine add_field(this)
    class(scratch_registry_t), intent(inout) :: this

    if (this%n_fields() == size(this%fields)) then
       call this%expand()
    end if

    this%n = this%n + 1

    ! initialize the field at the appropraite index
    call field_init(this%fields(this%n), this%dof, 'wrk')

  end subroutine add_field

  !> Get the number of fields stored in the registry
  pure function n_fields(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%n  
  end function n_fields

  !> Get the size of the fields array
  pure function get_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = size(this%fields)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Get a field from the registry by assigning it to a pointer
  subroutine get_field(this, f, index, i)
    class(scratch_registry_t), target, intent(inout) :: this
    type(field_t), pointer, intent(inout) :: f
    integer, intent(inout) :: index !< The index of the field in the inuse array
    integer :: i
    
    do i=1,this%n
       if (this%inuse(i) .eqv. .false.) then
         f => this%fields(i)
         index = i
         this%inuse(i) = .true.
         return
       end if
    end do
    call add_field(this)
    index = this%n
    this%inuse(this%n) = .true.
    f => this%fields(this%n)
  end subroutine get_field
  
  !> Relinquish the use of a field in the registry
  subroutine relinquish_field(this, index)
    class(scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the field to free
    
    this%inuse(index) = .false.
  end subroutine relinquish_field

end module scratch_registry
