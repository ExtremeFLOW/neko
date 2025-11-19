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
!> Defines a registry for storing and requesting temporary vectors
!! This can be used when you have a function that will be called
!! often and you don't want to create temporary vectors (work arrays) inside
!! it on each call.
module vector_scratch_registry
  use vector, only : vector_t, vector_ptr_t
  use vector_math, only : vector_rzero
  implicit none
  private

  type, public :: vector_scratch_registry_t
     !> list of scratch vectors
     type(vector_ptr_t), private, allocatable :: vectors(:)
     !> Tracks which vectors are used
     logical, private, allocatable :: inuse(:)
     !> number of registered vectors
     integer, private :: nvectors = 0
     !> number of vectors in use
     integer, private :: nvectors_inuse = 0
     !> the size the vectors array is increased by upon reallocation
     integer, private :: expansion_size = 10
   contains
     procedure, private, pass(this) :: expand
     !> Constructor
     procedure, pass(this) :: init => scratch_registry_init
     !> Destructor
     procedure, pass(this) :: free => scratch_registry_free
     !> Getter for nvectors
     procedure, pass(this) :: get_nvectors
     !> Getter for nvectors_inuse
     procedure, pass(this) :: get_nvectors_inuse
     !> Getter for expansion_size
     procedure, pass(this) :: get_expansion_size
     !> Return size of allocated vectors
     procedure, pass(this) :: get_size
     !> Get value of inuse for a given index
     procedure, pass(this) :: get_inuse
     !> Get a new scratch vector
     procedure, pass(this) :: request_vector
     procedure, pass(this) :: relinquish_vector_single
     procedure, pass(this) :: relinquish_vector_multiple
     !> Free a vector for later reuse
     generic :: relinquish_vector => relinquish_vector_single, &
          relinquish_vector_multiple
  end type vector_scratch_registry_t

  !> Global scratch registry
  type(vector_scratch_registry_t), public, target :: &
       neko_vector_scratch_registry

contains

  !> Constructor, optionally taking initial registry and expansion
  !! size as argument
  subroutine scratch_registry_init(this, size, expansion_size)
    class(vector_scratch_registry_t), intent(inout) :: this
    integer, optional, intent(in) :: size
    integer, optional, intent(in) :: expansion_size
    integer :: i, size_

    call this%free()

    ! Set sizes
    if (present(size)) then
       size_ = size
    else
       size_ = 10
    end if

    if (present(expansion_size)) then
       this%expansion_size = expansion_size
    else
       this%expansion_size = 10
    end if

    ! Initial allocation
    allocate (this%vectors(size_))
    do i = 1, size_
       allocate(this%vectors(i)%ptr)
    end do
    allocate (this%inuse(size_))
    this%inuse(:) = .false.

  end subroutine scratch_registry_init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(vector_scratch_registry_t), intent(inout):: this
    integer :: i

    if (allocated(this%vectors)) then
       do i = 1, this%nvectors
          call this%vectors(i)%ptr%free()
          deallocate(this%vectors(i)%ptr)
       end do

       deallocate(this%vectors)
       deallocate(this%inuse)
    end if

    this%nvectors = 0
    this%nvectors_inuse = 0
    this%expansion_size = 10
  end subroutine scratch_registry_free


  !> Get the number of vectors stored in the registry
  pure function get_nvectors(this) result(n)
    class(vector_scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%nvectors
  end function get_nvectors

  pure function get_nvectors_inuse(this) result(n)
    class(vector_scratch_registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1,this%get_size()
       if (this%inuse(i)) n = n + 1
    end do
  end function get_nvectors_inuse

  !> Get the size of the vectors array
  pure function get_size(this) result(n)
    class(vector_scratch_registry_t), intent(in) :: this
    integer :: n

    n = size(this%vectors)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(vector_scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Expand the registry by the expansion size
  subroutine expand(this)
    class(vector_scratch_registry_t), intent(inout) :: this
    type(vector_ptr_t), allocatable :: temp(:)
    logical, allocatable :: temp2(:)
    integer :: i

    allocate(temp(this%get_size() + this%expansion_size))
    temp(1:this%nvectors) = this%vectors(1:this%nvectors)

    do i = this%nvectors + 1, size(temp)
       allocate(temp(i)%ptr)
    enddo

    call move_alloc(temp, this%vectors)

    allocate(temp2(this%get_size() + this%expansion_size))
    temp2(1:this%nvectors) = this%inuse(1:this%nvectors)
    temp2(this%nvectors+1:) = .false.
    this%inuse = temp2
  end subroutine expand

  !> Get a vector from the registry by assigning it to a pointer.
  !! @param n Size of the requested vector.
  !! @param v Pointer to the requested vector.
  !! @param index Index of the vector in the registry (for relinquishing later).
  subroutine request_vector(this, n, v, index)
    class(vector_scratch_registry_t), target, intent(inout) :: this
    integer, intent(in) :: n
    type(vector_t), pointer, intent(inout) :: v
    integer, intent(inout) :: index

    associate(nvectors => this%nvectors, nvectors_inuse => this%nvectors_inuse)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. allocated(this%vectors(index)%ptr%x)) then
               call this%vectors(index)%ptr%init(n)
               nvectors = nvectors + 1
            else if (this%vectors(index)%ptr%size() .ne. n) then
               cycle
            end if

            v => this%vectors(index)%ptr
            call vector_rzero(v)
            this%inuse(index) = .true.
            this%nvectors_inuse = this%nvectors_inuse + 1
            return
         end if
      end do

      ! all existing vectors in use, we need to expand to add a new one
      index = nvectors + 1
      call this%expand()
      nvectors = nvectors + 1
      nvectors_inuse = nvectors_inuse + 1
      this%inuse(nvectors) = .true.
      call this%vectors(nvectors)%ptr%init(n)
      v => this%vectors(nvectors)%ptr

    end associate
  end subroutine request_vector

  !> Relinquish the use of a vector in the registry
  subroutine relinquish_vector_single(this, index)
    class(vector_scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the vector to free

    this%inuse(index) = .false.
    this%nvectors_inuse = this%nvectors_inuse - 1
  end subroutine relinquish_vector_single

  subroutine relinquish_vector_multiple(this, indices)
    class(vector_scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:) !< The indices of the vector to free
    integer :: i

    do i = 1, size(indices)
       this%inuse(indices(i)) = .false.
    end do
    this%nvectors_inuse = this%nvectors_inuse - size(indices)
  end subroutine relinquish_vector_multiple

  !> Get whether a vector is in use
  pure logical function get_inuse(this, index)
    class(vector_scratch_registry_t), intent(in) :: this
    integer, intent(in) :: index

    get_inuse = this%inuse(index)
  end function get_inuse

end module vector_scratch_registry
