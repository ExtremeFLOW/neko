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
!> Defines a registry for storing and requesting temporary matrices
!! This can be used when you have a function that will be called
!! often and you don't want to create temporary matrices (work arrays) inside
!! it on each call.
module matrix_scratch_registry
  use matrix, only : matrix_t, matrix_ptr_t
  use matrix_math, only : matrix_rzero
  implicit none
  private

  type, public :: matrix_scratch_registry_t
     !> list of scratch matrices
     type(matrix_ptr_t), private, allocatable :: matrices(:)
     !> Tracks which matrices are used
     logical, private, allocatable :: inuse(:)
     !> number of registered matrices
     integer, private :: nmatrices = 0
     !> number of matrices in use
     integer, private :: nmatrices_inuse = 0
     !> the size the matrices array is increased by upon reallocation
     integer, private :: expansion_size = 10
   contains
     procedure, private, pass(this) :: expand
     !> Constructor
     procedure, pass(this) :: init => scratch_registry_init
     !> Destructor
     procedure, pass(this) :: free => scratch_registry_free
     !> Getter for nmatrices
     procedure, pass(this) :: get_nmatrices
     !> Getter for nmatrices_inuse
     procedure, pass(this) :: get_nmatrices_inuse
     !> Getter for expansion_size
     procedure, pass(this) :: get_expansion_size
     !> Return size of allocated matrices
     procedure, pass(this) :: get_size
     !> Get value of inuse for a given index
     procedure, pass(this) :: get_inuse
     !> Get a new scratch matrix
     procedure, pass(this) :: request_matrix
     procedure, pass(this) :: relinquish_matrix_single
     procedure, pass(this) :: relinquish_matrix_multiple
     !> Free a matrix for later reuse
     generic :: relinquish_matrix => relinquish_matrix_single, &
          relinquish_matrix_multiple
  end type matrix_scratch_registry_t

contains

  !> Constructor, optionally taking initial registry and expansion
  !! size as argument
  subroutine scratch_registry_init(this, size, expansion_size)
    class(matrix_scratch_registry_t), intent(inout) :: this
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
    allocate(this%matrices(size_))
    do i = 1, size_
       allocate(this%matrices(i)%ptr)
    end do
    allocate(this%inuse(size_))
    this%inuse(:) = .false.

  end subroutine scratch_registry_init

  !> Destructor
  subroutine scratch_registry_free(this)
    class(matrix_scratch_registry_t), intent(inout):: this
    integer :: i

    if (allocated(this%matrices)) then
       do i = 1, this%nmatrices
          call this%matrices(i)%ptr%free()
          deallocate(this%matrices(i)%ptr)
       end do

       deallocate(this%matrices)
       deallocate(this%inuse)
    end if

    this%nmatrices = 0
    this%nmatrices_inuse = 0
    this%expansion_size = 10
  end subroutine scratch_registry_free


  !> Get the number of matrices stored in the registry
  pure function get_nmatrices(this) result(n)
    class(matrix_scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%nmatrices
  end function get_nmatrices

  pure function get_nmatrices_inuse(this) result(n)
    class(matrix_scratch_registry_t), intent(in) :: this
    integer :: n, i

    n = 0
    do i = 1,this%get_size()
       if (this%inuse(i)) n = n + 1
    end do
  end function get_nmatrices_inuse

  !> Get the size of the matrices array
  pure function get_size(this) result(n)
    class(matrix_scratch_registry_t), intent(in) :: this
    integer :: n

    n = size(this%matrices)
  end function get_size

  !> Get the expansion size
  pure function get_expansion_size(this) result(n)
    class(matrix_scratch_registry_t), intent(in) :: this
    integer :: n

    n = this%expansion_size
  end function get_expansion_size

  !> Expand the registry by the expansion size
  subroutine expand(this)
    class(matrix_scratch_registry_t), intent(inout) :: this
    type(matrix_ptr_t), allocatable :: temp(:)
    logical, allocatable :: temp2(:)
    integer :: i

    allocate(temp(this%get_size() + this%expansion_size))
    temp(1:this%nmatrices) = this%matrices(1:this%nmatrices)

    do i = this%nmatrices + 1, size(temp)
       allocate(temp(i)%ptr)
    enddo

    call move_alloc(temp, this%matrices)

    allocate(temp2(this%get_size() + this%expansion_size))
    temp2(1:this%nmatrices) = this%inuse(1:this%nmatrices)
    temp2(this%nmatrices+1:) = .false.
    this%inuse = temp2
  end subroutine expand

  !> Get a matrix from the registry by assigning it to a pointer.
  !! @param nrows Number of rows of the requested matrix.
  !! @param ncols Number of columns of the requested matrix.
  !! @param m Pointer to the requested matrix.
  !! @param index Index of the matrix in the registry (for relinquishing later).
  subroutine request_matrix(this, nrows, ncols, m, index)
    class(matrix_scratch_registry_t), target, intent(inout) :: this
    integer, intent(in) :: nrows, ncols
    type(matrix_t), pointer, intent(inout) :: m
    integer, intent(inout) :: index

    associate(nmatrices => this%nmatrices, nmatrices_inuse => this%nmatrices_inuse)

      do index = 1, this%get_size()
         if (.not. this%inuse(index)) then

            if (.not. allocated(this%matrices(index)%ptr%x)) then
               call this%matrices(index)%ptr%init(nrows, ncols)
               nmatrices = nmatrices + 1
            else if (this%matrices(index)%ptr%get_nrows() .ne. nrows .and. &
                 this%matrices(index)%ptr%get_ncols() .ne. ncols) then
               cycle
            end if

            m => this%matrices(index)%ptr
            call matrix_rzero(m)
            this%inuse(index) = .true.
            this%nmatrices_inuse = this%nmatrices_inuse + 1
            return
         end if
      end do

      ! all existing matrices in use, we need to expand to add a new one
      index = nmatrices + 1
      call this%expand()
      nmatrices = nmatrices + 1
      nmatrices_inuse = nmatrices_inuse + 1
      this%inuse(nmatrices) = .true.
      call this%matrices(nmatrices)%ptr%init(nrows, ncols)
      m => this%matrices(nmatrices)%ptr

    end associate
  end subroutine request_matrix

  !> Relinquish the use of a matrix in the registry
  subroutine relinquish_matrix_single(this, index)
    class(matrix_scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: index !< The index of the matrix to free

    this%inuse(index) = .false.
    this%nmatrices_inuse = this%nmatrices_inuse - 1
  end subroutine relinquish_matrix_single

  subroutine relinquish_matrix_multiple(this, indices)
    class(matrix_scratch_registry_t), target, intent(inout) :: this
    integer, intent(inout) :: indices(:) !< The indices of the matrix to free
    integer :: i

    do i = 1, size(indices)
       this%inuse(indices(i)) = .false.
    end do
    this%nmatrices_inuse = this%nmatrices_inuse - size(indices)
  end subroutine relinquish_matrix_multiple

  !> Get whether a matrix is in use
  pure logical function get_inuse(this, index)
    class(matrix_scratch_registry_t), intent(in) :: this
    integer, intent(in) :: index

    get_inuse = this%inuse(index)
  end function get_inuse

end module matrix_scratch_registry
