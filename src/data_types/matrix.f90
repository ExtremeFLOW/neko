! Copyright (c) 2024-2026, The Neko Authors
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
!> Defines a matrix
module matrix
  use array, only : array_t
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use utils, only : neko_error, neko_warning
  use cpu_matrix_math, only : cpu_matrix_inverse
  implicit none
  private

  type, public, extends(array_t) :: matrix_t
     !> Matrix shaped pointer.
     real(kind=rp), pointer, dimension(:,:) :: x => null()
     !> Number of matrix rows.
     integer, dimension(2), private :: dims = [0, 0]

   contains
     !> Generic init method. Calls the matrix specific init method.
     generic :: init => init_matrix
     !> Initialise a matrix of size `nrows*ncols`.
     procedure, pass(this), private :: init_matrix => matrix_init
     !> Deallocate a matrix.
     procedure, pass(this) :: free => matrix_free

     !> Returns the dimensions of the matrix.
     procedure, pass(this) :: get_dims => matrix_dims
     !> Return the specific dimension of the matrix.
     procedure, pass(this) :: get_dim => matrix_dim
     !> Returns the number of rows in the matrix.
     procedure, pass(this) :: get_nrows => matrix_nrows
     !> Returns the number of columns in the matrix.
     procedure, pass(this) :: get_ncols => matrix_ncols

     !> Invert a matrix.
     procedure, pass(this) :: inverse => matrix_inverse
     !> Invert a matrix on the host.
     procedure, pass(this) :: inverse_on_host => matrix_inverse_on_host

  end type matrix_t

  type, public :: matrix_ptr_t
     type(matrix_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => matrix_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => matrix_ptr_free
  end type matrix_ptr_t

contains

  !> Initialise a matrix of size `nrows*ncols`.
  !! @param nrows Number of rows.
  !! @param ncols Number of columns.
  !! @param name Optional name of the object.
  subroutine matrix_init(this, nrows, ncols, name)
    class(matrix_t), intent(inout), target :: this
    integer, intent(in) :: nrows, ncols
    character(len=*), intent(in), optional :: name

    call this%init_base(nrows*ncols, name)

    this%dims = [nrows, ncols]
    this%x(1:nrows, 1:ncols) => this%data(:)

  end subroutine matrix_init

  !> Deallocate a matrix.
  subroutine matrix_free(this)
    class(matrix_t), intent(inout) :: this

    call this%free_base()
    this%dims = [0, 0]
    nullify(this%x)

  end subroutine matrix_free

  !> Returns the dimensions of the matrix.
  pure function matrix_dims(this) result(dims)
    class(matrix_t), intent(in) :: this
    integer :: dims(2)
    dims = this%dims
  end function matrix_dims

  !> Return the specific dimension of the matrix.
  function matrix_dim(this, dim) result(d)
    class(matrix_t), intent(in) :: this
    integer, intent(in) :: dim
    integer :: d
    d = this%dims(dim)
  end function matrix_dim

  !> Returns the number of rows in the matrix.
  pure function matrix_nrows(this) result(nr)
    class(matrix_t), intent(in) :: this
    integer :: nr
    nr = this%dims(1)
  end function matrix_nrows

  !> Returns the number of columns in the matrix.
  pure function matrix_ncols(this) result(nc)
    class(matrix_t), intent(in) :: this
    integer :: nc
    nc = this%dims(2)
  end function matrix_ncols

  ! ========================================================================== !
  ! matrix pointer type subroutines

  subroutine matrix_ptr_init(this, ptr)
    class(matrix_ptr_t), intent(inout) :: this
    type(matrix_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine matrix_ptr_init

  subroutine matrix_ptr_free(this)
    class(matrix_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine matrix_ptr_free

  ! ========================================================================== !
  ! Matrix Specific extensions

  !> Compute the inverse of the matrix.
  !! @details Invert the current matrix. If the optional argument `bcknd` is
  !! provided, then we expect the user have supplied a matrix which is ready to
  !! execute on the specified backend. Otherwise, we will call the appropriate
  !! backend based on the current configuration.
  subroutine matrix_inverse(this, bcknd)
    class(matrix_t), intent(inout) :: this
    integer, optional :: bcknd

    if (present(bcknd)) then
       if (bcknd .eq. 1) then
          call neko_error("matrix_inverse: GPU backend not implemented yet.")
       else if (bcknd .eq. 0) then
          call cpu_matrix_inverse(this%x, this%dims(1), this%dims(2))
       else
          call neko_error("matrix_inverse: Invalid backend specified. " // &
               "Use 0 for CPU or 1 for GPU.")
       end if
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_warning("matrix_inverse: GPU backend not implemented yet. " // &
            "Falling back to CPU.")
       call this%copy_from(DEVICE_TO_HOST, .true.)
       call cpu_matrix_inverse(this%x, this%dims(1), this%dims(2))
       call this%copy_from(HOST_TO_DEVICE, .true.)
    else
       call cpu_matrix_inverse(this%x, this%dims(1), this%dims(2))
    end if

  end subroutine matrix_inverse

  !> Compute the inverse of the matrix on the host.
  subroutine matrix_inverse_on_host(this)
    class(matrix_t), intent(inout) :: this
    call cpu_matrix_inverse(this%x, this%dims(1), this%dims(2))
  end subroutine matrix_inverse_on_host

end module matrix
