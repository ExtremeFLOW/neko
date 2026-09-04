! Copyright (c) 2026, The Neko Authors
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
!> Defines a tensor3
module tensor3
  use array, only : array_t
  use num_types, only : rp
  implicit none
  private

  type, public, extends(array_t) :: tensor3_t
     !> tensor3 shaped pointer.
     real(kind=rp), pointer, dimension(:,:,:) :: x => null()
     !> Number of tensor3 dimensions.
     integer, dimension(3), private :: dims = [0, 0, 0]

   contains
     !> Generic init method. Calls the tensor3 specific init method.
     generic :: init => init_tensor3
     !> Initialise a tensor3 of size `n1*n2*n3`.
     procedure, pass(this), private :: init_tensor3 => tensor3_init
     !> Deallocate a tensor3.
     procedure, pass(this) :: free => tensor3_free

     !> Returns the dimensions of the tensor3.
     procedure, pass(this) :: get_dims => tensor3_dims
     !> Return the specific dimension of the tensor3.
     procedure, pass(this) :: get_dim => tensor3_dim

  end type tensor3_t

  type, public :: tensor3_ptr_t
     type(tensor3_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => tensor3_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => tensor3_ptr_free
  end type tensor3_ptr_t

contains

  !> Initialise a tensor of size `n1*n2*n3`.
  !! @param t Tensor to initialise.
  !! @param n1 Size of the first dimension.
  !! @param n2 Size of the second dimension.
  !! @param n3 Size of the third dimension.
  !! @param name Optional name of the tensor.
  subroutine tensor3_init(this, n1, n2, n3, name)
    class(tensor3_t), intent(inout), target :: this
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    character(len=*), intent(in), optional :: name

    call this%init_base(n1*n2*n3, name)

    this%dims = [n1, n2, n3]
    this%x(1:n1, 1:n2, 1:n3) => this%data(:)

  end subroutine tensor3_init

  !> Deallocate a tensor3.
  subroutine tensor3_free(this)
    class(tensor3_t), intent(inout) :: this

    call this%free_base()
    this%dims = [0, 0, 0]
    nullify(this%x)

  end subroutine tensor3_free

  !> Returns the dimensions of the tensor3.
  pure function tensor3_dims(this) result(dims)
    class(tensor3_t), intent(in) :: this
    integer :: dims(3)
    dims = this%dims
  end function tensor3_dims

  !> Return the specific dimension of the tensor3.
  function tensor3_dim(this, dim) result(d)
    class(tensor3_t), intent(in) :: this
    integer, intent(in) :: dim
    integer :: d
    d = this%dims(dim)
  end function tensor3_dim

  ! ========================================================================== !
  ! tensor3 pointer type subroutines

  !> Constructor. Just assigns the pointer.
  !! @param this Pointer wrapper to initialise.
  !! @param ptr Tensor to point to.
  subroutine tensor3_ptr_init(this, ptr)
    class(tensor3_ptr_t), intent(inout) :: this
    type(tensor3_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine tensor3_ptr_init

  !> Destructor. Just nullifies the pointer.
  !! @param this Pointer wrapper to destroy.
  subroutine tensor3_ptr_free(this)
    class(tensor3_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine tensor3_ptr_free

end module tensor3
