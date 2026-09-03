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
!> Defines a tensor4
module tensor4
  use array, only : array_t
  use num_types, only : rp, xp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use utils, only : neko_warning, neko_error
  implicit none
  private

  type, public, extends(array_t) :: tensor4_t
     !> tensor4 shaped pointer.
     real(kind=rp), pointer, dimension(:,:,:,:) :: x => null()
     !> Number of tensor4 dimensions.
     integer, dimension(4), private :: dims = [0, 0, 0, 0]

   contains
     generic :: init => init_tensor4
     !> Initialise a tensor4 of size `n1*n2*n3*n4`.
     procedure, pass(this) :: init_tensor4 => tensor4_init
     !> Deallocate a tensor4.
     procedure, pass(this) :: free => tensor4_free

     !> Returns the dimensions of the tensor4.
     procedure, pass(this) :: get_dims => tensor4_dims
     !> Return the specific dimension of the tensor4.
     procedure, pass(this) :: get_dim => tensor4_dim

  end type tensor4_t

  type, public :: tensor4_ptr_t
     type(tensor4_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => tensor4_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => tensor4_ptr_free
  end type tensor4_ptr_t

contains

  !> Initialise a tensor of size `n1*n2*n3*n4`.
  !! @param this Tensor to initialise.
  !! @param n1 Size of the first dimension.
  !! @param n2 Size of the second dimension.
  !! @param n3 Size of the third dimension.
  !! @param n4 Size of the fourth dimension.
  !! @param name Optional name of the tensor.
  subroutine tensor4_init(this, n1, n2, n3, n4, name)
    class(tensor4_t), intent(inout), target :: this
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: n4
    character(len=*), intent(in), optional :: name

    call this%init_base(n1*n2*n3*n4, name)
    this%dims = [n1, n2, n3, n4]
    this%x(1:n1, 1:n2, 1:n3, 1:n4) => this%data(:)

  end subroutine tensor4_init

  !> Deallocate a tensor4.
  subroutine tensor4_free(this)
    class(tensor4_t), intent(inout) :: this

    call this%free_base()
    this%dims = [0, 0, 0, 0]
    nullify(this%x)

  end subroutine tensor4_free

  !> Returns the dimensions of the tensor4.
  pure function tensor4_dims(this) result(dims)
    class(tensor4_t), intent(in) :: this
    integer :: dims(4)
    dims = this%dims
  end function tensor4_dims

  !> Return the specific dimension of the tensor4.
  function tensor4_dim(this, dim) result(d)
    class(tensor4_t), intent(in) :: this
    integer, intent(in) :: dim
    integer :: d
    d = this%dims(dim)
  end function tensor4_dim

  ! ========================================================================== !
  ! tensor4 pointer type subroutines

  !> Constructor. Just assigns the pointer.
  !! @param this Pointer wrapper to initialise.
  !! @param ptr Tensor to point to.
  subroutine tensor4_ptr_init(this, ptr)
    class(tensor4_ptr_t), intent(inout) :: this
    type(tensor4_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine tensor4_ptr_init

  !> Destructor. Just nullifies the pointer.
  !! @param this Pointer wrapper to destroy.
  subroutine tensor4_ptr_free(this)
    class(tensor4_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine tensor4_ptr_free

end module tensor4
