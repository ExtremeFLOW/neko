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
!> Scaffolding for future global boundary resolution.
module bc_resolver
  use bc, only : bc_t
  use mask, only : mask_t
  use math, only : cfill_mask
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use utils, only : neko_error
  use device, only : device_get_ptr
  use device_math, only : device_cfill_mask
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, abstract :: bc_resolver_t
   contains
     procedure, pass(this) :: init => bc_resolver_init
     procedure, pass(this) :: free => bc_resolver_free
  end type bc_resolver_t

  type, public, extends(bc_resolver_t) :: scalar_bc_resolver_t
     type(mask_t) :: dof_mask
   contains
     procedure, pass(this) :: free => scalar_bc_resolver_free
     procedure, pass(this) :: mark_bc => scalar_bc_resolver_mark_bc
     procedure, pass(this) :: apply_scalar => scalar_bc_resolver_apply_scalar
  end type scalar_bc_resolver_t

  type, public, extends(bc_resolver_t) :: vector_bc_resolver_t
     type(scalar_bc_resolver_t) :: x
     type(scalar_bc_resolver_t) :: y
     type(scalar_bc_resolver_t) :: z
   contains
     procedure, pass(this) :: free => vector_bc_resolver_free
  end type vector_bc_resolver_t

  type, public, extends(bc_resolver_t) :: coupled_vector_bc_resolver_t
  end type coupled_vector_bc_resolver_t

contains

  !> Initialize the boundary resolver.
  subroutine bc_resolver_init(this)
    class(bc_resolver_t), intent(inout) :: this
  end subroutine bc_resolver_init

  !> Free the boundary resolver.
  subroutine bc_resolver_free(this)
    class(bc_resolver_t), intent(inout) :: this
  end subroutine bc_resolver_free

  !> Free a scalar boundary resolver.
  subroutine scalar_bc_resolver_free(this)
    class(scalar_bc_resolver_t), intent(inout) :: this

    call this%dof_mask%free()
  end subroutine scalar_bc_resolver_free

  !> Add the constrained dofs from a scalar boundary condition.
  subroutine scalar_bc_resolver_mark_bc(this, bc)
    class(scalar_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    integer, pointer :: current_mask(:)
    integer, allocatable :: merged_mask(:)
    integer :: current_size
    integer :: incoming_size
    integer :: merged_size
    integer :: i

    if (.not. allocated(bc%msk)) then
      call neko_error("Attempting to mark resolver from an unfinalized BC.")
    end if

    incoming_size = bc%msk(0)
    if (incoming_size .eq. 0) return

    if (.not. this%dof_mask%is_set()) then
       call this%dof_mask%init(bc%msk(1:incoming_size), incoming_size)
       return
    end if

    current_size = this%dof_mask%size()
    current_mask => this%dof_mask%get()
    allocate(merged_mask(current_size + incoming_size))

    merged_mask(1:current_size) = current_mask(1:current_size)
    merged_size = current_size

    do i = 1, incoming_size
       if (.not. any(merged_mask(1:merged_size) .eq. bc%msk(i))) then
          merged_size = merged_size + 1
          merged_mask(merged_size) = bc%msk(i)
       end if
    end do

    call this%dof_mask%set(merged_mask(1:merged_size), merged_size)
  end subroutine scalar_bc_resolver_mark_bc

  !> Apply the scalar boundary constraints by zeroing constrained dofs.
  subroutine scalar_bc_resolver_apply_scalar(this, x, n)
    class(scalar_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    type(c_ptr) :: x_d

    if (.not. this%dof_mask%is_set()) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       call device_cfill_mask(x_d, 0.0_rp, n, this%dof_mask%get_d(), &
            this%dof_mask%size())
    else
       call cfill_mask(x, 0.0_rp, n, this%dof_mask%get(), this%dof_mask%size())
    end if
  end subroutine scalar_bc_resolver_apply_scalar

  !> Free a vector boundary resolver.
  subroutine vector_bc_resolver_free(this)
    class(vector_bc_resolver_t), intent(inout) :: this

    call this%x%free()
    call this%y%free()
    call this%z%free()
  end subroutine vector_bc_resolver_free

end module bc_resolver
