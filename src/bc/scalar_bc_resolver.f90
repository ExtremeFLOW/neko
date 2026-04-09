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
!> Implements `scalar_resolver_t`.
module scalar_bc_resolver
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use mask, only : mask_t
  use field, only : field_t
  use fld_file, only : fld_file_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : cfill_mask, sort
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_get_ptr, DEVICE_TO_HOST
  use device_math, only : device_cfill_mask
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Resolver for scalar boundary conditions.
  !! @details Collects the constrained scalar degrees of freedom from one or
  !! more boundary conditions and applies them as a zeroing mask. For the scalar
  !! case this is it, i.e. there no actual resolution being done. Neumann
  !! conditions should simply never be marked at all. This type is therefore
  !! used purely for constraining the dirichlet boundary values to 0 in KSP
  !! solvers, the residual, etc.
  type, public :: scalar_bc_resolver_t
     !> Union of constrained scalar degrees of freedom.
     type(mask_t) :: dof_mask
   contains
     !> Free the scalar boundary-condition resolver.
     procedure, pass(this) :: free => scalar_bc_resolver_free
     !> Add the constrained dofs from a scalar boundary condition.
     procedure, pass(this) :: mark_bc => scalar_bc_resolver_mark_bc
     !> Add the constrained dofs from all boundary conditions in a list.
     procedure, pass(this) :: mark_bc_list => scalar_bc_resolver_mark_bc_list
     generic :: mark => mark_bc, mark_bc_list
     !> Apply the scalar boundary constraints by zeroing constrained dofs.
     procedure, pass(this) :: apply => scalar_bc_resolver_apply
     !> Write a field showing the resolved scalar BC mask.
     procedure, pass(this) :: debug_output => scalar_bc_resolver_debug_output
  end type scalar_bc_resolver_t

contains

  !> Free the scalar boundary-condition resolver.
  subroutine scalar_bc_resolver_free(this)
    class(scalar_bc_resolver_t), intent(inout) :: this

    call this%dof_mask%free()
  end subroutine scalar_bc_resolver_free

  !> Add the constrained dofs from a scalar boundary condition.
  !! @details Appends the bc mask to `this%dof_masks`, sorts, and takes the
  !! unique set.
  !! @param[in] bc Finalized scalar boundary condition contributing dofs to the
  !! resolver mask.
  subroutine scalar_bc_resolver_mark_bc(this, bc)
    class(scalar_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(in) :: bc

    integer, pointer :: current_mask(:)
    integer, allocatable :: merged_mask(:)
    integer, allocatable :: perm(:)
    integer :: current_size
    integer :: incoming_size
    integer :: merged_size
    integer :: i

    if (.not. allocated(bc%msk)) then
       call neko_error("Attempting to mark resolver from an unfinalized BC.")
    end if

    incoming_size = bc%msk(0)
    if (incoming_size .eq. 0) return

    ! If this is the first bc to be marked
    if (.not. this%dof_mask%is_set()) then
       allocate(merged_mask(incoming_size), perm(incoming_size))
       merged_mask = bc%msk(1:incoming_size)
       call sort(merged_mask, perm, incoming_size)

       merged_size = 1
       do i = 2, incoming_size
          if (merged_mask(i) .ne. merged_mask(merged_size)) then
             merged_size = merged_size + 1
             merged_mask(merged_size) = merged_mask(i)
          end if
       end do

       call this%dof_mask%init(merged_mask(1:merged_size), merged_size)
       return
    end if

    ! Append the incoming mask to the current mask
    current_size = this%dof_mask%size()
    current_mask => this%dof_mask%get()
    allocate(merged_mask(current_size + incoming_size))
    allocate(perm(current_size + incoming_size))

    merged_mask(1:current_size) = current_mask(1:current_size)
    merged_mask(current_size + 1:current_size + incoming_size) = &
         bc%msk(1:incoming_size)

    ! Sort the merged mask
    call sort(merged_mask, perm, current_size + incoming_size)

    ! Take out unique elements
    merged_size = 1
    do i = 2, current_size + incoming_size
       if (merged_mask(i) .ne. merged_mask(merged_size)) then
          merged_size = merged_size + 1
          merged_mask(merged_size) = merged_mask(i)
       end if
    end do

    call this%dof_mask%set(merged_mask(1:merged_size), merged_size)
  end subroutine scalar_bc_resolver_mark_bc

  !> Add the constrained dofs from all boundary conditions in a list.
  !! @param[in] bclst List of scalar boundary conditions to merge into the
  !! resolver mask.
  subroutine scalar_bc_resolver_mark_bc_list(this, bclst)
    class(scalar_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    integer :: i

    do i = 1, bclst%size()
       call this%mark_bc(bclst%get(i))
    end do
  end subroutine scalar_bc_resolver_mark_bc_list

  !> Apply the scalar boundary constraints by zeroing constrained dofs.
  !! @param[inout] x Scalar field values.
  !! @param[in] n Number of entries in `x`.
  !! @param[inout] strm Optional backend stream/queue used on device backends.
  subroutine scalar_bc_resolver_apply(this, x, n, strm)
    class(scalar_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    type(c_ptr), intent(inout), optional :: strm
    type(c_ptr) :: x_d

    if (.not. this%dof_mask%is_set()) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       call device_cfill_mask(x_d, 0.0_rp, n, this%dof_mask%get_d(), &
            this%dof_mask%size(), strm = strm)
    else
       call cfill_mask(x, 0.0_rp, n, this%dof_mask%get(), this%dof_mask%size())
    end if
  end subroutine scalar_bc_resolver_apply

  !> Write a field showing the resolved scalar BC mask.
  !! @param[in] field_name Optional base name for the output file. The .fld is
  !! appended automatically.
  subroutine scalar_bc_resolver_debug_output(this, field_name)
    class(scalar_bc_resolver_t), intent(inout) :: this
    character(len=*), intent(in), optional :: field_name
    type(field_t), pointer :: mask_field
    type(fld_file_t) :: mask_file
    integer :: scratch_idx
    character(len=:), allocatable :: field_name_

    if (present(field_name)) then
       field_name_ = trim(field_name)
    else
       field_name_ = 'scalar_bc_resolver_mask'
    end if

    call neko_scratch_registry%request_field(mask_field, scratch_idx, .true.)

    if (this%dof_mask%is_set()) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_cfill_mask(mask_field%x_d, 1.0_rp, mask_field%size(), &
               this%dof_mask%get_d(), this%dof_mask%size())
          call mask_field%copy_from(DEVICE_TO_HOST, .true.)
       else
          call cfill_mask(mask_field%x, 1.0_rp, mask_field%size(), &
               this%dof_mask%get(), this%dof_mask%size())
       end if
    end if

    call mask_file%init(field_name_ // '.fld')
    call mask_file%write(mask_field)
    call neko_scratch_registry%relinquish_field(scratch_idx)
  end subroutine scalar_bc_resolver_debug_output

end module scalar_bc_resolver
