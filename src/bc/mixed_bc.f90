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
!   * Redistributions in binary form must reproduce the above copyright
!     notice, this list of conditions and the following disclaimer in the
!     documentation and/or other materials provided with the distribution.
!
!   * Neither the name of the authors nor the names of its contributors may be
!     used to endorse or promote products derived from this software without
!     specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!> Implements `mixed_bc_t`.
module mixed_bc
  use bc, only : bc_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : DEVICE_TO_HOST
  use mask, only : mask_t
  use matrix, only : matrix_t
  use field, only : field_t
  use field_list, only : field_list_t
  use fld_file, only : fld_file_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : cfill_mask, masked_scatter_copy
  use device_math, only : device_cfill_mask
  use field_math, only : field_rzero
  implicit none
  private

  !> Base type for mixed boundary conditions that need resolver-provided
  !! local-basis data on the physical field.
  type, public, abstract, extends(bc_t) :: mixed_bc_t
     !> Resolved subset of local dofs where this mixed bc remains active after
     !! global boundary-condition resolution. Populated by a
     !! `coupled_vector_bc_resolver_t`.
     !! Note that the resolved mask is only used for contraining vector
     !! components. For the Neumann contributions of the bc, the face_msk is
     !! used.
     type(mask_t) :: resolved_msk
     !> Local basis vectors on resolved_msk.
     type(matrix_t) :: n
     type(matrix_t) :: t1
     type(matrix_t) :: t2
   contains
     !> Destructor for mixed-bc-specific resolved state.
     procedure, pass(this) :: free_mixed => mixed_bc_free
     !> Write debug fields for the resolved support and normals.
     procedure, pass(this) :: debug_output => mixed_bc_debug_output
  end type mixed_bc_t

contains

  !> Destructor for the mixed-bc extension state.
  subroutine mixed_bc_free(this)
    class(mixed_bc_t), intent(inout) :: this

    call this%free_base()
    call this%resolved_msk%free()
    call this%n%free()
    call this%t1%free()
    call this%t2%free()
  end subroutine mixed_bc_free

  !> Write debug fields for a mixed bc.
  !! @details Writes one field list with the resolved mask marked by `1`
  !! and the three components of the resolved normals.
  !! @param[in] field_name Base name for the output files.
  subroutine mixed_bc_debug_output(this, field_name)
    class(mixed_bc_t), intent(inout) :: this
    character(len=*), intent(in), optional :: field_name
    type(field_t), pointer :: mask_field, nx_field, ny_field, nz_field
    type(field_list_t) :: debug_fields
    type(fld_file_t) :: dump_file
    integer :: scratch_idx(4)
    character(len=:), allocatable :: field_name_

    if (present(field_name)) then
       field_name_ = trim(field_name)
    else if (allocated(this%name) .and. len_trim(this%name) .gt. 0) then
       field_name_ = trim(this%name)
    else
       field_name_ = 'mixed_bc'
    end if

    call neko_scratch_registry%request_field(mask_field, scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(nx_field, scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(ny_field, scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(nz_field, scratch_idx(4), .true.)

    call field_rzero(mask_field)
    call field_rzero(nx_field)
    call field_rzero(ny_field)
    call field_rzero(nz_field)

    if (this%resolved_msk%is_set()) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_cfill_mask(mask_field%x_d, 1.0_rp, mask_field%size(), &
               this%resolved_msk%get_d(), this%resolved_msk%size())
          call mask_field%copy_from(DEVICE_TO_HOST, .true.)
       else
          call cfill_mask(mask_field%x, 1.0_rp, mask_field%size(), &
               this%resolved_msk%get(), this%resolved_msk%size())
       end if

       if (this%n%size() .gt. 0) then
          call this%n%copy_from(DEVICE_TO_HOST, .true.)

          call masked_scatter_copy(nx_field%x(:,1,1,1), this%n%x(1,:), &
               this%resolved_msk%get(), nx_field%size(), &
               this%resolved_msk%size())
          call masked_scatter_copy(ny_field%x(:,1,1,1), this%n%x(2,:), &
               this%resolved_msk%get(), ny_field%size(), &
               this%resolved_msk%size())
          call masked_scatter_copy(nz_field%x(:,1,1,1), this%n%x(3,:), &
               this%resolved_msk%get(), nz_field%size(), &
               this%resolved_msk%size())
       end if
    end if

    call debug_fields%init(4)
    call debug_fields%assign(1, mask_field)
    call debug_fields%assign(2, nx_field)
    call debug_fields%assign(3, ny_field)
    call debug_fields%assign(4, nz_field)

    call dump_file%init(field_name_ // '.fld')
    call dump_file%write(debug_fields)
    call debug_fields%free()

    call neko_scratch_registry%relinquish_field(scratch_idx)
  end subroutine mixed_bc_debug_output
end module mixed_bc
