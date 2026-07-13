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
!> Implements sampling from GLL nodes at prescribed off-wall indices.
module wall_gll_sampler
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use wall_sampler, only : wall_sampler_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error, nonlinear_index, linear_index
  use math, only : masked_gather_copy
  use device_math, only : device_masked_gather_copy_aligned
  use device, only : device_map, device_unmap, device_memcpy, HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  !> Wall sampler implementation for GLL nodes at prescribed off-wall indices.
  type, public, extends(wall_sampler_t) :: wall_gll_sampler_t
     !> Off-wall indices for each sampling point.
     integer, allocatable :: indices(:)
     !> Linear indices into the sampled fields.
     integer, allocatable :: sample_idx(:)
     !> Device pointer to linear indices into the sampled fields.
     type(c_ptr) :: sample_idx_d = C_NULL_PTR
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => wall_gll_sampler_init
     !> Constructor from off-wall indices.
     procedure, pass(this) :: init_from_indices => &
          wall_gll_sampler_init_from_indices
     procedure, pass(this) :: finalize => wall_gll_sampler_finalize
     !> Sample the solution field at the sampling points.
     procedure, pass(this) :: sample => wall_gll_sampler_sample
     !> Destructor.
     procedure, pass(this) :: free => wall_gll_sampler_free
  end type wall_gll_sampler_t

contains

  !> Constructor from JSON.
  subroutine wall_gll_sampler_init(this, json)
    class(wall_gll_sampler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, allocatable :: indices(:)
    integer :: index, var_type
    logical :: found

    call json%info('value', found = found, var_type = var_type)
    if (.not. found) then
       call neko_error('Wall GLL sampler requires a value')
    end if

    select case (var_type)
    case (3) ! json_array
       call json_get(json, 'value', indices)
       call this%init_from_indices(indices)
    case (5) ! json_integer
       call json_get(json, 'value', index)
       call this%init_from_indices([index])
    case default
       call neko_error('Wall GLL sampler value must be an integer or array')
    end select
  end subroutine wall_gll_sampler_init

  subroutine wall_gll_sampler_init_from_indices(this, indices)
    class(wall_gll_sampler_t), intent(inout) :: this
    integer, intent(in) :: indices(:)

    if (any(indices < 1)) then 
       call neko_error('Wall GLL sampler indices must be positive')
    end if

    this%indices = indices
    this%n_samples = size(indices)
  end subroutine wall_gll_sampler_init_from_indices

  subroutine wall_gll_sampler_finalize(this, coef, msk, facet, n_x, n_y, n_z)
    class(wall_gll_sampler_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(0:)
    integer, intent(in) :: facet(0:)
    type(vector_t), intent(in) :: n_x, n_y, n_z
    integer :: i, j, p, fid, idx(4), sample(4), lx, ly, lz
    real(kind=rp) :: xw, yw, zw, dx, dy, dz

    if (.not. allocated(this%indices)) then
       call neko_error('Wall GLL sampler has not been initialized')
    end if

    call this%h%free()
    this%n_nodes = msk(0)
    this%n_samples = size(this%indices)
    allocate(this%sample_idx(this%n_nodes * this%n_samples))
    call this%h%init(size(this%sample_idx))

    lx = coef%Xh%lx
    ly = coef%Xh%ly
    lz = coef%Xh%lz
    do i = 1, this%n_nodes
       idx = nonlinear_index(msk(i), lx, ly, lz)
       xw = coef%dof%x(idx(1), idx(2), idx(3), idx(4))
       yw = coef%dof%y(idx(1), idx(2), idx(3), idx(4))
       zw = coef%dof%z(idx(1), idx(2), idx(3), idx(4))
       fid = facet(i)
       do j = 1, this%n_samples
          p = (i - 1) * this%n_samples + j
          sample = idx
          select case (fid)
          case (1)
             sample(1) = sample(1) + this%indices(j)
          case (2)
             sample(1) = sample(1) - this%indices(j)
          case (3)
             sample(2) = sample(2) + this%indices(j)
          case (4)
             sample(2) = sample(2) - this%indices(j)
          case (5)
             sample(3) = sample(3) + this%indices(j)
          case (6)
             sample(3) = sample(3) - this%indices(j)
          case default
             call neko_error('Invalid facet in wall GLL sampler')
          end select

          if (sample(1) < 1 .or. sample(1) > lx .or. &
               sample(2) < 1 .or. sample(2) > ly .or. &
               sample(3) < 1 .or. sample(3) > lz) then
             call neko_error('Wall GLL sampling index lies outside element')
          end if

          this%sample_idx(p) = linear_index(sample(1), sample(2), sample(3), &
               sample(4), lx, ly, lz)

          dx = coef%dof%x(sample(1), sample(2), sample(3), sample(4)) - xw
          dy = coef%dof%y(sample(1), sample(2), sample(3), sample(4)) - yw
          dz = coef%dof%z(sample(1), sample(2), sample(3), sample(4)) - zw
          this%h%x(p) = -(dx*n_x%x(i) + dy*n_y%x(i) + dz*n_z%x(i))
       end do
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%sample_idx, this%sample_idx_d, size(this%sample_idx))
       call device_memcpy(this%sample_idx, this%sample_idx_d, &
            size(this%sample_idx), HOST_TO_DEVICE, sync = .false.)
       call this%h%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine wall_gll_sampler_finalize

  !> Sample the solution field at the sampling points.
  !! @param field Field to be sampled.
  !! @param values Output sampled values in sampler ordering.
  subroutine wall_gll_sampler_sample(this, field, values)
    class(wall_gll_sampler_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%n_nodes * this%n_samples

    if (values%size() .ne. n) then
       call neko_error('Wall sampler destination has an invalid size')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_aligned(values%x_d, field%x_d, &
            this%sample_idx_d, field%size(), n)
    else
       call masked_gather_copy(values%x, field%x, this%sample_idx, &
            field%size(), n)
    end if
  end subroutine wall_gll_sampler_sample

  !> Destructor
  subroutine wall_gll_sampler_free(this)
    class(wall_gll_sampler_t), intent(inout) :: this

    if (allocated(this%sample_idx)) then
       if (c_associated(this%sample_idx_d)) then
          call device_unmap(this%sample_idx, this%sample_idx_d)
       end if
       deallocate(this%sample_idx)
    end if

    if (allocated(this%indices)) deallocate(this%indices)

    call this%h%free()

    this%sample_idx_d = C_NULL_PTR
    this%n_nodes = 0
    this%n_samples = 0
  end subroutine wall_gll_sampler_free

end module wall_gll_sampler
