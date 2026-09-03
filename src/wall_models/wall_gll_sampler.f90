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
  use json_utils, only : json_get, json_get_or_default
  use wall_sampler, only : wall_sampler_t
  use user_intf, only : user_t
  use mask, only : mask_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error, nonlinear_index, linear_index
  use math, only : masked_gather_copy
  use device_math, only : device_masked_gather_copy_aligned
  use device, only : HOST_TO_DEVICE
  implicit none
  private

  !> Wall sampler implementation for GLL nodes at prescribed off-wall indices.
  type, public, extends(wall_sampler_t) :: wall_gll_sampler_t
     !> Off-wall GLL indices for each sampling point per wall node.
     integer, allocatable :: indices(:)
     !> Linear indices into the sampled fields.
     type(mask_t) :: sample_idx
   contains
     !> Partial constructor from JSON.
     procedure, pass(this) :: init => wall_gll_sampler_init
     !> Partial constructor from off-wall GLL indices.
     procedure, pass(this) :: init_from_indices => &
          wall_gll_sampler_init_from_indices
     !> Partial constructor for user-provided GLL indices.
     procedure, private, pass(this) :: init_from_user => &
          wall_gll_sampler_init_from_user
     !> Fully construct a sampler from its resolved components.
     procedure, pass(this) :: init_from_components => &
          wall_gll_sampler_init_from_components
     !> Construct field indices and wall-normal distances for sampling.
     procedure, pass(this) :: finalize => wall_gll_sampler_finalize
     !> Sample the solution field at the sampling points.
     procedure, pass(this) :: sample => wall_gll_sampler_sample
     !> Destructor.
     procedure, pass(this) :: free => wall_gll_sampler_free
  end type wall_gll_sampler_t

contains

  !> Partial constructor from JSON.
  !! @param json Sampler configuration data.
  subroutine wall_gll_sampler_init(this, json)
    class(wall_gll_sampler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, allocatable :: indices(:)
    integer :: index, n_samples, var_type
    logical :: found, output_h
    character(len=:), allocatable :: value

    call json_get_or_default(json, 'output_h', output_h, .true.)

    call json%info('value', found = found, var_type = var_type)
    if (.not. found) then
       call neko_error('Wall GLL sampler requires a value')
    end if

    select case (var_type)
    case (3) ! json_array
       call json_get(json, 'value', indices)
       call this%init_from_indices(indices, output_h)
    case (5) ! json_integer
       call json_get(json, 'value', index)
       call this%init_from_indices([index], output_h)
    case (7) ! json_string
       call json_get(json, 'value', value)
       if (trim(value) /= 'user') then
          call neko_error('Wall GLL sampler value must be an integer, ' // &
               'array, or user')
       end if
       call json_get(json, 'n_samples', n_samples)
       call this%init_from_user(n_samples, output_h)
    case default
       call neko_error('Wall GLL sampler value must be an integer, ' // &
            'array, or user')
    end select
  end subroutine wall_gll_sampler_init

  !> Partial constructor from off-wall GLL indices.
  !! @param indices Positive GLL indices, one for each sampling point per
  !! wall node.
  !! @param output_h Whether to write the sampling-distance diagnostic.
  subroutine wall_gll_sampler_init_from_indices(this, indices, output_h)
    class(wall_gll_sampler_t), intent(inout) :: this
    integer, intent(in) :: indices(:)
    logical, intent(in) :: output_h

    if (any(indices < 1)) then
       call neko_error('Wall GLL sampler indices must be positive')
    end if

    this%indices = indices
    this%n_samples = size(indices)
    this%user_values = .false.
    this%output_h_enabled = output_h
  end subroutine wall_gll_sampler_init_from_indices

  !> Partial constructor for user-provided GLL sampling indices.
  !! @param n_samples Number of samples per wall node.
  !! @param output_h Whether to write the sampling-distance diagnostic.
  subroutine wall_gll_sampler_init_from_user(this, n_samples, output_h)
    class(wall_gll_sampler_t), intent(inout) :: this
    integer, intent(in) :: n_samples
    logical, intent(in) :: output_h

    if (n_samples < 1) then
       call neko_error('Wall GLL sampler n_samples must be positive')
    end if

    this%n_samples = n_samples
    this%user_values = .true.
    this%output_h_enabled = output_h
  end subroutine wall_gll_sampler_init_from_user

  !> Fully construct a GLL sampler from its resolved components.
  !! @note No compatibility checks are performed, input data is assumed to be
  !! valid. Checking would require geometric information.
  !! @param indices Off-wall GLL indices for each sample per wall node.
  !! @param sample_idx Linear sampled-field indices in sampler ordering.
  !! @param h Wall-normal distances in sampler ordering.
  !! @param output_h Sampling-distance diagnostic output flag.
  subroutine wall_gll_sampler_init_from_components(this, indices, sample_idx, &
       h, output_h)
    class(wall_gll_sampler_t), intent(inout) :: this
    integer, intent(in) :: indices(:)
    integer, intent(in) :: sample_idx(:)
    type(vector_t), intent(in) :: h
    logical, intent(in) :: output_h
    integer :: n_nodes

    if (size(indices) < 1 .or. any(indices < 1)) then
       call neko_error('Wall GLL sampler indices must be positive')
    end if
    if (size(sample_idx) < 1 .or. any(sample_idx < 1)) then
       call neko_error('Wall GLL sampler sample indices must be positive')
    end if
    if (mod(size(sample_idx), size(indices)) /= 0) then
       call neko_error('Wall GLL sampler sample indices have an invalid size')
    end if

    n_nodes = size(sample_idx) / size(indices)
    call this%free()
    call this%init_base(n_nodes, size(indices), h, output_h)
    this%indices = indices
    call this%sample_idx%init(sample_idx, size(sample_idx))
  end subroutine wall_gll_sampler_init_from_components

  !> Construct field indices and wall-normal distances for sampling. Completes
  !! initialization.
  !! @param coef SEM coefficients.
  !! @param msk Mask selecting local wall nodes.
  !! @param facet Facet index for each selected wall node.
  !! @param n_x X-component of wall-normal vectors at wall nodes.
  !! @param n_y Y-component of wall-normal vectors at wall nodes.
  !! @param n_z Z-component of wall-normal vectors at wall nodes.
  !! @param bc_name Name of the owning boundary condition for user values.
  !! @param user User interface providing the GLL-sampling callback.
  subroutine wall_gll_sampler_finalize(this, coef, msk, facet, n_x, n_y, &
       n_z, bc_name, user)
    class(wall_gll_sampler_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(0:)
    integer, intent(in) :: facet(0:)
    type(vector_t), intent(in) :: n_x, n_y, n_z
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user
    integer :: i, j, p, fid, idx(4), sample(4), lx, ly, lz, offwall
    real(kind=rp) :: xw, yw, zw, dx, dy, dz
    integer, allocatable :: user_indices(:,:), sample_idx(:)

    if (.not. this%user_values .and. .not. allocated(this%indices)) then
       call neko_error('Wall GLL sampler has not been initialized')
    end if

    call this%h%free()
    this%n_nodes = msk(0)
    if (.not. this%user_values) this%n_samples = size(this%indices)

    if (this%user_values) then
       if (.not. present(bc_name) .or. .not. present(user)) then
          call neko_error('Wall GLL user sampler has no configured callback')
       end if

       allocate(user_indices(this%n_samples, this%n_nodes))

       call user%wall_sampling_gll(bc_name, msk(1:this%n_nodes), user_indices)

       if (any(user_indices < 1)) then
          call neko_error('Wall GLL sampler indices must be positive')
       end if
    end if

    allocate(sample_idx(this%n_nodes * this%n_samples))
    call this%h%init(size(sample_idx))

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
          if (this%user_values) then
             offwall = user_indices(j,i)
          else
             offwall = this%indices(j)
          end if
          sample = idx
          select case (fid)
          case (1)
             sample(1) = sample(1) + offwall
          case (2)
             sample(1) = sample(1) - offwall
          case (3)
             sample(2) = sample(2) + offwall
          case (4)
             sample(2) = sample(2) - offwall
          case (5)
             sample(3) = sample(3) + offwall
          case (6)
             sample(3) = sample(3) - offwall
          case default
             call neko_error('Invalid facet in wall GLL sampler')
          end select

          if (sample(1) < 1 .or. sample(1) > lx .or. &
               sample(2) < 1 .or. sample(2) > ly .or. &
               sample(3) < 1 .or. sample(3) > lz) then
             call neko_error('Wall GLL sampling index lies outside element')
          end if

          sample_idx(p) = linear_index(sample(1), sample(2), sample(3), &
               sample(4), lx, ly, lz)

          dx = coef%dof%x(sample(1), sample(2), sample(3), sample(4)) - xw
          dy = coef%dof%y(sample(1), sample(2), sample(3), sample(4)) - yw
          dz = coef%dof%z(sample(1), sample(2), sample(3), sample(4)) - zw
          this%h%x(p) = -(dx*n_x%x(i) + dy*n_y%x(i) + dz*n_z%x(i))
       end do
    end do

    call this%sample_idx%init(sample_idx, size(sample_idx))
    if (NEKO_BCKND_DEVICE .eq. 1) &
         call this%h%copy_from(HOST_TO_DEVICE, sync = .true.)

    if (this%output_h_enabled .and. present(bc_name)) then
       call this%output_h(coef, msk, bc_name)
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
            this%sample_idx%get_d(), field%size(), n)
    else
       call masked_gather_copy(values%x, field%x, this%sample_idx%get(), &
            field%size(), n)
    end if
  end subroutine wall_gll_sampler_sample

  !> Destructor
  subroutine wall_gll_sampler_free(this)
    class(wall_gll_sampler_t), intent(inout) :: this

    call this%sample_idx%free()

    if (allocated(this%indices)) deallocate(this%indices)

    call this%h%free()

    this%n_nodes = 0
    this%n_samples = 0
    this%user_values = .false.
    this%output_h_enabled = .true.
  end subroutine wall_gll_sampler_free

end module wall_gll_sampler
