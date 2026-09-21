!> @file normal_vec_bcs.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!! Dirichlet condition applied in the facet normal direction
module normal_vec_bcs
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use vector, only : vector_t
  use coefs, only : coef_t
  use bc, only : bc_t
  use utils, only : neko_error, nonlinear_index
  use json_module, only : json_file
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated
  use htable, only : htable_i4_t
  use device, only : device_map, device_memcpy, device_unmap, &
       HOST_TO_DEVICE
  use time_state, only : time_state_t
  implicit none
  private

  !> Dirichlet condition in facet normal direction
  type, public, extends(bc_t) :: normal_vec_bcs_t
     integer, allocatable :: unique_mask(:)
     type(c_ptr) :: unique_mask_d = c_null_ptr
     type(vector_t) :: nx, ny, nz, work
   contains
     procedure, pass(this) :: apply_scalar => normal_vec_bcs_apply_scalar
     procedure, pass(this) :: apply_scalar_dev => &
          normal_vec_bcs_apply_scalar_dev
     procedure, pass(this) :: apply_vector => normal_vec_bcs_apply_vector
     procedure, pass(this) :: apply_vector_dev => &
          normal_vec_bcs_apply_vector_dev
     procedure, pass(this) :: apply_n_dot => normal_vec_bcs_apply_n_dot
     procedure, pass(this) :: apply_n_cross => normal_vec_bcs_apply_n_cross
     !> Constructor.
     procedure, pass(this) :: init => normal_vec_bcs_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          normal_vec_bcs_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => normal_vec_bcs_free
     !> Finalize.
     procedure, pass(this) :: finalize => normal_vec_bcs_finalize
  end type normal_vec_bcs_t

contains

  !> Constructor.
  !! @param[inout] this The boundary condition object.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine normal_vec_bcs_init(this, coef, json)
    class(normal_vec_bcs_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_from_components(coef)
  end subroutine normal_vec_bcs_init

  !> Constructor from components.
  !! @param[inout] this The boundary condition object.
  !! @param[in] coef The SEM coefficients.
  subroutine normal_vec_bcs_init_from_components(this, coef)
    class(normal_vec_bcs_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef

    call this%init_base(coef)
  end subroutine normal_vec_bcs_init_from_components

  !> No-op scalar apply.
  !! @param[inout] this The boundary condition object.
  !! @param[inout] x Scalar field values.
  !! @param[in] n Number of scalar entries.
  !! @param[in] time Time state (unused).
  !! @param[in] strong Strong enforcement flag (unused).
  subroutine normal_vec_bcs_apply_scalar(this, x, n, time, strong)
    class(normal_vec_bcs_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine normal_vec_bcs_apply_scalar

  !> No-op scalar apply on device.
  !! @param[inout] this The boundary condition object.
  !! @param[inout] x_d Device pointer to scalar field values.
  !! @param[in] time Time state (unused).
  !! @param[in] strong Strong enforcement flag (unused).
  !! @param[inout] strm Device stream handle.
  subroutine normal_vec_bcs_apply_scalar_dev(this, x_d, time, strong, strm)
    class(normal_vec_bcs_t), intent(inout), target :: this
    type(c_ptr),intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

  end subroutine normal_vec_bcs_apply_scalar_dev

  !> No-op vector apply on device.
  !! @param[inout] this The boundary condition object.
  !! @param[inout] x_d Device pointer to x-component.
  !! @param[inout] y_d Device pointer to y-component.
  !! @param[inout] z_d Device pointer to z-component.
  !! @param[in] time Time state (unused).
  !! @param[in] strong Strong enforcement flag (unused).
  !! @param[inout] strm Device stream handle.
  subroutine normal_vec_bcs_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(normal_vec_bcs_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

  end subroutine normal_vec_bcs_apply_vector_dev

  !> No-op vector apply.
  !! @param[inout] this The boundary condition object.
  !! @param[inout] x Vector x-component.
  !! @param[inout] y Vector y-component.
  !! @param[inout] z Vector z-component.
  !! @param[in] n Number of entries.
  !! @param[in] time Time state (unused).
  !! @param[in] strong Strong enforcement flag (unused).
  subroutine normal_vec_bcs_apply_vector(this, x, y, z, n, time, strong)
    class(normal_vec_bcs_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine normal_vec_bcs_apply_vector

  !> Apply the dot product of a vector field with facet normal.
  !! @param[in] this The boundary condition object.
  !! @param[inout] f Output array for the dot product.
  !! @param[in] u Vector x-component.
  !! @param[in] v Vector y-component.
  !! @param[in] w Vector z-component.
  !! @param[in] n Number of entries.
  !! @param[in] time Time state (unused).
  subroutine normal_vec_bcs_apply_n_dot(this, f, u, v, w, n, time)
    class(normal_vec_bcs_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: f
    real(kind=rp), intent(inout), dimension(n) :: u
    real(kind=rp), intent(inout), dimension(n) :: v
    real(kind=rp), intent(inout), dimension(n) :: w
    type(time_state_t), intent(in), optional :: time
    integer :: i, m, k

    m = this%unique_mask(0)

    ! This should be implemented n dot u weighted by the 2D mass matrix, which
    ! I believe is the area, and it appears that n is weighted by the area.
    do i = 1, m
       k = this%unique_mask(i)
       f(k) = u(k) * this%nx%x(i) + v(k) * this%ny%x(i) + w(k) * this%nz%x(i)
    end do

  end subroutine normal_vec_bcs_apply_n_dot

  !> Apply n x u.
  !! @param[in] this The boundary condition object.
  !! @param[inout] x Output x-component of n x u.
  !! @param[inout] y Output y-component of n x u.
  !! @param[inout] z Output z-component of n x u.
  !! @param[in] u Vector x-component.
  !! @param[in] v Vector y-component.
  !! @param[in] w Vector z-component.
  !! @param[in] n Number of entries.
  !! @param[in] time Time state (unused).
  subroutine normal_vec_bcs_apply_n_cross(this, x, y, z, u, v, w, n, time)
    class(normal_vec_bcs_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(inout), dimension(n) :: u
    real(kind=rp), intent(inout), dimension(n) :: v
    real(kind=rp), intent(inout), dimension(n) :: w
    type(time_state_t), intent(in), optional :: time
    integer :: i, m, k

    m = this%unique_mask(0)

    do i = 1, m
       k = this%unique_mask(i)
       x(k) = this%ny%x(i) * w(k) - this%nz%x(i) * v(k)
       y(k) = this%nz%x(i) * u(k) - this%nx%x(i) * w(k)
       z(k) = this%nx%x(i) * v(k) - this%ny%x(i) * u(k)
    end do

  end subroutine normal_vec_bcs_apply_n_cross

  !> Destructor.
  !! @param[inout] this The boundary condition object.
  subroutine normal_vec_bcs_free(this)
    class(normal_vec_bcs_t), target, intent(inout) :: this

    call this%free_base()
    if (allocated(this%unique_mask)) then
       if (c_associated(this%unique_mask_d)) then
          call device_unmap(this%unique_mask, this%unique_mask_d)
       end if
       deallocate(this%unique_mask)
    end if

    call this%nx%free()
    call this%ny%free()
    call this%nz%free()
    call this%work%free()

  end subroutine normal_vec_bcs_free

  !> Finalize.
  !! @param[inout] this The boundary condition object.
  subroutine normal_vec_bcs_finalize(this)
    class(normal_vec_bcs_t), target, intent(inout) :: this
    type(htable_i4_t) :: unique_point_idx
    integer :: htable_data, rcode, i, j, idx(4), facet
    real(kind=rp) :: area, normal(3)

    call this%finalize_base()
    ! This part is purely needed to ensure that contributions
    ! for all faces a point is on is properly summed up.
    ! If one simply uses the original mask, if a point is on a corner
    ! where both faces are on the boundary
    ! one will only get the contribution from one face, not both
    ! We solve this by adding up the normals of both faces for these points
    ! and storing this sum in this%nx, this%ny, this%nz.
    ! As both contrbutions are added already,
    ! we also ensure that we only visit each point once
    ! and create a new mask with only unique points (this%unique_mask).
    if (allocated(this%unique_mask)) then
       if (c_associated(this%unique_mask_d)) then
          call device_unmap(this%unique_mask, this%unique_mask_d)
       end if
       deallocate(this%unique_mask)
    end if

    call unique_point_idx%init(this%facet_node_msk(0), htable_data)
    j = 0
    do i = 1, this%facet_node_msk(0)
       if (unique_point_idx%get(this%facet_node_msk(i), htable_data) .ne. 0) then
          j = j + 1
          htable_data = j
          call unique_point_idx%set(this%facet_node_msk(i), j)
       end if
    end do

    ! Only allocate work vectors if size is non-zero
    if (unique_point_idx%num_entries() .gt. 0 ) then
       call this%nx%init(unique_point_idx%num_entries())
       call this%ny%init(unique_point_idx%num_entries())
       call this%nz%init(unique_point_idx%num_entries())
       call this%work%init(unique_point_idx%num_entries())
    end if
    allocate(this%unique_mask(0:unique_point_idx%num_entries()))

    this%unique_mask(0) = unique_point_idx%num_entries()
    do i = 1, this%unique_mask(0)
       this%unique_mask(i) = 0
    end do


    do i = 1, this%facet_node_msk(0)
       rcode = unique_point_idx%get(this%facet_node_msk(i), htable_data)
       if (rcode .ne. 0) call neko_error("Facet normal: htable get failed.")
       this%unique_mask(htable_data) = this%facet_node_msk(i)
       facet = this%facet(i)

       idx = nonlinear_index(this%facet_node_msk(i), this%Xh%lx, this%Xh%lx, &
            this%Xh%lx)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), facet)
       area = this%coef%get_area(idx(1), idx(2), idx(3), idx(4), facet)
       normal = normal * area !Scale normal by area
       this%nx%x(htable_data) = this%nx%x(htable_data) + normal(1)
       this%ny%x(htable_data) = this%ny%x(htable_data) + normal(2)
       this%nz%x(htable_data) = this%nz%x(htable_data) + normal(3)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1 .and. &
         (unique_point_idx%num_entries() .gt. 0 )) then
       call device_map(this%unique_mask, this%unique_mask_d, &
            size(this%unique_mask))
       call device_memcpy(this%unique_mask, this%unique_mask_d, &
            size(this%unique_mask), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%nx%x, this%nx%x_d, &
            this%nx%size(), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%ny%x, this%ny%x_d, &
            this%ny%size(), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%nz%x, this%nz%x_d, &
            this%nz%size(), HOST_TO_DEVICE, sync = .true.)
    end if

    call unique_point_idx%free()

  end subroutine normal_vec_bcs_finalize

end module normal_vec_bcs
