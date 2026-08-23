! Copyright (c) 2020-2026, The Neko Authors
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
!> Slip boundary condition for arbitrarily oriented stationary walls.
module slip
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use facet_normal, only : facet_normal_t
  use vector, only : vector_t
  use coefs, only : coef_t
  use json_module, only : json_file
  use time_state, only : time_state_t
  use utils, only : neko_error
  use device, only : device_memcpy, HOST_TO_DEVICE
  use device_math, only : device_col3, device_addcol3, device_subcol3, &
       device_masked_gather_copy_0, device_masked_scatter_copy_0
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Impermeable, free-slip wall.
  !! Removes the local wall-normal component of velocity while preserving
  !! both tangential components.
  type, public, extends(facet_normal_t) :: slip_t
     type(vector_t) :: work_x, work_y, work_z
   contains
     procedure, pass(this) :: init => slip_init
     procedure, pass(this) :: finalize => slip_finalize
     procedure, pass(this) :: free => slip_free
     procedure, pass(this) :: recompute_normals => slip_recompute_normals
     procedure, pass(this) :: apply_vector => slip_apply_vector
     procedure, pass(this) :: apply_vector_dev => slip_apply_vector_dev
  end type slip_t

contains

  !> Constructor.
  subroutine slip_init(this, coef, json)
    class(slip_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%facet_normal_t%init(coef, json)

    ! Slip is a mixed condition: normal velocity is prescribed, whereas
    ! tangential velocity is unconstrained.
    this%strong = .false.
  end subroutine slip_init

  !> Finalize the mask and construct unit normals.
  subroutine slip_finalize(this, only_facets)
    class(slip_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    integer :: m

    if (present(only_facets)) then
       call this%facet_normal_t%finalize(only_facets)
    else
       call this%facet_normal_t%finalize()
    end if
    call slip_normalize_normals(this)

    m = this%unique_mask(0)
    if (m .gt. 0) then
       call this%work_x%init(m)
       call this%work_y%init(m)
       call this%work_z%init(m)
    end if
  end subroutine slip_finalize

  !> Recompute and normalize wall normals after a mesh update.
  subroutine slip_recompute_normals(this)
    class(slip_t), target, intent(inout) :: this

    call this%facet_normal_t%recompute_normals()
    call slip_normalize_normals(this)
  end subroutine slip_recompute_normals

  !> Normalize the area-weighted nodal normals constructed by facet_normal_t.
  subroutine slip_normalize_normals(this)
    class(slip_t), target, intent(inout) :: this
    integer :: i, m
    real(kind=rp) :: normal_magnitude

    if (.not. allocated(this%unique_mask)) return

    m = this%unique_mask(0)
    do i = 1, m
       normal_magnitude = sqrt(this%nx%x(i)**2 + this%ny%x(i)**2 + &
            this%nz%x(i)**2)
       if (normal_magnitude .le. tiny(normal_magnitude)) then
          call neko_error("Slip boundary condition encountered a zero normal.")
       end if
       this%nx%x(i) = this%nx%x(i) / normal_magnitude
       this%ny%x(i) = this%ny%x(i) / normal_magnitude
       this%nz%x(i) = this%nz%x(i) / normal_magnitude
    end do

    if (NEKO_BCKND_DEVICE .eq. 1 .and. m .gt. 0) then
       call device_memcpy(this%nx%x, this%nx%x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%ny%x, this%ny%x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%nz%x, this%nz%x_d, m, HOST_TO_DEVICE, &
            sync = .true.)
    end if
  end subroutine slip_normalize_normals

  !> Project velocity onto the local tangent plane.
  subroutine slip_apply_vector(this, x, y, z, n, time, strong)
    class(slip_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x, y, z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i, k, m
    real(kind=rp) :: normal_velocity

    strong_ = .true.
    if (present(strong)) strong_ = strong
    if (.not. strong_) return

    m = this%unique_mask(0)
    !$omp parallel do private(k, normal_velocity)
    do i = 1, m
       k = this%unique_mask(i)
       normal_velocity = x(k) * this%nx%x(i) + &
            y(k) * this%ny%x(i) + z(k) * this%nz%x(i)
       x(k) = x(k) - normal_velocity * this%nx%x(i)
       y(k) = y(k) - normal_velocity * this%ny%x(i)
       z(k) = z(k) - normal_velocity * this%nz%x(i)
    end do
    !$omp end parallel do
  end subroutine slip_apply_vector

  !> Project velocity onto the local tangent plane on a device.
  subroutine slip_apply_vector_dev(this, x_d, y_d, z_d, time, strong, strm)
    class(slip_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d, y_d, z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_
    integer :: m, n

    strong_ = .true.
    if (present(strong)) strong_ = strong
    if (.not. strong_) return

    m = this%unique_mask(0)
    if (m .eq. 0) return
    n = this%coef%dof%size()

    call device_masked_gather_copy_0(this%work_x%x_d, x_d, &
         this%unique_mask_d, n, m, strm)
    call device_masked_gather_copy_0(this%work_y%x_d, y_d, &
         this%unique_mask_d, n, m, strm)
    call device_masked_gather_copy_0(this%work_z%x_d, z_d, &
         this%unique_mask_d, n, m, strm)

    call device_col3(this%work%x_d, this%work_x%x_d, this%nx%x_d, m, strm)
    call device_addcol3(this%work%x_d, this%work_y%x_d, this%ny%x_d, &
         m, strm)
    call device_addcol3(this%work%x_d, this%work_z%x_d, this%nz%x_d, &
         m, strm)
    call device_subcol3(this%work_x%x_d, this%work%x_d, this%nx%x_d, &
         m, strm)
    call device_subcol3(this%work_y%x_d, this%work%x_d, this%ny%x_d, &
         m, strm)
    call device_subcol3(this%work_z%x_d, this%work%x_d, this%nz%x_d, &
         m, strm)

    call device_masked_scatter_copy_0(x_d, this%work_x%x_d, &
         this%unique_mask_d, n, m, strm)
    call device_masked_scatter_copy_0(y_d, this%work_y%x_d, &
         this%unique_mask_d, n, m, strm)
    call device_masked_scatter_copy_0(z_d, this%work_z%x_d, &
         this%unique_mask_d, n, m, strm)
  end subroutine slip_apply_vector_dev

  !> Destructor.
  subroutine slip_free(this)
    class(slip_t), target, intent(inout) :: this

    call this%work_x%free()
    call this%work_y%free()
    call this%work_z%free()
    call this%facet_normal_t%free()
  end subroutine slip_free

end module slip
