! Copyright (c) 2020-2025, The Neko Authors
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
!> Implements `non_normal_aligned_t`.
module non_normal_aligned
  use json_module, only : json_file
  use bc, only : bc_t, BC_TYPES
  use dirichlet, only : dirichlet_t
  use device_inhom_dirichlet, only : device_inhom_dirichlet_apply_scalar
  use import_field_utils, only : import_fields
  use num_types, only : rp, dp
  use tuple, only : tuple_i4_t
  use coefs, only : coef_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup
  use utils, only : neko_error
  use vector, only : vector_t
  use vector_math, only : vector_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry
  use field, only : field_t
  use global_interpolation, only : GLOB_INTERP_TOL, GLOB_INTERP_PAD
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Axis-aligned mixed Dirichlet condition in the non-normal direction.
  !! @warning Only works for axis-aligned plane boundaries.
  !! @details Since `dirichlet_t`currently only supports constant values,
  !! we store the values separately in vectors, and use the nested `dirichlet_t`
  !! only for its mask.
  type, public, extends(bc_t) :: non_normal_aligned_t
     !> Nested dirichlet bcs for each component. Used only for their masks.
     type(dirichlet_t) :: bc_x
     type(dirichlet_t) :: bc_y
     type(dirichlet_t) :: bc_z
     !> The applied values for each component. Allocated at finalization.
     type(vector_t) :: value_x
     type(vector_t) :: value_y
     type(vector_t) :: value_z
     !> Helper to parse uniform input from the case file.
     real(kind=rp), private :: constant_value(3) = 0.0_rp

     !> Bookkeeping for initialization. Necessary because we can only init
     !! in finalize()
     logical, private :: use_constant_value = .true.
     logical, private :: read_values_from_field = .false.
     logical, private :: field_interpolate = .false.
     character(len=:), allocatable, private :: field_file_name
     character(len=:), allocatable, private :: field_mesh_file_name
     real(kind=dp), private :: field_interp_tolerance = GLOB_INTERP_TOL
     real(kind=dp), private :: field_interp_padding = GLOB_INTERP_PAD
   contains
     !> No-op scalar application.
     procedure, pass(this) :: apply_scalar => non_normal_aligned_apply_scalar
     !> Apply the tangential components of the prescribed vector on the CPU.
     procedure, pass(this) :: apply_vector => non_normal_aligned_apply_vector
     !> No-op scalar application on the device.
     procedure, pass(this) :: apply_scalar_dev => &
          non_normal_aligned_apply_scalar_dev
     !> Apply the tangential components of the prescribed vector on the device.
     procedure, pass(this) :: apply_vector_dev => &
          non_normal_aligned_apply_vector_dev
     !> Construct the boundary condition from JSON.
     procedure, pass(this) :: init => non_normal_aligned_init
     !> Construct the boundary condition from a uniform global vector.
     procedure, pass(this) :: init_from_components => &
          non_normal_aligned_init_from_components
     !> Set inhomogeneous prescribed global vector components.
     procedure, pass(this) :: set_values => non_normal_aligned_set_values
     !> Free the boundary condition and its nested storage.
     procedure, pass(this) :: free => non_normal_aligned_free
     !> Finalize the boundary condition.
     procedure, pass(this) :: finalize => non_normal_aligned_finalize
     !> Estimate which global axis is normal to a marked facet.
     procedure, pass(this) :: get_normal_axis => &
          non_normal_aligned_get_normal_axis
  end type non_normal_aligned_t

contains

  !> Construct the boundary condition from JSON.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine non_normal_aligned_init(this, coef, json)
    class(non_normal_aligned_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: value(:)
    real(kind=rp) :: value_3(3)
    logical :: found_file_name, found_value

    call json%info("file_name", found = found_file_name)
    call json%info("value", found = found_value)

    if (found_file_name .and. found_value) then
       call neko_error("non_normal_aligned accepts either 'file_name' or " // &
            "'value', but not both.")
    end if

    if (found_file_name) then
       call this%free()
       call this%init_base(coef)
       call this%bc_x%init_from_components(coef, 0.0_rp)
       call this%bc_y%init_from_components(coef, 0.0_rp)
       call this%bc_z%init_from_components(coef, 0.0_rp)
       call json_get(json, "file_name", this%field_file_name)
       call json_get_or_default(json, "interpolate", this%field_interpolate, &
            .false.)
       call json_get_or_default(json, "mesh_file_name", &
            this%field_mesh_file_name, "none")
       call json_get_or_default(json, "interpolation.tolerance", &
            this%field_interp_tolerance, GLOB_INTERP_TOL)
       call json_get_or_default(json, "interpolation.padding", &
            this%field_interp_padding, GLOB_INTERP_PAD)
       this%read_values_from_field = .true.
       this%use_constant_value = .false.
       this%constraints = [.false., .true., .true.]
       this%bc_type = BC_TYPES%MIXED_CONSTRAINS_TANGENT
    else
       value_3 = 0.0_rp
       call json_get_or_lookup(json, "value", value)
       if (size(value) .ne. 3) then
          call neko_error("The non_normal boundary condition requires a " // &
               "3-component value vector.")
       end if
       value_3 = value

       call this%init_from_components(coef, value_3)
    end if
  end subroutine non_normal_aligned_init

  !> Construct the boundary condition from a uniform global vector.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] value Global vector whose tangential components are enforced.
  subroutine non_normal_aligned_init_from_components(this, coef, value)
    class(non_normal_aligned_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    real(kind=rp), intent(in) :: value(3)

    call this%free()
    call this%init_base(coef)
    call this%bc_x%init_from_components(coef, 0.0_rp)
    call this%bc_y%init_from_components(coef, 0.0_rp)
    call this%bc_z%init_from_components(coef, 0.0_rp)
    this%constant_value = value
    this%use_constant_value = .true.
    this%constraints = [.false., .true., .true.]
    this%bc_type = BC_TYPES%MIXED_CONSTRAINS_TANGENT
  end subroutine non_normal_aligned_init_from_components

  !> Set inhomogeneous prescribed global vector components.
  !! @details The vectors are stored in compact form and must match the
  !! component-wise mask sizes of `bc_x`, `bc_y` and `bc_z`. Therefore this
  !! routine must be called only after the nested masks have been finalized.
  !! @param[in] value_x Prescribed x values ordered like `bc_x%msk(1:)`.
  !! @param[in] value_y Prescribed y values ordered like `bc_y%msk(1:)`.
  !! @param[in] value_z Prescribed z values ordered like `bc_z%msk(1:)`.
  subroutine non_normal_aligned_set_values(this, value_x, value_y, value_z)
    class(non_normal_aligned_t), intent(inout) :: this
    type(vector_t), intent(in) :: value_x
    type(vector_t), intent(in) :: value_y
    type(vector_t), intent(in) :: value_z

    if (.not. allocated(this%bc_x%msk) .or. .not. allocated(this%bc_y%msk) &
         .or. .not. allocated(this%bc_z%msk)) then
       call neko_error("non_normal_aligned_set_values requires finalized " // &
            "nested boundary-condition masks.")
    end if

    if (value_x%size() .ne. this%bc_x%msk(0) .or. &
         value_y%size() .ne. this%bc_y%msk(0) .or. &
         value_z%size() .ne. this%bc_z%msk(0)) then
       call neko_error("non_normal_aligned_set_values requires compact " // &
            "component vectors matching bc_x, bc_y and bc_z mask sizes.")
    end if

    call this%value_x%free()
    call this%value_y%free()
    call this%value_z%free()
    this%value_x = value_x
    this%value_y = value_y
    this%value_z = value_z
    this%read_values_from_field = .false.
    this%use_constant_value = .false.
  end subroutine non_normal_aligned_set_values

  !> No-op scalar application.
  !! @param x Scalar field values.
  !! @param n Number of entries in `x`.
  !! @param time Current time state.
  !! @param strong Whether to apply the strong form.
  subroutine non_normal_aligned_apply_scalar(this, x, n, time, strong)
    class(non_normal_aligned_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine non_normal_aligned_apply_scalar

  !> Apply the tangential components of the prescribed vector on the CPU.
  !! @details Each nested Dirichlet boundary condition contributes only its
  !! compact mask. The corresponding component values are stored in the same
  !! compact ordering as that mask.
  !! @param x x-component of the field.
  !! @param y y-component of the field.
  !! @param z z-component of the field.
  !! @param n Number of entries in each component array.
  !! @param time Current time state.
  !! @param strong Whether to apply the strong form.
  subroutine non_normal_aligned_apply_vector(this, x, y, z, n, time, strong)
    class(non_normal_aligned_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    ! The local coordinates are aligned with the global axes, so the nested
    ! Dirichlet conditions only provide the masks. The prescribed values are
    ! stored compactly in the same order as each nested mask.
    if (strong_) then
       do i = 1, this%bc_x%msk(0)
          x(this%bc_x%msk(i)) = this%value_x%x(i)
       end do
       do i = 1, this%bc_y%msk(0)
          y(this%bc_y%msk(i)) = this%value_y%x(i)
       end do
       do i = 1, this%bc_z%msk(0)
          z(this%bc_z%msk(i)) = this%value_z%x(i)
       end do
    end if
  end subroutine non_normal_aligned_apply_vector

  !> No-op scalar application on the device.
  !! @param x_d Device pointer to the scalar field.
  !! @param time Current time state.
  !! @param strong Whether to apply the strong form.
  !! @param strm Device stream.
  subroutine non_normal_aligned_apply_scalar_dev(this, x_d, time, strong, strm)
    class(non_normal_aligned_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine non_normal_aligned_apply_scalar_dev

  !> Apply the tangential components of the prescribed vector on the device.
  !! @details Uses the compact per-component storage together with the masks of
  !! the nested Dirichlet boundary conditions.
  !! @param x_d Device pointer to the x-component.
  !! @param y_d Device pointer to the y-component.
  !! @param z_d Device pointer to the z-component.
  !! @param time Current time state.
  !! @param strong Whether to apply the strong form.
  !! @param strm Device stream.
  subroutine non_normal_aligned_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(non_normal_aligned_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (this%bc_x%msk(0) .gt. 0) then
          call device_inhom_dirichlet_apply_scalar(this%bc_x%msk_d, x_d, &
               this%value_x%x_d, this%bc_x%msk(0), strm)
       end if
       if (this%bc_y%msk(0) .gt. 0) then
          call device_inhom_dirichlet_apply_scalar(this%bc_y%msk_d, y_d, &
               this%value_y%x_d, this%bc_y%msk(0), strm)
       end if
       if (this%bc_z%msk(0) .gt. 0) then
          call device_inhom_dirichlet_apply_scalar(this%bc_z%msk_d, z_d, &
               this%value_z%x_d, this%bc_z%msk(0), strm)
       end if
    end if
  end subroutine non_normal_aligned_apply_vector_dev

  !> Estimate which global axis is normal to a marked facet.
  !! @param[out] sx Deviation from an x-aligned normal.
  !! @param[out] sy Deviation from a y-aligned normal.
  !! @param[out] sz Deviation from a z-aligned normal.
  !! @param[in] facet Facet id on the reference hex.
  !! @param[in] el Element id.
  subroutine non_normal_aligned_get_normal_axis(this, sx, sy, sz, facet, el)
    class(non_normal_aligned_t), target, intent(inout) :: this
    real(kind=rp), intent(out) :: sx, sy, sz
    integer, intent(in) :: facet
    integer, intent(in) :: el
    integer :: j, l

    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
         nz => this%coef%nz)
      sx = 0.0_rp
      sy = 0.0_rp
      sz = 0.0_rp
      select case (facet)
      case (1, 2)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx -1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      case (3, 4)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx - 1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      case (5, 6)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx - 1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      end select
      sx = sx / (c%Xh%lx - 2)**2
      sy = sy / (c%Xh%lx - 2)**2
      sz = sz / (c%Xh%lx - 2)**2
    end associate
  end subroutine non_normal_aligned_get_normal_axis

  !> Finalize the boundary condition.
  !! @details Marks the two tangential component conditions for each facet,
  !! finalizes the nested Dirichlet masks, allocates compact component-value
  !! storage, and populates it for the uniform-value case.
  subroutine non_normal_aligned_finalize(this)
    class(non_normal_aligned_t), target, intent(inout) :: this
    integer :: i
    integer :: scratch_idx(3)
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el
    type(field_t), pointer :: value_x_field, value_y_field, value_z_field
    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
         nz => this%coef%nz)
      bfp => this%marked_facet%array()
      do i = 1, this%marked_facet%size()
         bc_facet = bfp(i)
         facet = bc_facet%x(1)
         el = bc_facet%x(2)
         call this%get_normal_axis(sx, sy, sz, facet, el)

         if (sx .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sy .lt. TOL) then
            call this%bc_x%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sz .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
            call this%bc_x%mark_facet(facet, el)
         end if
      end do
    end associate
    call this%bc_x%finalize()
    call this%bc_y%finalize()
    call this%bc_z%finalize()

    call this%value_x%init(this%bc_x%msk(0), 'non_normal_aligned_x')
    call this%value_y%init(this%bc_y%msk(0), 'non_normal_aligned_y')
    call this%value_z%init(this%bc_z%msk(0), 'non_normal_aligned_z')

    if (this%use_constant_value) then
      this%value_x = this%constant_value(1)
      this%value_y = this%constant_value(2)
      this%value_z = this%constant_value(3)
    else if (this%read_values_from_field) then
      call neko_scratch_registry%request_field(value_x_field, scratch_idx(1), &
           .true.)
      call neko_scratch_registry%request_field(value_y_field, scratch_idx(2), &
           .true.)
      call neko_scratch_registry%request_field(value_z_field, scratch_idx(3), &
           .true.)

      call import_fields(this%field_file_name, this%field_mesh_file_name, &
           u = value_x_field, v = value_y_field, w = value_z_field, &
           interpolate = this%field_interpolate, &
           tolerance = this%field_interp_tolerance, &
           padding = this%field_interp_padding)

      call vector_masked_gather_copy_0(this%value_x, &
           value_x_field%x(:,1,1,1), this%bc_x%msk, value_x_field%size(), &
           this%bc_x%msk(0))
      call vector_masked_gather_copy_0(this%value_y, &
           value_y_field%x(:,1,1,1), this%bc_y%msk, value_y_field%size(), &
           this%bc_y%msk(0))
      call vector_masked_gather_copy_0(this%value_z, &
           value_z_field%x(:,1,1,1), this%bc_z%msk, value_z_field%size(), &
           this%bc_z%msk(0))
      call neko_scratch_registry%relinquish_field(scratch_idx)
    end if

    call this%finalize_base()
  end subroutine non_normal_aligned_finalize

  !> Free the boundary condition and its nested storage.
  subroutine non_normal_aligned_free(this)
    class(non_normal_aligned_t), target, intent(inout) :: this

    call this%bc_x%free()
    call this%bc_y%free()
    call this%bc_z%free()
    call this%value_x%free()
    call this%value_y%free()
    call this%value_z%free()

    if (allocated(this%field_file_name)) then
       deallocate(this%field_file_name)
    end if
    if (allocated(this%field_mesh_file_name)) then
       deallocate(this%field_mesh_file_name)
    end if

    this%field_interpolate = .false.
    this%field_interp_tolerance = GLOB_INTERP_TOL
    this%field_interp_padding = GLOB_INTERP_PAD
    this%read_values_from_field = .false.
    call this%free_base()
  end subroutine non_normal_aligned_free
end module non_normal_aligned
