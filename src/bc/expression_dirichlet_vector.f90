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
!> Defines a vector valued Dirichlet condition prescribed by mathematical
!! expressions
module expression_dirichlet_vector
  use num_types, only : rp
  use bc, only : bc_t
  use coefs, only : coef_t
  use expression, only : expression_t, expression_check_finite, NEKO_EXPR_LEN
  use expression_dirichlet, only : expression_mask_coords
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_unmap, device_memcpy, HOST_TO_DEVICE
  use device_inhom_dirichlet, only : device_inhom_dirichlet_apply_vector
  use json_module, only : json_file
  use json_utils, only : json_get
  use time_state, only : time_state_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Vector valued Dirichlet condition, with one mathematical expression per
  !! component, given in the case file.
  !! @details The vector valued counterpart of `expression_dirichlet_t`, see
  !! that type for how the expressions are evaluated and cached.
  type, public, extends(bc_t) :: expression_dirichlet_vector_t
     !> The compiled expressions, one per component.
     type(expression_t) :: expr(3)
     !> The values in each point of the mask, one array per component.
     real(kind=rp), allocatable :: gx(:)
     real(kind=rp), allocatable :: gy(:)
     real(kind=rp), allocatable :: gz(:)
     !> Coordinates of the points of the mask.
     real(kind=rp), allocatable :: xm(:)
     real(kind=rp), allocatable :: ym(:)
     real(kind=rp), allocatable :: zm(:)
     !> Device pointers for `gx`, `gy` and `gz`.
     type(c_ptr) :: gx_d = C_NULL_PTR
     type(c_ptr) :: gy_d = C_NULL_PTR
     type(c_ptr) :: gz_d = C_NULL_PTR
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => expression_dirichlet_vector_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          expression_dirichlet_vector_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => expression_dirichlet_vector_free
     !> Finalize.
     procedure, pass(this) :: finalize => expression_dirichlet_vector_finalize
     !> (No-op) Apply scalar.
     procedure, pass(this) :: apply_scalar => &
          expression_dirichlet_vector_apply_scalar
     !> Apply the condition to a vector field.
     procedure, pass(this) :: apply_vector => &
          expression_dirichlet_vector_apply_vector
     !> (No-op) Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          expression_dirichlet_vector_apply_scalar_dev
     !> Apply the condition to a vector field (device).
     procedure, pass(this) :: apply_vector_dev => &
          expression_dirichlet_vector_apply_vector_dev
     !> Bring the values up to date with the current time.
     procedure, pass(this) :: update => expression_dirichlet_vector_update
  end type expression_dirichlet_vector_t

contains

  !> Constructor from JSON.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine expression_dirichlet_vector_init(this, coef, json)
    class(expression_dirichlet_vector_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=NEKO_EXPR_LEN), allocatable :: str(:)

    call json_get(json, "value", str, filler = '')

    if (size(str) .ne. 3) then
       call neko_error("An expression velocity boundary condition takes " // &
            "exactly three expressions, one per component")
    end if

    call this%init_from_components(coef, str(1), str(2), str(3))
    if (allocated(str)) deallocate(str)

  end subroutine expression_dirichlet_vector_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] str_x The expression for the x-component.
  !! @param[in] str_y The expression for the y-component.
  !! @param[in] str_z The expression for the z-component.
  subroutine expression_dirichlet_vector_init_from_components(this, coef, &
       str_x, str_y, str_z)
    class(expression_dirichlet_vector_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: str_x
    character(len=*), intent(in) :: str_y
    character(len=*), intent(in) :: str_z

    call this%free()
    call this%init_base(coef)

    if (len_trim(str_x) .eq. 0 .or. len_trim(str_y) .eq. 0 .or. &
         len_trim(str_z) .eq. 0) then
       call neko_error("An expression boundary condition needs a non-empty " &
            // "expression for every component")
    end if

    call this%expr(1)%init(str_x)
    call this%expr(2)%init(str_y)
    call this%expr(3)%init(str_z)

  end subroutine expression_dirichlet_vector_init_from_components

  !> Destructor.
  subroutine expression_dirichlet_vector_free(this)
    class(expression_dirichlet_vector_t), target, intent(inout) :: this
    integer :: i

    call this%free_base()

    do i = 1, 3
       call this%expr(i)%free()
    end do

    if (allocated(this%gx)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%gx, this%gx_d)
       end if
       deallocate(this%gx)
    end if

    if (allocated(this%gy)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%gy, this%gy_d)
       end if
       deallocate(this%gy)
    end if

    if (allocated(this%gz)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%gz, this%gz_d)
       end if
       deallocate(this%gz)
    end if

    if (allocated(this%xm)) deallocate(this%xm)
    if (allocated(this%ym)) deallocate(this%ym)
    if (allocated(this%zm)) deallocate(this%zm)

  end subroutine expression_dirichlet_vector_free

  !> Finalize.
  !! @details Tabulates the coordinates of the points of the mask and
  !! evaluates whichever of the three expressions do not depend on time.
  !! @param[in] only_facets Whether to only mark the facets of the mask.
  subroutine expression_dirichlet_vector_finalize(this, only_facets)
    class(expression_dirichlet_vector_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_
    integer :: m

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = .false.
    end if

    call this%finalize_base(only_facets_)

    m = this%msk(0)
    if (m .eq. 0) return

    allocate(this%xm(m), this%ym(m), this%zm(m))
    allocate(this%gx(m), this%gy(m), this%gz(m))
    call expression_mask_coords(this, this%xm, this%ym, this%zm)
    this%gx = 0.0_rp
    this%gy = 0.0_rp
    this%gz = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%gx, this%gx_d, m)
       call device_map(this%gy, this%gy_d, m)
       call device_map(this%gz, this%gz_d, m)
    end if

    if (.not. this%expr(1)%time_dependent) then
       call this%expr(1)%eval(this%gx, m, this%xm, this%ym, this%zm)
       call expression_check_finite(this%expr(1)%src, this%gx, m, &
            "boundary condition")
    end if
    if (.not. this%expr(2)%time_dependent) then
       call this%expr(2)%eval(this%gy, m, this%xm, this%ym, this%zm)
       call expression_check_finite(this%expr(2)%src, this%gy, m, &
            "boundary condition")
    end if
    if (.not. this%expr(3)%time_dependent) then
       call this%expr(3)%eval(this%gz, m, this%xm, this%ym, this%zm)
       call expression_check_finite(this%expr(3)%src, this%gz, m, &
            "boundary condition")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%gx, this%gx_d, m, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%gy, this%gy_d, m, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%gz, this%gz_d, m, HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine expression_dirichlet_vector_finalize

  !> Bring the values up to date with the current time.
  !! @details A no-op unless at least one of the expressions depends on time
  !! and the values have not already been updated in this timestep. Only the
  !! time dependent components are re-evaluated, but all three are copied to
  !! the device, since they share one transfer anyway.
  !! @param[in] time The current time state.
  !! @param[inout] strm The device stream to issue the copies on.
  subroutine expression_dirichlet_vector_update(this, time, strm)
    class(expression_dirichlet_vector_t), intent(inout) :: this
    type(time_state_t), intent(in), optional :: time
    type(c_ptr), intent(inout), optional :: strm
    integer :: m

    if (.not. (this%expr(1)%time_dependent .or. &
         this%expr(2)%time_dependent .or. &
         this%expr(3)%time_dependent)) return
    if (this%updated) return

    if (.not. present(time)) then
       call neko_error("A boundary condition expression depends on time, " // &
            "but the solver did not provide a time state")
    end if

    m = this%msk(0)

    if (this%expr(1)%time_dependent) then
       call this%expr(1)%eval(this%gx, m, this%xm, this%ym, this%zm, &
            time%t, time%dt)
       call expression_check_finite(this%expr(1)%src, this%gx, m, &
            "boundary condition")
    end if
    if (this%expr(2)%time_dependent) then
       call this%expr(2)%eval(this%gy, m, this%xm, this%ym, this%zm, &
            time%t, time%dt)
       call expression_check_finite(this%expr(2)%src, this%gy, m, &
            "boundary condition")
    end if
    if (this%expr(3)%time_dependent) then
       call this%expr(3)%eval(this%gz, m, this%xm, this%ym, this%zm, &
            time%t, time%dt)
       call expression_check_finite(this%expr(3)%src, this%gz, m, &
            "boundary condition")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! The three copies go on the same stream, so waiting for the last one
       ! covers all of them. Synchronous on purpose, see
       ! expression_dirichlet_update.
       call device_memcpy(this%gx, this%gx_d, m, HOST_TO_DEVICE, &
            sync = .false., strm = strm)
       call device_memcpy(this%gy, this%gy_d, m, HOST_TO_DEVICE, &
            sync = .false., strm = strm)
       call device_memcpy(this%gz, this%gz_d, m, HOST_TO_DEVICE, &
            sync = .true., strm = strm)
    end if

    this%updated = .true.

  end subroutine expression_dirichlet_vector_update

  !> (No-op) Apply scalar.
  !! @param[inout] x The field onto which to apply the values.
  !! @param[in] n The size of `x`.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  subroutine expression_dirichlet_vector_apply_scalar(this, x, n, time, strong)
    class(expression_dirichlet_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine expression_dirichlet_vector_apply_scalar

  !> (No-op) Apply scalar (device version).
  !! @param[inout] x_d Device pointer to the field.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  !! @param[inout] strm The device stream to issue the work on.
  subroutine expression_dirichlet_vector_apply_scalar_dev(this, x_d, time, &
       strong, strm)
    class(expression_dirichlet_vector_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine expression_dirichlet_vector_apply_scalar_dev

  !> Apply the condition to a vector field.
  !! @param[inout] x The x-component of the field.
  !! @param[inout] y The y-component of the field.
  !! @param[inout] z The z-component of the field.
  !! @param[in] n The size of `x`, `y` and `z`.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  subroutine expression_dirichlet_vector_apply_vector(this, x, y, z, n, &
       time, strong)
    class(expression_dirichlet_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    m = this%msk(0)
    if (.not. strong_ .or. m .eq. 0) return

    ! The bc list applies the host conditions from inside an OpenMP parallel
    ! region, and evaluating an expression mutates the shared evaluation stack
    ! of `expr`, so only one thread may run the update. The implicit barrier of
    ! `single` also makes the values visible to every thread before the loop
    ! below.
    !$omp single
    call this%update(time)
    !$omp end single

    !$omp do
    do i = 1, m
       k = this%msk(i)
       x(k) = this%gx(i)
       y(k) = this%gy(i)
       z(k) = this%gz(i)
    end do
    !$omp end do

  end subroutine expression_dirichlet_vector_apply_vector

  !> Apply the condition to a vector field (device version).
  !! @param[inout] x_d Device pointer to the x-component of the field.
  !! @param[inout] y_d Device pointer to the y-component of the field.
  !! @param[inout] z_d Device pointer to the z-component of the field.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  !! @param[inout] strm The device stream to issue the work on.
  subroutine expression_dirichlet_vector_apply_vector_dev(this, x_d, y_d, &
       z_d, time, strong, strm)
    class(expression_dirichlet_vector_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    integer :: m
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    m = this%msk(0)
    if (.not. strong_ .or. m .eq. 0) return

    call this%update(time, strm)

    call device_inhom_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
         this%gx_d, this%gy_d, this%gz_d, m, strm)

  end subroutine expression_dirichlet_vector_apply_vector_dev

end module expression_dirichlet_vector
