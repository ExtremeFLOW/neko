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
!> Defines a Dirichlet condition prescribed by a mathematical expression
module expression_dirichlet
  use num_types, only : rp
  use bc, only : bc_t
  use coefs, only : coef_t
  use expression, only : expression_t, expression_check_finite
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_unmap, device_memcpy, HOST_TO_DEVICE
  use device_inhom_dirichlet, only : device_inhom_dirichlet_apply_scalar
  use json_module, only : json_file
  use json_utils, only : json_get
  use time_state, only : time_state_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Dirichlet condition \f$ x = g(x, y, z, t) \f$ on \f$\partial \Omega\f$,
  !! where \f$ g \f$ is a mathematical expression given in the case file.
  !! @details The expression is compiled once, at setup, and evaluated on the
  !! host in the points of the mask only. The resulting values are scattered
  !! into the solution with the same kernel that the Blasius condition uses,
  !! so no expression evaluation happens on the device.
  !!
  !! An expression that does not depend on time is evaluated once, in
  !! `finalize`, and from then on the condition costs exactly what a constant
  !! Dirichlet condition costs. A time dependent one is re-evaluated once per
  !! timestep, gated on the `updated` flag of the base type.
  type, public, extends(bc_t) :: expression_dirichlet_t
     !> The compiled expression prescribing the boundary value.
     type(expression_t) :: expr
     !> The value in each point of the mask.
     real(kind=rp), allocatable :: g(:)
     !> Coordinates of the points of the mask.
     real(kind=rp), allocatable :: xm(:)
     real(kind=rp), allocatable :: ym(:)
     real(kind=rp), allocatable :: zm(:)
     !> Device pointer for `g`.
     type(c_ptr) :: g_d = C_NULL_PTR
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => expression_dirichlet_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          expression_dirichlet_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => expression_dirichlet_free
     !> Finalize.
     procedure, pass(this) :: finalize => expression_dirichlet_finalize
     !> Apply the condition to a scalar field.
     procedure, pass(this) :: apply_scalar => expression_dirichlet_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => expression_dirichlet_apply_vector
     !> Apply the condition to a scalar field (device).
     procedure, pass(this) :: apply_scalar_dev => &
          expression_dirichlet_apply_scalar_dev
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          expression_dirichlet_apply_vector_dev
     !> Bring `g` up to date with the current time.
     procedure, pass(this) :: update => expression_dirichlet_update
  end type expression_dirichlet_t

  public :: expression_mask_coords

contains

  !> Constructor from JSON.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine expression_dirichlet_init(this, coef, json)
    class(expression_dirichlet_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: str

    call json_get(json, "value", str)
    call this%init_from_components(coef, str)
    if (allocated(str)) deallocate(str)

  end subroutine expression_dirichlet_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] str The expression prescribing the boundary value.
  subroutine expression_dirichlet_init_from_components(this, coef, str)
    class(expression_dirichlet_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: str

    call this%free()
    call this%init_base(coef)

    if (len_trim(str) .eq. 0) then
       call neko_error("An expression boundary condition needs a non-empty " &
            // "expression under the value keyword")
    end if

    call this%expr%init(str)

  end subroutine expression_dirichlet_init_from_components

  !> Destructor.
  subroutine expression_dirichlet_free(this)
    class(expression_dirichlet_t), target, intent(inout) :: this

    call this%free_base()
    call this%expr%free()

    if (allocated(this%g)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%g, this%g_d)
       end if
       deallocate(this%g)
    end if

    if (allocated(this%xm)) deallocate(this%xm)
    if (allocated(this%ym)) deallocate(this%ym)
    if (allocated(this%zm)) deallocate(this%zm)

  end subroutine expression_dirichlet_free

  !> Finalize.
  !! @details Tabulates the coordinates of the points of the mask, which is
  !! only known once the mask has been built, and evaluates the expression
  !! right away if it does not depend on time.
  subroutine expression_dirichlet_finalize(this)
    class(expression_dirichlet_t), target, intent(inout) :: this
    integer :: m

    call this%finalize_base()

    m = this%msk(0)
    if (m .eq. 0) return

    allocate(this%xm(m), this%ym(m), this%zm(m), this%g(m))
    call expression_mask_coords(this, this%xm, this%ym, this%zm)
    this%g = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%g, this%g_d, m)
    end if

    if (.not. this%expr%time_dependent) then
       call this%expr%eval(this%g, m, this%xm, this%ym, this%zm)
       call expression_check_finite(this%expr%src, this%g, m, &
            "boundary condition")
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%g, this%g_d, m, HOST_TO_DEVICE, &
               sync = .true.)
       end if
    end if

  end subroutine expression_dirichlet_finalize

  !> Bring `g` up to date with the current time.
  !! @details A no-op unless the expression depends on time and the values
  !! have not already been updated in this timestep.
  !! @param[in] time The current time state.
  !! @param[inout] strm The device stream to issue the copy on.
  subroutine expression_dirichlet_update(this, time, strm)
    class(expression_dirichlet_t), intent(inout) :: this
    type(time_state_t), intent(in), optional :: time
    type(c_ptr), intent(inout), optional :: strm
    integer :: m

    if (.not. this%expr%time_dependent) return
    if (this%updated) return

    if (.not. present(time)) then
       call neko_error("The boundary condition expression '" // &
            this%expr%src // "' depends on time, but the solver did not " // &
            "provide a time state")
    end if

    m = this%msk(0)
    call this%expr%eval(this%g, m, this%xm, this%ym, this%zm, time%t, time%dt)
    call expression_check_finite(this%expr%src, this%g, m, &
         "boundary condition")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! Synchronous on purpose: the host buffer is reused every timestep,
       ! and the copy is only as large as the mask.
       call device_memcpy(this%g, this%g_d, m, HOST_TO_DEVICE, &
            sync = .true., strm = strm)
    end if

    this%updated = .true.

  end subroutine expression_dirichlet_update

  !> Apply the condition to a scalar field.
  !! @param[inout] x The field onto which to apply the values.
  !! @param[in] n The size of `x`.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  subroutine expression_dirichlet_apply_scalar(this, x, n, time, strong)
    class(expression_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m
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
    ! `single` also makes `g` visible to every thread before the loop below.
    !$omp single
    call this%update(time)
    !$omp end single

    !$omp do
    do i = 1, m
       x(this%msk(i)) = this%g(i)
    end do
    !$omp end do

  end subroutine expression_dirichlet_apply_scalar

  !> (No-op) Apply vector.
  !! @param[inout] x The x-component of the field.
  !! @param[inout] y The y-component of the field.
  !! @param[inout] z The z-component of the field.
  !! @param[in] n The size of `x`, `y` and `z`.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  subroutine expression_dirichlet_apply_vector(this, x, y, z, n, time, strong)
    class(expression_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine expression_dirichlet_apply_vector

  !> Apply the condition to a scalar field (device version).
  !! @param[inout] x_d Device pointer to the field.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  !! @param[inout] strm The device stream to issue the work on.
  subroutine expression_dirichlet_apply_scalar_dev(this, x_d, time, strong, &
       strm)
    class(expression_dirichlet_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
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

    call device_inhom_dirichlet_apply_scalar(this%msk_d, x_d, this%g_d, m, strm)

  end subroutine expression_dirichlet_apply_scalar_dev

  !> (No-op) Apply vector (device version).
  !! @param[inout] x_d Device pointer to the x-component of the field.
  !! @param[inout] y_d Device pointer to the y-component of the field.
  !! @param[inout] z_d Device pointer to the z-component of the field.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether the condition is applied strongly.
  !! @param[inout] strm The device stream to issue the work on.
  subroutine expression_dirichlet_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(expression_dirichlet_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine expression_dirichlet_apply_vector_dev

  !> Tabulate the coordinates of the points of the mask of a boundary
  !! condition.
  !! @param[in] bc The boundary condition, which must have been finalized.
  !! @param[inout] xm The x-coordinates of the points of the mask.
  !! @param[inout] ym The y-coordinates of the points of the mask.
  !! @param[inout] zm The z-coordinates of the points of the mask.
  subroutine expression_mask_coords(bc, xm, ym, zm)
    class(bc_t), intent(in) :: bc
    real(kind=rp), intent(inout) :: xm(:)
    real(kind=rp), intent(inout) :: ym(:)
    real(kind=rp), intent(inout) :: zm(:)

    call gather_coords(bc%msk, bc%msk(0), bc%dof%x, bc%dof%y, bc%dof%z, &
         size(bc%dof%x), xm, ym, zm)

  end subroutine expression_mask_coords

  !> Gather the coordinates of the masked points into contiguous arrays.
  !! @details Split out from `expression_mask_coords` so that the rank 4
  !! coordinate arrays of the dofmap are linearised by the call.
  !! @param[in] msk The mask, where `msk(0)` is its length.
  !! @param[in] m The number of masked points.
  !! @param[in] x The x-coordinates of every point of the dofmap.
  !! @param[in] y The y-coordinates of every point of the dofmap.
  !! @param[in] z The z-coordinates of every point of the dofmap.
  !! @param[in] n The number of points of the dofmap.
  !! @param[inout] xm The x-coordinates of the masked points.
  !! @param[inout] ym The y-coordinates of the masked points.
  !! @param[inout] zm The z-coordinates of the masked points.
  subroutine gather_coords(msk, m, x, y, z, n, xm, ym, zm)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: msk(0:m)
    real(kind=rp), intent(in) :: x(n)
    real(kind=rp), intent(in) :: y(n)
    real(kind=rp), intent(in) :: z(n)
    real(kind=rp), intent(inout) :: xm(m)
    real(kind=rp), intent(inout) :: ym(m)
    real(kind=rp), intent(inout) :: zm(m)
    integer :: i, k

    do i = 1, m
       k = msk(i)
       xm(i) = x(k)
       ym(i) = y(k)
       zm(i) = z(k)
    end do

  end subroutine gather_coords

end module expression_dirichlet
