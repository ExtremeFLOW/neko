! Copyright (c) 2020-2021, The Neko Authors
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
!> Mixed Dirichlet condition constraining the tangential vector components.
module non_normal
  use json_module, only : json_file
  use bc, only : BC_TYPES
  use mixed_bc, only : mixed_bc_t
  use num_types, only : rp
  use device_constrain_mixed_bc, only : device_constrain_mixed_bc_set_const
  use coefs, only : coef_t
  use json_utils, only : json_get_or_lookup
  use utils, only : neko_error
  use time_state, only : time_state_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

!> Mixed Dirichlet condition constraining the tangential vector components.
  type, public, extends(mixed_bc_t) :: non_normal_t
     real(kind=rp) :: value(3) = 0.0_rp
   contains
     procedure, pass(this) :: apply_scalar => non_normal_apply_scalar
     procedure, pass(this) :: apply_vector => non_normal_apply_vector
     procedure, pass(this) :: apply_scalar_dev => non_normal_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => non_normal_apply_vector_dev
     procedure, pass(this) :: init => non_normal_init
     procedure, pass(this) :: init_from_components => &
          non_normal_init_from_components
     procedure, pass(this) :: free => non_normal_free
     procedure, pass(this) :: finalize => non_normal_finalize
  end type non_normal_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine non_normal_init(this, coef, json)
    class(non_normal_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json
    real(kind=rp), allocatable :: value(:)
    real(kind=rp) :: value_3(3)
    logical :: found
    integer :: var_type

    value_3 = 0.0_rp
    call json_get_or_lookup(json, "value", value)
    if (size(value) .ne. 3) then
       call neko_error("The non_normal boundary condition requires a " // &
            "3-component value vector.")
    end if
    value_3 = value
    call this%init_from_components(coef, value_3)
  end subroutine non_normal_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] value The tangential value in global coordinates.
  subroutine non_normal_init_from_components(this, coef, value)
    class(non_normal_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    real(kind=rp), intent(in) :: value(3)

    call this%free()
    call this%init_base(coef)
    this%constraints = [.false., .true., .true.]
    this%bc_type = BC_TYPES%MIXED_CONSTRAINS_TANGENT
    this%value = value
  end subroutine non_normal_init_from_components

  subroutine non_normal_apply_scalar(this, x, n, time, strong)
    class(non_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine non_normal_apply_scalar

  !> Strong application preserving the normal component while enforcing the
  !! tangential projections of the configured global value.
  subroutine non_normal_apply_vector(this, x, y, z, n, time, strong)
    class(non_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i, m, k
    real(kind=rp) :: normal(3), t1(3), t2(3)
    real(kind=rp) :: u_n, g_t1, g_t2

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. strong_) return

    m = this%resolved_msk%size()
    do i = 1, m
       k = this%resolved_msk%get(i)
       normal = this%n%x(:, i)
       t1 = this%t1%x(:, i)
       t2 = this%t2%x(:, i)

       ! Normal component, will be reserved
       u_n = x(k) * normal(1) + y(k) * normal(2) + z(k) * normal(3)
       ! Project global input onto local tangential directions
       g_t1 = this%value(1) * t1(1) + this%value(2) * t1(2) + &
            this%value(3) * t1(3)
       g_t2 = this%value(1) * t2(1) + this%value(2) * t2(2) + &
            this%value(3) * t2(3)

       ! Reconstruct in global coordinates
       x(k) = u_n * normal(1) + g_t1 * t1(1) + g_t2 * t2(1)
       y(k) = u_n * normal(2) + g_t1 * t1(2) + g_t2 * t2(2)
       z(k) = u_n * normal(3) + g_t1 * t1(3) + g_t2 * t2(3)
    end do
  end subroutine non_normal_apply_vector

  subroutine non_normal_apply_scalar_dev(this, x_d, time, strong, strm)
    class(non_normal_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine non_normal_apply_scalar_dev

  subroutine non_normal_apply_vector_dev(this, x_d, y_d, z_d, time, strong, &
       strm)
    class(non_normal_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_
    integer :: m

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       m = this%resolved_msk%size()
       if (m .gt. 0) then
          call device_constrain_mixed_bc_set_const( &
               this%resolved_msk%get_d(), x_d, y_d, z_d, 0, 1, 1, &
               this%n%x_d, this%t1%x_d, this%t2%x_d, this%value(1), &
               this%value(2), this%value(3), m, strm)
       end if
    end if
  end subroutine non_normal_apply_vector_dev

  !> Finalize generic mixed non-normal bc.
  subroutine non_normal_finalize(this)
    class(non_normal_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine non_normal_finalize

  !> Destructor for generic mixed non-normal bc.
  subroutine non_normal_free(this)
    class(non_normal_t), target, intent(inout) :: this

    call this%free_mixed()
    call this%free_base()
  end subroutine non_normal_free
end module non_normal
