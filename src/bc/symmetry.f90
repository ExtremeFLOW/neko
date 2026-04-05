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
!> Generic mixed Dirichlet-Neumann symmetry plane condition.
module symmetry
  use num_types, only : rp
  use bc, only : mixed_bc_t, BC_TYPES
  use device_constrain_mixed_bc, only : device_constrain_mixed_bc_zero
  use coefs, only : coef_t
  use json_module, only : json_file
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  type, public, extends(mixed_bc_t) :: symmetry_t
   contains
     procedure, pass(this) :: apply_scalar => symmetry_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_apply_vector
     procedure, pass(this) :: apply_scalar_dev => symmetry_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => symmetry_apply_vector_dev
     procedure, pass(this) :: init => symmetry_init
     procedure, pass(this) :: init_from_components => &
          symmetry_init_from_components
     procedure, pass(this) :: free => symmetry_free
     procedure, pass(this) :: finalize => symmetry_finalize
  end type symmetry_t

contains

  subroutine symmetry_init(this, coef, json)
    class(symmetry_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%init_from_components(coef)
  end subroutine symmetry_init

  subroutine symmetry_init_from_components(this, coef)
    class(symmetry_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef

    call this%free()

    call this%init_base(coef)
    this%constraints = (/ .true., .false., .false. /)
    this%bc_type = BC_TYPES%MIXED_CONSTRAINS_NORMAL
  end subroutine symmetry_init_from_components

  subroutine symmetry_finalize(this)
    class(symmetry_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine symmetry_finalize

  subroutine symmetry_apply_scalar(this, x, n, time, strong)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine symmetry_apply_scalar

  subroutine symmetry_apply_vector(this, x, y, z, n, time, strong)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i, m, k
    real(kind=rp) :: normal(3), u_n

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       m = this%resolved_msk%size()

       do i = 1, m
          k = this%resolved_msk%get(i)
          normal = this%n%x(:,i)
          u_n = x(k) * normal(1) + y(k) * normal(2) + z(k) * normal(3)

          x(k) = x(k) - u_n * normal(1)
          y(k) = y(k) - u_n * normal(2)
          z(k) = z(k) - u_n * normal(3)
       end do
    end if
  end subroutine symmetry_apply_vector

  subroutine symmetry_apply_scalar_dev(this, x_d, time, strong, strm)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine symmetry_apply_scalar_dev

  subroutine symmetry_apply_vector_dev(this, x_d, y_d, z_d, time, strong, strm)
    class(symmetry_t), intent(inout), target :: this
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
          call device_constrain_mixed_bc_zero(this%resolved_msk%get_d(), &
               x_d, y_d, z_d, 1, 0, 0, this%n%x_d, this%t1%x_d, &
               this%t2%x_d, m, strm)
       end if
    end if
  end subroutine symmetry_apply_vector_dev

  subroutine symmetry_free(this)
    class(symmetry_t), target, intent(inout) :: this

    call this%free_mixed()
    call this%free_base()
  end subroutine symmetry_free

end module symmetry
