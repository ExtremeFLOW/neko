! Copyright (c) 2020-2024, The Neko Authors
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
!> Defines a zero-valued Dirichlet boundary condition.
module zero_dirichlet
  use device_wall
  use num_types, only : rp
  use bc, only : bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use coefs, only : coef_t
  use json_module, only : json_file
  implicit none
  private

  !> Zero-valued Dirichlet boundary condition.
  !! Used for no-slip walls, but also for various auxillary conditions,
  !! such as for residuals.
  type, public, extends(bc_t) :: zero_dirichlet_t
   contains
     procedure, pass(this) :: apply_scalar => no_slip_wall_apply_scalar
     procedure, pass(this) :: apply_vector => no_slip_wall_apply_vector
     procedure, pass(this) :: apply_scalar_dev => no_slip_wall_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => no_slip_wall_apply_vector_dev
     !> Constructor.
     procedure, pass(this) :: init => no_slip_wall_init
     !> Destructor.
     procedure, pass(this) :: free => no_slip_wall_free
  end type zero_dirichlet_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine no_slip_wall_init(this, coef, json)
    class(zero_dirichlet_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_base(coef)
  end subroutine no_slip_wall_init

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x
  subroutine no_slip_wall_apply_scalar(this, x, n, t, tstep)
    class(zero_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
    end do

  end subroutine no_slip_wall_apply_scalar

  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z
  subroutine no_slip_wall_apply_vector(this, x, y, z, n, t, tstep)
    class(zero_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
       y(k) = 0d0
       z(k) = 0d0
    end do

  end subroutine no_slip_wall_apply_vector

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x (device version)
  subroutine no_slip_wall_apply_scalar_dev(this, x_d, t, tstep)
    class(zero_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_no_slip_wall_apply_scalar(this%msk_d, x_d, size(this%msk))

  end subroutine no_slip_wall_apply_scalar_dev

  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine no_slip_wall_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(zero_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_no_slip_wall_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                          size(this%msk))

  end subroutine no_slip_wall_apply_vector_dev

  !> Destructor
  subroutine no_slip_wall_free(this)
    class(zero_dirichlet_t), target, intent(inout) :: this

    call this%free_base()

  end subroutine no_slip_wall_free

end module zero_dirichlet
