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
!> Defines inflow dirichlet conditions
module inflow
  use device_inflow
  use num_types, only : rp
  use bc, only : bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_loc
  use coefs, only : coef_t
  use json_module, only : json_file
  use json_utils, only : json_get
  implicit none
  private

  !> Dirichlet condition for inlet (vector valued)
  type, public, extends(bc_t) :: inflow_t
     real(kind=rp), dimension(3) :: x = (/0d0, 0d0, 0d0 /)
   contains
     procedure, pass(this) :: apply_scalar => inflow_apply_scalar
     procedure, pass(this) :: apply_vector => inflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => inflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => inflow_apply_vector_dev
     !> Constructor
     procedure, pass(this) :: init => inflow_init
     !> Destructor.
     procedure, pass(this) :: free => inflow_free
     !> Finalize.
     procedure, pass(this) :: finalize => inflow_finalize
  end type inflow_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine inflow_init(this, coef, json)
    class(inflow_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) ::json
    real(kind=rp), allocatable :: x(:)

    call this%init_base(coef)
    !call json_get(json, 'case.fluid.inflow_condition.value', x)
    call json_get(json, 'value', x)
    this%x = x
  end subroutine inflow_init

  !> No-op scalar apply
  subroutine inflow_apply_scalar(this, x, n, t, tstep)
    class(inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine inflow_apply_scalar

  !> No-op scalar apply (device version)
  subroutine inflow_apply_scalar_dev(this, x_d, t, tstep)
    class(inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine inflow_apply_scalar_dev

  !> Apply inflow conditions (vector valued)
  subroutine inflow_apply_vector(this, x, y, z, n, t, tstep)
    class(inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k


    m = this%msk(0)

    write(*,*) "Applying inflow", m

    do i = 1, m
       k = this%msk(i)
       x(k) = this%x(1)
       y(k) = this%x(2)
       z(k) = this%x(3)
    end do
  end subroutine inflow_apply_vector

  !> Apply inflow conditions (vector valued) (device version)
  subroutine inflow_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_inflow_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                    c_loc(this%x), this%msk(0))

  end subroutine inflow_apply_vector_dev

  !> Destructor
  subroutine inflow_free(this)
    class(inflow_t), target, intent(inout) :: this

    call this%free_base()
  end subroutine inflow_free

  !> Finalize
  subroutine inflow_finalize(this)
    class(inflow_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine inflow_finalize


end module inflow
