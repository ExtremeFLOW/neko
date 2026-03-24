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
!> Defines no-slip boundary condition (extends zero_dirichlet)
module no_slip
  use num_types, only : rp
  use coefs, only : coef_t
  use field, only : field_t
  use zero_dirichlet, only : zero_dirichlet_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use time_state, only : time_state_t
  use math, only : masked_copy_0
  use device_math, only : device_masked_copy_0
  use utils, only : neko_error
  use registry, only : neko_registry
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  type, public, extends(zero_dirichlet_t) :: no_slip_t
     ! mesh velocity fields.
     type(field_t), pointer :: wx => null()
     type(field_t), pointer :: wy => null()
     type(field_t), pointer :: wz => null()
     logical :: is_moving = .false.
   contains
     procedure, pass(this) :: init => no_slip_init
     procedure, pass(this) :: apply_vector => no_slip_apply_vector
     procedure, pass(this) :: apply_vector_dev => no_slip_apply_vector_dev
     procedure, pass(this) :: free => no_slip_free
  end type no_slip_t

contains

  subroutine no_slip_init(this, coef, json)
    class(no_slip_t), intent(inout), target :: this
    type(coef_t), intent(in), target :: coef
    type(json_file), intent(inout) :: json

    ! Normal init (zero_dirichlet)
    call this%zero_dirichlet_t%init(coef, json)
    call json_get_or_default(json, "moving", this%is_moving, .false.)

    if (neko_registry%field_exists('wm_x')) then
       this%wx => neko_registry%get_field('wm_x')
       this%wy => neko_registry%get_field('wm_y')
       this%wz => neko_registry%get_field('wm_z')
    end if

    ! If moving is requested, mesh velocitiy fileds should be linked already
    ! in ale_manager.
    if (this%is_moving) then
       if (.not. associated(this%wx) .or. .not. associated(this%wy) .or. &
            .not. associated(this%wz)) then
          call neko_error("Velocity BC 'no_slip' with moving: true is &
          &found, but ALE is not activated in case file.")
       end if
    end if
  end subroutine no_slip_init


  subroutine no_slip_apply_vector(this, x, y, z, n, time, strong)
    class(no_slip_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x, y, z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    strong_ = .true.
    if (present(strong)) strong_ = strong
    if (.not. strong_) return

    if (this%is_moving) then
       ! moving wall: u_wall = w_mesh
       call masked_copy_0(x, this%wx%x, this%msk, n, this%msk(0))
       call masked_copy_0(y, this%wy%x, this%msk, n, this%msk(0))
       call masked_copy_0(z, this%wz%x, this%msk, n, this%msk(0))
    else
       call this%zero_dirichlet_t%apply_vector(x, y, z, n, time, strong_)
    end if
  end subroutine no_slip_apply_vector


  subroutine no_slip_apply_vector_dev(this, x_d, y_d, z_d, time, strong, strm)
    class(no_slip_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d, y_d, z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    strong_ = .true.
    if (present(strong)) strong_ = strong
    if (.not. strong_) return

    if (this%is_moving) then
       call device_masked_copy_0(x_d, this%wx%x_d, this%msk_d, &
            this%wx%dof%size(), this%msk(0), strm)
       call device_masked_copy_0(y_d, this%wy%x_d, this%msk_d, &
            this%wy%dof%size(), this%msk(0), strm)
       call device_masked_copy_0(z_d, this%wz%x_d, this%msk_d, &
            this%wz%dof%size(), this%msk(0), strm)
    else
       call this%zero_dirichlet_t%apply_vector_dev(x_d, y_d, z_d, time, &
            strong_, strm)
    end if
  end subroutine no_slip_apply_vector_dev


  subroutine no_slip_free(this)
    class(no_slip_t), intent(inout), target :: this

    call this%zero_dirichlet_t%free()
    nullify(this%wx)
    nullify(this%wy)
    nullify(this%wz)
  end subroutine no_slip_free

end module no_slip
