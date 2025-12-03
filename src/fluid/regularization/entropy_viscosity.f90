! Copyright (c) 2025, The Neko Authors
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
module entropy_viscosity
  use num_types, only : rp
  use regularization, only: regularization_t
  use json_module, only : json_file
  use json_utils, only: json_get_or_default
  use field, only: field_t
  use field_math, only: field_cfill, field_glsum, field_cadd, field_copy
  use field_series, only: field_series_t
  use coefs, only: coef_t
  use dofmap, only: dofmap_t
  use time_state, only: time_state_t
  use bdf_time_scheme, only: bdf_time_scheme_t
  use operators, only: div
  use math, only: glmax, absval
  use scratch_registry, only: neko_scratch_registry
  use mesh, only: mesh_t
  use space, only: space_t
  use gather_scatter, only: gs_t
  use gs_ops, only: GS_OP_ADD, GS_OP_MAX
  implicit none
  private

  type, public, extends(regularization_t) :: entropy_viscosity_t
     real(kind=rp) :: c_entropy
     real(kind=rp) :: c_max
     type(field_t) :: entropy_residual
     type(field_series_t) :: S_lag
     type(field_t), pointer :: S => null()
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(field_t), pointer :: h => null()
     type(field_t), pointer :: max_wave_speed => null()
     type(mesh_t), pointer :: msh => null()
     type(space_t), pointer :: Xh => null()
     type(gs_t), pointer :: gs => null()
   contains
     procedure, pass(this) :: init => entropy_viscosity_init
     procedure, pass(this) :: free => entropy_viscosity_free
     procedure, pass(this) :: compute => entropy_viscosity_compute
     procedure, pass(this) :: update_lag => entropy_viscosity_update_lag
    procedure, pass(this), private :: compute_residual => entropy_viscosity_compute_residual
    procedure, pass(this), private :: compute_viscosity => entropy_viscosity_compute_viscosity
    procedure, pass(this), private :: smooth_viscosity => entropy_viscosity_smooth_viscosity
    procedure, pass(this), private :: apply_element_max => entropy_viscosity_apply_element_max
    procedure, pass(this), private :: low_order_viscosity => entropy_viscosity_low_order
  end type entropy_viscosity_t

  public :: entropy_viscosity_set_fields

contains

  subroutine entropy_viscosity_init(this, json, coef, dof, reg_coeff)
    class(entropy_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof
    type(field_t), intent(in), target :: reg_coeff

    call this%init_base(json, coef, dof, reg_coeff)

    call json_get_or_default(json, 'c_max', this%c_max, 1.0_rp)
    call json_get_or_default(json, 'c_entropy', this%c_entropy, 1.0_rp)

    call this%entropy_residual%init(dof, 'entropy_residual')

    nullify(this%S)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%h)
    nullify(this%max_wave_speed)
    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%gs)

  end subroutine entropy_viscosity_init

  subroutine entropy_viscosity_free(this)
    class(entropy_viscosity_t), intent(inout) :: this

    call this%free_base()
    call this%entropy_residual%free()
    call this%S_lag%free()

    nullify(this%S)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%h)
    nullify(this%max_wave_speed)
    nullify(this%msh)
    nullify(this%Xh)
    nullify(this%gs)

  end subroutine entropy_viscosity_free

  subroutine entropy_viscosity_compute(this, time, tstep, dt)
    class(entropy_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: tstep
    real(kind=rp), intent(in) :: dt
    integer :: i, n

    n = this%dof%size()

    if (this%c_entropy >= 1.0e10_rp) then
       do i = 1, n
          this%reg_coeff%x(i,1,1,1) = this%low_order_viscosity(i)
       end do
       return
    end if

    call this%compute_residual(tstep, dt, time%dtlag)
    call this%compute_viscosity(tstep)

  end subroutine entropy_viscosity_compute

  subroutine entropy_viscosity_compute_residual(this, tstep, dt, dt_lag)
    class(entropy_viscosity_t), intent(inout) :: this
    integer, intent(in) :: tstep
    real(kind=rp), intent(in) :: dt
    real(kind=rp), intent(in) :: dt_lag(10)
    integer :: i, n
    type(field_t), pointer :: us_field, vs_field, ws_field, div_field
    integer :: temp_indices(4)
    real(kind=rp) :: bdf_coeffs(4)
    type(bdf_time_scheme_t) :: bdf_scheme
    real(kind=rp) :: dt_local(10)

    if (tstep .le. 3) then
       return
    end if

    n = this%dof%size()
    call field_cfill(this%entropy_residual, 0.0_rp, n)

    bdf_coeffs = 0.0_rp
    dt_local = dt_lag
    
    call bdf_scheme%compute_coeffs(bdf_coeffs, dt_local, 3)

    do i = 1, n
       this%entropy_residual%x(i,1,1,1) = (bdf_coeffs(1) * this%S%x(i,1,1,1) &
                                           - bdf_coeffs(2) * this%S_lag%lf(1)%x(i,1,1,1) &
                                           - bdf_coeffs(3) * this%S_lag%lf(2)%x(i,1,1,1) &
                                           - bdf_coeffs(4) * this%S_lag%lf(3)%x(i,1,1,1)) / dt
    end do

    call neko_scratch_registry%request_field(us_field, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(vs_field, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ws_field, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(div_field, temp_indices(4), .false.)

    do i = 1, n
       us_field%x(i,1,1,1) = this%u%x(i,1,1,1) * this%S%x(i,1,1,1)
       vs_field%x(i,1,1,1) = this%v%x(i,1,1,1) * this%S%x(i,1,1,1)  
       ws_field%x(i,1,1,1) = this%w%x(i,1,1,1) * this%S%x(i,1,1,1)
    end do

    call div(div_field%x, us_field%x, vs_field%x, ws_field%x, this%coef)

    do i = 1, n
       this%entropy_residual%x(i,1,1,1) = abs(this%entropy_residual%x(i,1,1,1) + div_field%x(i,1,1,1))
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine entropy_viscosity_compute_residual

  subroutine entropy_viscosity_compute_viscosity(this, tstep)
    class(entropy_viscosity_t), intent(inout) :: this
    integer, intent(in) :: tstep
    integer :: i, n, temp_indices(1)
    real(kind=rp) :: S_mean, n_S
    type(field_t), pointer :: temp_field

    n = this%dof%size()

    if (tstep .le. 3) then
       call field_cfill(this%reg_coeff, 0.0_rp, n)
       return
    end if
    
    call neko_scratch_registry%request_field(temp_field, temp_indices(1), .false.)
    
    call field_cfill(temp_field, 1.0_rp, n)
    S_mean = field_glsum(this%S, n) / field_glsum(temp_field, n)
    
    call field_copy(temp_field, this%S, n)
    call field_cadd(temp_field, -S_mean, n)
    
    call absval(temp_field%x, n)
    
    n_S = glmax(temp_field%x, n)
    
    call neko_scratch_registry%relinquish_field(temp_indices)
    
    if (n_S < 1.0e-12_rp) then
       n_S = 1.0e-12_rp
    end if
    
    do i = 1, n
       this%reg_coeff%x(i,1,1,1) = this%c_entropy * this%h%x(i,1,1,1)**2 * &
                                   this%entropy_residual%x(i,1,1,1) / n_S
    end do

    ! Apply element max to make uniform within elements
    call this%apply_element_max()

    ! Clamp to low-order viscosity
    do i = 1, n
       this%reg_coeff%x(i,1,1,1) = min(this%reg_coeff%x(i,1,1,1), &
                                       this%low_order_viscosity(i))
    end do

    ! Apply Gaussian smoothing for smooth final viscosity
    call this%smooth_viscosity()

  end subroutine entropy_viscosity_compute_viscosity

  !> Cross-element smoothing via gather-scatter averaging.
  !! Averages viscosity values at shared nodes between elements.
  subroutine entropy_viscosity_smooth_viscosity(this)
    class(entropy_viscosity_t), intent(inout) :: this
    integer :: i, n
    type(field_t), pointer :: temp_field, mult_field
    integer :: temp_indices(2)

    n = this%dof%size()

    call neko_scratch_registry%request_field(temp_field, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(mult_field, temp_indices(2), .false.)

    ! Sum values at shared nodes
    call field_copy(temp_field, this%reg_coeff, n)
    call this%gs%op(temp_field, GS_OP_ADD)

    ! Compute multiplicity (number of elements sharing each node)
    call field_cfill(mult_field, 1.0_rp, n)
    call this%gs%op(mult_field, GS_OP_ADD)

    ! Average: divide sum by multiplicity
    do i = 1, n
       temp_field%x(i,1,1,1) = temp_field%x(i,1,1,1) / mult_field%x(i,1,1,1)
    end do

    call field_copy(this%reg_coeff, temp_field, n)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine entropy_viscosity_smooth_viscosity

  subroutine entropy_viscosity_apply_element_max(this)
    class(entropy_viscosity_t), intent(inout) :: this
    integer :: i, j, k, el, lx
    real(kind=rp) :: max_visc_el

    lx = this%Xh%lx
    
    do el = 1, this%msh%nelv
       max_visc_el = 0.0_rp
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                max_visc_el = max(max_visc_el, this%reg_coeff%x(i,j,k,el))
             end do
          end do
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                this%reg_coeff%x(i,j,k,el) = max_visc_el
             end do
          end do
       end do
    end do

  end subroutine entropy_viscosity_apply_element_max

  subroutine entropy_viscosity_set_fields(this, S, u, v, w, h, max_wave_speed, &
       msh, Xh, gs)
    class(entropy_viscosity_t), intent(inout) :: this
    type(field_t), target, intent(inout) :: S
    type(field_t), target, intent(in) :: u, v, w, h, max_wave_speed
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: Xh
    type(gs_t), target, intent(in) :: gs

    this%S => S
    this%u => u
    this%v => v
    this%w => w
    this%h => h
    this%max_wave_speed => max_wave_speed
    this%msh => msh
    this%Xh => Xh
    this%gs => gs

    call this%S_lag%init(S, 3)

  end subroutine entropy_viscosity_set_fields

  subroutine entropy_viscosity_update_lag(this)
    class(entropy_viscosity_t), intent(inout) :: this

    call this%S_lag%update()

  end subroutine entropy_viscosity_update_lag

  !> Compute low-order viscosity at point i: c_max * h * max_wave_speed
  pure function entropy_viscosity_low_order(this, i) result(visc)
    class(entropy_viscosity_t), intent(in) :: this
    integer, intent(in) :: i
    real(kind=rp) :: visc

    visc = this%c_max * this%h%x(i,1,1,1) * this%max_wave_speed%x(i,1,1,1)

  end function entropy_viscosity_low_order

end module entropy_viscosity
