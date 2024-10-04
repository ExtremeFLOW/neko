! Copyright (c) 2024, The Neko Authors
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
!> Implements the `coriolis_source_term_t` type.
!! Maintainer: Timofey Mukha

module coriolis_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use coriolis_source_term_cpu, only : coriolis_source_term_compute_cpu
  implicit none
  private

  !> This source term adds the Coriolis force.
  type, public, extends(source_term_t) :: coriolis_source_term_t
     !> The rotation vector.
     real(kind=rp) :: omega(3)
     !> The geostrophic wind.
     real(kind=rp) :: u_geo(3) = 0
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => coriolis_source_term_init_from_json
     !> The costrucructor from type components.
     procedure, pass(this) :: init_from_compenents => &
       coriolis_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => coriolis_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => coriolis_source_term_compute
  end type coriolis_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine coriolis_source_term_init_from_json(this, json, fields, coef)
    class(coriolis_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    ! Rotation vector and geostrophic wind
    real(kind=rp), allocatable :: rotation_vec(:), u_geo(:)
    ! Alternative parameters to set the rotation vector
    real(kind=rp) :: omega, phi, f, pi
    real(kind=rp) :: start_time, end_time

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    if (json%valid_path("geostrophic_wind")) then
       call json_get(json, "geostrophic_wind", u_geo)
    else
       allocate(u_geo(3))
       u_geo = 0.0_rp
    end if

    if (json%valid_path("rotation_vector")) then
       call json_get(json, "rotation_vector", rotation_vec)
    else if (json%valid_path("omega") .and. json%valid_path("phi")) then
       call json_get(json, "phi", phi)
       call json_get(json, "omega", omega)

       allocate(rotation_vec(3))
       pi = 4 * atan(1.0_rp)
       rotation_vec(1) = 0.0_rp
       rotation_vec(2) = omega * cos(phi * pi / 180 )
       rotation_vec(3) = omega * sin(phi * pi / 180)
    else if (json%valid_path("f")) then
       call json_get(json, "f", phi)

       allocate(rotation_vec(3))
       rotation_vec(1) = 0.0_rp
       rotation_vec(2) = 0.0_rp
       rotation_vec(3) = 0.5_rp * f
    else
       call neko_error("Specify either rotation_vector, phi and omega, or f &
             & for the Coriolis source term.")
    end if



    call coriolis_source_term_init_from_components(this, fields, rotation_vec, &
          u_geo, coef, start_time, end_time)

  end subroutine coriolis_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param omega The rotation vector.
  !! @param u_geo The geostrophic wind.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine coriolis_source_term_init_from_components(this, fields, omega, &
       u_geo, coef, start_time, end_time)
    class(coriolis_source_term_t), intent(inout) :: this
    class(field_list_t), intent(inout), target :: fields
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: u_geo(3)
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (fields%size() .ne. 3) then
       call neko_error("Number of fields for the Coriolis force must be 3.")
    end if

    this%omega = omega
    this%u_geo = u_geo
  end subroutine coriolis_source_term_init_from_components

  !> Destructor.
  subroutine coriolis_source_term_free(this)
    class(coriolis_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine coriolis_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine coriolis_source_term_compute(this, t, tstep)
    class(coriolis_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("The Coriolis force is only implemented on the CPU")
    else
       call coriolis_source_term_compute_cpu(this%fields, this%omega, &
            this%u_geo)
    end if
  end subroutine coriolis_source_term_compute

end module coriolis_source_term
