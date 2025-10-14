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
!> Implements the `centrifugal_source_term_t` type.
!! Maintainer: Adam Peplinski.

module centrifugal_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use centrifugal_source_term_cpu, only : centrifugal_source_term_compute_cpu
  use centrifugal_source_term_device, only : &
       centrifugal_source_term_compute_device
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use time_state, only : time_state_t
  implicit none
  private

  !> This source term adds the centrifugal force.
  !! @details This forcing can be used to perform simulation in a rotating
  !! reference frame and adds a source term in the form
  !! \f$ - \Omega \times (\Omega \times r) \f$, where \f$ \Omega \f$ is the
  !! rotation vector and \f$ r \f$ is the position relative to the reference
  !! point, which is any point lying on the rotation axis. To perform simulation
  !! in the rotating reference frame one needs to use Coriolis source term as
  !! well, taking care of a consistent definition of \f$ \Omega \f$.
  type, public, extends(source_term_t) :: centrifugal_source_term_t
     !> The rotation vector.
     real(kind=rp) :: omega(3)
     !> The reference point
     real(kind=rp) :: ref_point(3)
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => centrifugal_source_term_init_from_json
     !> The costrucructor from type components.
     procedure, pass(this) :: init_from_compenents => &
          centrifugal_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => centrifugal_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => centrifugal_source_term_compute
  end type centrifugal_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable for this source term.
  subroutine centrifugal_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(centrifugal_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    ! Rotation vector and reference point
    real(kind=rp), allocatable :: rotation_vec(:), ref_point(:)
    real(kind=rp) :: start_time, end_time

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    if (json%valid_path("rotation_vector")) then
       call json_get(json, "rotation_vector", rotation_vec)
    else
       call neko_error("Specify rotation_vector &
       &for the centrifugal source term.")
    end if

    if (json%valid_path("reference_point")) then
       call json_get(json, "reference_point", ref_point)
    else
       call neko_error("Specify reference_point &
       &for the centrifugal source term.")
    end if

    call centrifugal_source_term_init_from_components(this, fields, &
         rotation_vec, ref_point, coef, start_time, end_time)

  end subroutine centrifugal_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param omega The rotation vector.
  !! @param ref_point The reference point.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine centrifugal_source_term_init_from_components(this, fields, omega, &
       ref_point, coef, start_time, end_time)
    class(centrifugal_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: ref_point(3)
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (fields%size() .ne. 3) then
       call neko_error("Number of fields for the centrifugal force must be 3.")
    end if

    this%omega = omega
    this%ref_point = ref_point
  end subroutine centrifugal_source_term_init_from_components

  !> Destructor.
  subroutine centrifugal_source_term_free(this)
    class(centrifugal_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine centrifugal_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine centrifugal_source_term_compute(this, time)
    class(centrifugal_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call centrifugal_source_term_compute_device(this%omega, this%ref_point, &
            this%fields)
    else
       call centrifugal_source_term_compute_cpu(this%omega, this%ref_point, &
            this%fields)
    end if
  end subroutine centrifugal_source_term_compute

end module centrifugal_source_term
