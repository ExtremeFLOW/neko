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
!> Implements the `translation_source_term_t` type.

module translation_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get_or_lookup, json_get_or_lookup_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use field, only : field_t
  use field_math, only : field_col3, field_add2
  use operators, only : div
  use device, only : device_sync
  use registry, only : neko_registry
  use time_state, only : time_state_t
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  !> This source term adds the translation force.
  type, public, extends(source_term_t) :: translation_source_term_t
     !> The domain velocity.
     real(kind=rp) :: domain_vel(3)
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => translation_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          translation_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => translation_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => translation_source_term_compute
  end type translation_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable for this source term.
  subroutine translation_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(translation_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp), allocatable :: domain_vel(:)
    ! Alternative parameters to set the rotation vector
    real(kind=rp) :: start_time, end_time

    call json_get_or_lookup_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_lookup_or_default(json, "end_time", end_time, huge(0.0_rp))

    if (json%valid_path("domain_velocity")) then
       call json_get_or_lookup(json, "domain_velocity", domain_vel)

       if (size(domain_vel) .ne. 3) then
          call neko_error("The domain velocity should have 3 components.")
       end if
    else
       allocate(domain_vel(3))
       domain_vel = 0.0_rp
    end if

    call translation_source_term_init_from_components(this, fields, &
         domain_vel, coef, start_time, end_time)

  end subroutine translation_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param domain_vel The domain velocity.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine translation_source_term_init_from_components(this, fields, domain_vel, &
       coef, start_time, end_time)
    class(translation_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    real(kind=rp), intent(in) :: domain_vel(3)
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (fields%size() .ne. 3) then
       call neko_error("Number of fields for the translation force must be 3.")
    end if

    this%domain_vel = domain_vel
  end subroutine translation_source_term_init_from_components

  !> Destructor.
  subroutine translation_source_term_free(this)
    class(translation_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine translation_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  !! @note In this case we get the convection term produced by
  !! the mesh velocity, which comes from using the relative velocity
  !! as advection velocity. The term comes from ((u-u_mesh) . div) u
  !! Here we calculate - (u_mesh . div) u. In conservative form so
  !! - div (u_mesh u). Then we put it in the RHS as +
  subroutine translation_source_term_compute(this, time)
    class(translation_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: work1, work2, work3
    type(field_t), pointer :: work4, work5, work6, work7
    integer :: tmp_index(7)
    type(field_t), pointer :: fx, fy, fz

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    fx => this%fields%get_by_index(1)
    fy => this%fields%get_by_index(2)
    fz => this%fields%get_by_index(3)


    call neko_scratch_registry%request_field(work1, tmp_index(1), .false.)
    call neko_scratch_registry%request_field(work2, tmp_index(2), .false.)
    call neko_scratch_registry%request_field(work3, tmp_index(3), .false.)
    call neko_scratch_registry%request_field(work4, tmp_index(4), .false.)
    call neko_scratch_registry%request_field(work5, tmp_index(5), .false.)
    call neko_scratch_registry%request_field(work6, tmp_index(6), .false.)
    call neko_scratch_registry%request_field(work7, tmp_index(7), .false.)


    ! Assign the velocities
    work1 = this%domain_vel(1)
    work2 = this%domain_vel(2)
    work3 = this%domain_vel(3)

    ! umesh_j * du/dx_j = div (umesh u) - u * div umesh but div umesh = 0 so we just calculate div (umesh u)
    call field_col3(work4, u, work1, u%size())
    call field_col3(work5, u, work2, u%size())
    call field_col3(work6, u, work3, u%size())
    call div(work7%x, work4%x, work5%x, work6%x, this%coef)
    call field_add2(fx, work7, work7%size())

    ! umesh_j * dv/dx_j
    call field_col3(work4, v, work1, v%size())
    call field_col3(work5, v, work2, v%size())
    call field_col3(work6, v, work3, v%size())
    call div(work7%x, work4%x, work5%x, work6%x, this%coef)
    call field_add2(fy, work7, work7%size())

    ! umesh_j * dw/dx_j
    call field_col3(work4, w, work1, w%size())
    call field_col3(work5, w, work2, w%size())
    call field_col3(work6, w, work3, w%size())
    call div(work7%x, work4%x, work5%x, work6%x, this%coef)
    call field_add2(fz, work7, work7%size())

    call device_sync()

    ! Release the scratch fields
    call neko_scratch_registry%relinquish_field(tmp_index)


  end subroutine translation_source_term_compute

end module translation_source_term
