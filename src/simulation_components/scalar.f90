! Copyright (c) 2023, The Neko Authors
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
!
!> Implements the `scalar_t` type.

module scalar
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use scalar_pnpn, only : scalar_pnpn_t
  use logger, only : neko_log, LOG_SIZE
  use scalar_ic, only : set_scalar_ic
  use json_utils, only : json_get
  use fld_file_output, only : fld_file_output_t
  use mpi_f08, only : MPI_WTIME
  implicit none
  private

  !> A simulation component that solves a scalar transport equation.
  type, public, extends(simulation_component_t) :: scalar_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> The scalar scheme.
     type(scalar_pnpn_t) :: scheme

     !> Output writer.
     type(fld_file_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => scalar_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        scalar_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => scalar_free
     !> Compute the scalar field.
     procedure, pass(this) :: compute_ => scalar_compute
  end type scalar_t

contains

  !> Constructor from json.
  subroutine scalar_init_from_json(this, json, case)
    class(scalar_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=:), allocatable :: ic_type

    call this%init_base(json, case)

    call json_get(json, 'initial_condition.type', ic_type)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call scalar_init_from_attributes(this, ic_type, filename, dp)
           else
              call scalar_init_from_attributes(this, ic_type, filename, sp)
           end if
       else
           call scalar_init_from_attributes(this, ic_type, filename)
       end if
    else
       call scalar_init_from_attributes(this, ic_type)
    end if

  end subroutine scalar_init_from_json

  !> Actual constructor.
  subroutine scalar_init_from_attributes(this, ic_type, filename, precision)
    class(scalar_t), target, intent(inout) :: this
    character(len=*), intent(in) :: ic_type
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    !allocate(this%scheme)

    call this%scheme%init(this%case%msh, this%case%fluid%c_Xh, &
       this%case%fluid%gs_Xh, this%case%params, this%case%usr, &
       this%case%material_properties)

    call this%scheme%set_user_bc(this%case%usr%scalar_user_bc)

    call this%scheme%slag%set(this%scheme%s)
    call this%scheme%validate()

    call this%case%fluid%chkp%add_scalar(this%scheme%s)
    this%case%fluid%chkp%abs1 => this%scheme%abx1
    this%case%fluid%chkp%abs2 => this%scheme%abx2
    this%case%fluid%chkp%slag => this%scheme%slag

    ! User boundary condition
    call this%scheme%set_user_bc(this%case%usr%scalar_user_bc)

    ! Initial condition
    if (trim(ic_type) .ne. 'user') then
       call set_scalar_ic(this%scheme%s, this%scheme%c_Xh, this%scheme%gs_Xh,&
                          ic_type, this%case%params)
    else
       call set_scalar_ic(this%scheme%s,  this%scheme%c_Xh, this%scheme%gs_Xh,&
                          this%case%usr%scalar_user_ic, this%case%params)
    end if

    ! Setup file output
    if (present(filename)) then
       if (present(precision)) then
          call this%output%init(precision, filename, 1)
       else
          call this%output%init(sp, filename, 1)
       end if
       this%output%fields%fields(1)%f => this%scheme%s
       call this%case%s%add(this%output, this%output_controller%control_value, &
                            this%output_controller%control_mode)
    else
       call this%case%f_out%fluid%append(this%scheme%s)
    end if
  end subroutine scalar_init_from_attributes

  !> Destructor.
  subroutine scalar_free(this)
    class(scalar_t), intent(inout) :: this
    call this%free_base()

    call this%scheme%free()
    !deallocate(this%scheme)
  end subroutine scalar_free

  !> Solve the scalar equation.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine scalar_compute(this, t, tstep)
    class(scalar_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=dp) :: start_time, end_time
    character(len=LOG_SIZE) :: log_buf


    start_time = MPI_WTIME()
    call neko_log%section('Scalar simcomp')
    call this%scheme%step(t, tstep, this%case%dt, this%case%ext_bdf)
    end_time = MPI_WTIME()
    write(log_buf, '(A,E15.7)') 'Step time:', end_time-start_time
    call neko_log%end_section(log_buf)

  end subroutine scalar_compute

end module scalar
