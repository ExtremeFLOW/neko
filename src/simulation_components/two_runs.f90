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

module two_runs
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fluid_pnpn, only : fluid_pnpn_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use flow_ic, only : set_flow_ic
  use json_utils, only : json_get, json_get_or_default
  use fld_file_output, only : fld_file_output_t
  use mpi_f08, only : MPI_WTIME
  use material_properties, only : material_properties_t
  use utils, only : neko_error
  use time_step_controller
  implicit none
  private

  !> A simulation component that solves a scalar transport equation.
  type, public, extends(simulation_component_t) :: two_runs_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> The scalar scheme.
     !type(fluid_pnpn_t), allocatable :: scheme(:)
     type(fluid_pnpn_t) :: scheme

     !> Output writer.
     type(fld_file_output_t), private :: output

     type(time_step_controller_t) :: dt_controller_pert

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => two_runs_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        two_runs_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => two_runs_free
     !> Compute the scalar field.
     procedure, pass(this) :: compute_ => two_runs_compute
  end type two_runs_t

contains

  !> Constructor from json.
  !! @note We cannot use only this constructor, because `this` is not a
  !! `target´, which is necessary.
  subroutine two_runs_init_from_json(this, json, case)
    class(two_runs_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=:), allocatable :: ic_type
    logical :: dealias
    integer :: p_order, dealias_order

    ! base for simulation component
    call this%init_base(json, case)
    ! initial condition
    call json_get(json, 'initial_condition.type', ic_type)

    ! Uses global case parameter for dealiasing.
    call json_get(this%case%params, 'case.numerics.dealias', dealias)
    call json_get(this%case%params, 'case.numerics.polynomial_order', p_order)
    call json_get_or_default(this%case%params,&
                             'case.numerics.dealiased_polynomial_order',&
                             dealias_order,&
                             int(3.0_rp / 2.0_rp * (p_order + 1)))
    ! Todo
    ! they should all have potentially their own numerics...

    ! Handle output options and call the other constructor.
    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call two_runs_init_from_attributes(this, json, ic_type, dealias,&
                                               dealias_order, filename, dp)
           else
              call two_runs_init_from_attributes(this, json, ic_type, dealias,&
                                               dealias_order, filename, sp)
           end if
       else
            call two_runs_init_from_attributes(this, json, ic_type, dealias,&
                                             dealias_order, filename)
       end if
    else
       call two_runs_init_from_attributes(this, json, ic_type, dealias,&
                                        dealias_order)
    end if


    ! do something with our timestep controller
    ! for now we copy them?
    call this%dt_controller_pert%init(case%params)
  end subroutine two_runs_init_from_json

  !> The second constructor.
  !! Admittedly, the split between the two constructors is quite arbitrary.
  subroutine two_runs_init_from_attributes(this, json, ic_type, dealias, &
                                         dealias_order, filename, precision)
    class(two_runs_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: ic_type
    logical, intent(in) :: dealias
    integer, intent(in) :: dealias_order
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    integer i_fluid
    character(len=10) :: filenizzle

    ! Init the scheme
    ! I feel like I should ONLY have to do this...
    i_fluid = 4
    call this%scheme%init(this%case%msh,  &
                          this%case%fluid%xh%lx, this%case%params, this%case%usr, &
                          this%case%material_properties,i_fluid)

    
    ! Setup file output
    write(filenizzle,'(A,I0.3)') 'bigdog', i_fluid
       if (present(precision)) then
          call this%output%init(precision, filenizzle, 3)
       else
          call this%output%init(sp, filenizzle, 3)
       end if
       this%output%fields%fields(1)%f => this%scheme%u
       this%output%fields%fields(2)%f => this%scheme%v
       this%output%fields%fields(3)%f => this%scheme%w
       call this%case%s%add(this%output, this%output_controller%control_value, &
                            this%output_controller%control_mode)

! HARRY
! do I need all of this, or can I just initialize?
!    ! Set lag arrays
!    ! HARRY
!    call this%scheme%ulag%set(this%scheme%u)
!    call this%scheme%vlag%set(this%scheme%v)
!    call this%scheme%wlag%set(this%scheme%w)
!
!    ! Validate scheme initialization
!    call this%scheme%validate()
!
!    call this%case%fluid%chkp%add_scalar(this%scheme%s)
!    this%case%fluid%chkp%abs1 => this%scheme%abx1
!    this%case%fluid%chkp%abs2 => this%scheme%abx2
!    this%case%fluid%chkp%slag => this%scheme%slag
!
!    ! User boundary condition
!    call this%scheme%set_user_bc(this%case%usr%scalar_user_bc)
!
!    ! Initial condition
!    if (trim(ic_type) .ne. 'user') then
!       call set_scalar_ic(this%scheme%s, this%scheme%c_Xh, this%scheme%gs_Xh,&
!                          ic_type, this%case%params)
!    else
!       call set_scalar_ic(this%scheme%s,  this%scheme%c_Xh, this%scheme%gs_Xh,&
!                          this%case%usr%scalar_user_ic, this%case%params)
!    end if
!
!
!    call this%scheme%validate
  end subroutine two_runs_init_from_attributes

  !> Destructor.
  subroutine two_runs_free(this)
    class(two_runs_t), intent(inout) :: this
    call this%free_base()
    call this%scheme%free()
  end subroutine two_runs_free

  !> Solve the scalar equation.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine two_runs_compute(this, t, tstep)
    class(two_runs_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    ! -------------------------------------
    real(kind=dp) :: start_time, end_time
    character(len=LOG_SIZE) :: log_buf
    ! Tim,
    ! it looks to me that fluid has t and tstep intent inout
    ! and simcomps intent in, implying simcomps shouldn't change 
    ! the underlying time and timstep?
    !
    ! I guess we may need to change the time for our work...
    ! so we should introduce a pseudo time?
    real(kind=rp) :: p_t
    integer :: p_tstep

	 p_t = t
	 p_tstep = tstep

    start_time = MPI_WTIME()
    call neko_log%section('Scalar simcomp')
    call this%scheme%step(p_t, p_tstep, this%case%dt, this%case%ext_bdf,this%dt_controller_pert)
    !call this%scheme%compute(t, tstep, this%case%dt, this%case%ext_bdf)
    end_time = MPI_WTIME()
    write(log_buf, '(A,E15.7)') 'Step time:', end_time-start_time
    call neko_log%end_section(log_buf)

  end subroutine two_runs_compute


  subroutine set_material_properties(properties, json)
    type(material_properties_t), intent(inout) :: properties
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf


!    ! Incorrect user input
!    if (json%valid_path('Pe') .and. (json%valid_path('lambda') .or. &
!                                     json%valid_path('cp'))) then
!       call neko_error("To set the material properties for the scalar,&
!       & either provide Pe OR lambda  and cp in the case file.")
!
!    ! Non-dimensional case
!    else if (json%valid_path('Pe')) then
!       write(log_buf, '(A)') 'Non-dimensional scalar material properties &
!       & input.'
!       call neko_log%message(log_buf, lvl=NEKO_LOG_VERBOSE)
!       write(log_buf, '(A)') 'Specific heat capacity will be set to 1, &
!       & conductivity to 1/Pe.'
!       call neko_log%message(log_buf, lvl=NEKO_LOG_VERBOSE)
!
!       ! Read Pe into lambda for further manipulation.
!       call json_get(json, 'Pe', properties%lambda)
!       write(log_buf, '(A,ES13.6)') 'Pe         :',  properties%lambda
!       call neko_log%message(log_buf)
!
!       ! Set cp and rho to 1 since the setup is non-dimensional.
!       properties%cp = 1.0_rp
!       properties%rho = 1.0_rp
!       ! Invert the Pe to get conductivity
!       properties%lambda = 1.0_rp/properties%lambda
!    ! Dimensional case
!    else
!          call json_get(json, 'lambda', properties%lambda)
!          call json_get(json, 'cp', properties%cp)
!    end if

  end subroutine set_material_properties

end module two_runs
