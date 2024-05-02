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
!> Implements the `power_iterations_t` type.

module power_iterations
  use num_types, only: rp, dp, sp
  use json_module, only: json_file
  use simulation_component, only: simulation_component_t
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use operators, only: curl
  use case, only: case_t
  use fld_file_output, only: fld_file_output_t
  use json_utils, only: json_get, json_get_or_default
  use fluid_scheme, only: fluid_scheme_t
  use file, only: file_t
  use vector, only: vector_t
  use math, only: vdot3, glsc2, cmult
  use device_math, only: device_cmult
  use comm, only: pe_rank
  use logger, only: neko_log, NEKO_LOG_DEBUG
  implicit none
  private

  !> A simulation component that computes the power_iterations field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: power_iterations_t
     private

     ! ======================================================================= !
     ! Addressable attributes

     real(kind=rp) :: norm_scaling !< Constant for the norm of the velocity field.
     real(kind=rp) :: norm_target !< Target norm for the velocity field.
     real(kind=rp) :: norm_tolerance !< Tolerance for when to rescale the flow.

     ! ======================================================================= !
     ! Definition of shorthands and local variables

     !> Size of the arrays
     integer :: n

     !> X velocity component.
     real(kind=rp), dimension(:,:,:,:), pointer :: u
     !> Y velocity component.
     real(kind=rp), dimension(:,:,:,:), pointer :: v
     !> Z velocity component.
     real(kind=rp), dimension(:,:,:,:), pointer :: w

     !> Volume of the mesh
     real(kind=rp) :: vol
     !> Mass matrix
     real(kind=rp), dimension(:,:,:,:), pointer :: B

     !> The previously used timestep
     real(kind=rp) :: t_old = 0.0_rp
     !> The norm of the velocity field at the previous timestep
     real(kind=rp) :: norm_l2_old = -1.0_rp
     !> The slope of the norm of the velocity field
     real(kind=rp) :: slope_value = 0.0_rp
     !> The number of times the slope has been computed
     real(kind=rp) :: slope_count = 0.0_rp

     !> Upper limit for the norm
     real(kind=rp) :: norm_l2_upper = 1.0_rp
     !> Lower limit for the norm
     real(kind=rp) :: norm_l2_lower = 0.0_rp

     !> Flag to indicate if the flow has been rescaled.
     !! Used to wait for the slope to stabilize.
     logical :: has_rescaled = .false.

     !> Output file
     type(file_t) :: file_output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, public, pass(this) :: init => power_iterations_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
       power_iterations_init_from_attributes
     !> Destructor.
     procedure, public, pass(this) :: free => power_iterations_free
     !> Compute the power_iterations field.
     procedure, public, pass(this) :: compute_ => power_iterations_compute
  end type power_iterations_t

contains

  !> Constructor from json.
  subroutine power_iterations_init_from_json(this, json, case)
    class(power_iterations_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: file_name
    character(len=256) :: header_line

    call this%free()
    call this%init_base(json, case)

    ! Read the json file
    call json_get_or_default(json, 'norm_scaling', &
                             this%norm_scaling, 0.5_rp)
    call json_get_or_default(json, 'norm_tolerance', &
                             this%norm_tolerance, 10.0_rp)
    call json_get_or_default(json, 'output_file', &
                             file_name, 'power_iterations.csv')

    ! Build the header
    this%file_output = file_t(trim(file_name))
    write(header_line, '(A)') 'Time, Norm, Slope, Eigenvalue, Scaling'
    call this%file_output%set_header(header_line)

    ! Initialize the fields
    call power_iterations_init_from_attributes(this)
  end subroutine power_iterations_init_from_json

  !> Actual constructor.
  subroutine power_iterations_init_from_attributes(this)
    class(power_iterations_t), intent(inout) :: this

    ! Temporary field pointers
    type(field_t), pointer :: u_tmp, v_tmp, w_tmp
    real(kind=rp) :: norm_l2

    u_tmp => neko_field_registry%get_field("u")
    v_tmp => neko_field_registry%get_field("v")
    w_tmp => neko_field_registry%get_field("w")

    ! Point the internal fields to the correct data
    this%n = u_tmp%size()
    this%u => u_tmp%x
    this%v => v_tmp%x
    this%w => w_tmp%x

    ! Allocate the mass matrix
    this%B => this%case%fluid%c_Xh%B
    this%vol = this%case%fluid%c_Xh%volume

    ! Setup the target norm for the velocity field
    this%norm_target = sqrt(this%norm_scaling) * norm(this%u, this%v, this%w, &
                                                      this%B, this%vol, this%n)

    this%norm_l2_old = norm_l2
    this%norm_l2_upper = this%norm_tolerance * this%norm_target
    this%norm_l2_lower = this%norm_target / this%norm_tolerance

  end subroutine power_iterations_init_from_attributes

  !> Destructor.
  subroutine power_iterations_free(this)
    class(power_iterations_t), intent(inout) :: this

    ! Deallocate fields
    call this%free_base()
  end subroutine power_iterations_free

  subroutine rescale_fluid(fluid_data, scale)
    use neko_config, only: NEKO_BCKND_DEVICE
    implicit none

    !> Fluid data
    class(fluid_scheme_t), intent(inout) :: fluid_data
    !> Scaling factor
    real(kind=rp), intent(in) :: scale

    ! Local variables
    integer :: i

    ! Scale the velocity fields
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%u%x_d, scale, fluid_data%u%size())
       call device_cmult(fluid_data%v%x_d, scale, fluid_data%v%size())
       call device_cmult(fluid_data%w%x_d, scale, fluid_data%w%size())
    else
       call cmult(fluid_data%u%x, scale, fluid_data%u%size())
       call cmult(fluid_data%v%x, scale, fluid_data%v%size())
       call cmult(fluid_data%w%x, scale, fluid_data%w%size())
    end if

    ! Scale the right hand sides
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%f_x%x_d, scale, fluid_data%f_x%size())
       call device_cmult(fluid_data%f_y%x_d, scale, fluid_data%f_y%size())
       call device_cmult(fluid_data%f_z%x_d, scale, fluid_data%f_z%size())
    else
       call cmult(fluid_data%f_x%x, scale, fluid_data%f_x%size())
       call cmult(fluid_data%f_y%x, scale, fluid_data%f_y%size())
       call cmult(fluid_data%f_z%x, scale, fluid_data%f_z%size())
    end if

    ! Scale the lag terms
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, fluid_data%ulag%size()
          call device_cmult(fluid_data%ulag%lf(i)%x_d, &
                            scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call device_cmult(fluid_data%vlag%lf(i)%x_d, &
                            scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call device_cmult(fluid_data%wlag%lf(i)%x_d, &
                            scale, fluid_data%wlag%lf(i)%size())
       end do
    else
       do i = 1, fluid_data%ulag%size()
          call cmult(fluid_data%ulag%lf(i)%x, &
                     scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call cmult(fluid_data%vlag%lf(i)%x, &
                     scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call cmult(fluid_data%wlag%lf(i)%x, &
                     scale, fluid_data%wlag%lf(i)%size())
       end do
    end if

  end subroutine rescale_fluid

  function norm(x, y, z, B, volume, n)
    real(kind=rp), dimension(n), intent(in) :: x, y, z
    real(kind=rp), dimension(n), intent(in) :: B
    real(kind=rp), intent(in) :: volume
    integer, intent(in) :: n

    real(kind=rp) :: norm
    real(kind=rp), dimension(n) :: tmp

    call vdot3(tmp, x, y, z, x, y, z, n)
    norm = sqrt(glsc2(tmp, B, n) / volume)

  end function norm

  ! function device_norm(x_d, y_d, z_d, B_d, volume, n)
  !   use neko_config, only: NEKO_BCKND_DEVICE
  !   implicit none

  !   type(c_ptr), intent(in) :: x_d, y_d, z_d
  !   type(c_ptr), intent(in) :: B_d
  !   real(kind=rp), intent(in) :: volume
  !   integer, intent(in) :: n

  !   real(kind=rp) :: norm
  !   real(kind=rp), dimension(n) :: tmp

  !   call vdot3(tmp, x_d, y_d, z_d, x_d, y_d, z_d, n)
  !   norm = sqrt(glsc2(tmp, B_d, n) / volume)

  ! end function device_norm

  !> Compute the power_iterations field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine power_iterations_compute(this, t, tstep)
    class(power_iterations_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Local variables
    real(kind=rp) :: scaling_factor
    real(kind=rp) :: norm_l2
    real(kind=rp) :: lambda
    real(kind=rp) :: dt
    character(len=256) :: log_message
    type(vector_t) :: data_line

    dt = t - this%t_old

    ! Compute the norm of the velocity field
    norm_l2 = sqrt(this%norm_scaling) * norm(this%u, this%v, this%w, &
                                             this%B, this%vol, this%n)
    lambda = (log(norm_l2) - log(this%norm_l2_old)) / dt
    scaling_factor = 1.0_rp

    if (norm_l2 .gt. this%norm_l2_upper) then
       this%has_rescaled = .true.

       scaling_factor = this%norm_l2_lower / norm_l2
       call rescale_fluid(this%case%fluid, scaling_factor)
       norm_l2 = this%norm_l2_lower

       this%t_old = t
       This%norm_l2_old = norm_l2
    else if ( norm_l2 .lt. this%norm_l2_lower) then
       this%has_rescaled = .true.

       scaling_factor = this%norm_l2_upper / norm_l2
       call rescale_fluid(this%case%fluid, scaling_factor)
       norm_l2 = this%norm_l2_upper

       this%t_old = t
       This%norm_l2_old = norm_l2
    else if (this%has_rescaled) then
       this%slope_count = this%slope_count + 1.0_rp

       this%slope_value = this%slope_value + &
         ( lambda - this%slope_value ) / (this%slope_count)

    end if

    ! Log the results
    call neko_log%section('Power Iterations')

    write (log_message, '(A7,E20.14)') 'Norm: ', norm_l2
    call neko_log%message(log_message, lvl=NEKO_LOG_DEBUG)
    write (log_message, '(A7,E20.14)') 'Slope: ', this%slope_value
    call neko_log%message(log_message, lvl=NEKO_LOG_DEBUG)
    write (log_message, '(A7,E20.14)') 'Eigen: ', abs(lambda)
    call neko_log%message(log_message, lvl=NEKO_LOG_DEBUG)
    write (log_message, '(A7,E20.14)') 'Scaling: ', scaling_factor

    ! Save to file
    call data_line%init(4)
    data_line%x = [norm_l2, this%slope_value, abs(lambda), scaling_factor]
    call this%file_output%write(data_line, t)

    call neko_log%end_section('Power Iterations')
  end subroutine power_iterations_compute

end module power_iterations

