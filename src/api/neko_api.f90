! Copyright (c) 2022-2025, The Neko Authors
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
!> Neko C API
module neko_api
  use neko
  use, intrinsic :: iso_c_binding
  implicit none
  private

contains

  !> Initialise Neko
  subroutine neko_api_init() bind(c, name="neko_init")

    call neko_init()

  end subroutine neko_api_init

  !> Finalize Neko
  subroutine neko_api_finalize() bind(c, name="neko_finalize")

    call neko_finalize()

  end subroutine neko_api_finalize

  !> Display job information
  subroutine neko_api_job_info() bind(c, name="neko_job_info")
    logical :: initialized

    call MPI_Initialized(initialized)

    if (.not.initialized) then
       call neko_warning('Neko has not been initialised')
    else
       call neko_job_info()
       call neko_log%newline()
    end if

  end subroutine neko_api_job_info

  !> Initalise a Neko case
  !! @param case_json Serialised JSON object describing the case
  !! @param case_iptr Opaque pointer for the Neko case
  subroutine neko_api_case_init(case_json, case_len, case_iptr) &
       bind(c, name="neko_case_init")
    type(c_ptr), intent(in) :: case_json
    integer(c_int), value :: case_len
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(json_file) :: json_case
    type(case_t), pointer :: C
    type(c_ptr) :: cp

    ! Convert passed in serialised JSON object into a Fortran
    ! character string and create a json_file object
    if (c_associated(case_json)) then
       block
         character(kind=c_char,len=case_len+1),pointer :: s
         character(len=:), allocatable :: fcase_json
         call c_f_pointer(case_json, s)
         fcase_json = s(1:case_len)
         call json_case%load_from_string(fcase_json)
         deallocate(fcase_json)
         nullify(s)
       end block
    end if

    allocate(C)

    !
    ! Create case
    !
    call case_init(C, json_case)

    !
    ! Create simulation components
    !
    call neko_simcomps%init(C)


    cp = c_loc(C)
    case_iptr = transfer(cp, 0_c_intptr_t)

  end subroutine neko_api_case_init


  !> Destroy a Neko case
  !! @param case_iptr Opaque pointer for the Neko case
  subroutine neko_api_case_free(case_iptr) bind(c, name="neko_case_free")
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cp

    cp = transfer(case_iptr, c_null_ptr)
    if (c_associated(cp)) then
       call c_f_pointer(cp, C)
       call case_free(c)
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_case_free

  !> Solve a neko case
  !! @param case_iptr Opaque pointer for the Neko case
  subroutine neko_api_solve(case_iptr) bind(c, name="neko_solve")
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cp

    cp = transfer(case_iptr, c_null_ptr)
    if (c_associated(cp)) then
       call c_f_pointer(cp, C)
       call neko_solve(C)
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_solve

  !> Compute a time-step for a neko case
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param t The time value
  !! @param tstep The current time-stepper iteration
  subroutine neko_api_step(case_iptr, t, tstep) bind(c, name="neko_step")
    use time_step_controller, only : time_step_controller_t
    integer(c_intptr_t), intent(inout) :: case_iptr
    real(kind=c_double), value :: t
    integer(c_int), value :: tstep
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: rho, mu, cp, lambda
    type(time_step_controller_t) :: dt_controller
    real(kind=rp) :: cfl_avrg = 0.0_rp
    real(kind=rp) :: cfl
    integer :: i

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
    else
       call neko_error('Invalid Neko case')
    end if

    call dt_controller%init(C%params)

    write(log_buf, '(A4,E15.7)') 't = ', t
    call neko_log%message(repeat('-', 64), NEKO_LOG_QUIET)
    call neko_log%message(log_buf, NEKO_LOG_QUIET)
    call neko_log%message(repeat('-', 64), NEKO_LOG_QUIET)

    call neko_log%begin()

    cfl = C%fluid%compute_cfl(C%dt)
    call dt_controller%set_dt(C%dt, cfl, cfl_avrg, tstep)

    do i = 10, 2, -1
       C%tlag(i) = C%tlag(i-1)
       C%dtlag(i) = C%dtlag(i-1)
    end do

    C%dtlag(i) = C%dt
    C%tlag(1) = t
    if (C%ext_bdf%ndiff .eq. 0 ) then
       C%dtlag(2) = C%dt
       C%tlag(2) = t
    end if

    t = t + C%dt

    call C%ext_bdf%set_coeffs(C%dtlag)

    ! Run the preprocessing
    call neko_log%section('Preprocessing')
    call neko_simcomps%preprocess(t, tstep)
    call neko_log%end_section()

    call neko_log%section('Fluid')
    call C%fluid%step(t, tstep, C%dt, C%ext_bdf, dt_controller)
    call neko_log%end_section()

    ! Scalar step
    if (allocated(C%scalar)) then
       call neko_log%section('Scalar')
       call C%scalar%step(t, tstep, C%dt, C%ext_bdf, dt_controller)
       call neko_log%end_section()

       !> @todo Temporary fix until we have reworked the material properties
       cp = C%scalar%cp
       lambda = C%scalar%lambda
    end if

    call neko_log%section('Postprocessing')
    ! Execute all simulation components
    call neko_simcomps%compute(t, tstep)
    call neko_log%end_section()


    !> @todo Temporary fix until we have reworked the material properties
    rho = C%fluid%rho
    mu = C%fluid%mu

    ! Update material properties
    call C%usr%material_properties(t, tstep, rho, mu, cp, lambda, C%params)

    !> @todo Temporary fix until we have reworked the material properties
    C%fluid%rho = rho
    C%fluid%mu = mu
    call C%fluid%update_material_properties()

    if (allocated(C%scalar)) then
       C%scalar%cp = cp
       C%scalar%lambda = lambda
       call C%scalar%update_material_properties()
    end if

    call neko_log%end()
    call neko_log%newline()

  end subroutine neko_api_step

  !> Execute the Case's output controller
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param t The time value
  !! @param tstep The current time-stepper iteration
  subroutine neko_api_output_ctrl_execute(case_iptr, t, tstep, force_output) &
       bind(c, name="neko_output_ctrl_execute")
    integer(c_intptr_t), intent(inout) :: case_iptr
    real(kind=c_double), value :: t
    integer(c_int), value :: tstep
    logical(kind=c_bool), value :: force_output
    real(kind=rp) :: f_t
    logical :: f_force_output
    type(case_t), pointer :: C
    type(c_ptr) :: cp

    cp = transfer(case_iptr, c_null_ptr)
    if (c_associated(cp)) then
       call c_f_pointer(cp, C)
    else
       call neko_error('Invalid Neko case')
    end if

    f_force_output = transfer(force_output, f_force_output)
    f_t = transfer(t, f_t)

    call C%output_controller%execute(f_t, tstep, f_force_output)

  end subroutine neko_api_output_ctrl_execute



end module neko_api
