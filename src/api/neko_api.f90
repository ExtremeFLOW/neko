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

  !> Retrive the current time of a case
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param time The case's current time
  function neko_api_case_time(case_iptr) result(time) &
       bind(c, name="neko_case_time")
    integer(c_intptr_t), intent(inout) :: case_iptr
    real(kind=c_double) :: time
    type(case_t), pointer :: C
    type(c_ptr) :: cptr

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       time = C%time%t
    else
       call neko_error('Invalid Neko case')
    end if

  end function neko_api_case_time

  !> Retrive the end time of a case
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param end_time The end time of a case
  function neko_api_case_end_time(case_iptr) result(end_time) &
       bind(c, name="neko_case_end_time")
    integer(c_intptr_t), intent(inout) :: case_iptr
    real(kind=c_double) :: end_time
    type(case_t), pointer :: C
    type(c_ptr) :: cptr

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       end_time = C%time%end_time
    else
       call neko_error('Invalid Neko case')
    end if

  end function neko_api_case_end_time

  !> Retrive the time-step of a case
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param tstep The current time-step of a case
  function neko_api_case_tstep(case_iptr) result(tstep) &
    bind(c, name="neko_case_tstep")
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    integer(c_int) :: tstep

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       tstep = C%time%tstep
    else
       call neko_error('Invalid Neko case')
    end if

  end function neko_api_case_tstep

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
  subroutine neko_api_step(case_iptr) &
       bind(c, name="neko_step")
    use time_step_controller, only : time_step_controller_t
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    type(time_step_controller_t), save, allocatable :: dt_controller
    real(kind=dp), save :: step_loop_start = 0.0_dp
    real(kind=rp), save :: cfl

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)

       if (.not. allocated(dt_controller)) then
          allocate(dt_controller)
          call dt_controller%init(C%params)
       end if

       if (C%time%tstep .eq. 0) then
          call simulation_init(C, dt_controller)

          cfl = C%fluid%compute_cfl(C%time%dt)
          step_loop_start = MPI_Wtime()
       end if

       if (C%time%t .lt. C%time%end_time) then
          call simulation_step(C, dt_controller, cfl, step_loop_start)
       end if

       if (C%time%t .ge. C%time%end_time) then
          call simulation_finalize(C)
          if (allocated(dt_controller)) then
             deallocate(dt_controller)
          end if
       end if

    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_step

  !> Execute the Case's output controller
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param t The time value
  !! @param tstep The current time-stepper iteration
  subroutine neko_api_output_ctrl_execute(case_iptr, force_output) &
       bind(c, name="neko_output_ctrl_execute")
    integer(c_intptr_t), intent(inout) :: case_iptr
    logical(kind=c_bool), value :: force_output
    logical :: f_force_output
    type(case_t), pointer :: C
    type(c_ptr) :: cp
    type(time_state_t) :: f_time

    cp = transfer(case_iptr, c_null_ptr)
    if (c_associated(cp)) then
       call c_f_pointer(cp, C)
    else
       call neko_error('Invalid Neko case')
    end if

    f_force_output = transfer(force_output, f_force_output)
    call C%output_controller%execute(C%time, f_force_output)

  end subroutine neko_api_output_ctrl_execute

  !> Retrive a pointer to a flow field
  !! @param field_name Field registry entry
  function neko_api_field(field_name) result(field_ptr) &
       bind(c, name='neko_field')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    character(len=8192) :: name
    type(field_t), pointer :: field
    type(c_ptr) :: field_ptr
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))

    field_ptr = c_loc(field%x)

  end function neko_api_field

  !> Retrive the order of a field
  !! @param field_name Field registry entry
  function neko_api_field_order(field_name) result(field_lx) &
       bind(c, name='neko_field_order')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    character(len=8192) :: name
    type(field_t), pointer :: field
    integer(c_int) :: field_lx
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))

    field_lx = field%Xh%lx

  end function neko_api_field_order

  !> Retrive the number of elements in a field
  !! @param field_name Field registry entry
  function neko_api_field_nelements(field_name) result(field_nelv) &
       bind(c, name='neko_field_nelements')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    character(len=8192) :: name
    type(field_t), pointer :: field
    integer(c_int) :: field_nelv
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))

    field_nelv = field%msh%nelv

  end function neko_api_field_nelements

  !> Retrive the total number of degrees of freedom of a field
  !! @param field_name Field registry entry
  function neko_api_field_size(field_name) result(field_size) &
       bind(c, name='neko_field_size')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    character(len=8192) :: name
    type(field_t), pointer :: field
    integer(c_int) :: field_size
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))

    field_size = field%dof%size()

  end function neko_api_field_size


end module neko_api
