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

  interface
     !> Register callbacks
     !! @param user User interface type
     !! @param initial_cb Initial condition callback
     !! @param preprocess_cb Pre timestep callback
     !! @param compute_cb End of timestep callback
     !! @param dirichlet_cb User boundary condition callback
     !! @param material_cb Material properties callback
     !! @param source_cb Source term callback
     module subroutine neko_api_user_cb_register(user, initial_cb, &
          preprocess_cb, compute_cb, dirichlet_cb, material_cb, source_cb)
       type(user_t), intent(inout) :: user
       type(c_funptr), value :: initial_cb, preprocess_cb, compute_cb
       type(c_funptr), value :: dirichlet_cb, material_cb, source_cb
     end subroutine neko_api_user_cb_register
  end interface

  interface
     !> Retrive a pointer to a field for the currently active callback
     !! @param field_name Field list entry
     module function neko_api_user_cb_get_field_by_name(field_name) result(f)
       character(len=*), intent(in) :: field_name
       type(field_t), pointer :: f
     end function neko_api_user_cb_get_field_by_name
  end interface

  interface
     !> Retrive a pointer to a field for the currently active callback
     !! @param field_idx Field index in the field list
     module function neko_api_user_cb_get_field_by_index(field_idx) result(f)
       integer, intent(in) :: field_idx
       type(field_t), pointer :: f
     end function neko_api_user_cb_get_field_by_index
  end interface

  interface neko_api_user_cb_get_field
     module procedure neko_api_user_cb_get_field_by_name, &
          neko_api_user_cb_get_field_by_index
  end interface neko_api_user_cb_get_field

contains

  !> Initialise Neko
  subroutine neko_api_init() bind(c, name="neko_init")

    call neko_init()

  end subroutine neko_api_init

  !> Finalize Neko
  subroutine neko_api_finalize() bind(c, name="neko_finalize")

    call neko_finalize()

  end subroutine neko_api_finalize

  !> Initialise Neko device layer
  subroutine neko_api_device_init() bind(c, name="neko_device_init")

    call device_init

  end subroutine neko_api_device_init

  !> Finalize Neko device layer
  subroutine neko_api_device_finalize() bind(c, name="neko_device_finalize")

    call device_finalize

  end subroutine neko_api_device_finalize

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

  !> Initialise a Neko field registry
  subroutine neko_api_field_registry_init() bind(c, name="neko_field_registry_init")

    call neko_field_registry%init()

  end subroutine neko_api_field_registry_init

  !> Destroy a Neko field registry
  subroutine neko_api_field_registry_free() bind(c, name="neko_field_registry_free")

    call neko_field_registry%free()

  end subroutine neko_api_field_registry_free

  !> Allocate memory for a Neko case
  !! @param case_iptr Opaque pointer for the created Neko case
  subroutine neko_api_case_allocate(case_iptr) &
       bind(c, name="neko_case_allocate")
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cp

    allocate(C)
    cp = c_loc(C)
    case_iptr = transfer(cp, 0_c_intptr_t)

  end subroutine neko_api_case_allocate

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

    ! Check if the case has already been allocated
    ! e.g. if a user callback has been injected
    cp = transfer(case_iptr, c_null_ptr)
    if (c_associated(cp)) then
       call c_f_pointer(cp, C)
    else
       allocate(C)
       cp = c_loc(C)
       case_iptr = transfer(cp, 0_c_intptr_t)
    end if

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

    !
    ! Create case
    !
    call case_init(C, json_case)

    !
    ! Create simulation components
    !
    call neko_simcomps%init(C)

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
    real(kind=c_rp) :: time
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
    real(kind=c_rp) :: end_time
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
    type(json_file) :: dt_params
    type(time_step_controller_t), save, allocatable :: dt_controller
    real(kind=dp), save :: step_loop_start = 0.0_dp

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)

       if (.not. allocated(dt_controller)) then
          allocate(dt_controller)
          call json_get(C%params, 'case.time', dt_params)
          call dt_controller%init(dt_params)
       end if

       if (C%time%tstep .eq. 0) then
          call simulation_init(C, dt_controller)

          step_loop_start = MPI_Wtime()
       end if

       if (.not. C%time%is_done()) then
          call simulation_step(C, dt_controller, step_loop_start)
       end if

       if (C%time%is_done()) then
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

  !> Retrive the dofmap associated with a field
  !! @param field_name Field registry entry
  !! @param dof_ptr Pointer to  unique degrees of freedom
  !! @param x_ptr Pointer to x-coordinates
  !! @param x_ptr Pointer to y-coordinates
  !! @param x_ptr Pointer to z-coordinates
  subroutine neko_api_field_dofmap(field_name, dof_ptr, x_ptr, y_ptr, z_ptr) &
       bind(c, name='neko_field_dofmap')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    type(c_ptr), intent(inout) :: dof_ptr, x_ptr, y_ptr, z_ptr
    character(len=8192) :: name
    type(field_t), pointer :: field
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))
    call neko_api_wrap_dofmap(field%dof, dof_ptr, x_ptr, y_ptr, z_ptr)

  end subroutine neko_api_field_dofmap

  !> Retrive the dofmap associated with a case's fluid solver
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param dof_ptr Pointer to  unique degrees of freedom
  !! @param x_ptr Pointer to x-coordinates
  !! @param x_ptr Pointer to y-coordinates
  !! @param x_ptr Pointer to z-coordinates
  !! @param size Number of dofs
  subroutine neko_api_case_fluid_dofmap(case_iptr, dof_ptr, &
       x_ptr, y_ptr, z_ptr, size) bind(c, name='neko_case_fluid_dofmap')
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    type(c_ptr), intent(inout) :: dof_ptr, x_ptr, y_ptr, z_ptr
    integer, intent(inout) :: size

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       call neko_api_wrap_dofmap(C%fluid%dm_Xh, dof_ptr, x_ptr, y_ptr, z_ptr)
       size = C%fluid%dm_Xh%size()
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_case_fluid_dofmap

  !> Helper function to assign pointers to a dofmap's data
  subroutine neko_api_wrap_dofmap(dm, dof_ptr, x_ptr, y_ptr, z_ptr)
    type(dofmap_t), target, intent(in) :: dm
    type(c_ptr), intent(inout) :: dof_ptr, x_ptr, y_ptr, z_ptr

    dof_ptr = c_loc(dm%dof)
    x_ptr = c_loc(dm%x)
    y_ptr = c_loc(dm%y)
    z_ptr = c_loc(dm%z)

  end subroutine neko_api_wrap_dofmap

  !> Retrive the space associated with a field
  !! @param field_name Field registry entry
  !! @param lx Polynomial dimension in each direction
  !! @param zg Pointer to quadrature points
  !! @param dr_inv Pointer to 1/dist quadrature points
  !! @param ds_inv Pointer to 1/dist quadrature points
  !! @param dt_inv Pointer to 1/dist quadrature points
  !! @param wx Pointer to quadrature weights
  !! @param wy Pointer to quadrature weights
  !! @param wz Pointer to quadrature weights
  !! @param dx Pointer to derivative operator \f$ D_1 \f$
  !! @param dy Pointer to derivative operator \f$ D_2 \f$
  !! @param dz Pointer to derivative operator \f$ D_3 \f$
  subroutine neko_api_field_space(field_name, lx, zg, &
       dr_inv, ds_inv, dt_inv, wx, wy, wz, dx, dy, dz) &
       bind(c, name='neko_field_space')
    character(kind=c_char), dimension(*), intent(in) :: field_name
    integer, intent(inout) :: lx
    type(c_ptr), intent(inout) :: zg, dr_inv, ds_inv, dt_inv
    type(c_ptr), intent(inout) :: wx, wy, wz, dx, dy, dz
    character(len=8192) :: name
    type(field_t), pointer :: field
    integer :: len

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    field => neko_field_registry%get_field(trim(name(1:len)))
    call neko_api_wrap_space(field%Xh, lx, zg, dr_inv, ds_inv, dt_inv, &
         wx, wy, wz, dx, dy, dz)

  end subroutine neko_api_field_space

  !> Retrive the space associated with a case's fluid solver
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param lx Polynomial dimension in each direction
  !! @param zg Pointer to quadrature points
  !! @param dr_inv Pointer to 1/dist quadrature points
  !! @param ds_inv Pointer to 1/dist quadrature points
  !! @param dt_inv Pointer to 1/dist quadrature points
  !! @param wx Pointer to quadrature weights
  !! @param wy Pointer to quadrature weights
  !! @param wz Pointer to quadrature weights
  !! @param dx Pointer to derivative operator \f$ D_1 \f$
  !! @param dy Pointer to derivative operator \f$ D_2 \f$
  !! @param dz Pointer to derivative operator \f$ D_3 \f$
  subroutine neko_api_case_fluid_space(case_iptr, lx, zg, &
       dr_inv, ds_inv, dt_inv, wx, wy, wz, dx, dy, dz) &
       bind(c, name='neko_case_fluid_space')
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    integer, intent(inout) :: lx
    type(c_ptr), intent(inout) :: zg, dr_inv, ds_inv, dt_inv
    type(c_ptr), intent(inout) :: wx, wy, wz, dx, dy, dz


    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       call neko_api_wrap_space(C%fluid%Xh, lx, zg, dr_inv, ds_inv, dt_inv, &
            wx, wy, wz, dx, dy, dz)
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_case_fluid_space

  !> Helper function to assign pointers to a space's data
  subroutine neko_api_wrap_space(Xh, lx, zg, dr_inv, ds_inv, dt_inv, &
       wx, wy, wz, dx, dy, dz)
    type(space_t), target, intent(in) :: Xh
    integer, intent(inout) :: lx
    type(c_ptr), intent(inout) :: zg, dr_inv, ds_inv, dt_inv
    type(c_ptr), intent(inout) :: wx, wy, wz, dx, dy, dz

    lx = Xh%lx
    zg = c_loc(Xh%zg)
    dr_inv = c_loc(Xh%dr_inv)
    ds_inv = c_loc(Xh%ds_inv)
    dt_inv = c_loc(Xh%dt_inv)
    wx = c_loc(Xh%wx)
    wy = c_loc(Xh%wy)
    wz = c_loc(Xh%wz)
    dx = c_loc(Xh%dx)
    dy = c_loc(Xh%dy)
    dz = c_loc(Xh%dz)

  end subroutine neko_api_wrap_space

  !> Retrive the coefficient associated with a case's fluid solver
  !! @param case_iptr Opaque pointer for the Neko case
  subroutine neko_api_case_fluid_coef(case_iptr, G11, G22, G33, G12, G13, G23, &
       mult, dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, &
       drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jac, B, area, nx, ny, nz) bind(c, name='neko_case_fluid_coef')
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(case_t), pointer :: C
    type(c_ptr) :: cptr
    type(c_ptr), intent(inout) :: G11, G22, G33, G12, G13, G23
    type(c_ptr), intent(inout) :: mult
    type(c_ptr), intent(inout) :: dxdr, dydr, dzdr
    type(c_ptr), intent(inout) :: dxds, dyds, dzds
    type(c_ptr), intent(inout) :: dxdt, dydt, dzdt
    type(c_ptr), intent(inout) :: drdx, drdy, drdz
    type(c_ptr), intent(inout) :: dsdx, dsdy, dsdz
    type(c_ptr), intent(inout) :: dtdx, dtdy, dtdz
    type(c_ptr), intent(inout) :: jac, B, area, nx, ny, nz


    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       G11 = c_loc(C%fluid%c_Xh%G11)
       G22 = c_loc(C%fluid%c_Xh%G22)
       G33 = c_loc(C%fluid%c_Xh%G33)
       G12 = c_loc(C%fluid%c_Xh%G12)
       G13 = c_loc(C%fluid%c_Xh%G13)
       G23 = c_loc(C%fluid%c_Xh%G23)
       mult = c_loc(C%fluid%c_Xh%mult)
       dxdr = c_loc(C%fluid%c_Xh%dxdr)
       dydr = c_loc(C%fluid%c_Xh%dydr)
       dzdr = c_loc(C%fluid%c_Xh%dzdr)
       dxds = c_loc(C%fluid%c_Xh%dxds)
       dyds = c_loc(C%fluid%c_Xh%dyds)
       dzds = c_loc(C%fluid%c_Xh%dzds)
       dxdt = c_loc(C%fluid%c_Xh%dxdt)
       dydt = c_loc(C%fluid%c_Xh%dydt)
       dzdt = c_loc(C%fluid%c_Xh%dzdt)
       drdx = c_loc(C%fluid%c_Xh%drdx)
       drdy = c_loc(C%fluid%c_Xh%drdy)
       drdz = c_loc(C%fluid%c_Xh%drdz)
       dsdx = c_loc(C%fluid%c_Xh%dsdx)
       dsdy = c_loc(C%fluid%c_Xh%dsdy)
       dsdz = c_loc(C%fluid%c_Xh%dsdz)
       dtdx = c_loc(C%fluid%c_Xh%dtdx)
       dtdy = c_loc(C%fluid%c_Xh%dtdy)
       dtdz = c_loc(C%fluid%c_Xh%dtdz)
       jac = c_loc(C%fluid%c_Xh%jac)
       B = c_loc(C%fluid%c_Xh%B)
       area = c_loc(C%fluid%c_Xh%area)
       nx = c_loc(C%fluid%c_Xh%nx)
       ny = c_loc(C%fluid%c_Xh%ny)
       nz = c_loc(C%fluid%c_Xh%nz)
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_case_fluid_coef

  !> Setup user-provided callbacks
  !! @param case_iptr Opaque pointer for the Neko case
  !! @param initial_cb Initial condition callback
  !! @param preprocess_cb Pre timestep callback
  !! @param compute_cb End of timestep callback
  !! @param dirichlet_cb User boundary condition callback
  !! @param material_cb Material properties callback
  !! @param source_cb Source term callback
  subroutine neko_api_user_setup(case_iptr, initial_cb, preprocess_cb, &
       compute_cb, dirichlet_cb, material_cb, source_cb) &
       bind(c, name='neko_user_setup')
    integer(c_intptr_t), intent(inout) :: case_iptr
    type(c_funptr), value :: initial_cb, preprocess_cb, compute_cb
    type(c_funptr), value :: dirichlet_cb, material_cb, source_cb
    type(case_t), pointer :: C
    type(c_ptr) :: cptr

    cptr = transfer(case_iptr, c_null_ptr)
    if (c_associated(cptr)) then
       call c_f_pointer(cptr, C)
       call neko_api_user_cb_register(C%user, initial_cb, preprocess_cb, &
            compute_cb, dirichlet_cb, material_cb, source_cb)
    else
       call neko_error('Invalid Neko case')
    end if

  end subroutine neko_api_user_setup

  !> Retrive a pointer to a user callback field
  !! @param field_name Field list entry
  function neko_api_user_cb_field_by_name(field_name) result(field_ptr) &
       bind(c, name='neko_cb_field_by_name')
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

    field => neko_api_user_cb_get_field(trim(name(1:len)))

    field_ptr = c_loc(field%x)

  end function neko_api_user_cb_field_by_name

  !> Retrive a pointer to a user callback field
  !! @param field_idx Field index in the field list
  function neko_api_user_cb_field_by_index(field_idx) &
       result(field_ptr) bind(c, name='neko_cb_field_by_index')
    integer, intent(in) :: field_idx
    type(field_t), pointer :: field
    type(c_ptr) :: field_ptr

    field => neko_api_user_cb_get_field(field_idx)

    field_ptr = c_loc(field%x)

  end function neko_api_user_cb_field_by_index

  !> Check if the user callback field at a given index has a given name
  !! @param field_idx Field index in the field list
  !! @param field_name Field name to compare against
  function neko_api_user_cb_field_name_at_index(field_idx, field_name) &
       result(same_name) bind(c, name='neko_cb_field_name_at_index')
    integer, intent(in) :: field_idx
    character(kind=c_char), dimension(*), intent(in) :: field_name
    character(len=8192) :: name
    type(field_t), pointer :: f1, f2
    type(c_ptr) :: field_ptr
    integer :: len
    logical(c_bool) :: same_name

    len = 0
    do
       if (field_name(len+1) .eq. C_NULL_CHAR) exit
       len = len + 1
       name(len:len) = field_name(len)
    end do

    f1 => neko_api_user_cb_get_field(field_idx)
    f2 => neko_api_user_cb_get_field(trim(name(1:len)))

    same_name = trim(f1%name) .eq. trim(f2%name)

  end function neko_api_user_cb_field_name_at_index


end module neko_api
