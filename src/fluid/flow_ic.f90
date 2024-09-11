! Copyright (c) 2021, The Neko Authors
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
!> Initial flow condition
module flow_ic
  use num_types, only : rp
  use logger, only: neko_log, LOG_SIZE
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use flow_profile, only : blasius_profile, blasius_linear, blasius_cubic, &
    blasius_quadratic, blasius_quartic, blasius_sin
  use device, only: device_memcpy, HOST_TO_DEVICE
  use field, only : field_t
  use utils, only : neko_error, filename_suffix, filename_chsuffix, &
       neko_warning, NEKO_FNAME_LEN
  use coefs, only : coef_t
  use math, only : col2, cfill, cfill_mask, copy
  use device_math, only : device_col2, device_cfill, device_cfill_mask
  use user_intf, only : useric
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use fld_file_data, only: fld_file_data_t
  use fld_file, only: fld_file_t
  use checkpoint, only: chkp_t
  use file, only: file_t
  use global_interpolation, only: global_interpolation_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GLL
  implicit none
  private

  interface set_flow_ic
     module procedure set_flow_ic_int, set_flow_ic_usr
  end interface set_flow_ic

  public :: set_flow_ic

contains

  !> Set initial flow condition (builtin)
  subroutine set_flow_ic_int(u, v, w, p, coef, gs, type, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file), intent(inout) :: params
    real(kind=rp) :: delta, tol
    real(kind=rp), allocatable :: uinf(:)
    real(kind=rp), allocatable :: zone_value(:)
    character(len=:), allocatable :: read_str
    character(len=80) :: suffix
    character(len=NEKO_FNAME_LEN) :: fname
    integer :: sample_idx, fpos, sample_mesh_idx, default_sample_mesh_idx
    logical :: found_previous_mesh, interpolate
    character(len=LOG_SIZE) :: log_buf

    !
    ! Uniform (Uinf, Vinf, Winf)
    !
    if (trim(type) .eq. 'uniform') then

       call json_get(params, 'case.fluid.initial_condition.value', uinf)
       call set_flow_ic_uniform(u, v, w, uinf)
       write (log_buf, '(A,"[",2(F10.6,","),F10.6,"]")') "Value: ", uinf(1), &
            uinf(2), uinf(3)
       call neko_log%message(log_buf)

    !
    ! Blasius boundary layer
    !
    else if (trim(type) .eq. 'blasius') then

       call json_get(params, 'case.fluid.blasius.delta', delta)
       call json_get(params, 'case.fluid.blasius.approximation', &
                     read_str)
       call json_get(params, 'case.fluid.blasius.freestream_velocity', uinf)

       write (log_buf, '(A,F10.6)') "delta       : ", delta
       call neko_log%message(log_buf)
       call neko_log%message("Approximation: " // trim(read_str))
       write (log_buf, '(A,"[",2(F10.6,","),F10.6,"]")') "Value: ", uinf(1), &
            uinf(2), uinf(3)
       call neko_log%message(log_buf)

       call set_flow_ic_blasius(u, v, w, delta, uinf, read_str)

    !
    ! Point zone initial condition
    !
    else if (trim(type) .eq. 'point_zone') then

       call json_get(params, 'case.fluid.initial_condition.base_value', uinf)
       call json_get(params, 'case.fluid.initial_condition.zone_name', &
                     read_str)
       call json_get(params, 'case.fluid.initial_condition.zone_value', &
            zone_value)

       write (log_buf, '(A,F10.6)') "Base value: ", uinf
       call neko_log%message(log_buf)
       call neko_log%message("Zone name : " // trim(read_str))
       write (log_buf, '(A,"[",2(F10.6,","),F10.6,"]")') "Value: ", &
            zone_value(1), zone_value(2), zone_value(3)
       call neko_log%message(log_buf)

       call set_flow_ic_point_zone(u, v, w, uinf, read_str, zone_value)

    !
    ! Field initial condition (from fld file)
    !
    else if (trim(type) .eq. 'field') then

       call json_get(params, 'case.fluid.initial_condition.file_name', &
            read_str)
       fname = trim(read_str)
       call filename_suffix(fname, suffix)
       call neko_log%message("File name: " // trim(fname))

       !
       ! Get the index of the file to sample
       !
       call json_get_or_default(params, &
            'case.fluid.initial_condition.sample_index', sample_idx, -1)

       ! In case we interpolate, the default index of the file in which to
       ! look for the coordinates
       default_sample_mesh_idx = 0

       if (sample_idx .ne. -1) then
          write (log_buf, '(A,I5)') "Sample index: ", sample_idx
          call neko_log%message(log_buf)
       end if

       ! If it's not fld assume it's either .nek5000 or .f00*
       if (trim(suffix) .ne. "fld") then

          ! Check if the suffix is of type "f000*", and if so extract
          ! the index e.g. "f00035" --> 35
          ! NOTE: overwrites whatever is in sampled_index
          fpos = scan(suffix, 'f')
          if (fpos .eq. 1) then
             if (sample_idx .ne. -1) &
                  call neko_warning("Overwriting sample index")
             read (suffix(2:), "(I5.5)") sample_idx

             ! In case we interpolate we assume the mesh is also in this file
             default_sample_mesh_idx = sample_idx
          end if

          call filename_chsuffix(fname, fname, 'fld')
       end if

       !
       ! Check if we want to interpolate from the field file, default is no
       !
       call json_get_or_default(params, &
            'case.fluid.initial_condition.interpolate', interpolate, &
            .false.)

       !
       ! Get the tolerance for potential interpolation, defaults to 1e-6
       !
       call json_get_or_default(params, &
            'case.fluid.initial_condition.tolerance', tol, 0.000001_rp)

       ! If no mesh file is provided, try to get the index of the
       ! file that contains the mesh.
       call json_get_or_default(params, &
            'case.fluid.initial_condition.sample_mesh_index', &
            sample_mesh_idx, default_sample_mesh_idx)

       if (interpolate) then
          call neko_log%message("Interpolation   : yes")
          write (log_buf, '(A,E15.7)') "Tolerance       : ", tol
          call neko_log%message(log_buf)
          write (log_buf, '(A,I5.5)') "Coordinates at sample index : ", &
               sample_mesh_idx
          call neko_log%message(log_buf)
       end if

       call set_flow_ic_fld(u, v, w, p, fname, sample_idx, interpolate, &
            tolerance = tol, sample_mesh_idx = sample_mesh_idx)

    else
       call neko_error('Invalid initial condition')
    end if

    call set_flow_ic_common(u, v, w, p, coef, gs)

  end subroutine set_flow_ic_int

  !> Set intial flow condition (user defined)
  subroutine set_flow_ic_usr(u, v, w, p, coef, gs, usr_ic, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(useric) :: usr_ic
    type(json_file), intent(inout) :: params

    call usr_ic(u, v, w, p, params)

    call set_flow_ic_common(u, v, w, p, coef, gs)

  end subroutine set_flow_ic_usr

  subroutine set_flow_ic_common(u, v, w, p, coef, gs)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(v%x, v%x_d, n, &
                          HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(w%x, w%x_d, n, &
                          HOST_TO_DEVICE, sync = .false.)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(u%x, u%dof%size(), GS_OP_ADD)
    call gs%op(v%x, v%dof%size(), GS_OP_ADD)
    call gs%op(w%x, w%dof%size(), GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(u%x_d, coef%mult_d, u%dof%size())
       call device_col2(v%x_d, coef%mult_d, v%dof%size())
       call device_col2(w%x_d, coef%mult_d, w%dof%size())
    else
       call col2(u%x, coef%mult, u%dof%size())
       call col2(v%x, coef%mult, v%dof%size())
       call col2(w%x, coef%mult, w%dof%size())
    end if

  end subroutine set_flow_ic_common

  !> Uniform initial condition
  subroutine set_flow_ic_uniform(u, v, w, uinf)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in) :: uinf(3)
    integer :: n
    u = uinf(1)
    v = uinf(2)
    w = uinf(3)
    n = u%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(u%x, uinf(1), n)
       call cfill(v%x, uinf(2), n)
       call cfill(w%x, uinf(3), n)
    end if

  end subroutine set_flow_ic_uniform

  !> Set a Blasius profile as initial condition
  !! @note currently limited to axis aligned flow
  subroutine set_flow_ic_blasius(u, v, w, delta, uinf, type)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in) :: delta
    real(kind=rp), intent(in) :: uinf(3)
    character(len=*), intent(in) :: type
    procedure(blasius_profile), pointer :: bla => null()
    integer :: i

    select case (trim(type))
    case ('linear')
       bla => blasius_linear
    case ('quadratic')
       bla => blasius_quadratic
    case ('cubic')
       bla => blasius_cubic
    case ('quartic')
       bla => blasius_quartic
    case ('sin')
       bla => blasius_sin
    case default
       call neko_error('Invalid Blasius approximation')
    end select

    if ((uinf(1) .gt. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
       .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = bla(u%dof%z(i,1,1,1), delta, uinf(1))
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .gt. 0.0_rp) &
            .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = bla(u%dof%x(i,1,1,1), delta, uinf(2))
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
            .and. (uinf(3) .gt. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = bla(u%dof%y(i,1,1,1), delta, uinf(3))
       end do
    end if

  end subroutine set_flow_ic_blasius

  !> Set the initial condition of the flow based on a point zone.
  !! @details The initial condition is set to the base value and then the
  !! zone is filled with the zone value.
  !! @param u The x-component of the velocity field.
  !! @param v The y-component of the velocity field.
  !! @param w The z-component of the velocity field.
  !! @param base_value The base value of the initial condition.
  !! @param zone_name The name of the point zone.
  !! @param zone_value The value of the point zone.
  subroutine set_flow_ic_point_zone(u, v, w, base_value, zone_name, zone_value)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in), dimension(3) :: base_value
    character(len=*), intent(in) :: zone_name
    real(kind=rp), intent(in) :: zone_value(:)

    ! Internal variables
    class(point_zone_t), pointer :: zone
    integer :: size

    call set_flow_ic_uniform(u, v, w, base_value)
    size = u%dof%size()

    zone => neko_point_zone_registry%get_point_zone(trim(zone_name))

    call cfill_mask(u%x, zone_value(1), size, zone%mask, zone%size)
    call cfill_mask(v%x, zone_value(2), size, zone%mask, zone%size)
    call cfill_mask(w%x, zone_value(3), size, zone%mask, zone%size)

  end subroutine set_flow_ic_point_zone

  !> Set the initial condition of the flow based on a point zone.
  !! @details The initial condition is set to the base value and then the
  !! zone is filled with the zone value.
  !! @param u The x-component of the velocity field.
  !! @param v The y-component of the velocity field.
  !! @param w The z-component of the velocity field.
  !! @param p The pressure field.
  !! @param file_name The name of the "fld" file series.
  !! @param sample_idx index of the field file .f000* to read, default is
  !! -1..
  subroutine set_flow_ic_fld(u, v, w, p, file_name, sample_idx, &
       interpolate, tolerance, sample_mesh_idx)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: sample_idx
    logical, intent(in) :: interpolate
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: sample_mesh_idx

    integer :: last_index
    type(fld_file_data_t) :: fld_data
    type(file_t) :: f

    ! ---- For the mesh to mesh interpolation
    type(global_interpolation_t) :: global_interp
    ! -----

    ! ---- For space to space interpolation
    type(space_t) :: prev_Xh
    type(interpolator_t) :: space_interp
    ! ----

    call fld_data%init
    f = file_t(trim(file_name))

    ! If no sample is specified, we read as a classic fld series
    ! with the mesh located at index 0
    if (sample_idx .eq. -1) then

       ! Read one time to get all information about the fld series
       ! and read the first file in the series.
       call f%read(fld_data)

       last_index = fld_data%meta_nsamples + fld_data%meta_start_counter - 1

       if (interpolate) then
          ! Read the coordinates if they are neither at index 0 or at the
          ! last index
          if (sample_mesh_idx .ne. 0 .and. &
               sample_mesh_idx .ne. last_index) then
             call f%set_counter(sample_mesh_idx)
             call f%read(fld_data)
          end if
       end if ! interpolate

       ! If there is more than one file, read the last one in the series.
       if (fld_data%meta_nsamples .gt. 1) then
          call f%set_counter(last_index)
          call f%read(fld_data)
       end if

    else ! In this case we have a specific sample to read in the fld series

       if (interpolate) then
          ! Read the mesh coordinates if they are not in our fld file
          if (sample_mesh_idx .ne. sample_idx) then
             call f%set_counter(sample_mesh_idx)
             call f%read(fld_data)
          end if

          ! Read the actual fld file
          call f%read(fld_data)
          f = file_t(trim(file_name))
       end if ! interpolate

       call f%set_counter(sample_idx)
       call f%read(fld_data)

    end if ! sample_idx .eq. -1

    !
    ! Check if the data in the fld file matches the current case.
    ! Note that this is a safeguard and there are corner cases where
    ! two different meshes have the same dimension and same # of elements
    ! but this should be enough to cover the most obvious cases.
    !
    if ( fld_data%glb_nelv .ne. u%msh%glb_nelv .or. &
         fld_data%gdim .ne. u%msh%gdim .and. &
         .not. interpolate) then
       call neko_error("The fld file must match the current mesh! &
&Use 'interpolate': 'true' to enable interpolation.")
    else if (interpolate) then
       call neko_log%warning("You have activated interpolation but you might &
&still be using the same mesh.")
    end if

    ! Mesh interpolation if specified
    if (interpolate) then

       ! Issue a warning if the mesh is in single precision
       select type (ft => f%file_type)
       type is (fld_file_t)
          if (.not. ft%dp_precision) then
             call neko_warning("The coordinates read from the field file are &
&in single precision.")
             call neko_log%message("It is recommended to use a mesh in double &
                  &precision for better interpolation results.")
             call neko_log%message("Reduce the tolerance if the interpolation &
                  &does not work.")
          end if
       class default
       end select

       global_interp = fld_data%generate_interpolator(u%dof, u%msh, &
            tolerance)

       ! Evaluate velocities and pressure
       call global_interp%evaluate(u%x, fld_data%u%x)
       call global_interp%evaluate(v%x, fld_data%v%x)
       call global_interp%evaluate(w%x, fld_data%w%x)
       call global_interp%evaluate(p%x, fld_data%p%x)

       call global_interp%free

    else ! No interpolation, but potentially just from different spaces

       ! Build a space_t object from the data in the fld file
       call prev_Xh%init(GLL, fld_data%lx, fld_data%ly, fld_data%lz)
       call space_interp%init(u%Xh, prev_Xh)

       call space_interp%map_host(u%x, fld_data%u%x, fld_data%nelv, u%Xh)
       call space_interp%map_host(v%x, fld_data%v%x, fld_data%nelv, u%Xh)
       call space_interp%map_host(w%x, fld_data%w%x, fld_data%nelv, u%Xh)
       call space_interp%map_host(p%x, fld_data%p%x, fld_data%nelv, u%Xh)

       call space_interp%free

    end if

    ! NOTE: we do not copy (u,v,w) on the GPU since `set_flow_ic_common` does
    ! the copy for us
    if (NEKO_BCKND_DEVICE .eq. 1) call device_memcpy(p%x, p%x_d, p%dof%size(), &
         HOST_TO_DEVICE, sync = .false.)

    call fld_data%free

  end subroutine set_flow_ic_fld

end module flow_ic
