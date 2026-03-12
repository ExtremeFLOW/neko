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
!> Implements the `force_torque_t` type.

module force_torque
  use num_types, only : rp, dp, sp
  use time_based_controller, only : time_based_controller_t
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use coefs, only : coef_t
  use operators, only : strain_rate
  use vector, only : vector_t
  use dirichlet, only : dirichlet_t
  use drag_torque, only : calc_force_array, device_calc_force_array, &
       setup_normals
  use logger, only : LOG_SIZE, neko_log
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : masked_gather_copy_0, cadd, glsum, vcross
  use device_math, only : device_masked_gather_copy_0, device_cadd, &
       device_glsum, device_vcross
  use mpi_f08, only : MPI_INTEGER, MPI_SUM, MPI_Allreduce
  use comm, only : NEKO_COMM
  use device, only : device_memcpy, HOST_TO_DEVICE
  use ale_manager, only : neko_ale
  use ale_rigid_kinematics, only : pivot_state_t
  use utils, only : neko_error

  implicit none
  private

  !> A simulation component that computes the force and torque on a given
  !! boundary zone.
  type, public, extends(simulation_component_t) :: force_torque_t
     !> X velocity component.
     type(field_t), pointer :: u => null()
     !> Y velocity component.
     type(field_t), pointer :: v => null()
     !> Z velocity component.
     type(field_t), pointer :: w => null()
     !> Pressure.
     type(field_t), pointer :: p => null()
     !> Total dynamic viscosity.
     type(field_t), pointer :: mu => null()

     ! Masked working arrays
     type(vector_t) :: n1, n2, n3
     type(vector_t) :: r1, r2, r3
     type(vector_t) :: force1, force2, force3
     type(vector_t) :: force4, force5, force6
     type(vector_t) :: pmsk
     type(vector_t) :: mu_msk
     type(vector_t) :: s11msk, s22msk, s33msk, s12msk, s13msk, s23msk
     real(kind=rp) :: center(3) = 0.0_rp
     real(kind=rp) :: scale
     integer :: zone_id
     character(len=20) :: zone_name
     type(coef_t), pointer :: coef => null()
     type(dirichlet_t) :: bc
     character(len=80) :: print_format
     ! Pointer to the live pivot state inside ale_manager
     type(pivot_state_t), pointer :: pivot_link => null()
     logical :: moving_center = .false.
     logical :: update_normals = .false.
     character(len=64) :: linked_body_name = 'NOT_LINKED'
     ! Stores the Time=0 offset from the initial pivot
     real(kind=rp) :: local_offset(3) = 0.0_rp
     ! Current Pivot Position
     real(kind=rp), pointer :: body_P(:) => null()
     ! Current Rotation Matrix
     real(kind=rp), pointer :: body_R(:,:) => null()

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => force_torque_init_from_json
     !> Common part of constructors.
     procedure, private, pass(this) :: init_common => force_torque_init_common
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          force_torque_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          force_torque_init_from_controllers_properties
     !> Destructor.
     procedure, pass(this) :: free => force_torque_free
     !> Compute the force_torque field.
     procedure, pass(this) :: compute_ => force_torque_compute
     !> Routine to setup ALE links
     procedure, private, pass(this) :: ale_link => setup_ale_link
  end type force_torque_t

contains

  !> Constructor from json.
  subroutine force_torque_init_from_json(this, json, case)
    class(force_torque_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    integer :: zone_id
    real(kind=rp), allocatable :: center(:)
    character(len=:), allocatable :: zone_name, fluid_name, center_type
    character(len=:), allocatable :: effective_center_type
    character(len=:), allocatable :: name
    real(kind=rp) :: scale
    logical :: long_print

    call json_get_or_default(json, "name", name, "force_torque")
    call this%init_base(json, case)

    call json_get_or_default(json, 'fluid_name', fluid_name, 'fluid')
    call json_get_or_lookup(json, 'zone_id', zone_id)
    call json_get_or_default(json, 'zone_name', zone_name, ' ')
    call json_get_or_lookup_or_default(json, 'scale', scale, 1.0_rp)
    call json_get_or_default(json, 'long_print', long_print, .false.)
    call json_get_or_lookup(json, 'center', center)
    call json_get_or_default(json, 'center_type', center_type, 'fixed')
    if (trim(center_type) /= 'fixed' .and. &
         trim(center_type) /= 'pivot' .and. &
         trim(center_type) /= 'body_attached') then
       call neko_error("force_torque: center_type must be 'fixed'"// &
            ", 'pivot', or 'body_attached'.")
    end if

    if (allocated(center)) this%center = center

    call this%init_common(name, fluid_name, zone_id, zone_name, this%center, &
         scale, case%fluid%c_xh, long_print, center_type=center_type)
  end subroutine force_torque_init_from_json

  !> Constructor from components, passing controllers.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param fluid_name The name of the fluid solver.
  !! @param zone_id The id of the boundary zone.
  !! @param zone_name The name of the boundary zone, to use in the log.
  !! @param center The center of the torque calculation.
  !! @param scale Normalization factor.
  !! @param coef The SEM coefficients.
  !! @param long_print If true, use a more precise print format.
  subroutine force_torque_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       fluid_name, zone_id, zone_name, center, scale, coef, long_print)
    class(force_torque_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*), intent(in) :: fluid_name
    character(len=*), intent(in) :: zone_name
    integer, intent(in) :: zone_id
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: scale
    type(coef_t), target, intent(in) :: coef
    logical, intent(in) :: long_print

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%init_common(name, fluid_name, zone_id, zone_name, center, scale, &
         coef, long_print)

  end subroutine force_torque_init_from_controllers

  !> Constructor from components, passing properties to the
  !! time_based_controller` components in the base type.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param fluid_name The name of the fluid solver.
  !! @param zone_id The id of the boundary zone.
  !! @param zone_name The name of the boundary zone, to use in the log.
  !! @param center The center of the torque calculation.
  !! @param scale Normalization factor.
  !! @param coef The SEM coefficients.
  !! @param long_print If true, use a more precise print format.
  subroutine force_torque_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, fluid_name, zone_name, &
       zone_id, center, scale, coef, long_print)
    class(force_torque_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*), intent(in) :: fluid_name
    character(len=*), intent(in) :: zone_name
    integer, intent(in) :: zone_id
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: scale
    type(coef_t), target, intent(in) :: coef
    logical, intent(in) :: long_print

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%init_common(name, fluid_name, zone_id, zone_name, center, scale, &
         coef, long_print)

  end subroutine force_torque_init_from_controllers_properties

  !> Common part of constructors.
  !! @param name The unique name of the simcomp.
  !! @param fluid_name The name of the fluid solver.
  !! @param zone_id The id of the boundary zone.
  !! @param zone_name The name of the boundary zone, to use in the log.
  !! @param center The center of the torque calculation.
  !! @param center_type The type of center used.
  !! @param scale Normalization factor.
  !! @param coef The SEM coefficients.
  !! @param long_print If true, use a more precise print format.
  subroutine force_torque_init_common(this, name, fluid_name, zone_id, &
       zone_name, center, scale, coef, long_print, center_type)
    class(force_torque_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: center(3)
    real(kind=rp), intent(in) :: scale
    character(len=*), intent(in) :: fluid_name
    character(len=*), intent(in) :: zone_name
    integer, intent(in) :: zone_id
    type(coef_t), target, intent(in) :: coef
    logical, intent(in) :: long_print
    character(len=*), intent(in), optional :: center_type
    character(len=:), allocatable :: ctype_str
    integer :: n_pts, glb_n_pts, ierr
    real(kind=rp) :: avg_r(3)
    character(len=1000) :: log_buf
    this%name = name
    this%coef => coef
    this%zone_id = zone_id
    this%scale = scale
    this%zone_name = zone_name

    if (present(center_type)) then
       ctype_str = center_type
    else
       ctype_str = 'fixed' ! Default behavior
    end if

    call this%ale_link(zone_id, ctype_str, center)

    if (ctype_str /= 'fixed' .and. .not. this%moving_center) then
       ctype_str = 'fixed (reverted from ' // ctype_str // ')'
    end if

    ! Set fixed center if not linked to an ALE body
    if (.not. this%moving_center) then
       this%center = center
    end if


    if (long_print) then
       this%print_format = '(I7,E20.10,E20.10,E20.10,E20.10,A)'
    else
       this%print_format = '(I7,E13.5,E13.5,E13.5,E13.5,A)'
    end if

    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")
    this%p => neko_registry%get_field_by_name("p")
    this%mu => neko_registry%get_field_by_name(fluid_name // '_mu_tot')


    call this%bc%init_base(this%coef)
    call this%bc%mark_zone(this%case%msh%labeled_zones(this%zone_id))
    call this%bc%finalize()
    n_pts = this%bc%msk(0)
    if (n_pts .gt. 0) then
       call this%n1%init(n_pts)
       call this%n2%init(n_pts)
       call this%n3%init(n_pts)
       call this%r1%init(n_pts)
       call this%r2%init(n_pts)
       call this%r3%init(n_pts)
       call this%force1%init(n_pts)
       call this%force2%init(n_pts)
       call this%force3%init(n_pts)
       call this%force4%init(n_pts)
       call this%force5%init(n_pts)
       call this%force6%init(n_pts)
       call this%s11msk%init(n_pts)
       call this%s22msk%init(n_pts)
       call this%s33msk%init(n_pts)
       call this%s12msk%init(n_pts)
       call this%s13msk%init(n_pts)
       call this%s23msk%init(n_pts)
       call this%pmsk%init(n_pts)
       call this%mu_msk%init(n_pts)
    end if

    call setup_normals(this%coef, this%bc%msk, this%bc%facet, &
         this%n1%x, this%n2%x, this%n3%x, n_pts)
    call masked_gather_copy_0(this%r1%x, this%coef%dof%x, this%bc%msk, &
         this%u%size(), n_pts)
    call masked_gather_copy_0(this%r2%x, this%coef%dof%y, this%bc%msk, &
         this%u%size(), n_pts)
    call masked_gather_copy_0(this%r3%x, this%coef%dof%z, this%bc%msk, &
         this%u%size(), n_pts)

    call MPI_Allreduce(n_pts, glb_n_pts, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    ! Calculatig avg pos here
    avg_r(1) = glsum(this%r1%x, n_pts)/glb_n_pts
    avg_r(2) = glsum(this%r2%x, n_pts)/glb_n_pts
    avg_r(3) = glsum(this%r3%x, n_pts)/glb_n_pts
    ! Print some information
    call neko_log%section('Force/torque calculation')
    write(log_buf, '(A,I4,A,A)') 'Zone ', zone_id, '  ', trim(zone_name)
    call neko_log%message(log_buf)

    write(log_buf, '(A,I6, I6)') 'Global number of GLL points in zone: ', &
         glb_n_pts
    call neko_log%message(log_buf)

    write(log_buf, '(A,A)') 'Center Type: ', trim(ctype_str)
    call neko_log%message(log_buf)

    if (trim(ctype_str) == 'pivot' .or. &
    trim(ctype_str) == 'body_attached') then
       write(log_buf, '(A,A)') 'Linked to ALE Body movement: ', &
            trim(this%linked_body_name)
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7,E15.7,E15.7)') &
            'Initial center for torque calculation: ', this%center
       call neko_log%message(log_buf)
    else
       write(log_buf, '(A,E15.7,E15.7,E15.7)') &
            'Fixed center for torque calculation: ', this%center
       call neko_log%message(log_buf)
    end if

    write(log_buf, '(A,E15.7,E15.7,E15.7)') &
         'Average of zone''s coordinates: ', avg_r
    call neko_log%message(log_buf)

    write(log_buf, '(A,E15.7)') 'Scale: ', scale
    call neko_log%message(log_buf)
    call neko_log%end_section()


    call cadd(this%r1%x, -this%center(1), n_pts)
    call cadd(this%r2%x, -this%center(2), n_pts)
    call cadd(this%r3%x, -this%center(3), n_pts)
    if (NEKO_BCKND_DEVICE .eq. 1 .and. n_pts .gt. 0) then
       call device_memcpy(this%n1%x, this%n1%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)
       call device_memcpy(this%n2%x, this%n2%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)
       call device_memcpy(this%n3%x, this%n3%x_d, n_pts, HOST_TO_DEVICE, &
            .true.)
       call device_memcpy(this%r1%x, this%r1%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)
       call device_memcpy(this%r2%x, this%r2%x_d, n_pts, HOST_TO_DEVICE, &
            .false.)
       call device_memcpy(this%r3%x, this%r3%x_d, n_pts, HOST_TO_DEVICE, &
            .true.)
    end if

  end subroutine force_torque_init_common

  !> Destructor.
  subroutine force_torque_free(this)
    class(force_torque_t), intent(inout) :: this
    call this%free_base()

    call this%n1%free()
    call this%n2%free()
    call this%n3%free()

    call this%r1%free()
    call this%r2%free()
    call this%r3%free()

    call this%force1%free()
    call this%force2%free()
    call this%force3%free()

    call this%force4%free()
    call this%force5%free()
    call this%force6%free()

    call this%pmsk%free()
    call this%mu_msk%free()
    call this%s11msk%free()
    call this%s22msk%free()
    call this%s33msk%free()
    call this%s12msk%free()
    call this%s13msk%free()
    call this%s23msk%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%coef)
    nullify(this%mu)
    nullify(this%pivot_link)
  end subroutine force_torque_free

  !> Compute the force_torque field.
  !! @param time The current time state.
  subroutine force_torque_compute(this, time)
    class(force_torque_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    real(kind=rp) :: dgtq(12) = 0.0_rp
    integer :: n_pts, temp_indices(6)
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    character(len=1000) :: log_buf
    real(kind=rp) :: rot_offset(3)
    n_pts = this%bc%msk(0)


    ! body_attached
    if (this%moving_center .and. associated(this%body_P)) then
       if (associated(this%body_R)) then
          ! R * Offset
          rot_offset(1) = this%body_R(1,1)*this%local_offset(1) + &
                          this%body_R(1,2)*this%local_offset(2) + &
                          this%body_R(1,3)*this%local_offset(3)
          rot_offset(2) = this%body_R(2,1)*this%local_offset(1) + &
                          this%body_R(2,2)*this%local_offset(2) + &
                          this%body_R(2,3)*this%local_offset(3)
          rot_offset(3) = this%body_R(3,1)*this%local_offset(1) + &
                          this%body_R(3,2)*this%local_offset(2) + &
                          this%body_R(3,3)*this%local_offset(3)

          this%center = this%body_P + rot_offset
       end if
    end if

    if (this%update_normals) then

       call setup_normals(this%coef, this%bc%msk, this%bc%facet, &
            this%n1%x, this%n2%x, this%n3%x, n_pts)

       call masked_gather_copy_0(this%r1%x, this%coef%dof%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%r2%x, this%coef%dof%y, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%r3%x, this%coef%dof%z, this%bc%msk, &
            this%u%size(), n_pts)

       call cadd(this%r1%x, -this%center(1), n_pts)
       call cadd(this%r2%x, -this%center(2), n_pts)
       call cadd(this%r3%x, -this%center(3), n_pts)

       if (NEKO_BCKND_DEVICE .eq. 1 .and. n_pts .gt. 0) then
          call device_memcpy(this%n1%x, this%n1%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n2%x, this%n2%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n3%x, this%n3%x_d, n_pts, &
               HOST_TO_DEVICE, .true.)
          call device_memcpy(this%r1%x, this%r1%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%r2%x, this%r2%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%r3%x, this%r3%x_d, n_pts, &
               HOST_TO_DEVICE, .true.)
       end if
    end if

    call neko_scratch_registry%request_field(s11, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(s12, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(s13, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(s22, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(s23, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(s33, temp_indices(6), .false.)

    call strain_rate(s11%x, s22%x, s33%x, s12%x, &
         s13%x, s23%x, this%u, this%v, this%w, this%coef)

    ! On the CPU we can actually just use the original subroutines...
    if (NEKO_BCKND_DEVICE .eq. 0) then
       call masked_gather_copy_0(this%s11msk%x, s11%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%s22msk%x, s22%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%s33msk%x, s33%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%s12msk%x, s12%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%s13msk%x, s13%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%s23msk%x, s23%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%pmsk%x, this%p%x, this%bc%msk, &
            this%u%size(), n_pts)
       call masked_gather_copy_0(this%mu_msk%x, this%mu%x, this%bc%msk, &
            this%u%size(), n_pts)
       call calc_force_array(this%force1%x, this%force2%x, this%force3%x, &
            this%force4%x, this%force5%x, this%force6%x, &
            this%s11msk%x, &
            this%s22msk%x, &
            this%s33msk%x, &
            this%s12msk%x, &
            this%s13msk%x, &
            this%s23msk%x, &
            this%pmsk%x, &
            this%n1%x, &
            this%n2%x, &
            this%n3%x, &
            this%mu_msk%x, &
            n_pts)
       dgtq(1) = glsum(this%force1%x, n_pts)
       dgtq(2) = glsum(this%force2%x, n_pts)
       dgtq(3) = glsum(this%force3%x, n_pts)
       dgtq(4) = glsum(this%force4%x, n_pts)
       dgtq(5) = glsum(this%force5%x, n_pts)
       dgtq(6) = glsum(this%force6%x, n_pts)
       call vcross(this%s11msk%x, this%s22msk%x, this%s33msk%x, &
            this%r1%x, this%r2%x, this%r3%x, &
            this%force1%x, this%force2%x, this%force3%x, n_pts)

       dgtq(7) = glsum(this%s11msk%x, n_pts)
       dgtq(8) = glsum(this%s22msk%x, n_pts)
       dgtq(9) = glsum(this%s33msk%x, n_pts)
       call vcross(this%s11msk%x, this%s22msk%x, this%s33msk%x, &
            this%r1%x, this%r2%x, this%r3%x, &
            this%force4%x, this%force5%x, this%force6%x, n_pts)
       dgtq(10) = glsum(this%s11msk%x, n_pts)
       dgtq(11) = glsum(this%s22msk%x, n_pts)
       dgtq(12) = glsum(this%s33msk%x, n_pts)
    else
       if (n_pts .gt. 0) then
          call device_masked_gather_copy_0(this%s11msk%x_d, s11%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%s22msk%x_d, s22%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%s33msk%x_d, s33%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%s12msk%x_d, s12%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%s13msk%x_d, s13%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%s23msk%x_d, s23%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%pmsk%x_d, this%p%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)
          call device_masked_gather_copy_0(this%mu_msk%x_d, this%mu%x_d, &
               this%bc%msk_d, this%u%size(), n_pts)

          call device_calc_force_array(this%force1%x_d, this%force2%x_d, &
               this%force3%x_d, &
               this%force4%x_d, &
               this%force5%x_d, &
               this%force6%x_d, &
               this%s11msk%x_d, &
               this%s22msk%x_d, &
               this%s33msk%x_d, &
               this%s12msk%x_d, &
               this%s13msk%x_d, &
               this%s23msk%x_d, &
               this%pmsk%x_d, &
               this%n1%x_d, &
               this%n2%x_d, &
               this%n3%x_d, &
               this%mu%x_d, &
               n_pts)
          ! Overwriting masked s11, s22, s33 as they are no longer needed
          call device_vcross(this%s11msk%x_d, this%s22msk%x_d, &
               this%s33msk%x_d, &
               this%r1%x_d, this%r2%x_d, this%r3%x_d, &
               this%force1%x_d, this%force2%x_d, &
               this%force3%x_d, n_pts)
          call device_vcross(this%s12msk%x_d, this%s13msk%x_d, this%s23msk%x_d,&
               this%r1%x_d, this%r2%x_d, this%r3%x_d, &
               this%force4%x_d, this%force5%x_d, this%force6%x_d, n_pts)
       end if
       dgtq(1) = device_glsum(this%force1%x_d, n_pts)
       dgtq(2) = device_glsum(this%force2%x_d, n_pts)
       dgtq(3) = device_glsum(this%force3%x_d, n_pts)
       dgtq(4) = device_glsum(this%force4%x_d, n_pts)
       dgtq(5) = device_glsum(this%force5%x_d, n_pts)
       dgtq(6) = device_glsum(this%force6%x_d, n_pts)
       dgtq(7) = device_glsum(this%s11msk%x_d, n_pts)
       dgtq(8) = device_glsum(this%s22msk%x_d, n_pts)
       dgtq(9) = device_glsum(this%s33msk%x_d, n_pts)
       dgtq(10) = device_glsum(this%s12msk%x_d, n_pts)
       dgtq(11) = device_glsum(this%s13msk%x_d, n_pts)
       dgtq(12) = device_glsum(this%s23msk%x_d, n_pts)
    end if
    dgtq = this%scale*dgtq
    write(log_buf, '(A, I4, A, A)') 'Force and torque on zone ', &
         this%zone_id, '  ', this%zone_name
    call neko_log%message(log_buf)
    write(log_buf, '(A)') &
         'Time step, time, total force/torque, pressure, viscous, direction'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(1) + dgtq(4), dgtq(1), dgtq(4), ', forcex'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(2) + dgtq(5), dgtq(2), dgtq(5), ', forcey'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(3) + dgtq(6), dgtq(3), dgtq(6), ', forcez'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(7) + dgtq(10), dgtq(7), dgtq(10), ', torquex'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(8) + dgtq(11), dgtq(8), dgtq(11), ', torquey'
    call neko_log%message(log_buf)
    write(log_buf, this%print_format) &
         time%tstep, time%t, dgtq(9) + dgtq(12), dgtq(9), dgtq(12), ', torquez'
    call neko_log%message(log_buf)
    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine force_torque_compute

  !> Routine to configure ALE connectivity for force/torque module
  subroutine setup_ale_link(this, zone_id, center_type, center_in)
    class(force_torque_t), intent(inout) :: this
    integer, intent(in) :: zone_id
    character(len=*), intent(in) :: center_type
    real(kind=rp), intent(in) :: center_in(3)

    character(len=512) :: log_buf
    integer :: i, j, nbodies, nindices, body_id
    logical :: body_found
    logical :: ale_active = .false.

    this%moving_center = .false.
    this%update_normals = .false.
    this%local_offset = 0.0_rp
    nullify(this%body_P)
    nullify(this%body_R)

    ! Check Global Pointer for ALE Activity
    if (associated(neko_ale)) then
       if (neko_ale%active) then
          ale_active = .true.
          this%update_normals = .true.

          if (trim(center_type) == 'pivot' .or. &
               trim(center_type) == 'body_attached') then

             body_found = .false.
             nbodies = neko_ale%config%nbodies
             i = 1
             do while (i <= nbodies .and. .not. body_found)
                if (allocated(neko_ale%config%bodies(i)%zone_indices)) then
                   nindices = size(neko_ale%config%bodies(i)%zone_indices)
                   j = 1
                   do while (j <= nindices .and. .not. body_found)
                      if (neko_ale%config%bodies(i)%zone_indices(j) == zone_id) then

                         ! Body found
                         this%moving_center = .true.
                         this%pivot_link => neko_ale%ale_pivot(i)
                         this%linked_body_name = neko_ale%config%bodies(i)%name
                         body_id = neko_ale%config%bodies(i)%id
                         body_found = .true.

                         ! Point to Live Data
                         this%body_P => neko_ale%ale_pivot(i)%pos
                         this%body_R => neko_ale%body_rot_matrices(:, :, i)

                         ! Calculate local offset (Using Time=0 data)
                         if (trim(center_type) == 'pivot') then
                            ! Attached to center -> Offset is 0
                            this%local_offset = 0.0_rp
                            this%center = this%body_P ! Initial Position
                         else
                            ! Detached from pivot -->
                            ! offset = JSON_Input - Initial_Pivot
                            ! Note: we assume center_in is the Time=0 global coord
                            this%local_offset = center_in - &
                                 neko_ale%config%bodies(i)%rot_center
                            ! Set initial position for init_common
                            this%center = center_in
                         end if

                      end if
                      j = j + 1
                   end do
                end if
                i = i + 1
             end do

             if (.not. body_found) then
                call neko_log%message(' ')
                write(log_buf, '(A,I0,A)') 'Warning: Zone ', zone_id, &
                     ' requested "' // trim(center_type) // '" center, but is not'// &
                     ' registered as an ALE body.'
                call neko_log%message(log_buf)
                call neko_log%message('Reverting to FIXED center'// &
                     ' using JSON coordinates.')

                this%moving_center = .false.
             end if

          end if
       end if
    end if

    if ((.not. ale_active) .and. (trim(center_type) == 'pivot' .or. &
                                  trim(center_type) == 'body_attached')) then
       call neko_log%message(' ')
       write(log_buf, '(A,I0,A)') "Warning: Zone ", zone_id, &
            " requested '" // trim(center_type) // "' center, but ALE is not active."
       call neko_log%message(log_buf)
       call neko_log%message("pivot and body_attached work only for ALE simulations.")
       call neko_log%message("Reverting to 'fixed' center using JSON coordinates.")
    end if

  end subroutine setup_ale_link

end module force_torque
