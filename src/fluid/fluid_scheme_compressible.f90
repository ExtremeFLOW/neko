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
module fluid_scheme_compressible
  use field, only : field_t
  use field_math, only : field_col2, field_col3, &
       field_cmult2, field_cmult, field_addcol3, field_add2, &
       field_cfill

  use registry, only : neko_registry
  use bc, only : bc_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_module, only : json_file
  use num_types, only : rp
  use mesh, only : mesh_t
  use scratch_registry, only : neko_scratch_registry
  use space, only : GLL
  use euler_idp, only : euler_idp_config_t, euler_idp_diagnostics_t
  use user_intf, only : user_t, user_material_properties_intf, &
       dummy_user_material_properties
  use json_utils, only : json_get_or_default, json_get_or_lookup_or_default
  use mpi_f08
  use operators, only : cfl_compressible
  use device, only : device_memcpy, HOST_TO_DEVICE
  use compressible_ops_cpu, only : &
       compressible_ops_cpu_compute_max_wave_speed, &
       compressible_ops_cpu_compute_entropy, &
       compressible_ops_cpu_conserved_to_primitive, EULER_STATE_OK
  use compressible_ops_device, only : &
       compressible_ops_device_compute_max_wave_speed, &
       compressible_ops_device_compute_entropy
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  use time_state, only : time_state_t
  use logger, only : neko_log, LOG_SIZE
  use math, only : glsum
  use num_types, only : i8
  use comm, only : NEKO_COMM, pe_rank
  use utils, only : neko_error
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  public :: fluid_scheme_compressible_validate
  public :: fluid_scheme_compressible_compute_cfl

  !> Base type of compressible fluid formulations.
  type, public, abstract, extends(fluid_scheme_base_t) :: &
       fluid_scheme_compressible_t
     !> The momentum field
     type(field_t), pointer :: m_x => null() !< x-component of Momentum
     type(field_t), pointer :: m_y => null() !< y-component of Momentum
     type(field_t), pointer :: m_z => null() !< z-component of Momentum
     type(field_t), pointer :: E => null() !< Total energy
     type(field_t), pointer :: temperature => null() !< Temperature field
     type(field_t), pointer :: max_wave_speed => null() !< Maximum wave speed field
     type(field_t), pointer :: S => null() !< Entropy field
     type(field_t), pointer :: artificial_visc => null() !< Artificial viscosity field (without physical)
     type(field_t), pointer :: kappa => null() !< Thermal conductivity

     real(kind=rp) :: gamma

     !> Configuration of the opt-in Euler IDP path
     type(euler_idp_config_t) :: euler_idp
     !> Diagnostics produced by the Euler IDP path
     type(euler_idp_diagnostics_t) :: euler_idp_diagnostics

     !> Global number of GLL points for the fluid (not unique)
     integer(kind=i8) :: glb_n_points
     !> Global number of GLL points for the fluid (unique)
     integer(kind=i8) :: glb_unique_points

   contains
     !> Constructors
     procedure, pass(this) :: scheme_init => fluid_scheme_compressible_init
     !> Destructor for the base type
     procedure, pass(this) :: scheme_free => fluid_scheme_compressible_free

     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => fluid_scheme_compressible_validate
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl &
          => fluid_scheme_compressible_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: set_material_properties => &
          fluid_scheme_compressible_set_material_properties
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          fluid_scheme_compressible_update_material_properties
     !> Compute entropy field
     procedure, pass(this) :: compute_entropy => &
          fluid_scheme_compressible_compute_entropy
     !> Compute maximum wave speed
     procedure, pass(this) :: compute_max_wave_speed => &
          fluid_scheme_compressible_compute_max_wave_speed
     !> Log solver information
     procedure, pass(this) :: log_solver_info => &
          fluid_scheme_compressible_log_solver_info

  end type fluid_scheme_compressible_t

contains
  !> Initialize common data for compressible fluid scheme
  !> @param this The compressible fluid scheme object
  !> @param msh Mesh data structure
  !> @param lx Polynomial order in x-direction
  !> @param params JSON configuration parameters
  !> @param scheme Name of the numerical scheme
  !> @param user User-defined parameters and functions
  subroutine fluid_scheme_compressible_init(this, msh, lx, params, scheme, user)
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user

    !
    ! SEM simulation fundamentals
    !

    this%msh => msh

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Assign Dofmap to scratch registry
    call neko_scratch_registry%set_dofmap(this%dm_Xh)

    ! Case parameters
    this%params => params

    ! Assign a name
    call json_get_or_default(params, 'case.fluid.name', this%name, "fluid")

    ! Material properties will be set up via set_material_properties
    call neko_registry%add_field(this%dm_Xh, this%name // "_rho")
    this%rho => neko_registry%get_field(this%name // "_rho")

    ! Assign momentum fields
    call neko_registry%add_field(this%dm_Xh, "m_x")
    call neko_registry%add_field(this%dm_Xh, "m_y")
    call neko_registry%add_field(this%dm_Xh, "m_z")
    this%m_x => neko_registry%get_field("m_x")
    this%m_y => neko_registry%get_field("m_y")
    this%m_z => neko_registry%get_field("m_z")
    call this%m_x%init(this%dm_Xh, "m_x")
    call this%m_y%init(this%dm_Xh, "m_y")
    call this%m_z%init(this%dm_Xh, "m_z")

    ! Assign energy field
    call neko_registry%add_field(this%dm_Xh, "E")
    this%E => neko_registry%get_field("E")
    call this%E%init(this%dm_Xh, "E")

    ! Assign temperature field
    call neko_registry%add_field(this%dm_Xh, "temperature")
    this%temperature => neko_registry%get_field("temperature")
    call this%temperature%init(this%dm_Xh, "temperature")

    ! Assign maximum wave speed field
    call neko_registry%add_field(this%dm_Xh, "max_wave_speed")
    this%max_wave_speed => neko_registry%get_field("max_wave_speed")
    call this%max_wave_speed%init(this%dm_Xh, "max_wave_speed")

    ! Assign entropy field
    call neko_registry%add_field(this%dm_Xh, "S")
    this%S => neko_registry%get_field("S")
    call this%S%init(this%dm_Xh, "S")

    ! Assign artificial viscosity field (without physical viscosity)
    call neko_registry%add_field(this%dm_Xh, "artificial_visc")
    this%artificial_visc => neko_registry%get_field("artificial_visc")
    call this%artificial_visc%init(this%dm_Xh, "artificial_visc")

    ! ! Assign velocity fields
    call neko_registry%add_field(this%dm_Xh, "u")
    call neko_registry%add_field(this%dm_Xh, "v")
    call neko_registry%add_field(this%dm_Xh, "w")
    this%u => neko_registry%get_field("u")
    this%v => neko_registry%get_field("v")
    this%w => neko_registry%get_field("w")
    call this%u%init(this%dm_Xh, "u")
    call this%v%init(this%dm_Xh, "v")
    call this%w%init(this%dm_Xh, "w")
    call neko_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_registry%get_field('p')
    call this%p%init(this%dm_Xh, "p")

    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_x)
    allocate(this%f_y)
    allocate(this%f_z)
    call this%f_x%init(this%dm_Xh, fld_name = "fluid_rhs_x")
    call this%f_y%init(this%dm_Xh, fld_name = "fluid_rhs_y")
    call this%f_z%init(this%dm_Xh, fld_name = "fluid_rhs_z")

    ! Material properties
    call this%set_material_properties(params, user)

    ! Compressible parameters
    call json_get_or_default(params, 'case.fluid.gamma', this%gamma, 1.4_rp)

    call this%euler_idp%init(params)
    call fluid_scheme_compressible_validate_idp_setup(this, params, user)

    ! Calculate global points for logging
    this%glb_n_points = int(this%msh%glb_nelv, i8)*int(this%Xh%lxyz, i8)
    this%glb_unique_points = int(glsum(this%c_Xh%mult, this%dm_Xh%size()), i8)

    !
    ! Log solver information
    !
    call this%log_solver_info(params, scheme, lx)
  end subroutine fluid_scheme_compressible_init

  !> Free allocated memory and cleanup resources
  !> @param this The compressible fluid scheme object to destroy
  subroutine fluid_scheme_compressible_free(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    class(bc_t), pointer :: bc
    integer :: i

    do i = 1, this%bcs_vel%size()
       bc => this%bcs_vel%get(i)
       if (associated(bc)) then
          call bc%free()
          deallocate(bc)
       end if
    end do
    call this%bcs_vel%free()

    do i = 1, this%bcs_prs%size()
       bc => this%bcs_prs%get(i)
       if (associated(bc)) then
          call bc%free()
          deallocate(bc)
       end if
    end do
    call this%bcs_prs%free()

    call this%dm_Xh%free()
    call this%gs_Xh%free()
    call this%c_Xh%free()
    call this%Xh%free()

    if (associated(this%m_x)) then
       call this%m_x%free()
    end if

    if (associated(this%m_y)) then
       call this%m_y%free()
    end if

    if (associated(this%m_z)) then
       call this%m_z%free()
    end if

    if (associated(this%E)) then
       call this%E%free()
    end if

    if (associated(this%temperature)) then
       call this%temperature%free()
    end if

    if (associated(this%max_wave_speed)) then
       call this%max_wave_speed%free()
    end if

    if (associated(this%S)) then
       call this%S%free()
    end if

    if (associated(this%f_x)) then
       call this%f_x%free()
       deallocate(this%f_x)
    end if

    if (associated(this%f_y)) then
       call this%f_y%free()
       deallocate(this%f_y)
    end if

    if (associated(this%f_z)) then
       call this%f_z%free()
       deallocate(this%f_z)
    end if

    nullify(this%f_x)
    nullify(this%f_y)
    nullify(this%f_z)

    nullify(this%m_x)
    nullify(this%m_y)
    nullify(this%m_z)
    nullify(this%E)
    nullify(this%temperature)
    nullify(this%max_wave_speed)
    nullify(this%S)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%rho)
    nullify(this%mu)
    nullify(this%kappa)

    call this%material_properties%free()

  end subroutine fluid_scheme_compressible_free

  !> Validate field initialization and compute derived quantities
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_validate(this)
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    integer :: n, i, first_invalid
    type(field_t), pointer :: temp
    integer :: temp_indices(1)
    integer, allocatable :: state_status(:)
    real(kind=rp), allocatable :: rho_primitive(:), sound_speed(:)
    real(kind=rp), allocatable :: internal_energy(:)
    character(len=LOG_SIZE) :: log_buf

    n = this%dm_Xh%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1), .false.)

    if (this%euler_idp%enabled) then
       do i = 1, n
          if (.not. ieee_is_finite(this%rho%x(i,1,1,1)) .or. &
               this%rho%x(i,1,1,1) .le. 0.0_rp) then
             write(log_buf, '(A,I0,A,I0)') &
                  'Euler IDP initial density is invalid on rank ', pe_rank, &
                  ' at local node ', i
             call neko_error(trim(log_buf))
          end if
          if (.not. ieee_is_finite(this%p%x(i,1,1,1)) .or. &
               this%p%x(i,1,1,1) / (this%gamma - 1.0_rp) .lt. &
               this%euler_idp%internal_energy_floor) then
             write(log_buf, '(A,I0,A,I0)') &
                  'Euler IDP initial internal energy is invalid on rank ', &
                  pe_rank, ' at local node ', i
             call neko_error(trim(log_buf))
          end if
       end do
    end if

    !> Initialize the momentum field
    call field_col3(this%m_x, this%u, this%rho)
    call field_col3(this%m_y, this%v, this%rho)
    call field_col3(this%m_z, this%w, this%rho)

    !> Initialize total energy
    !> Total energy E := p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2)
    call field_cmult2(this%E, this%p, 1.0_rp/(this%gamma - 1.0_rp), n)
    call field_col3(temp, this%u, this%u, n)
    call field_addcol3(temp, this%v, this%v, n)
    call field_addcol3(temp, this%w, this%w, n)
    call field_col2(temp, this%rho, n)
    call field_cmult(temp, 0.5_rp, n)
    call field_add2(this%E, temp, n)

    !> Initialize temperature T = p / (rho * (gamma - 1))
    do i = 1, n
       this%temperature%x(i,1,1,1) = this%p%x(i,1,1,1) / &
            (this%rho%x(i,1,1,1) * (this%gamma - 1.0_rp))
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

    if (this%euler_idp%enabled) then
       allocate(rho_primitive(n), sound_speed(n), internal_energy(n))
       allocate(state_status(n))
       call compressible_ops_cpu_conserved_to_primitive(this%rho%x, &
            this%m_x%x, this%m_y%x, this%m_z%x, this%E%x, this%gamma, &
            this%euler_idp%internal_energy_floor, rho_primitive, &
            this%u%x, this%v%x, this%w%x, this%p%x, sound_speed, &
            internal_energy, state_status, n)
       if (any(state_status .ne. EULER_STATE_OK)) then
          first_invalid = minloc(state_status, dim = 1, &
               mask = state_status .ne. EULER_STATE_OK)
          write(log_buf, '(A,I0,A,I0,A,I0)') &
               'Euler IDP initial conserved state is invalid on rank ', &
               pe_rank, ' at local node ', first_invalid, ', status ', &
               state_status(first_invalid)
          call neko_error(trim(log_buf))
       end if
       deallocate(rho_primitive, sound_speed, internal_energy, state_status)
    end if

    !> Compute initial entropy and maximum wave speed from initial conditions
    call this%compute_entropy()
    call this%compute_max_wave_speed()

  end subroutine fluid_scheme_compressible_validate

  !> Compute CFL number
  !> @param this The compressible fluid scheme object
  !> @param dt Current timestep size
  !> @return Computed CFL number
  function fluid_scheme_compressible_compute_cfl(this, dt) result(c)
    class(fluid_scheme_compressible_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c
    integer :: n

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         rho => this%rho, Xh => this%Xh, c_Xh => this%c_Xh, &
         msh => this%msh, gamma => this%gamma, &
         max_wave_speed => this%max_wave_speed)

      n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

      ! Use the compressible CFL function with precomputed maximum wave speed
      c = cfl_compressible(dt, max_wave_speed%x, Xh, c_Xh, msh%nelv, msh%gdim)
    end associate

  end function fluid_scheme_compressible_compute_cfl

  !> Set material properties mu and rho
  !> @param this The compressible fluid scheme object
  !> @param params The case parameter file
  !> @param user The user interface
  subroutine fluid_scheme_compressible_set_material_properties(this, &
       params, user)
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    procedure(user_material_properties_intf), pointer :: dummy_mp_ptr
    type(time_state_t) :: dummy_time_state
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: const_mu, const_kappa

    dummy_mp_ptr => dummy_user_material_properties

    call neko_registry%add_field(this%dm_Xh, this%name // "_mu")
    this%mu => neko_registry%get_field(this%name // "_mu")

    call neko_registry%add_field(this%dm_Xh, this%name // "_kappa")
    this%kappa => neko_registry%get_field(this%name // "_kappa")

    call this%material_properties%init(3)
    call this%material_properties%assign(1, this%rho)
    call this%material_properties%assign(2, this%mu)
    call this%material_properties%assign(3, this%kappa)

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then
       this%user_material_properties => user%material_properties
       call user%material_properties(this%name, this%material_properties, &
            dummy_time_state)
    else
       this%user_material_properties => dummy_user_material_properties
       call json_get_or_lookup_or_default(params, 'case.fluid.mu', const_mu, &
            0.0_rp)
       call json_get_or_lookup_or_default(params, 'case.fluid.kappa', &
            const_kappa, 0.0_rp)

       call field_cfill(this%mu, const_mu)
       call field_cfill(this%kappa, const_kappa)

       write(log_buf, '(A,ES13.6)') 'mu         :', const_mu
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'kappa      :', const_kappa
       call neko_log%message(log_buf)
    end if

  end subroutine fluid_scheme_compressible_set_material_properties

  !> Update variable material properties
  !> @param this The compressible fluid scheme object
  !> @param time The time state
  subroutine fluid_scheme_compressible_update_material_properties(this, time)
    use device, only : device_memcpy, HOST_TO_DEVICE
    class(fluid_scheme_compressible_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%user_material_properties(this%name, this%material_properties, &
         time)

  end subroutine fluid_scheme_compressible_update_material_properties

  !> Compute entropy field S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_compute_entropy(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    integer :: n

    n = this%S%dof%size()

    !> TODO: Add support for SX
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call compressible_ops_device_compute_entropy(this%S, this%p, this%rho, &
            this%gamma, n)
    else
       call compressible_ops_cpu_compute_entropy(this%S%x, this%p%x, &
            this%rho%x, this%gamma, n)
    end if

  end subroutine fluid_scheme_compressible_compute_entropy

  !> Compute maximum wave speed for compressible flows
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_compute_max_wave_speed(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    integer :: n

    n = this%u%dof%size()

    !> TODO: Add support for SX
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call compressible_ops_device_compute_max_wave_speed( &
            this%max_wave_speed, this%u, this%v, this%w, this%gamma, this%p, &
            this%rho, n)
    else
       call compressible_ops_cpu_compute_max_wave_speed(this%max_wave_speed%x, &
            this%u%x, this%v%x, this%w%x, this%gamma, this%p%x, this%rho%x, n)
    end if

  end subroutine fluid_scheme_compressible_compute_max_wave_speed

  !> Log comprehensive solver information
  !> @param this The compressible fluid scheme object
  !> @param params JSON configuration parameters
  !> @param scheme Name of the numerical scheme
  !> @param lx Polynomial order in x-direction
  subroutine fluid_scheme_compressible_log_solver_info(this, params, scheme, lx)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    character(len=*), intent(in) :: scheme
    integer, intent(in) :: lx
    character(len=LOG_SIZE) :: log_buf
    logical :: logical_val
    real(kind=rp) :: real_val
    integer :: integer_val

    call neko_log%section('Fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)

    ! Polynomial order
    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'Poly order : ', lx-1
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'Poly order : ', lx-1
    else
       write(log_buf, '(A, I3)') 'Poly order : ', lx-1
    end if
    call neko_log%message(log_buf)

    ! Global points information
    write(log_buf, '(A, I0)') 'GLL points : ', this%glb_n_points
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') 'Unique pts.: ', this%glb_unique_points
    call neko_log%message(log_buf)

    ! Material properties
    write(log_buf, '(A,ES13.6)') 'gamma      :', this%gamma
    call neko_log%message(log_buf)

    ! Compressible-specific parameters
    call json_get_or_default(params, 'case.numerics.c_avisc_low', real_val, &
         0.5_rp)
    write(log_buf, '(A,ES13.6)') 'c_avisc_low:', real_val
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.numerics.time_order', integer_val, 4)
    write(log_buf, '(A, I0)') 'RK order   : ', integer_val
    call neko_log%message(log_buf)
    if (this%euler_idp%enabled) then
       write(log_buf, '(A)') 'Euler IDP  : enabled'
       call neko_log%message(log_buf)
    end if
    call neko_log%end_section()

  end subroutine fluid_scheme_compressible_log_solver_info

  !> Validate the restrictions of the opt-in Euler IDP path.
  subroutine fluid_scheme_compressible_validate_idp_setup(this, params, user)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), intent(in) :: user
    procedure(user_material_properties_intf), pointer :: dummy_mp_ptr
    character(len=LOG_SIZE) :: message
    logical :: dealias, oifs, cyclic, ale_enabled
    logical :: local_valid, global_valid
    integer :: time_order, ierr, local_periodic, global_periodic
    integer :: e, f, i, n_facets
    real(kind=rp) :: metric_tolerance

    if (.not. this%euler_idp%enabled) return

    if (.not. this%euler_idp%valid(message)) then
       call neko_error('Invalid Euler IDP configuration: ' // trim(message))
    end if
    if (NEKO_BCKND_DEVICE .eq. 1 .or. NEKO_BCKND_SX .eq. 1) then
       call neko_error('Euler IDP currently requires the CPU backend')
    end if

    call json_get_or_default(params, 'case.numerics.time_order', &
         time_order, 4)
    if (time_order .ne. 1 .and. time_order .ne. 3) then
       call neko_error('Euler IDP requires time_order 1 or 3')
    end if
    call json_get_or_default(params, 'case.numerics.dealias', &
         dealias, .false.)
    if (dealias) then
       call neko_error('Euler IDP does not support dealiasing')
    end if
    call json_get_or_default(params, 'case.numerics.oifs', oifs, .false.)
    if (oifs) then
       call neko_error('Euler IDP does not support OIFS')
    end if
    if (this%Xh%t .ne. GLL) then
       call neko_error('Euler IDP requires collocated GLL nodes')
    end if

    cyclic = .false.
    if (params%valid_path('case.fluid.cyclic')) then
       call json_get_or_default(params, 'case.fluid.cyclic', cyclic, .false.)
    end if
    if (cyclic .or. this%c_Xh%cyclic) then
       call neko_error('Euler IDP does not support rotational periodicity')
    end if
    if (params%valid_path('case.fluid.boundary_conditions')) then
       call neko_error('Euler IDP requires a fully periodic domain')
    end if
    if (params%valid_path('case.fluid.source_terms')) then
       call neko_error('Euler IDP does not support fluid source terms')
    end if

    ale_enabled = .false.
    if (params%valid_path('case.fluid.ale.enabled')) then
       call json_get_or_default(params, 'case.fluid.ale.enabled', &
            ale_enabled, .false.)
    end if
    if (ale_enabled) then
       call neko_error('Euler IDP does not support ALE')
    end if

    dummy_mp_ptr => dummy_user_material_properties
    if (.not. associated(user%material_properties, dummy_mp_ptr)) then
       call neko_error( &
            'Euler IDP does not support variable material properties')
    end if
    if (any(this%mu%x .ne. 0.0_rp) .or. &
         any(this%kappa%x .ne. 0.0_rp)) then
       call neko_error( &
            'Euler IDP requires zero physical viscosity and conductivity')
    end if
    if (.not. ieee_is_finite(this%gamma) .or. this%gamma .le. 1.0_rp) then
       call neko_error('Euler IDP requires a finite gamma greater than one')
    end if
    if (this%gamma .gt. 5.0_rp / 3.0_rp) then
       call neko_error( &
            'Euler IDP guaranteed wave speed requires gamma <= 5/3')
    end if

    local_valid = .true.
    do i = 1, size(this%msh%labeled_zones)
       if (this%msh%labeled_zones(i)%size .gt. 0) local_valid = .false.
    end do
    if (this%msh%gdim .eq. 2) then
       n_facets = 4
    else if (this%msh%gdim .eq. 3) then
       n_facets = 6
    else
       local_valid = .false.
       n_facets = 0
    end if
    do e = 1, this%msh%nelv
       do f = 1, n_facets
          if (this%msh%facet_neigh(f,e) .eq. 0 .and. &
               .not. fluid_scheme_compressible_is_periodic_facet( &
               this%msh, f, e)) local_valid = .false.
       end do
    end do
    local_periodic = this%msh%periodic%size
    call MPI_Allreduce(local_periodic, global_periodic, 1, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM, ierr)
    if (global_periodic .eq. 0) local_valid = .false.
    call MPI_Allreduce(local_valid, global_valid, 1, MPI_LOGICAL, MPI_LAND, &
         NEKO_COMM, ierr)
    if (.not. global_valid) then
       call neko_error( &
            'Euler IDP requires complete translational periodic pairing')
    end if

    metric_tolerance = 100.0_rp * epsilon(1.0_rp) * &
         real(max(1, this%Xh%lxyz), rp)
    local_valid = .true.
    do e = 1, this%msh%nelv
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dxdr, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dydr, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dzdr, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dxds, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dyds, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dzds, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dxdt, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dydt, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%dzdt, e, &
            metric_tolerance)
       local_valid = local_valid .and. &
            fluid_scheme_compressible_constant_metric(this%c_Xh%jac, e, &
            metric_tolerance)
    end do
    call MPI_Allreduce(local_valid, global_valid, 1, MPI_LOGICAL, MPI_LAND, &
         NEKO_COMM, ierr)
    if (.not. global_valid) then
       call neko_error('Euler IDP requires affine element metrics')
    end if
  end subroutine fluid_scheme_compressible_validate_idp_setup

  !> Return whether a local mesh facet has periodic metadata.
  logical function fluid_scheme_compressible_is_periodic_facet(msh, facet, &
       element) result(is_periodic)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: facet, element
    integer :: i

    is_periodic = .false.
    do i = 1, msh%periodic%size
       if (msh%periodic%facet_el(i)%x(1) .eq. facet .and. &
            msh%periodic%facet_el(i)%x(2) .eq. element) then
          is_periodic = .true.
          return
       end if
    end do
  end function fluid_scheme_compressible_is_periodic_facet

  !> Return whether a geometric metric is constant on one element.
  logical function fluid_scheme_compressible_constant_metric(metric, &
       element, tolerance) result(is_constant)
    real(kind=rp), intent(in) :: metric(:,:,:,:)
    integer, intent(in) :: element
    real(kind=rp), intent(in) :: tolerance
    real(kind=rp) :: reference, scale

    reference = metric(1,1,1,element)
    scale = max(1.0_rp, maxval(abs(metric(:,:,:,element))))
    is_constant = maxval(abs(metric(:,:,:,element) - reference)) .le. &
         tolerance * scale
  end function fluid_scheme_compressible_constant_metric

end module fluid_scheme_compressible
