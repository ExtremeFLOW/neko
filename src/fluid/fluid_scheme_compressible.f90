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
  use field_math, only : field_cfill, field_col2, field_col3, &
       field_cmult2, field_cmult, field_addcol3, field_add2

  use field_registry, only : neko_field_registry
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_module, only : json_file
  use num_types, only : rp
  use mesh, only : mesh_t
  use scratch_registry, only : scratch_registry_t
  use space, only : GLL
  use user_intf, only : user_t
  use json_utils, only : json_get_or_default
  use mpi_f08
  use operators, only : cfl_compressible
  use compressible_ops_cpu, only : compressible_ops_cpu_compute_max_wave_speed, &
       compressible_ops_cpu_compute_entropy
  use compressible_ops_device, only : compressible_ops_device_compute_max_wave_speed, &
       compressible_ops_device_compute_entropy
  use neko_config, only : NEKO_BCKND_DEVICE
  use time_state, only : time_state_t
  use logger, only : neko_log, LOG_SIZE
  use math, only : glsum
  use num_types, only : i8
  implicit none
  private

  !> Base type of compressible fluid formulations.
  type, public, abstract, extends(fluid_scheme_base_t) :: fluid_scheme_compressible_t
     !> The momentum field
     type(field_t), pointer :: m_x => null() !< x-component of Momentum
     type(field_t), pointer :: m_y => null() !< y-component of Momentum
     type(field_t), pointer :: m_z => null() !< z-component of Momentum
     type(field_t), pointer :: E => null() !< Total energy
     type(field_t), pointer :: max_wave_speed => null() !< Maximum wave speed field
     type(field_t), pointer :: S => null() !< Entropy field

     real(kind=rp) :: gamma

     !> Global number of GLL points for the fluid (not unique)
     integer(kind=i8) :: glb_n_points
     !> Global number of GLL points for the fluid (unique)
     integer(kind=i8) :: glb_unique_points

     type(scratch_registry_t) :: scratch !< Manager for temporary fields

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

    ! Local scratch registry
    call this%scratch%init(this%dm_Xh, 10, 2)

    ! Case parameters
    this%params => params

    ! Assign a name
    call json_get_or_default(params, 'case.fluid.name', this%name, "fluid")

    ! Fill mu and rho field with the physical value
    call neko_field_registry%add_field(this%dm_Xh, this%name // "_mu")
    call neko_field_registry%add_field(this%dm_Xh, this%name // "_rho")
    this%mu => neko_field_registry%get_field(this%name // "_mu")
    this%rho => neko_field_registry%get_field(this%name // "_rho")
    call field_cfill(this%mu, 0.0_rp, this%mu%size())

    ! Assign momentum fields
    call neko_field_registry%add_field(this%dm_Xh, "m_x")
    call neko_field_registry%add_field(this%dm_Xh, "m_y")
    call neko_field_registry%add_field(this%dm_Xh, "m_z")
    this%m_x => neko_field_registry%get_field("m_x")
    this%m_y => neko_field_registry%get_field("m_y")
    this%m_z => neko_field_registry%get_field("m_z")
    call this%m_x%init(this%dm_Xh, "m_x")
    call this%m_y%init(this%dm_Xh, "m_y")
    call this%m_z%init(this%dm_Xh, "m_z")

    ! Assign energy field
    call neko_field_registry%add_field(this%dm_Xh, "E")
    this%E => neko_field_registry%get_field("E")
    call this%E%init(this%dm_Xh, "E")

    ! Assign maximum wave speed field
    call neko_field_registry%add_field(this%dm_Xh, "max_wave_speed")
    this%max_wave_speed => neko_field_registry%get_field("max_wave_speed")
    call this%max_wave_speed%init(this%dm_Xh, "max_wave_speed")

    ! Assign entropy field
    call neko_field_registry%add_field(this%dm_Xh, "S")
    this%S => neko_field_registry%get_field("S")
    call this%S%init(this%dm_Xh, "S")

    ! ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, "u")
    call neko_field_registry%add_field(this%dm_Xh, "v")
    call neko_field_registry%add_field(this%dm_Xh, "w")
    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    call this%u%init(this%dm_Xh, "u")
    call this%v%init(this%dm_Xh, "v")
    call this%w%init(this%dm_Xh, "w")
    call neko_field_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_field_registry%get_field('p')
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

    ! Compressible parameters
    call json_get_or_default(params, 'case.fluid.gamma', this%gamma, 1.4_rp)

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
    nullify(this%max_wave_speed)
    nullify(this%S)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%rho)
    nullify(this%mu)

  end subroutine fluid_scheme_compressible_free

  !> Validate field initialization and compute derived quantities
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_validate(this)
    class(fluid_scheme_compressible_t), target, intent(inout) :: this
    integer :: n
    type(field_t), pointer :: temp
    integer :: temp_indices(1)

    n = this%dm_Xh%size()
    call this%scratch%request_field(temp, temp_indices(1))

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

    call this%scratch%relinquish_field(temp_indices)

    !> Compute initial maximum wave speed from initial conditions
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

  !> Set rho and mu
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_update_material_properties(this, time)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
  end subroutine fluid_scheme_compressible_update_material_properties

  !> Compute entropy field S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
  !> @param this The compressible fluid scheme object
  subroutine fluid_scheme_compressible_compute_entropy(this)
    class(fluid_scheme_compressible_t), intent(inout) :: this
    integer :: n

    n = this%S%dof%size()

    !> TODO: Add support for SX
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call compressible_ops_device_compute_entropy(this%S, this%p, this%rho, this%gamma, n)
    else
       call compressible_ops_cpu_compute_entropy(this%S%x, this%p%x, this%rho%x, this%gamma, n)
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
       call compressible_ops_device_compute_max_wave_speed(this%max_wave_speed, &
            this%u, this%v, this%w, this%gamma, this%p, this%rho, n)
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
    call json_get_or_default(params, 'case.numerics.c_avisc_low', real_val, 0.5_rp)
    write(log_buf, '(A,ES13.6)') 'c_avisc_low:', real_val
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.numerics.time_order', integer_val, 4)
    write(log_buf, '(A, I0)') 'RK order   : ', integer_val
    call neko_log%message(log_buf)

  end subroutine fluid_scheme_compressible_log_solver_info

end module fluid_scheme_compressible
