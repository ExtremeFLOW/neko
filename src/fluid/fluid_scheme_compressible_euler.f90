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
module fluid_scheme_compressible_euler
  use comm, only : NEKO_COMM
  use advection, only : advection_t
  use device, only : device_memcpy, HOST_TO_DEVICE
  use field_math, only : field_add2, field_cfill, field_cmult, &
       field_copy, field_col2, field_col3, &
       field_addcol3, field_sub2, field_invcol2
  use math, only : col2
  use device_math, only : device_col2
  use field, only : field_t
  use fluid_scheme_compressible, only: fluid_scheme_compressible_t
  use gs_ops, only : GS_OP_ADD
  use num_types, only : rp
  use mesh, only : mesh_t
  use checkpoint, only : chkp_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use profiler, only : profiler_start_region, profiler_end_region
  use user_intf, only : user_t
  use time_step_controller, only : time_step_controller_t
  use ax_product, only : ax_t, ax_helm_factory
  use coefs, only: coef_t
  use euler_residual, only: euler_rhs_t, euler_rhs_factory
  use neko_config, only : NEKO_BCKND_DEVICE
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use bc_list, only: bc_list_t
  use bc, only : bc_t
  use utils, only : neko_error, neko_type_error
  use logger, only : LOG_SIZE
  use time_state, only : time_state_t
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) &
       :: fluid_scheme_compressible_euler_t
     type(field_t) :: rho_res, m_x_res, m_y_res, m_z_res, m_E_res
     type(field_t) :: drho, dm_x, dm_y, dm_z, dE
     type(field_t) :: h
     real(kind=rp) :: c_avisc_low
     class(advection_t), allocatable :: adv
     class(ax_t), allocatable :: Ax
     class(euler_rhs_t), allocatable :: euler_rhs
     type(runge_kutta_time_scheme_t) :: rk_scheme

     ! List of boundary conditions for velocity
     type(bc_list_t) :: bcs_density
   contains
     procedure, pass(this) :: init => fluid_scheme_compressible_euler_init
     procedure, pass(this) :: free => fluid_scheme_compressible_euler_free
     procedure, pass(this) :: step => fluid_scheme_compressible_euler_step
     procedure, pass(this) :: restart => fluid_scheme_compressible_euler_restart
     !> Set up boundary conditions.
     procedure, pass(this) :: setup_bcs &
          => fluid_scheme_compressible_euler_setup_bcs
     procedure, pass(this) :: compute_h
  end type fluid_scheme_compressible_euler_t

  interface
     !> Boundary condition factory for density.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_euler_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine density_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_euler_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine density_bc_factory
  end interface

  interface
     !> Boundary condition factory for pressure.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_euler_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine pressure_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_euler_t), intent(inout) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine pressure_bc_factory
  end interface

  interface
     !> Boundary condition factory for velocity
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_euler_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine velocity_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_euler_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine velocity_bc_factory
  end interface

contains
  !> Initialize the compressible Euler fluid scheme
  !! @param this The fluid scheme object
  !! @param msh Mesh data structure
  !! @param lx Polynomial order in x-direction
  !! @param params JSON configuration parameters
  !! @param user User-defined parameters and functions
  !! @param chkp Checkpoint to write to
  subroutine fluid_scheme_compressible_euler_init(this, msh, lx, params, user, &
       chkp)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    character(len=12), parameter :: scheme = 'compressible'
    integer :: rk_order

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, scheme, user)

    call euler_rhs_factory(this%euler_rhs)

    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%drho%init(dm_Xh, 'drho')
      call this%dm_x%init(dm_Xh, 'dm_x')
      call this%dm_y%init(dm_Xh, 'dm_y')
      call this%dm_z%init(dm_Xh, 'dm_z')
      call this%dE%init(dm_Xh, 'dE')
      call this%h%init(dm_Xh, 'h')

    end associate

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(p => this%p, rho => this%rho, &
            u => this%u, v => this%v, w => this%w, &
            m_x => this%m_x, m_y => this%m_y, m_z => this%m_z)
         call device_memcpy(p%x, p%x_d, p%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(rho%x, rho%x_d, rho%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(u%x, u%x_d, u%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(v%x, v%x_d, v%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(w%x, w%x_d, w%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_x%x, m_x%x_d, m_x%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_y%x, m_y%x_d, m_y%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_z%x, m_z%x_d, m_z%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
       end associate
    end if

    ! Initialize the diffusion operator
    call ax_helm_factory(this%Ax, full_formulation = .false.)

    ! Compute h
    call this%compute_h()
    call json_get_or_default(params, 'case.numerics.c_avisc_low', &
         this%c_avisc_low, 0.5_rp)

    ! Initialize Runge-Kutta scheme
    call json_get_or_default(params, 'case.numerics.time_order', rk_order, 4)
    call this%rk_scheme%init(rk_order)

    ! Set up boundary conditions
    call this%setup_bcs(user, params)

  end subroutine fluid_scheme_compressible_euler_init

  !> Free allocated memory and cleanup
  !> @param this The fluid scheme object to destroy
  subroutine fluid_scheme_compressible_euler_free(this)
    class(fluid_scheme_compressible_euler_t), intent(inout) :: this

    call this%scheme_free()

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    call this%drho%free()
    call this%dm_x%free()
    call this%dm_y%free()
    call this%dm_z%free()
    call this%dE%free()

    ! call this%scheme_free()
  end subroutine fluid_scheme_compressible_euler_free

  !> Advance the fluid simulation one timestep
  !> @param this The fluid scheme object
  !> @param time Current simulation time state
  !> @param ext_bdf Time integration controller
  !> @param dt_controller Timestep size controller
  subroutine fluid_scheme_compressible_euler_step(this, time, dt_controller)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    type(field_t), pointer :: temp
    integer :: temp_indices(1)
    ! number of degrees of freedom
    integer :: n
    integer :: i
    class(bc_t), pointer :: b

    n = this%dm_Xh%size()
    call this%scratch%request_field(temp, temp_indices(1))
    b => null()

    call profiler_start_region('Fluid compressible', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         m_x=> this%m_x, m_y => this%m_y, m_z => this%m_z, &
         Xh => this%Xh, msh => this%msh, Ax => this%Ax, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         E => this%E, rho => this%rho, mu => this%mu, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
         drho => this%drho, dm_x => this%dm_x, dm_y => this%dm_y, &
         dm_z => this%dm_z, dE => this%dE, &
         euler_rhs => this%euler_rhs, h => this%h, &
         t => time%t, tstep => time%tstep, dt => time%dt, &
         ext_bdf => this%ext_bdf, &
         c_avisc_low => this%c_avisc_low, rk_scheme => this%rk_scheme)

      call euler_rhs%step(rho, m_x, m_y, m_z, E, &
           p, u, v, w, Ax, &
           c_Xh, gs_Xh, h, c_avisc_low, &
           rk_scheme, dt)

      !> Apply density boundary conditions
      call this%bcs_density%apply(rho, time)

      !> Update variables
      ! Update u, v, w
      call field_copy(u, m_x, n)
      call field_invcol2(u, rho, n)
      call field_copy(v, m_y, n)
      call field_invcol2(v, rho, n)
      call field_copy(w, m_z, n)
      call field_invcol2(w, rho, n)

      !> Apply velocity boundary conditions
      call this%bcs_vel%apply_vector(u%x, v%x, w%x, &
           dm_Xh%size(), time, strong = .true.)
      call field_copy(m_x, u, n)
      call field_col2(m_x, rho, n)
      call field_copy(m_y, v, n)
      call field_col2(m_y, rho, n)
      call field_copy(m_z, w, n)
      call field_col2(m_z, rho, n)

      !> Update p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2 + w^2))
      call field_col3(temp, u, u, n)
      call field_addcol3(temp, v, v, n)
      call field_addcol3(temp, w, w, n)
      call field_col2(temp, rho, n)
      call field_cmult(temp, 0.5_rp, n)
      call field_copy(p, E, n)
      call field_sub2(p, temp, n)
      call field_cmult(p, this%gamma - 1.0_rp, n)

      !> Apply pressure boundary conditions
      call this%bcs_prs%apply(p, time)
      ! TODO: Make sure pressure is positive
      ! E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2)
      call field_copy(E, p, n)
      call field_cmult(E, 1.0_rp / (this%gamma - 1.0_rp), n)
      ! temp = 0.5 * rho * (u^2 + v^2 + w^2)
      call field_add2(E, temp, n)

      !> Compute entropy S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
      call this%compute_entropy()

      !> Update maximum wave speed for CFL computation
      call this%compute_max_wave_speed()

      do i = 1, this%bcs_vel%size()
         b => this%bcs_vel%get(i)
         b%updated = .false.
      end do

      do i = 1, this%bcs_prs%size()
         b => this%bcs_prs%get(i)
         b%updated = .false.
      end do

      do i = 1, this%bcs_density%size()
         b => this%bcs_density%get(i)
         b%updated = .false.
      end do
      nullify(b)

    end associate
    call profiler_end_region('Fluid compressible', 1)

    call this%scratch%relinquish_field(temp_indices)

  end subroutine fluid_scheme_compressible_euler_step

  !> Set up boundary conditions for the fluid scheme
  !> @param this The fluid scheme object
  !> @param user User-defined boundary conditions
  !> @param params Configuration parameters
  subroutine fluid_scheme_compressible_euler_setup_bcs(this, user, params)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: i, n_bcs, zone_index, j, zone_size, global_zone_size, ierr
    class(bc_t), pointer :: bc_i
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found
    integer, allocatable :: zone_indices(:)
    character(len=LOG_SIZE) :: log_buf

    ! Process boundary conditions
    if (params%valid_path('case.fluid.boundary_conditions')) then
       call params%info('case.fluid.boundary_conditions', n_children = n_bcs)
       call params%get_core(core)
       call params%get('case.fluid.boundary_conditions', bc_object, found)

       !
       ! Velocity bcs
       !
       call this%bcs_vel%init(n_bcs)

       do i = 1, n_bcs
          ! Extract BC configuration
          call json_extract_item(core, bc_object, i, bc_subdict)
          call json_get(bc_subdict, "zone_indices", zone_indices)

          ! Validate zones
          do j = 1, size(zone_indices)
             zone_size = this%msh%labeled_zones(zone_indices(j))%size
             call MPI_Allreduce(zone_size, global_zone_size, 1, &
                  MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

             if (global_zone_size .eq. 0) then
                write(log_buf,'(A,I0,A)') "Error: Zone ", zone_indices(j), &
                     " has zero size"
                call neko_error(log_buf)
             end if
          end do

          ! Create BC
          bc_i => null()
          call velocity_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Add to appropriate lists
          if (associated(bc_i)) then
             call this%bcs_vel%append(bc_i)
          end if
       end do

       !
       ! Pressure bcs
       !
       call this%bcs_prs%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call pressure_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_prs%append(bc_i)
          end if
       end do

       !
       ! Density bcs
       !
       call this%bcs_density%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call density_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_density%append(bc_i)
          end if
       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             call neko_error("No boundary_conditions entry in the case file!")
          end if
       end do

       ! For a pure periodic case, we still need to initilise the bc lists
       ! to a zero size to avoid issues with apply() in step()
       call this%bcs_prs%init()
       call this%bcs_vel%init()
       call this%bcs_density%init()

    end if
  end subroutine fluid_scheme_compressible_euler_setup_bcs

  !> Copied from les_model_compute_delta in les_model.f90
  !> TODO: move to a separate module
  !> Compute characteristic mesh size h
  !> @param this The fluid scheme object
  subroutine compute_h(this)
    class(fluid_scheme_compressible_euler_t), intent(inout) :: this
    integer :: e, i, j, k
    integer :: im, ip, jm, jp, km, kp
    real(kind=rp) :: di, dj, dk, ndim_inv
    integer :: lx_half, ly_half, lz_half

    lx_half = this%c_Xh%Xh%lx / 2
    ly_half = this%c_Xh%Xh%ly / 2
    lz_half = this%c_Xh%Xh%lz / 2

    do e = 1, this%c_Xh%msh%nelv
       do k = 1, this%c_Xh%Xh%lz
          km = max(1, k-1)
          kp = min(this%c_Xh%Xh%lz, k+1)

          do j = 1, this%c_Xh%Xh%ly
             jm = max(1, j-1)
             jp = min(this%c_Xh%Xh%ly, j+1)

             do i = 1, this%c_Xh%Xh%lx
                im = max(1, i-1)
                ip = min(this%c_Xh%Xh%lx, i+1)

                di = (this%c_Xh%dof%x(ip, j, k, e) - &
                     this%c_Xh%dof%x(im, j, k, e))**2 &
                     + (this%c_Xh%dof%y(ip, j, k, e) - &
                     this%c_Xh%dof%y(im, j, k, e))**2 &
                     + (this%c_Xh%dof%z(ip, j, k, e) - &
                     this%c_Xh%dof%z(im, j, k, e))**2

                dj = (this%c_Xh%dof%x(i, jp, k, e) - &
                     this%c_Xh%dof%x(i, jm, k, e))**2 &
                     + (this%c_Xh%dof%y(i, jp, k, e) - &
                     this%c_Xh%dof%y(i, jm, k, e))**2 &
                     + (this%c_Xh%dof%z(i, jp, k, e) - &
                     this%c_Xh%dof%z(i, jm, k, e))**2

                dk = (this%c_Xh%dof%x(i, j, kp, e) - &
                     this%c_Xh%dof%x(i, j, km, e))**2 &
                     + (this%c_Xh%dof%y(i, j, kp, e) - &
                     this%c_Xh%dof%y(i, j, km, e))**2 &
                     + (this%c_Xh%dof%z(i, j, kp, e) - &
                     this%c_Xh%dof%z(i, j, km, e))**2

                di = sqrt(di) / (ip - im)
                dj = sqrt(dj) / (jp - jm)
                dk = sqrt(dk) / (kp - km)
                this%h%x(i,j,k,e) = (di * dj * dk)**(1.0_rp / 3.0_rp)

             end do
          end do
       end do
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%h%x, this%h%x_d, this%h%dof%size(),&
            HOST_TO_DEVICE, sync = .false.)
       call this%gs_Xh%op(this%h, GS_OP_ADD)
       call device_col2(this%h%x_d, this%c_Xh%mult_d, this%h%dof%size())
    else
       call this%gs_Xh%op(this%h, GS_OP_ADD)
       call col2(this%h%x, this%c_Xh%mult, this%h%dof%size())
    end if

  end subroutine compute_h

  !> Restart the simulation from saved state
  !! @param this The fluid scheme object
  !! @param dtlag Previous timestep sizes
  !! @param tlag Previous time values
  subroutine fluid_scheme_compressible_euler_restart(this, chkp)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
  end subroutine fluid_scheme_compressible_euler_restart

end module fluid_scheme_compressible_euler
