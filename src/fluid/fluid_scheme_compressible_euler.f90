module fluid_scheme_compressible_euler
  use advection, only : advection_t, advection_factory
  use device, only : device_memcpy, HOST_TO_DEVICE
  use dofmap, only : dofmap_t
  use field_math, only : field_add2, field_cfill, field_cmult, field_cadd, field_copy, field_col2, &
                         field_col3, field_addcol3, field_sub2, field_invcol2
  use field, only : field_t
  use fluid_scheme_compressible, only: fluid_scheme_compressible_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use mesh, only : mesh_t
  use operators, only: div, grad
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use profiler, only : profiler_start_region, profiler_end_region
  use user_intf, only : user_t
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  use ax_product, only : ax_t, ax_helm_factory
  use field_list, only : field_list_t
  use coefs, only: coef_t
  use space, only : space_t
  use euler_residual, only: euler_rhs_t, euler_rhs_factory
  use neko_config, only : NEKO_BCKND_DEVICE
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) :: fluid_scheme_compressible_euler_t
     type(field_t) :: rho_res, m_x_res, m_y_res, m_z_res, m_E_res
     type(field_t) :: drho, dm_x, dm_y, dm_z, dE
     class(advection_t), allocatable :: adv
     class(ax_t), allocatable :: Ax
     class(euler_rhs_t), allocatable :: euler_rhs
   contains
     procedure, pass(this) :: init => fluid_scheme_compressible_euler_init
     procedure, pass(this) :: free => fluid_scheme_compressible_euler_free
     procedure, pass(this) :: step => fluid_scheme_compressible_euler_step
     procedure, pass(this) :: restart => fluid_scheme_compressible_euler_restart
     procedure, pass(this) :: rk4
  end type fluid_scheme_compressible_euler_t

contains
  subroutine fluid_scheme_compressible_euler_init(this, msh, lx, params, user, time_scheme)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    character(len=12), parameter :: scheme = 'compressible'
    logical :: advection

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

    end associate

    if (NEKO_BCKND_DEVICE .eq. 1) then
      associate(p => this%p, rho_field => this%rho_field, &
           u => this%u, v => this%v, w => this%w, &
           m_x => this%m_x, m_y => this%m_y, m_z => this%m_z)
        call device_memcpy(p%x, p%x_d, p%dof%size(), &
                           HOST_TO_DEVICE, sync = .false.)
        call device_memcpy(rho_field%x, rho_field%x_d, rho_field%dof%size(), &
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

  end subroutine fluid_scheme_compressible_euler_init

  subroutine fluid_scheme_compressible_euler_free(this)
    class(fluid_scheme_compressible_euler_t), intent(inout) :: this

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

  !> Advance fluid simulation in time.
  !! @param t The time value.
  !! @param tstep The current interation.
  !! @param dt The timestep
  !! @param ext_bdf Time integration logic.
  !! @param dt_controller timestep controller
  subroutine fluid_scheme_compressible_euler_step(this, t, tstep, dt, ext_bdf, dt_controller)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    type(field_t), pointer :: temp
    integer :: temp_indices(1)
    ! number of degrees of freedom
    integer :: n
    type(field_list_t) :: rhs_fields

    n = this%dm_Xh%size()
    call this%scratch%request_field(temp, temp_indices(1))

    call profiler_start_region('Fluid compressible', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
      m_x=> this%m_x, m_y => this%m_y, m_z => this%m_z, &
      Xh => this%Xh, msh => this%msh, Ax => this%Ax, &
      c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
      rho => this%rho, mu => this%mu, E => this%E, &
      rho_field => this%rho_field, mu_field => this%mu_field, &
      ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
      f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
      drho => this%drho, dm_x => this%dm_x, dm_y => this%dm_y, &
      dm_z => this%dm_z, dE => this%dE, &
      euler_rhs => this%euler_rhs)

      ! WIP: debugging setting m to (rho, 0, 0)
      ! call field_copy(m_x, rho_field, n)
      ! call field_col2(m_x, m_x, n) ! burgers
      ! call field_cfill(m_y, 0.0_rp, n)
      call field_cfill(m_z, 0.0_rp, n)

      call rhs_fields%init(5)
      call rhs_fields%assign(1, rho_field)
      call rhs_fields%assign(2, m_x)
      call rhs_fields%assign(3, m_y)
      call rhs_fields%assign(4, m_z)
      call rhs_fields%assign(5, E)

      ! call euler_rhs%compute(drho, dm_x, dm_y, dm_z, dE, &
      !     rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
      !     c_Xh, gs_Xh)

      call this%rk4(dt, rho_field, m_x, m_y, m_z, E, drho, dm_x, dm_y, dm_z, dE)

      !> TODO: apply boundary conditions

      ! Update variables
      ! Update u, v, w
      call field_copy(u, m_x, n)
      call field_invcol2(u, rho_field, n)
      call field_copy(v, m_y, n)
      call field_invcol2(v, rho_field, n)
      call field_copy(w, m_z, n)
      call field_invcol2(w, rho_field, n)

      ! Update p = (gamma - 1) * (E / rho - 0.5 * (u^2 + v^2 + w^2))
      call field_col3(temp, u, u, n)
      call field_addcol3(temp, v, v, n)
      call field_addcol3(temp, w, w, n)
      call field_cmult(temp, 0.5_rp, n)
      call field_copy(p, E, n)
      call field_invcol2(p, rho_field, n)
      call field_add2(p, temp, n)
      call field_cmult(p, this%gamma - 1.0_rp, n)

      !> TODO: Update maximum wave speed

      ! WIP: debugging visualizing rho
      ! call field_copy(p, rho_field, n);
      call field_copy(w, rho_field, n)

    end associate
    call profiler_end_region('Fluid compressible', 1)

    call this%scratch%relinquish_field(temp_indices)

  end subroutine fluid_scheme_compressible_euler_step

  subroutine rk4(this, dt, rho_field, m_x, m_y, m_z, E, rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: dt
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(inout) :: rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    integer :: n

    call this%euler_rhs%compute(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
          rho_field, m_x, m_y, m_z, E, this%p, this%u, this%v, this%w, this%Ax, &
          this%c_Xh, this%gs_Xh)

    n = this%dm_Xh%size()

    !> WIP: starting with forward Euler first

    !> rho = rho - dt * div(m)
    call field_cmult(rhs_rho_field, dt, n)
    call field_sub2(rho_field, rhs_rho_field, n)
    call this%gs_Xh%op(rho_field, GS_OP_MIN)
    !> m = m - dt * div(rho * u * u^T + p*I)
    call field_cmult(rhs_m_x, dt, n)
    call field_sub2(m_x, rhs_m_x, n)
    call this%gs_Xh%op(m_x, GS_OP_MIN)
    call field_cmult(rhs_m_y, dt, n)
    call field_sub2(m_y, rhs_m_y, n)
    call this%gs_Xh%op(m_y, GS_OP_MIN)
    call field_cmult(rhs_m_z, dt, n)
    call field_sub2(m_z, rhs_m_z, n)
    call this%gs_Xh%op(m_z, GS_OP_MIN)
    !> E = E - dt * div(u * (E + p))
    call field_cmult(rhs_E, dt, n)
    call field_sub2(E, rhs_E, n)
    call this%gs_Xh%op(E, GS_OP_MIN)

  end subroutine rk4

  subroutine fluid_scheme_compressible_euler_restart(this, dtlag, tlag)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
  end subroutine fluid_scheme_compressible_euler_restart

end module fluid_scheme_compressible_euler