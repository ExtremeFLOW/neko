module fluid_scheme_compressible_euler
  use advection, only : advection_t, advection_factory
  use field_math, only : field_add2, field_cfill, field_cmult, field_cadd, field_copy, field_col2, field_col3, field_addcol3
  use field, only : field_t
  use fluid_scheme_compressible, only: fluid_scheme_compressible_t
  use gs_ops, only : GS_OP_ADD
  use num_types, only : rp
  use math, only: subcol3, copy, sub2, add2, add3, col2, col3, addcol3, cmult, cfill, invcol3
  use mesh, only : mesh_t
  use operators, only: div, grad
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use profiler, only : profiler_start_region, profiler_end_region
  use user_intf, only : user_t
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  use ax_product, only : ax_t, ax_helm_factory
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) :: fluid_scheme_compressible_euler_t
     type(field_t) :: rho_res, m_x_res, m_y_res, m_z_res, m_E_res
     type(field_t) :: drho, dm_x, dm_y, dm_z, dE
     class(advection_t), allocatable :: adv
     class(ax_t), allocatable :: Ax
   contains
     procedure, pass(this) :: init => fluid_scheme_compressible_euler_init
     procedure, pass(this) :: free => fluid_scheme_compressible_euler_free
     procedure, pass(this) :: step => fluid_scheme_compressible_euler_step
     procedure, pass(this) :: restart => fluid_scheme_compressible_euler_restart
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

    ! Initialize the advection factory
    call json_get_or_default(params, 'case.fluid.advection', advection, .true.)
    call advection_factory(this%adv, params, this%c_Xh, &
                           this%ulag, this%vlag, this%wlag, &
                           this%chkp%dtlag, this%chkp%tlag, time_scheme, &
                           .not. advection)

    ! Initialize the diffusion operator
    call ax_helm_factory(this%Ax, full_formulation = .false.)

  end subroutine fluid_scheme_compressible_euler_init

  subroutine fluid_scheme_compressible_euler_free(this)
    class(fluid_scheme_compressible_euler_t), intent(inout) :: this

    if (allocated(this%Ax)) then
      deallocate(this%Ax)
    end if

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
    real(kind=rp), allocatable :: temp(:)
    ! number of degrees of freedom
    integer :: n
    real(kind=rp) :: h, c_avisc

    h = 0.03_rp / 5.0_rp ! grid size / polynomial degree
    c_avisc = 1.0_rp*h

    n = this%dm_Xh%size()
    allocate(temp(n))

    call profiler_start_region('Fluid compressible', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
      m_x=> this%m_x, m_y => this%m_y, m_z => this%m_z, &
      Xh => this%Xh, msh => this%msh, Ax => this%Ax, &
      c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
      rho => this%rho, mu => this%mu, E => this%E, &
      rho_field => this%rho_field, mu_field => this%mu_field, &
      f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
      ulag => this%ulag, vlag => this%vlag, wlag => this%wlag)

      ! WIP: debugging setting m to (rho, 0, 0)
      call field_copy(m_x, rho_field, n)
      call field_col2(m_x, m_x, n) ! burgers
      call field_cfill(m_y, 0.0_rp, n)
      call field_cfill(m_z, 0.0_rp, n)
      
      !> rho = rho - dt * div(m)
      call div(temp, m_x%x, m_y%x, m_z%x, c_Xh)
      call cmult(temp, dt, n)
      call sub2(rho_field%x, temp, n)
      ! artificial diffusion for rho
      call Ax%compute(temp, rho_field%x, c_Xh, msh, Xh)
      call gs_Xh%op(temp, n, GS_OP_ADD)
      call col2(temp, c_Xh%Binv, n)
      call cmult(temp, c_avisc, n) ! first-order viscosity
      call cmult(temp, dt, n)
      call sub2(rho_field%x, temp, n)

      ! m = m - dt * div(rho * u * u^T + p*I)
      !> m_x
      call copy(f_x%x, p%x, n)
      call addcol3(f_x%x, m_x%x, u%x, n)
      call col3(f_y%x, m_x%x, v%x, n)
      call col3(f_z%x, m_x%x, w%x, n)
      call div(temp, f_x%x, f_y%x, f_z%x, this%c_Xh)
      call cmult(temp, dt, n)
      call sub2(m_x%x, temp, n)
      !> m_y
      call col3(f_x%x, m_y%x, u%x, n)
      call copy(f_y%x, p%x, n)
      call addcol3(f_y%x, m_y%x, v%x, n)
      call col3(f_z%x, m_y%x, w%x, n)
      call div(temp, f_x%x, f_y%x, f_z%x, this%c_Xh)
      call cmult(temp, dt, n)
      call sub2(m_y%x, temp, n)
      ! m_z
      call col3(f_x%x, m_z%x, u%x, n)
      call col3(f_y%x, m_z%x, v%x, n)
      call copy(f_z%x, p%x, n)
      call addcol3(f_z%x, m_z%x, w%x, n)
      call div(temp, f_x%x, f_y%x, f_z%x, this%c_Xh)
      call cmult(temp, dt, n)
      call sub2(m_z%x, temp, n)

      ! E = E - dt * div(u * (E + p))
      call add3(temp, E%x, p%x, n)
      call col3(f_x%x, u%x, temp, n)
      call col3(f_y%x, v%x, temp, n)
      call col3(f_z%x, w%x, temp, n)
      call div(temp, f_x%x, f_y%x, f_z%x, this%c_Xh)
      call cmult(temp, dt, n)
      call sub2(E%x, temp, n)

      !> TODO: apply boundary conditions

      ! Update variables
      ! Update u, v, w
      call invcol3(u%x, m_x%x, rho_field%x, n)
      call invcol3(v%x, m_y%x, rho_field%x, n)
      call invcol3(w%x, m_z%x, rho_field%x, n)

      ! Update p = (gamma - 1) * (E / rho - 0.5 * (u^2 + v^2 + w^2))
      call col3(temp, u%x, u%x, n)
      call addcol3(temp, v%x, v%x, n)
      call addcol3(temp, w%x, w%x, n)
      call cmult(temp, 0.5_rp, n)
      call invcol3(p%x, E%x, rho_field%x, n)
      call add2(p%x, temp, n)
      call cmult(p%x, this%gamma - 1.0_rp, n)

      !> TODO: Update maximum wave speed

      ! WIP: debugging visualizing rho
      call field_copy(p, rho_field, n);

    end associate
    call profiler_end_region('Fluid compressible', 1)

  end subroutine fluid_scheme_compressible_euler_step

  subroutine fluid_scheme_compressible_euler_restart(this, dtlag, tlag)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
  end subroutine fluid_scheme_compressible_euler_restart

end module fluid_scheme_compressible_euler