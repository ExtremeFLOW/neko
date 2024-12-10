module fluid_scheme_compressible_euler
  use advection, only : advection_t, advection_factory
  use field_math, only : field_add2, field_cfill, field_cmult, field_cadd, field_copy
  use field, only : field_t
  use fluid_scheme_compressible, only: fluid_scheme_compressible_t
  use num_types, only : rp
  use math, only: subcol3, copy, sub2
  use mesh, only : mesh_t
  use operators, only: conv1
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use profiler, only : profiler_start_region, profiler_end_region
  use user_intf, only : user_t
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) :: fluid_scheme_compressible_euler_t
     type(field_t) :: rho_res, m_x_res, m_y_res, m_z_res, m_E_res
     type(field_t) :: drho, dm_x, dm_y, dm_z, dE
     class(advection_t), allocatable :: adv
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

  end subroutine fluid_scheme_compressible_euler_init

  subroutine fluid_scheme_compressible_euler_free(this)
    class(fluid_scheme_compressible_euler_t), intent(inout) :: this

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
    ! number of degrees of freedom
    integer :: n

    n = this%dm_Xh%size()

    call profiler_start_region('Fluid compressible', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
      m_x=> this%m_x, m_y => this%m_y, m_z => this%m_z, &
      Xh => this%Xh, msh => this%msh, &
      c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
      rho => this%rho, mu => this%mu, &
      rho_field => this%rho_field, mu_field => this%mu_field, &
      f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
      ulag => this%ulag, vlag => this%vlag, wlag => this%wlag)

      ! WIP: debugging setting m_x to 1.0
      call field_cfill(m_x, 1.0_rp, n)
      
      ! rho = rho - dt * div(m)
      call conv1(f_x%x, rho_field%x, m_x%x, m_y%x, m_z%x, Xh, this%c_Xh)
      call field_cmult(f_x, dt, n)
      call sub2(rho_field%x, f_x%x, n)

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