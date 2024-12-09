module fluid_scheme_compressible_euler
  use fluid_scheme_compressible, only: fluid_scheme_compressible_t
  use num_types, only : rp
  use mesh, only : mesh_t
  use json_module, only : json_file
  use user_intf, only : user_t
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) :: fluid_scheme_compressible_euler_t
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

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, scheme, user)

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
  end subroutine fluid_scheme_compressible_euler_step

  subroutine fluid_scheme_compressible_euler_restart(this, dtlag, tlag)
    class(fluid_scheme_compressible_euler_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
  end subroutine fluid_scheme_compressible_euler_restart

end module fluid_scheme_compressible_euler