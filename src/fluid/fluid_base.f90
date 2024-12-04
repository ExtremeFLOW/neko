module fluid_base
  use bc, only : bc_t, bc_list_t
  use checkpoint, only : chkp_t
  use coefs, only: coef_t
  use dirichlet, only : dirichlet_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_series, only : field_series_t
  use gather_scatter, only : gs_t
  use json_module, only : json_file
  use logger, only : LOG_SIZE
  use num_types, only : rp
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBL_LEN
  use space, only : space_t, GLL
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t
  use wall, only : no_slip_wall_t
  implicit none

  !> Base type of all fluid formulations.
  type, abstract :: fluid_base_t
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$

     !> The velocity field
     type(field_t), pointer :: u => null()    !< x-component of Velocity
     type(field_t), pointer :: v => null()    !< y-component of Velocity
     type(field_t), pointer :: w => null()    !< z-component of Velocity
     type(field_t), pointer :: p => null()    !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)

     !> The variable density field
     type(field_t) :: rho_field

     !> Boundary conditions
     type(field_t) :: bdry                     !< Boundary markings
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     class(bc_t), allocatable :: bc_inflow     !< Dirichlet inflow for velocity
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_vel_neumann      !< List of neumann velocity conditions

     type(json_file), pointer :: params        !< Parameters
     type(mesh_t), pointer :: msh => null()    !< Mesh
     type(chkp_t) :: chkp                      !< Checkpoint

     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)

     !> Dynamic viscosity
     real(kind=rp) :: mu

     !> The variable mu field
     type(field_t) :: mu_field
     
   contains
     !> Constructor
     procedure(fluid_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor
     procedure(fluid_scheme_free_intrf), pass(this), deferred :: free
     !> Advance one step in time
     procedure(fluid_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint
     procedure(fluid_scheme_restart_intrf), pass(this), deferred :: restart

     !> Compute the CFL number
     procedure(fluid_base_compute_cfl_intrf), pass(this), deferred :: compute_cfl
  end type fluid_base_t

  !> Initialize all fields
  abstract interface
     subroutine fluid_base_init_all_intrf(this, msh, lx, params, kspv_init, &
                                          kspp_init, scheme, user)
       import fluid_base_t
       import mesh_t
       import json_file
       import user_t
       import rp
       import LOG_SIZE
       class(fluid_base_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       logical :: kspv_init
       logical :: kspp_init
       character(len=*), intent(in) :: scheme
       real(kind=rp) :: abs_tol
       integer :: integer_val, ierr
       logical :: logical_val
       character(len=:), allocatable :: solver_type, precon_type
       character(len=LOG_SIZE) :: log_buf
       real(kind=rp) :: GJP_param_a, GJP_param_b
     end subroutine fluid_base_init_all_intrf
  end interface

  !> Initialize common data for the current scheme
  abstract interface
     subroutine fluid_base_init_common_intrf(this, msh, lx, params, scheme, user, &
                                             kspv_init)
       import fluid_base_t
       import mesh_t
       import json_file
       import user_t
       import dirichlet_t
       import LOG_SIZE
       import rp
       class(fluid_base_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       character(len=*), intent(in) :: scheme
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       logical, intent(in) :: kspv_init
       type(dirichlet_t) :: bdry_mask
       character(len=LOG_SIZE) :: log_buf
       real(kind=rp), allocatable :: real_vec(:)
       real(kind=rp) :: real_val, kappa, B, z0
       logical :: logical_val
       integer :: integer_val, ierr
       type(json_file) :: wm_json
       character(len=:), allocatable :: string_val1, string_val2
     end subroutine fluid_base_init_common_intrf
  end interface

  !> Deallocate a fluid formulation
  abstract interface
     subroutine fluid_base_free_intrf(this)
       import fluid_base_t
       class(fluid_base_t), intent(inout) :: this
     end subroutine fluid_base_free_intrf
  end interface

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_scheme_init_intrf(this, msh, lx, params, user, &
          time_scheme)
       import fluid_base_t
       import json_file
       import mesh_t
       import user_t
       import time_scheme_controller_t
       class(fluid_base_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(inout) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
     end subroutine fluid_scheme_init_intrf
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine fluid_scheme_free_intrf(this)
       import fluid_base_t
       class(fluid_base_t), intent(inout) :: this
     end subroutine fluid_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine fluid_scheme_step_intrf(this, t, tstep, dt, ext_bdf, &
                                        dt_controller)
       import fluid_base_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(fluid_base_t), target, intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(inout) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine fluid_scheme_step_intrf
  end interface

  !> Abstract interface to restart a fluid scheme
  abstract interface
     subroutine fluid_scheme_restart_intrf(this, dtlag, tlag)
       import fluid_base_t
       import rp
       class(fluid_base_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)

     end subroutine fluid_scheme_restart_intrf
  end interface

  !> Compute the CFL number
  abstract interface
     function fluid_base_compute_cfl_intrf(this, dt) result(c)
      import fluid_base_t
      import rp
      class(fluid_base_t), intent(in) :: this
      real(kind=rp), intent(in) :: dt
      real(kind=rp) :: c
     end function fluid_base_compute_cfl_intrf
  end interface

  public :: fluid_base_t
end module fluid_base