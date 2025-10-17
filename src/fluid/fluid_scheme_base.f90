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
module fluid_scheme_base
  use coefs, only: coef_t
  use dirichlet, only : dirichlet_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_series, only : field_series_t
  use gather_scatter, only : gs_t
  use json_module, only : json_file
  use logger, only : LOG_SIZE
  use num_types, only : rp
  use checkpoint, only : chkp_t
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBL_LEN
  use space, only : space_t
  use time_scheme_controller, only : time_scheme_controller_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t, user_material_properties_intf
  use utils, only : neko_error
  use bc_list, only : bc_list_t
  use field_list, only : field_list_t
  use time_state, only: time_state_t
  implicit none
  private
  public :: fluid_scheme_base_t, fluid_scheme_base_factory

  !> Base type of all fluid formulations.
  type, abstract :: fluid_scheme_base_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name

     type(space_t) :: Xh !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh !< Coefficients associated with \f$ X_h \f$

     type(time_scheme_controller_t), allocatable :: ext_bdf

     !> The velocity field
     type(field_t), pointer :: u => null() !< x-component of Velocity
     type(field_t), pointer :: v => null() !< y-component of Velocity
     type(field_t), pointer :: w => null() !< z-component of Velocity
     type(field_t), pointer :: p => null() !< Pressure
     type(field_series_t) :: ulag, vlag, wlag !< fluid field (lag)

     !> Checkpoint
     type(chkp_t), pointer :: chkp => null()

     !> X-component of the right-hand side.
     type(field_t), pointer :: f_x => null()
     !> Y-component of the right-hand side.
     type(field_t), pointer :: f_y => null()
     !> Z-component of the right-hand side.
     type(field_t), pointer :: f_z => null()

     !> Boundary conditions
     ! List of boundary conditions for pressure
     type(bc_list_t) :: bcs_prs
     ! List of boundary conditions for velocity
     type(bc_list_t) :: bcs_vel

     type(json_file), pointer :: params !< Parameters
     type(mesh_t), pointer :: msh => null() !< Mesh

     !> Boundary condition labels (if any)
     character(len=NEKO_MSH_MAX_ZLBL_LEN), allocatable :: bc_labels(:)

     !> Density field
     type(field_t), pointer :: rho => null()

     !> The dynamic viscosity
     type(field_t), pointer :: mu => null()

     !> A helper that packs material properties to pass to the user routine.
     type(field_list_t) :: material_properties

     !> Is the fluid frozen at the moment
     logical :: freeze = .false.

     !> User material properties routine
     procedure(user_material_properties_intf), nopass, pointer :: &
          user_material_properties => null()

   contains
     !> Constructor
     procedure(fluid_scheme_base_init_intrf), pass(this), deferred :: init
     !> Destructor
     procedure(fluid_scheme_base_free_intrf), pass(this), deferred :: free
     !> Advance one step in time
     procedure(fluid_scheme_base_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint
     procedure(fluid_scheme_base_restart_intrf), pass(this), deferred :: restart
     ! Setup boundary conditions
     procedure(fluid_scheme_setup_bcs_intrf), pass(this), deferred :: setup_bcs

     !> Set the user inflow
     procedure(validate_intrf), pass(this), deferred :: validate
     !> Compute the CFL number
     procedure(fluid_scheme_base_compute_cfl_intrf), pass(this), deferred :: compute_cfl
     !> Set rho and mu
     procedure(update_material_properties), pass(this), deferred:: update_material_properties
  end type fluid_scheme_base_t

  !> Initialize all fields
  abstract interface
     subroutine fluid_base_init_all_intrf(this, msh, lx, params, kspv_init, &
          kspp_init, scheme, user)
       import fluid_scheme_base_t
       import mesh_t
       import json_file
       import user_t
       import rp
       import LOG_SIZE
       class(fluid_scheme_base_t), target, intent(inout) :: this
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
     subroutine fluid_base_init_common_intrf(this, msh, lx, params, scheme, &
          user, kspv_init)
       import fluid_scheme_base_t
       import mesh_t
       import json_file
       import user_t
       import dirichlet_t
       import LOG_SIZE
       import rp
       class(fluid_scheme_base_t), target, intent(inout) :: this
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
       import fluid_scheme_base_t
       class(fluid_scheme_base_t), intent(inout) :: this
     end subroutine fluid_base_free_intrf
  end interface

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_scheme_base_init_intrf(this, msh, lx, params, user, chkp)
       import fluid_scheme_base_t
       import json_file
       import mesh_t
       import user_t
       import chkp_t
       import time_scheme_controller_t
       class(fluid_scheme_base_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh
       integer, intent(in) :: lx
       type(json_file), target, intent(inout) :: params
       type(user_t), target, intent(in) :: user
       type(chkp_t), target, intent(inout) :: chkp
     end subroutine fluid_scheme_base_init_intrf
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine fluid_scheme_base_free_intrf(this)
       import fluid_scheme_base_t
       class(fluid_scheme_base_t), intent(inout) :: this
     end subroutine fluid_scheme_base_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine fluid_scheme_base_step_intrf(this, time, dt_controller)
       import time_state_t
       import fluid_scheme_base_t
       import time_step_controller_t
       class(fluid_scheme_base_t), target, intent(inout) :: this
       type(time_state_t), intent(in) :: time
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine fluid_scheme_base_step_intrf
  end interface

  !> Abstract interface to restart a fluid scheme
  abstract interface
     subroutine fluid_scheme_base_restart_intrf(this, chkp)
       import fluid_scheme_base_t
       import chkp_t
       class(fluid_scheme_base_t), target, intent(inout) :: this
       type(chkp_t), intent(inout) :: chkp
     end subroutine fluid_scheme_base_restart_intrf
  end interface

  !> Abstract interface to setup boundary conditions
  abstract interface
     subroutine fluid_scheme_setup_bcs_intrf(this, user, params)
       import fluid_scheme_base_t, user_t, json_file
       class(fluid_scheme_base_t), target, intent(inout) :: this
       type(user_t), target, intent(in) :: user
       type(json_file), intent(inout) :: params
     end subroutine fluid_scheme_setup_bcs_intrf
  end interface

  !> Abstract interface to validate the user inflow
  abstract interface
     subroutine validate_intrf(this)
       import fluid_scheme_base_t
       class(fluid_scheme_base_t), target, intent(inout) :: this
     end subroutine validate_intrf
  end interface

  !> Abstract interface to sets rho and mu
  abstract interface
     subroutine update_material_properties(this, time)
       import fluid_scheme_base_t, time_state_t
       class(fluid_scheme_base_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine update_material_properties
  end interface

  !> Compute the CFL number
  abstract interface
     function fluid_scheme_base_compute_cfl_intrf(this, dt) result(c)
       import fluid_scheme_base_t
       import rp
       class(fluid_scheme_base_t), intent(in) :: this
       real(kind=rp), intent(in) :: dt
       real(kind=rp) :: c
     end function fluid_scheme_base_compute_cfl_intrf
  end interface

  interface
     !> Initialise a fluid scheme
     module subroutine fluid_scheme_base_factory(object, type_name)
       class(fluid_scheme_base_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine fluid_scheme_base_factory
  end interface
end module fluid_scheme_base
