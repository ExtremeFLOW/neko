! Copyright (c) 2020-2025, The Neko Authors
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
!> Interfaces for user interaction with NEKO
module user_intf
  use field, only : field_t
  use field_list, only : field_list_t
  use fluid_user_source_term, only : fluid_user_source_term_t, &
       fluid_source_compute_pointwise, fluid_source_compute_vector
  use scalar_user_source_term, only : scalar_user_source_term_t, &
       scalar_source_compute_pointwise, scalar_source_compute_vector
  use coefs, only : coef_t
  use bc_list, only : bc_list_t
  use mesh, only : mesh_t
  use usr_inflow, only : usr_inflow_t, usr_inflow_eval
  use usr_scalar, only : usr_scalar_t, usr_scalar_bc_eval
  use field_dirichlet, only: field_dirichlet_update
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_extract_item, json_get, json_get_or_default
  use utils, only : neko_error, neko_warning
  use logger, only : neko_log
  use bc, only : bc_t
  use field_dirichlet, only : field_dirichlet_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Abstract interface for a user start-up routine
  abstract interface
     subroutine user_startup(params)
       import json_file
       type(json_file), intent(inout) :: params
     end subroutine user_startup
  end interface

  !> Abstract interface for user defined initial conditions
  abstract interface
     subroutine useric(u, v, w, p, params)
       import field_t
       import json_file
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(json_file), intent(inout) :: params
     end subroutine useric
  end interface

  !> Abstract interface for user defined initial conditions
  abstract interface
     subroutine useric_compressible(rho, u, v, w, p, params)
       import field_t
       import json_file
       type(field_t), intent(inout) :: rho
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(json_file), intent(inout) :: params
     end subroutine useric_compressible
  end interface

  !> Abstract interface for user defined scalar initial conditions
  abstract interface
     subroutine useric_scalar(s, params)
       import field_t
       import json_file
       type(field_t), intent(inout) :: s
       type(json_file), intent(inout) :: params
     end subroutine useric_scalar
  end interface

  !> Abstract interface for initilialization of modules
  abstract interface
     subroutine user_initialize(time)
       import time_state_t
       type(time_state_t), intent(in) :: time
     end subroutine user_initialize
  end interface

  !> Abstract interface for user defined mesh deformation functions
  abstract interface
     subroutine user_mesh_setup(time, msh)
       import mesh_t
       import time_state_t
       type(time_state_t), intent(in) :: time
       type(mesh_t), intent(inout) :: msh
     end subroutine user_mesh_setup
  end interface

  !> Abstract interface for user defined check functions
  abstract interface
     subroutine user_compute(time)
       import time_state_t
       type(time_state_t), intent(in) :: time
     end subroutine user_compute
  end interface

  !> Abstract interface for finalizating user variables
  abstract interface
     subroutine user_finalize(time)
       import json_file
       import time_state_t
       type(time_state_t), intent(in) :: time
     end subroutine user_finalize
  end interface

  !> Abstract interface for setting material properties.
  !! @param t Time value.
  !! @param tstep Current time step.
  !! @param name The name of the solver calling the routine. By default
  !! "fluid" or "scalar"
  !! @param properties Array of properties, defined by convention for each
  !! scheme.
  !! @param params The JSON configuration of the scheme.
  abstract interface
     subroutine user_material_properties(t, tstep, name, properties)
       import rp, field_list_t
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       character(len=*), intent(in) :: name
       type(field_list_t), intent(inout) :: properties
     end subroutine user_material_properties
  end interface

  !> A type collecting all the overridable user routines and flag to suppress
  !! type injection from custom modules.
  type, public :: user_t
     !> Setting this to true in the user_setup routine in the user file will
     !! suppress custom modules registering their types in the factories. So you
     !! have to take care of type injection in `user_startup`. Use if you really
     !! want full control over type injection for some reason.
     logical :: suppress_type_injection = .false.
     !> Run as soon as the case file is read, with nothing else initialized.
     !! Use to manipulate the case file, and define custom parameters.
     procedure(user_startup), nopass, pointer :: startup => null()
     !> Run after the entire case is initialized and restarted, but before the
     !! time loop. Good place to create auxillary fields, etc.
     procedure(user_initialize), nopass, pointer :: initialize => null()
     !> Compute user initial conditions for the incompressible fluid.
     procedure(useric), nopass, pointer :: fluid_user_ic => null()
     !> Compute user initial conditions for the compressible fluid.
     procedure(useric_compressible), nopass, pointer :: &
          fluid_compressible_user_ic => null()
     !> Compute user initial conditions for the scalar.
     procedure(useric_scalar), nopass, pointer :: scalar_user_ic => null()
     !> Run right after reading the mesh and allows to manipulate it.
     procedure(user_mesh_setup), nopass, pointer :: mesh_setup => null()
     !> Run at the end of each time-step in the time loop, right before field
     !! output to disk.
     procedure(user_compute), nopass, pointer :: compute => null()
     !> Runs in the end of the simulation, after the last output. Mean as a
     !! place to run `free()` on user-allocated objects.
     procedure(user_finalize), nopass, pointer :: &
          finalize => null()
     !> User forcing for the fluid, pointwise interface.
     procedure(fluid_source_compute_pointwise), nopass, pointer :: &
          fluid_user_f => null()
     !> User forcing for the fluid, field (vector) interface.
     procedure(fluid_source_compute_vector), nopass, pointer :: &
          fluid_user_f_vector => null()
     !> User forcing for the scalar, pointwise interface.
     procedure(scalar_source_compute_pointwise), nopass, pointer :: &
          scalar_user_f => null()
     !> User forcing for the scalar, field (vector) interface.
     procedure(scalar_source_compute_vector), nopass, pointer :: &
          scalar_user_f_vector => null()
     !> User boundary condition for the fluid or the scalar, field interface
     !! (much more powerful than pointwise in terms of what can be done).
     procedure(field_dirichlet_update), nopass, pointer :: &
          user_dirichlet_update => null()
     !> Routine to set material properties.
     procedure(user_material_properties), nopass, pointer :: &
          material_properties => null()
   contains
     !> Constructor that points non-associated routines to dummy ones.
     !! Calling a dummy routine causes an error in most cases, but sometimes
     !! the dummy routine just does nothing. E.g., the dummmy `user_startup`
     !! just does nothing. The association of the actual routines defined by
     !! the user happens in the `user_setup` routine inside the user file.
     !! A call to this routine is injected into the neko executable by
     !! `makeneko`.
     procedure, pass(this) :: init => user_intf_init
  end type user_t

  public :: useric, useric_scalar, useric_compressible, &
       user_initialize, user_mesh_setup, dummy_user_material_properties, &
       user_material_properties, user_startup
contains

  !> Constructor.
  subroutine user_intf_init(this)
    class(user_t), intent(inout) :: this
    logical :: user_extended = .false.
    character(len=256), dimension(14) :: extensions
    integer :: i, n

    n = 0
    if (.not. associated(this%startup)) then
       this%startup => dummy_startup
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Startup'
    end if

    if (.not. associated(this%fluid_user_ic)) then
       this%fluid_user_ic => dummy_user_ic
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Fluid initial condition'
    end if

    if (.not. associated(this%scalar_user_ic)) then
       this%scalar_user_ic => dummy_user_ic_scalar
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Scalar initial condition'
    end if

    if (.not. associated(this%fluid_compressible_user_ic)) then
       this%fluid_compressible_user_ic => dummy_user_ic_compressible
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Compressible fluid initial condition'
    end if

    if (.not. associated(this%fluid_user_f)) then
       this%fluid_user_f => dummy_user_f
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Fluid source term'
    end if

    if (.not. associated(this%fluid_user_f_vector)) then
       this%fluid_user_f_vector => dummy_user_f_vector
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Fluid source term vector'
    end if

    if (.not. associated(this%scalar_user_f)) then
       this%scalar_user_f => dummy_scalar_user_f
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Scalar source term'
    end if

    if (.not. associated(this%scalar_user_f_vector)) then
       this%scalar_user_f_vector => dummy_user_scalar_f_vector
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Scalar source term vector'
    end if

    if (.not. associated(this%user_dirichlet_update)) then
       this%user_dirichlet_update => dirichlet_do_nothing
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Dirichlet boundary condition'
    end if

    if (.not. associated(this%mesh_setup)) then
       this%mesh_setup => dummy_user_mesh_setup
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Mesh setup'
    end if

    if (.not. associated(this%compute)) then
       this%compute => dummy_user_check
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- User check'
    end if

    if (.not. associated(this%initialize)) then
       this%initialize => dummy_initialize
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Initialize modules'
    end if

    if (.not. associated(this%finalize)) then
       this%finalize => dummy_user_finalize
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Finalize modules'
    end if

    if (.not. associated(this%material_properties)) then
       this%material_properties => dummy_user_material_properties
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Material properties'
    end if

    if (user_extended) then
       call neko_log%section('User defined extensions')

       do i = 1, n
          call neko_log%message(extensions(i))
       end do

       call neko_log%end_section()
    end if

  end subroutine user_intf_init


  !
  ! Below is the dummy user interface
  ! when running in pure turboNEKO mode
  !

  !> Dummy user startup
  subroutine dummy_startup(params)
    type(json_file), intent(inout) :: params
  end subroutine dummy_startup

  !> Dummy user initial condition
  subroutine dummy_user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    call neko_error('Dummy user defined initial condition set')
  end subroutine dummy_user_ic

  !> Dummy user initial condition
  subroutine dummy_user_ic_compressible(rho, u, v, w, p, params)
    type(field_t), intent(inout) :: rho
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    call neko_error('Dummy user defined initial condition set')
  end subroutine dummy_user_ic_compressible

  !> Dummy user initial condition for scalar field
  !! @param s Scalar field.
  !! @param params JSON parameters.
  subroutine dummy_user_ic_scalar(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    call neko_error('Dummy user defined scalar initial condition set')
  end subroutine dummy_user_ic_scalar

  subroutine dummy_user_ic_scalars(s, field_name, params)
    type(field_t), intent(inout) :: s
    character(len=*), intent(in) :: field_name
    type(json_file), intent(inout) :: params
    call neko_warning('Dummy multiple scalar initial condition called')
  end subroutine dummy_user_ic_scalars

  !> Dummy user (fluid) forcing
  subroutine dummy_user_f_vector(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined vector valued forcing set')
  end subroutine dummy_user_f_vector

  !> Dummy user (fluid) forcing
  subroutine dummy_user_f(u, v, w, j, k, l, e, t)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined forcing set')
  end subroutine dummy_user_f

  !> Dummy user (scalar) forcing
  subroutine dummy_user_scalar_f_vector(field_name, f, t)
    character(len=*), intent(in) :: field_name
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined vector valued forcing set')
  end subroutine dummy_user_scalar_f_vector

  !> Dummy user (scalar) forcing
  subroutine dummy_scalar_user_f(field_name, s, j, k, l, e, t)
    character(len=*), intent(in) :: field_name
    real(kind=rp), intent(inout) :: s
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined forcing set')
  end subroutine dummy_scalar_user_f

  !> Dummy user mesh apply
  subroutine dummy_user_mesh_setup(time, msh)
    type(time_state_t), intent(in) :: time
    type(mesh_t), intent(inout) :: msh
  end subroutine dummy_user_mesh_setup

  !> Dummy user compute
  subroutine dummy_user_check(time)
    type(time_state_t), intent(in) :: time
  end subroutine dummy_user_check

  subroutine dummy_initialize(time)
      type(time_state_t), intent(in) :: time
  end subroutine dummy_initialize

  subroutine dummy_user_init_no_simcomp(params)
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_init_no_simcomp

  subroutine dummy_user_finalize(time)
    type(time_state_t), intent(in) :: time
  end subroutine dummy_user_finalize

  subroutine dirichlet_do_nothing(dirichlet_field_list, dirichlet_bc, &
       coef, t, tstep)
    type(field_list_t), intent(inout) :: dirichlet_field_list
    type(field_dirichlet_t), intent(in) :: dirichlet_bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
  end subroutine dirichlet_do_nothing

  subroutine dummy_user_material_properties(t, tstep, name, properties)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    character(len=*), intent(in) :: name
    type(field_list_t), intent(inout) :: properties
  end subroutine dummy_user_material_properties


end module user_intf
