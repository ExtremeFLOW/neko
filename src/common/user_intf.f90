! Copyright (c) 2020-2024, The Neko Authors
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
  implicit none
  private

  !> Abstract interface for a user start-up routine
  abstract interface
     subroutine user_startup_intrf(params)
       import json_file
       type(json_file), intent(inout) :: params
     end subroutine user_startup_intrf
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
     subroutine user_initialize_modules(t, u, v, w, p, coef, params)
       import field_t
       import json_file
       import coef_t
       import rp
       real(kind=rp) :: t
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(coef_t), intent(inout) :: coef
       type(json_file), intent(inout) :: params
     end subroutine user_initialize_modules
  end interface

  !> Abstract interface for adding user defined simulation components
  abstract interface
     subroutine user_simcomp_init(params)
       import json_file
       type(json_file), intent(inout) :: params
     end subroutine user_simcomp_init
  end interface

  !> Abstract interface for user defined mesh deformation functions
  abstract interface
     subroutine usermsh(msh)
       import mesh_t
       type(mesh_t), intent(inout) :: msh
     end subroutine usermsh
  end interface

  !> Abstract interface for user defined check functions
  abstract interface
     subroutine usercheck(t, tstep, u, v, w, p, coef, param)
       import field_t
       import coef_t
       import json_file
       import rp
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(coef_t), intent(inout) :: coef
       type(json_file), intent(inout) :: param
     end subroutine usercheck
  end interface

  !> Abstract interface for finalizating user variables
  abstract interface
     subroutine user_final_modules(t, param)
       import json_file
       import rp
       real(kind=rp) :: t
       type(json_file), intent(inout) :: param
     end subroutine user_final_modules
  end interface

  !> Abstract interface for setting material properties.
  !! @param t Time value.
  !! @param tstep Current time step.
  !! @param rho Fluid density.
  !! @param mu Fluid dynamic viscosity.
  !! @param cp Scalar specific heat capacity.
  !! @param lambda Scalar thermal conductivity.
  abstract interface
     subroutine user_material_properties(t, tstep, rho, mu, cp, lambda, params)
       import rp
       import json_file
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       real(kind=rp), intent(inout) :: rho, mu, cp, lambda
       type(json_file), intent(inout) :: params
     end subroutine user_material_properties
  end interface

  !> A type collecting all the overridable user routines.
  type, public :: user_t
     !> Run as soon as the case file is read, with nothing else initialized.
     !! Use to manipulate the case file, and define custom parameters.
     procedure(user_startup_intrf), nopass, pointer :: &
          user_startup => null()
     !> Run after the entire case is initialized and restarted, but before the
     !! time loop. Good place to create auxillary fields, etc.
     procedure(user_initialize_modules), nopass, pointer :: &
          user_init_modules => null()
     !> Compute user initial conditions for the incompressible fluid.
     procedure(useric), nopass, pointer :: fluid_user_ic => null()
     !> Compute user initial conditions for the compressible fluid.
     procedure(useric_compressible), nopass, pointer :: &
          fluid_compressible_user_ic => null()
     !> Compute user initial conditions for the scalar.
     procedure(useric_scalar), nopass, pointer :: scalar_user_ic => null()
     !> Constructor for the user simcomp. Ran in the constructor of
     !! neko_simcomps.
     procedure(user_simcomp_init), nopass, pointer :: &
          init_user_simcomp => null()
     !> Run right after reading the mesh and allows to manipulate it.
     procedure(usermsh), nopass, pointer :: user_mesh_setup => null()
     !> Run at the end of each time-step in the time loop, right before field
     !! output to disk.
     procedure(usercheck), nopass, pointer :: user_check => null()
     !> Runs in the end of the simulation, after the last output. Mean as a
     !! place to run `free()` on user-allocated objects.
     procedure(user_final_modules), nopass, pointer :: &
          user_finalize_modules => null()
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
     !> User boundary condition for the fluid, pointwise interface.
     procedure(usr_inflow_eval), nopass, pointer :: fluid_user_if => null()
     !> User boundary condition for the scalar, pointwise interface.
     procedure(usr_scalar_bc_eval), nopass, pointer :: scalar_user_bc => null()
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
       user_initialize_modules, usermsh, dummy_user_material_properties, &
       user_material_properties, user_simcomp_init, &
       simulation_component_user_settings, user_startup_intrf
contains

  !> Constructor.
  subroutine user_intf_init(this)
    class(user_t), intent(inout) :: this
    logical :: user_extended = .false.
    character(len=256), dimension(14) :: extensions
    integer :: i, n

    n = 0
    if (.not. associated(this%user_startup)) then
       this%user_startup => dummy_user_startup
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

    if (.not. associated(this%scalar_user_bc)) then
       this%scalar_user_bc => dummy_scalar_user_bc
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Scalar boundary condition'
    end if

    if (.not. associated(this%user_dirichlet_update)) then
       this%user_dirichlet_update => dirichlet_do_nothing
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Dirichlet boundary condition'
    end if

    if (.not. associated(this%user_mesh_setup)) then
       this%user_mesh_setup => dummy_user_mesh_setup
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Mesh setup'
    end if

    if (.not. associated(this%user_check)) then
       this%user_check => dummy_user_check
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- User check'
    end if

    if (.not. associated(this%user_init_modules)) then
       this%user_init_modules => dummy_user_init_no_modules
    else
       user_extended = .true.
       n = n + 1
       write(extensions(n), '(A)') '- Initialize modules'
    end if

    if (.not. associated(this%init_user_simcomp)) then
       this%init_user_simcomp => dummy_user_init_no_simcomp
    end if

    if (.not. associated(this%user_finalize_modules)) then
       this%user_finalize_modules => dummy_user_final_no_modules
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
  subroutine dummy_user_startup(params)
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_startup

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
  subroutine dummy_user_scalar_f_vector(f, t)
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined vector valued forcing set')
  end subroutine dummy_user_scalar_f_vector

  !> Dummy user (scalar) forcing
  subroutine dummy_scalar_user_f(s, j, k, l, e, t)
    real(kind=rp), intent(inout) :: s
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: t
    call neko_error('Dummy user defined forcing set')
  end subroutine dummy_scalar_user_f

  !> Dummy user boundary condition for scalar
  subroutine dummy_scalar_user_bc(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, &
       tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    call neko_warning('Dummy scalar user bc set, applied on all' // &
         ' non-labeled zones')
  end subroutine dummy_scalar_user_bc

  !> Dummy user mesh apply
  subroutine dummy_user_mesh_setup(msh)
    type(mesh_t), intent(inout) :: msh
  end subroutine dummy_user_mesh_setup

  !> Dummy user check
  subroutine dummy_user_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_check

  subroutine dummy_user_init_no_modules(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_init_no_modules

  subroutine dummy_user_init_no_simcomp(params)
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_init_no_simcomp

  subroutine dummy_user_final_no_modules(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_final_no_modules

  subroutine dirichlet_do_nothing(dirichlet_field_list, dirichlet_bc, &
       coef, t, tstep)
    type(field_list_t), intent(inout) :: dirichlet_field_list
    type(field_dirichlet_t), intent(in) :: dirichlet_bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
  end subroutine dirichlet_do_nothing

  subroutine dummy_user_material_properties(t, tstep, rho, mu, cp, lambda,&
       params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_material_properties

  ! ========================================================================== !
  ! Helper functions for user defined interfaces

  !> JSON extraction helper function for simulation components
  !! @param name The name of the object to be created.
  !! @param params The JSON object containing the user-defined component.
  !! @return The JSON object for initializing the simulation component.
  function simulation_component_user_settings(name, params) result(comp_subdict)
    character(len=*), intent(in) :: name
    type(json_file), intent(inout) :: params
    type(json_file) :: comp_subdict

    character(len=:), allocatable :: current_type
    integer :: n_simcomps
    integer :: i
    logical :: found, is_user

    call params%info('', n_children = n_simcomps)

    found = .false.
    do i = 1, n_simcomps
       call json_extract_item(params, "", i, comp_subdict)
       call json_get_or_default(comp_subdict, "is_user", is_user, .false.)
       if (.not. is_user) cycle

       call json_get(comp_subdict, "type", current_type)
       if (trim(current_type) .eq. trim(name)) then
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       call neko_error("User-defined simulation component " &
            // trim(name) // " not found in case file.")
    end if

  end function simulation_component_user_settings


  !> @example simulation_components/user_simcomp.f90
  !! @brief User defined simulation components.
  !! @details
  !! Example of how to use the simcomp_executor to add a user defined
  !! simulation component to the list.
  !! @include simulation_components/user_simcomp.case
end module user_intf
