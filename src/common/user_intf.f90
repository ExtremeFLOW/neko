! Copyright (c) 2020-2023, The Neko Authors
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
  use field
  use fluid_user_source_term
  use source_scalar
  use coefs
  use mesh
  use usr_inflow
  use usr_scalar
  use num_types
  use json_module, only : json_file
  implicit none
  private

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

  type, public :: user_t
     procedure(useric), nopass, pointer :: fluid_user_ic => null()
     procedure(user_initialize_modules), nopass, pointer :: user_init_modules => null()
     procedure(usermsh), nopass, pointer :: user_mesh_setup => null()
     procedure(usercheck), nopass, pointer :: user_check => null()
     procedure(user_final_modules), nopass, pointer :: user_finalize_modules => null()
     procedure(fluid_source_compute_pointwise), nopass, pointer :: fluid_user_f => null()
     procedure(fluid_source_compute_vector), nopass, pointer :: fluid_user_f_vector => null()
     procedure(source_scalar_term_pw), nopass, pointer :: scalar_user_f => null()
     procedure(source_scalar_term), nopass, pointer :: scalar_user_f_vector => null()
     procedure(usr_inflow_eval), nopass, pointer :: fluid_user_if => null()
     procedure(usr_scalar_bc_eval), nopass, pointer :: scalar_user_bc => null()
     !> Routine to set material properties
     procedure(user_material_properties), nopass, pointer :: material_properties => null()
   contains
     procedure, pass(u) :: init => user_intf_init
  end type user_t

  public :: useric, user_initialize_modules, usermsh, &
            dummy_user_material_properties, user_material_properties
contains
  
  !> User interface initialization
  subroutine user_intf_init(u)
    class(user_t), intent(inout) :: u

    if (.not. associated(u%fluid_user_ic)) then
       u%fluid_user_ic => dummy_user_ic
    end if

    if (.not. associated(u%fluid_user_f)) then
       u%fluid_user_f => dummy_user_f
    end if
    
    if (.not. associated(u%fluid_user_f_vector)) then
       u%fluid_user_f_vector => dummy_user_f_vector
    end if

    if (.not. associated(u%scalar_user_f)) then
       u%scalar_user_f => dummy_scalar_user_f
    end if
    
    if (.not. associated(u%scalar_user_f_vector)) then
       u%scalar_user_f_vector => dummy_user_scalar_f_vector
    end if

    if (.not. associated(u%scalar_user_bc)) then
       u%scalar_user_bc => dummy_scalar_user_bc
    end if
    
    if (.not. associated(u%user_mesh_setup)) then
       u%user_mesh_setup => dummy_user_mesh_setup
    end if

    if (.not. associated(u%user_check)) then
       u%user_check => dummy_user_check
    end if

    if (.not. associated(u%user_init_modules)) then
       u%user_init_modules => dummy_user_init_no_modules
    end if

    if (.not. associated(u%user_finalize_modules)) then
       u%user_finalize_modules => dummy_user_final_no_modules
    end if

    if (.not. associated(u%material_properties)) then
       u%material_properties => dummy_user_material_properties
    end if
  end subroutine user_intf_init

  
  !
  ! Below is the dummy user interface
  ! when running in pure turboNEKO mode
  !

  !> Dummy user initial condition
  subroutine dummy_user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    call neko_error('Dummy user defined initial condition set')    
  end subroutine dummy_user_ic

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
    class(source_scalar_t), intent(inout) :: f
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
  subroutine dummy_scalar_user_bc(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
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
    call neko_warning('Dummy scalar user bc set, applied on all non-labeled zones')    
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

  subroutine dummy_user_final_no_modules(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_final_no_modules

  subroutine dummy_user_material_properties(t, tstep, rho, mu, cp, lambda,&
                                            params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
  end subroutine dummy_user_material_properties

end module user_intf
