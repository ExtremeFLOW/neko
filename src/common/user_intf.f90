! Copyright (c) 2020-2021, The Neko Authors
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
  use source
  use coefs
  use usr_inflow
  use parameters
  use num_types
  implicit none

  !> Abstract interface for user defined initial conditions
  abstract interface
     subroutine useric(u, v, w, p, params)
       import field_t
       import param_t
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(param_t), intent(inout) :: params
     end subroutine useric
  end interface

  !> Abstract interface for initilialization of modules
  abstract interface
     subroutine user_initialize_modules(t, u, v, w, p, c_Xh, params)
       import field_t
       import param_t
       import coef_t
       import rp
       real(kind=rp) :: t
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
       type(coef_t), intent(inout) :: c_Xh
       type(param_t), intent(inout) :: params
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
     subroutine usercheck(t, dt, tstep, u, v, w, p, coef)
       import field_t
       import coef_t
       import rp
       real(kind=rp), intent(in) :: t, dt
       integer, intent(in) :: tstep
       type(coef_t), intent(inout) :: coef
       type(field_t), intent(inout) :: u
       type(field_t), intent(inout) :: v
       type(field_t), intent(inout) :: w
       type(field_t), intent(inout) :: p
     end subroutine usercheck
  end interface

  type :: user_t
     procedure(useric), nopass, pointer :: fluid_usr_ic => null()
     procedure(user_initialize_modules), nopass, pointer :: user_init_modules => null()
     procedure(usermsh), nopass, pointer :: usr_msh_setup => null()
     procedure(usercheck), nopass, pointer :: usr_chk => null()
     procedure(source_term_pw), nopass, pointer :: fluid_usr_f => null()
     procedure(source_term), nopass, pointer :: fluid_usr_f_vector => null()
     procedure(usr_inflow_eval), nopass, pointer :: fluid_usr_if => null()
   contains
     procedure, pass(u) :: init => user_intf_init
  end type user_t
  
contains
  !> User interface initialization
  subroutine user_intf_init(u)
    class(user_t), intent(inout) :: u

    if (.not. associated(u%fluid_usr_ic)) then
       u%fluid_usr_ic => dummy_user_ic
    end if

    if (.not. associated(u%fluid_usr_f)) then
       u%fluid_usr_f => dummy_user_f
    end if
    if (.not. associated(u%fluid_usr_f_vector)) then
       u%fluid_usr_f_vector => dummy_user_f_vector
    end if

    if (.not. associated(u%usr_msh_setup)) then
       u%usr_msh_setup => dummy_user_mesh_setup
    end if

    if (.not. associated(u%usr_chk)) then
       u%usr_chk => dummy_user_check
    end if
    if (.not. associated(u%user_init_modules)) then
       u%user_init_modules => dummy_user_init_no_modules
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
    type(param_t), intent(inout) :: params
    call neko_error('Dummy user defined initial condition set')    
  end subroutine dummy_user_ic

  subroutine dummy_user_f_vector(f)
     class(source_t) :: f
    call neko_error('Dummy user defined vector valued forcing set')    
  end subroutine dummy_user_f_vector


  !> Dummy user forcing
  subroutine dummy_user_f(u, v, w, j, k, l, e)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    call neko_error('Dummy user defined forcing set')    
  end subroutine dummy_user_f
 
  !> Dummy user mesh apply
  subroutine dummy_user_mesh_setup(msh)
    type(mesh_t), intent(inout) :: msh
  end subroutine dummy_user_mesh_setup
  
  !> Dummy user check
  subroutine dummy_user_check(t, dt, tstep, u, v, w, p, coef)
    real(kind=rp), intent(in) :: t, dt
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
  end subroutine dummy_user_check

  subroutine dummy_user_init_no_modules(t, u, v, w, p, c_Xh, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: c_Xh
    type(param_t), intent(inout) :: params
  end subroutine dummy_user_init_no_modules

end module user_intf
