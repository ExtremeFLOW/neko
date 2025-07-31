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
!> Neko API User driver
submodule(neko_api) neko_api_user
  implicit none

  abstract interface
     subroutine api_ic_callback
       use, intrinsic :: iso_c_binding
     end subroutine api_ic_callback
  end interface

  abstract interface
     subroutine api_if_callback(u, v, w, x, y, z, nx, ny, nz, &
          ix, iy, iz, ie, t, tstep) bind(C)
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       real(kind=c_rp) :: u
       real(kind=c_rp) :: v
       real(kind=c_rp) :: w
       real(kind=c_rp) :: x
       real(kind=c_rp) :: y
       real(kind=c_rp) :: z
       real(kind=c_rp) :: nx
       real(kind=c_rp) :: ny
       real(kind=c_rp) :: nz
       integer(c_int) :: ix
       integer(c_int) :: iy
       integer(c_int) :: iz
       integer(c_int) :: ie
       real(kind=c_rp) :: t
       integer(c_int) :: tstep
     end subroutine api_if_callback
  end interface

  type api_user_cb
     procedure(api_ic_callback), nopass, pointer :: fluid_ic => null()
     procedure(api_ic_callback), nopass, pointer :: fluid_cns_ic => null()
     procedure(api_ic_callback), nopass, pointer :: scalar_ic => null()
     procedure(api_if_callback), nopass, pointer :: fluid_if => null()
  end type api_user_cb

  type(api_user_cb) :: neko_api_user_cb

contains

  module subroutine neko_api_user_cb_register(user, &
       fluid_ic_cb, fluid_ic_cns_cb, scalar_ic_cb, fluid_if_cb)
    type(user_t), intent(inout) :: user
    type(c_funptr), value :: fluid_ic_cb, fluid_ic_cns_cb, scalar_ic_cb
    type(c_funptr), value :: fluid_if_cb

    if (c_associated(fluid_ic_cb)) then
       user%fluid_user_ic => neko_api_user_fluid_ic
       call c_f_procpointer(fluid_ic_cb, neko_api_user_cb%fluid_ic)
    end if

    if (c_associated(fluid_ic_cns_cb)) then
       user%fluid_compressible_user_ic => neko_api_user_fluid_cns_ic
       call c_f_procpointer(fluid_ic_cns_cb, neko_api_user_cb%fluid_cns_ic)
    end if

    if (c_associated(scalar_ic_cb)) then
       user%scalar_user_ic => neko_api_user_scalar_ic
       call c_f_procpointer(scalar_ic_cb, neko_api_user_cb%scalar_ic)
    end if

    if (c_associated(fluid_if_cb)) then
       user%fluid_user_if => neko_api_user_fluid_if
       call c_f_procpointer(fluid_if_cb, neko_api_user_cb%fluid_if)
    end if

  end subroutine neko_api_user_cb_register

  subroutine neko_api_user_fluid_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params

    if (associated(neko_api_user_cb%fluid_ic)) then
       call neko_api_user_cb%fluid_ic
    end if

  end subroutine neko_api_user_fluid_ic

  subroutine neko_api_user_fluid_cns_ic(rho, u, v, w, p, params)
    type(field_t), intent(inout) :: rho
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params

    if (associated(neko_api_user_cb%fluid_cns_ic)) then
       call neko_api_user_cb%fluid_cns_ic
    end if

  end subroutine neko_api_user_fluid_cns_ic

  subroutine neko_api_user_scalar_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params

    if (associated(neko_api_user_cb%scalar_ic)) then
       call neko_api_user_cb%scalar_ic
    end if

  end subroutine neko_api_user_scalar_ic

  subroutine neko_api_user_fluid_if(u, v, w, x, y, z, nx, ny, nz, &
       ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
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

    if (associated(neko_api_user_cb%fluid_if)) then
       call neko_api_user_cb%fluid_if(u, v, w, x, y, z, &
            nx, ny, nz, ix, iy, iz, ie, t, tstep)
    end if

  end subroutine neko_api_user_fluid_if

end submodule neko_api_user
