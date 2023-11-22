! Copyright (c) 2021, The Neko Authors
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
!> Scalar initial condition
module scalar_ic
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device_math, only : device_col2
  use device, only : device_memcpy, HOST_TO_DEVICE
  use field, only : field_t
  use utils, only : neko_error
  use coefs, only : coef_t
  use math, only : col2, cfill
  use user_intf, only : useric_scalar
  use json_module, only : json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  interface set_scalar_ic
     module procedure set_scalar_ic_int, set_scalar_ic_usr
  end interface

  public :: set_scalar_ic

contains

  !> Set scalar initial condition (builtin)
  subroutine set_scalar_ic_int(s, coef, gs, type, params)
    type(field_t), intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file), intent(inout) :: params
    ! Variables for retrieving json parameters
    real(kind=rp) :: ic_value

    if (trim(type) .eq. 'uniform') then
       call json_get_or_default(params, 'case.scalar.initial_condition.value', &
                                ic_value, 0.0_rp)
       call set_scalar_ic_uniform(s, ic_value)
    else
       call neko_error('Invalid initial condition')
    end if

    call set_scalar_ic_common(s, coef, gs)

  end subroutine set_scalar_ic_int

  !> Set scalar intial condition (user defined)
  subroutine set_scalar_ic_usr(s, coef, gs, usr_ic, params)
    type(field_t), intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(useric_scalar) :: usr_ic
    type(json_file), intent(inout) :: params

    call usr_ic(s, params)

    call set_scalar_ic_common(s, coef, gs)

  end subroutine set_scalar_ic_usr

  subroutine set_scalar_ic_common(s, coef, gs)
    type(field_t), intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), HOST_TO_DEVICE)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(s%x, s%dof%size(), GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(s%x_d, coef%mult_d, s%dof%size())
    else
       call col2(s%x, coef%mult, s%dof%size())
    end if

  end subroutine set_scalar_ic_common

  !> Uniform initial condition
  subroutine set_scalar_ic_uniform(s, ic_value)
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: ic_value
    integer :: n
    s = ic_value
    n = s%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(s%x, ic_value, n)
    end if

  end subroutine set_scalar_ic_uniform

end module scalar_ic
