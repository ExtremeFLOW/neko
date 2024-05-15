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
!> Initial flow condition
module baseflow
  use num_types, only : rp
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use flow_profile, only : blasius_profile, blasius_linear, blasius_cubic, &
    blasius_quadratic, blasius_quartic, blasius_sin
  use field, only : field_t
  use utils, only : neko_error
  use coefs, only : coef_t
  use math, only : col2, cfill, copy
  use device_math, only : device_col2
  use user_intf, only : user_baseflow_init
  use json_module, only : json_file
  use json_utils, only: json_get
  use file, only: file_t, fld_file_data_t, file_free
  use device, only: device_memcpy, HOST_TO_DEVICE
  implicit none
  private

  interface set_baseflow
     module procedure set_baseflow_int, set_baseflow_usr
  end interface set_baseflow

  public :: set_baseflow

contains

  !> Set initial flow condition (builtin)
  subroutine set_baseflow_int(u_b, v_b, w_b, coef, gs, type, params)
    type(field_t), intent(inout) :: u_b
    type(field_t), intent(inout) :: v_b
    type(field_t), intent(inout) :: w_b
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file), intent(inout) :: params
    ! Variables for retrieving json parameters
    character(len=:), allocatable :: file_name

    if (trim(type) .eq. 'file') then
       ! Read the baseflow file
       call json_get(params, 'case.fluid.baseflow.file_name', file_name)
       call set_baseflow_file(u_b, v_b, w_b, file_name)
    else
       call neko_error('Invalid initial condition')
    end if

    call set_baseflow_common(u_b, v_b, w_b, coef, gs)

  end subroutine set_baseflow_int

  !> Set intial flow condition (user defined)
  subroutine set_baseflow_usr(u_b, v_b, w_b, coef, gs, usr_baseflow, params)
    type(field_t), intent(inout) :: u_b
    type(field_t), intent(inout) :: v_b
    type(field_t), intent(inout) :: w_b
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(user_baseflow_init) :: usr_baseflow
    type(json_file), intent(inout) :: params

    call usr_baseflow(u_b, v_b, w_b, params)

    call set_baseflow_common(u_b, v_b, w_b, coef, gs)

  end subroutine set_baseflow_usr

  subroutine set_baseflow_common(u_b, v_b, w_b, coef, gs)
    type(field_t), intent(inout) :: u_b
    type(field_t), intent(inout) :: v_b
    type(field_t), intent(inout) :: w_b
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u_b%x, u_b%x_d, u_b%dof%size(), &
                                                        HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v_b%x, v_b%x_d, v_b%dof%size(), &
                                                        HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w_b%x, w_b%x_d, w_b%dof%size(), &
                                                        HOST_TO_DEVICE, sync=.false.)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(u_b%x, u_b%dof%size(), GS_OP_ADD)
    call gs%op(v_b%x, v_b%dof%size(), GS_OP_ADD)
    call gs%op(w_b%x, w_b%dof%size(), GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(u_b%x_d, coef%mult_d, u_b%dof%size())
       call device_col2(v_b%x_d, coef%mult_d, v_b%dof%size())
       call device_col2(w_b%x_d, coef%mult_d, w_b%dof%size())
    else
       call col2(u_b%x, coef%mult, u_b%dof%size())
       call col2(v_b%x, coef%mult, v_b%dof%size())
       call col2(w_b%x, coef%mult, w_b%dof%size())
    end if

  end subroutine set_baseflow_common

  !> Read the baseflow from file or user input.
  !!
  !! This baseflow reader is currently capable of:
  !!   - loading a baseflow from a file
  !!   - Todo: Prescribing one from the user
  subroutine set_baseflow_file(u_b, v_b, w_b, file_name)
    type(field_t), intent(inout) :: u_b
    type(field_t), intent(inout) :: v_b
    type(field_t), intent(inout) :: w_b
    character(len=:), allocatable :: file_name

    type(file_t) :: field_file
    type(fld_file_data_t) :: field_data

    integer :: n
    real(kind=rp) :: norm_l2_base

    ! Allocation of the baseflow fields
    n = u_b%size()

    ! assume they're on the same mesh
    field_file = file_t(trim(file_name), precision = rp)
    call field_data%init
    call field_file%read(field_data)
    call copy(u_b%x, field_data%u%x, n)
    call copy(v_b%x, field_data%v%x, n)
    call copy(w_b%x, field_data%w%x, n)

    call file_free(field_file)
    call field_data%free()

  end subroutine set_baseflow_file
end module baseflow
