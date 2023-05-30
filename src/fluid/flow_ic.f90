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
module flow_ic
  use gather_scatter
  use neko_config
  use flow_profile
  use device_math
  use field
  use utils
  use coefs
  use math
  use user_intf, only : useric
  use json_module, only : json_file_t => json_file
  implicit none
  private

  interface set_flow_ic
     module procedure set_flow_ic_int, set_flow_ic_usr
  end interface

  public :: set_flow_ic
  
contains

  !> Set initial flow condition (builtin)
  subroutine set_flow_ic_int(u, v, w, p, coef, gs, type, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file_t), intent(inout) :: params
    ! Variables for retrieving json parameters
    logical :: found
    real(kind=rp) :: delta
    real(kind=rp), allocatable :: uinf(:)
    character(len=:), allocatable :: blasius_approximation

    if (trim(type) .eq. 'uniform') then       
       call params%get('case.fluid.initial_condition.value', uinf, found)
       if (.not. found) then
          call neko_error("Parameter fluid.initial_condition.value missing in &
                          &the case file")
       end if

       call set_flow_ic_uniform(u, v, w, uinf)
    else if (trim(type) .eq. 'blasius') then
       call params%get('case.fluid.blasius.delta', delta, found)
       if (.not. found) then
          call neko_error("Parameter fluid.blasius.delta missing in &
                          &the case file")
       end if

       call params%get('case.fluid.blasius.approximation', &
                       blasius_approximation, found)
       if (.not. found) then
          call neko_error("Parameter fluid.blasius.approximation missing in &
                          &the case file")
       end if

       call params%get('case.fluid.blasius.freestream_velocity', &
                       uinf, found)
       if (.not. found) then
          call neko_error("Parameter fluid.blasius.freeestream_velocity missing&
                          & in the case file")
       end if

       call set_flow_ic_blasius(u, v, w, delta, uinf, blasius_approximation)
    else
       call neko_error('Invalid initial condition')
    end if
    
    call set_flow_ic_common(u, v, w, p, coef, gs)
    
  end subroutine set_flow_ic_int

  !> Set intial flow condition (user defined)
  subroutine set_flow_ic_usr(u, v, w, p, coef, gs, usr_ic, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(useric) :: usr_ic
    type(json_file_t), intent(inout) :: params

    call usr_ic(u, v, w, p, params)
    
    call set_flow_ic_common(u, v, w, p, coef, gs)
    
  end subroutine set_flow_ic_usr

  subroutine set_flow_ic_common(u, v, w, p, coef, gs)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(u%x, u%x_d, u%dof%size(), HOST_TO_DEVICE)
       call device_memcpy(v%x, v%x_d, v%dof%size(), HOST_TO_DEVICE)
       call device_memcpy(w%x, w%x_d, w%dof%size(), HOST_TO_DEVICE)
    end if
    
    ! Ensure continuity across elements for initial conditions
    call gs_op(gs, u%x, u%dof%size(), GS_OP_ADD) 
    call gs_op(gs, v%x, v%dof%size(), GS_OP_ADD) 
    call gs_op(gs, w%x, w%dof%size(), GS_OP_ADD) 

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(u%x_d, coef%mult_d, u%dof%size())
       call device_col2(v%x_d, coef%mult_d, v%dof%size())
       call device_col2(w%x_d, coef%mult_d, w%dof%size())
    else
       call col2(u%x, coef%mult, u%dof%size())
       call col2(v%x, coef%mult, v%dof%size())
       call col2(w%x, coef%mult, w%dof%size())
    end if
    
  end subroutine set_flow_ic_common

  !> Uniform initial condition
  subroutine set_flow_ic_uniform(u, v, w, uinf)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in) :: uinf(3)
    integer :: n
    u = uinf(1)
    v = uinf(2)
    w = uinf(3)
    n = u%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(u%x, uinf(1), n)
       call cfill(v%x, uinf(2), n)
       call cfill(w%x, uinf(3), n)
    end if

  end subroutine set_flow_ic_uniform

  !> Set a Blasius profile as initial condition
  !! @note currently limited to axis aligned flow
  subroutine set_flow_ic_blasius(u, v, w, delta, uinf, type)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in) :: delta
    real(kind=rp), intent(in) :: uinf(3)
    character(len=*), intent(in) :: type
    procedure(blasius_profile), pointer :: bla => null()
    integer :: i

    select case(trim(type))
    case('linear')
       bla => blasius_linear
    case('quadratic')
       bla => blasius_quadratic
    case('cubic')
       bla => blasius_cubic
    case('quartic')
       bla => blasius_quartic
    case('sin')
       bla => blasius_sin
    case default
       call neko_error('Invalid Blasius approximation')
    end select
    
    if ((uinf(1) .gt. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
         .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = bla(u%dof%z(i,1,1,1), delta, uinf(1))
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .gt. 0.0_rp) &
         .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = bla(u%dof%x(i,1,1,1), delta, uinf(2))
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
         .and. (uinf(3) .gt. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = bla(u%dof%y(i,1,1,1), delta, uinf(3))
       end do       
    end if
    
  end subroutine set_flow_ic_blasius
  
end module flow_ic
