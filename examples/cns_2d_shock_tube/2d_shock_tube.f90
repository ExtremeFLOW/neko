! Copyright (c) 2026, The Neko Authors
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
! Reflected shock-boundary layer interaction in a 2D shock tube.
module user
  use neko
  implicit none

  real(kind=rp) :: Re, gamma, Pr

contains

  !> Register user callbacks for the shock-tube example.
  !! @param user User interface object.
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initial_conditions => initial_conditions
    user%material_properties => material_properties
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, 'case.fluid.Re', Re)
    call json_get(params, 'case.fluid.gamma', gamma)
    call json_get(params, 'case.fluid.Pr', Pr)

  end subroutine startup

  !> Set the Riemann initial condition with the diaphragm at x = 0.5.
  !! @param scheme_name Name of the fluid scheme.
  !! @param fields Simulation fields supplied by Neko.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: rho, u, v, w, p
    integer :: i
    real(kind=rp) :: x

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    do i = 1, rho%dof%size()
       x = rho%dof%x(i, 1, 1, 1)

       u%x(i, 1, 1, 1) = 0.0_rp
       v%x(i, 1, 1, 1) = 0.0_rp
       w%x(i, 1, 1, 1) = 0.0_rp

       if (x < 0.5_rp) then
          rho%x(i, 1, 1, 1) = 120.0_rp
          p%x(i, 1, 1, 1) = 120.0_rp / gamma
       else
          rho%x(i, 1, 1, 1) = 1.2_rp
          p%x(i, 1, 1, 1) = 1.2_rp / gamma
       end if
    end do
  end subroutine initial_conditions

  !> Set dynamic viscosity and thermal conductivity.
  !! @param scheme_name Name of the fluid scheme.
  !! @param properties Material property fields supplied by Neko.
  !! @param time Current simulation time state.
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: mu, kappa
    integer :: i

    ! Only set material properties at the first time step
    if (time%tstep .ne. 0) return

    !! Either this or we can also just add mu and kappa to
    !! the json file under case.fluid
    !!     "mu": 0.001,
    !!     "kappa": 0.001917808219178082,
    if (scheme_name .eq. "fluid") then
       mu => properties%get_by_name("fluid_mu")
       kappa => properties%get_by_name("fluid_kappa")

       call field_cfill(mu, 1.0_rp / Re)

       ! kappa = mu * gamma / Pr
       call field_cmult2(kappa, mu, gamma / Pr)

    end if
  end subroutine material_properties

end module user
