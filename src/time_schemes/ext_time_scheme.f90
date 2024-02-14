! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
!> Explicit extrapolation scheme for time integration.
module ext_time_scheme
  use neko_config
  use num_types, only : rp
  use time_scheme, only: time_scheme_t
  use math, only : rzero
  use utils, only : neko_error
  implicit none
  private

  !> Explicit extrapolation scheme for time integration.
  !! @details
  !! Compute the value at the current time-step by evaluating a polynomial
  !! built using the values on previous time-steps.
  !!
  !! For a constant time-step corresponds to the following schemes of order
  !! 1 to 3:
  !! - Order 1: \f$  u^{n+1} = u^n \f$
  !! - Order 2: \f$  u^{n+1} = 2u^n - u^{n-1} \f$, linear extrapolation
  !! - Order 3: \f$  u^{n+1} = 3u^n - 3u^{n-1} + u^{n-2} \f$
  !! An additional scheme is available via `compute_modified_coeffs`, which is
  !! gives improved stability when combined with BDF2 in a advection-diffusion
  !! equation.
  !! - Order 3: \f$  u^{n+1} = 8/3u^n - 7/3u^{n-1} + 2/3u^{n-2} \f$
  type, public, extends(time_scheme_t) :: ext_time_scheme_t
   contains
     !> Compute the scheme coefficients
     procedure, nopass :: compute_coeffs => ext_time_scheme_compute_coeffs
     !> Compute the coefficients for the modified EXT scheme
     procedure, nopass :: compute_modified_coeffs => &
         ext_time_scheme_compute_modified_coeffs
  end type ext_time_scheme_t

contains

  !> Compute the scheme coefficients
  !! @param t Timestep values, first element is the current timestep.
  !! @param order Order the scheme.
  subroutine ext_time_scheme_compute_coeffs(coeffs, dt, order)
    implicit none
    real(kind=rp), intent(out) :: coeffs(4)
    real(kind=rp), intent(in) :: dt(10)
    integer, intent(in) :: order

    call rzero(coeffs, 4)

    select case (order)
    case (1)
       coeffs(1) = 1.0_rp
    case (2)
       coeffs(2) = -dt(1) / dt(2)
       coeffs(1) =  1.0_rp - coeffs(2)
    case (3)
       coeffs(3) =  dt(1) / (dt(2) + dt(3)) * (dt(1) + dt(2)) / dt(3)
       coeffs(2) = - dt(1) / dt(2) * (1.0_rp + dt(2) / dt(3) + dt(1) / dt(3))
       coeffs(1) =  1.0_rp - coeffs(2) - coeffs(3)
    case default
       call neko_error("The order of the EXT time scheme must be 1 to 3.")
    end select
  end subroutine ext_time_scheme_compute_coeffs

  !> Compute the modified scheme coefficients
  !! @param t Timestep values, first element is the current timestep.
  subroutine ext_time_scheme_compute_modified_coeffs(coeffs, dt)
    implicit none
    real(kind=rp), intent(out) :: coeffs(4)
    real(kind=rp), intent(in) :: dt(10)
    real(kind=rp) dta, dtb, dtc, dtd, dte, dts

    call rzero(coeffs, 4)

    dts =  dt(2) + dt(3)
    dta =  dt(1) / dt(2)
    dtb =  dt(2) / dt(3)
    dtc =  dt(1) / dt(3)
    dtd =  dts / dt(2)
    dte =  dt(1) / dts
    coeffs(3) =  2.0_rp / 3.0_rp * dtc * (1.0_rp / dtd + dte)
    coeffs(2) = -dta - coeffs(3) * dtd
    coeffs(1) =  1.0_rp - coeffs(2) - coeffs(3)

  end subroutine ext_time_scheme_compute_modified_coeffs
end module ext_time_scheme
