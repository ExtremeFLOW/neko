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
!> Backward-differencing scheme for time integration.
module bdf_time_scheme
  use neko_config
  use num_types, only : rp
  use time_scheme, only: time_scheme_t
  use math, only : rzero
  use utils, only : neko_error
  implicit none
  private
  
  !> Implicit backward-differencing scheme for time integration.
  !! @details
  !! The explicit forumlas for the coefficients are taken from the following
  !! techincal note:
  !! "Derivation of BDF2/BDF3 for Variable Step Size" by Hiroaki Nishikawa,
  !! which can be found on ResearchGate.
  !!
  !! For a contant time-step this corresponds to the following schemes for
  !! order 1 to 3:
  !! - Order 1: \f$ \frac{1}{\Delta t} u^{n+1} - \frac{1}{\Delta t} u^{n} \f$ 
  !! - Order 2: \f$ \frac{3}{2\Delta t} u^{n+1} - \frac{4}{2\Delta t} u^{n}  + 
  !!                \frac{1}{2\Delta t} u^{n-1}\f$ 
  !! - Order 3: \f$ \frac{11}{6\Delta t} u^{n+1} - \frac{18}{6\Delta t} u^{n}  + 
  !!                \frac{9}{6\Delta t} u^{n-1} - \frac{2}{6\Delta t} u^{n-2}\f$ 
  !!
  !! It is assumed that all the coefficients but the first one premultiply terms
  !! that go to the right-hand side of the equation. 
  !! Accordingly, the signs of these coefficients are reversed in the `coeffs`
  !! array. This is taken into account, for example, in the implemeation of the
  !! `rhs_maker` class.
  !!
  !! Another important convention is that the coefficients are meant to be
  !! later divided by the current value of the timestep.
  !!
  !! In line with the above assumptions, the first order scheme always returns
  !! the array \f$[1, 1]\f$, and **not** \f$[1/\Delta t, -1/\Delta t]\f$, as one
  !! might expect. Similar for the second and third order.
  !!
  !! @remark 
  !! The current implementation can be easily extended to schemes of arbitrary
  !! order, by using the `fd_weights_full` subroutine to compute the
  !! coefficients. A demonstration of this is implemented in a test in
  !! `tests/ext_bdf_scheme/test_bdf.pf`
  type, public, extends(time_scheme_t) :: bdf_time_scheme_t 
     contains
       !> Compute the scheme coefficients
       procedure, nopass :: compute_coeffs => bdf_time_scheme_compute_coeffs
  end type bdf_time_scheme_t
  
  contains

  !> Compute the scheme coefficients
  !! @param t Timestep values, first element is the current timestep.
  !! @param order Order the scheme.
  subroutine bdf_time_scheme_compute_coeffs(coeffs, dt, order)
    implicit none
    real(kind=rp), intent(out) :: coeffs(4)
    real(kind=rp), intent(in) :: dt(10)
    integer, intent(in) :: order

    call rzero(coeffs, 4)
    
    ! Note, these are true coeffs, multiplied by dt(1)
    select case (order)
    case (1)
       coeffs(1) = 1.0_rp
       coeffs(2) = 1.0_rp
    case (2)
       coeffs(1) = (1 + dt(1)/(dt(1) + dt(2)))
       coeffs(3) = -dt(1)**2/dt(2)/(dt(1) + dt(2))
       coeffs(2) =  coeffs(1) - coeffs(3)
    case (3)
       coeffs(2) = (dt(1) + dt(2)) * (dt(1) + dt(2) + dt(3)) / &
                   (dt(1) * dt(2) * (dt(2) + dt(3)))
       coeffs(3) = -dt(1) * (dt(1) + dt(2) + dt(3)) / &
                   (dt(2) * dt(3) * (dt(1) + dt(2)))
       coeffs(4) = dt(1) * (dt(1) + dt(2)) / &
                   (dt(3) * (dt(2) + dt(3)) * (dt(1) + dt(2) + dt(3)))
       coeffs(1) = coeffs(2) + coeffs(3) + coeffs(4)
       coeffs = coeffs * dt(1)
    case default 
      call neko_error("The order of the BDF time scheme must be 1 to 3.")
    end select
    
  end subroutine bdf_time_scheme_compute_coeffs

end module bdf_time_scheme