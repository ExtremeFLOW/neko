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
!> Explicit and Backward Differentiation time-integration schemes
module ext_bdf_scheme
  use neko_config
  use num_types
  use math
  use utils
  use device
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Base abstract class for time integration schemes
  !! @details
  !! An important detail here is the handling of the first timesteps where a 
  !! high-order scheme cannot be constructed. The parameters `n`, which is
  !! initialized to 0, must be incremented by 1 in the beggining of the 
  !! `set_coeffs` routine to determine the current scheme order.
  !! When `n == time_order`, the incrementation should stop.
  type, abstract, public :: time_scheme_t
     !> The coefficients of the scheme
     real(kind=rp), dimension(4) :: coeffs 
     !> Controls the actual order of the scheme, e.g. 1 at the first time-step
     integer :: n = 0
     !> Order of the scheme, defaults to 3
     integer :: time_order
     !> Device pointer for `coeffs`
     type(c_ptr) :: coeffs_d = C_NULL_PTR 
   contains
     !> Controls current scheme order and computes the coefficients
     procedure(set_coeffs), deferred, pass(this) :: set_coeffs
     !> Constructor
     procedure, pass(this) :: init => time_scheme_init 
     !> Destructor
     procedure, pass(this) :: free => time_scheme_free
  end type time_scheme_t
  
  abstract interface
     !> Interface for setting the scheme coefficients
     !! @param t Timestep values, first element is the current timestep.
     subroutine set_coeffs(this, dt)
       import time_scheme_t
       import rp
       class(time_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout), dimension(10) :: dt
     end subroutine
  end interface
  
  type, public, extends(time_scheme_t) :: ext_time_scheme_t 
     contains
       procedure, pass(this) :: set_coeffs => ext_time_scheme_set_coeffs
  end type ext_time_scheme_t

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
       procedure, pass(this) :: set_coeffs => bdf_time_scheme_set_coeffs
  end type bdf_time_scheme_t
  
  type, public :: ext_bdf_scheme_t
     type(ext_time_scheme_t) :: ext
     type(bdf_time_scheme_t) :: bdf
     
     contains
       procedure, pass(this) :: init => ext_bdf_scheme_init
       procedure, pass(this) :: free => ext_bdf_scheme_free
  end type ext_bdf_scheme_t

contains
  !> Constructor
  subroutine time_scheme_init(this, torder)
    class(time_scheme_t), intent(inout) :: this
    integer, intent(in) :: torder !< Desired order of the scheme: 1, 2 or 3.

    if(torder .le. 3 .and. torder .gt. 0) then
       this%time_order = torder
    else
       this%time_order = 3
       call neko_warning('Invalid time order, defaulting to 3')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%coeffs, this%coeffs_d, 4)
    end if
  end subroutine time_scheme_init

  !> Destructor
  subroutine time_scheme_free(this)
    class(time_scheme_t), intent(inout) :: this

    if (c_associated(this%coeffs_d)) then
       call device_free(this%coeffs_d)
    end if
  end subroutine time_scheme_free

  !> Compute Adams-Bashforth coefficients (order NAB, less or equal to 3)
  !!   
  !! NBD .EQ. 1
  !! Standard Adams-Bashforth coefficients 
  !!
  !! NBD .GT. 1
  !! Modified Adams-Bashforth coefficients to be used in con-
  !! junction with Backward Differentiation schemes (order NBD)
  !!
  subroutine ext_time_scheme_set_coeffs(this, dt)
    class(ext_time_scheme_t), intent(inout)  :: this
    real(kind=rp), intent(inout), dimension(10) :: dt
    real(kind=rp) :: dt0, dt1, dt2, dts, dta, dtb, dtc, dtd, dte
    real(kind=rp), dimension(4) :: ab_old
    associate(nab => this%n, nbd => this%n, ext => this%coeffs, ext_d => this%coeffs_d)
      ab_old = ext
      nab = nab + 1
      nab = min(nab, this%time_order)
    
      dt0 = dt(1)
      dt1 = dt(2)
      dt2 = dt(3)
      call rzero(ext, 4)
      
      if (nab .eq. 1) then
         ext(1) = 1.0_rp
      else if (nab .eq. 2) then
         dta =  dt0 / dt1
         if (nbd .eq. 1) then
            ext(2) = -0.5_rp * dta
            ext(1) =  1.0_rp - ext(2)
         else if (nbd .eq. 2) then
            ext(2) = -dta
            ext(1) =  1.0_rp - ext(2)
         endif
      else if (nab .eq. 3) then
         dts =  dt1 + dt2
         dta =  dt0 / dt1
         dtb =  dt1 / dt2
         dtc =  dt0 / dt2
         dtd =  dts / dt1
         dte =  dt0 / dts
         if (nbd .eq. 1) then
            ext(3) =  dte*( 0.5d0*dtb + dtc/3d0 )
            ext(2) = -0.5_rp * dta - ext(3) * dtd
            ext(1) =  1.0_rp - ext(2) - ext(3)
         elseif (nbd .eq. 2) then
            ext(3) =  2.0_rp / 3.0_rp * dtc * (1.0_rp / dtd + dte)
            ext(2) = -dta - ext(3) * dtd
            ext(1) =  1.0_rp - ext(2) - ext(3)
         elseif (nbd .eq. 3) then
            ext(3) =  dte * (dtb + dtc)
            ext(2) = -dta * (1.0_rp + dtb + dtc)
            ext(1) =  1.0_rp - ext(2) - ext(3)
         endif
      endif

      if (c_associated(ext_d)) then
         if (maxval(abs(ext - ab_old)) .gt. 1e-10_rp) then
            call device_memcpy(ext, ext_d, 10, HOST_TO_DEVICE)
         end if
      end if
    end associate
    
  end subroutine ext_time_scheme_set_coeffs

  subroutine bdf_time_scheme_set_coeffs(this, dt)
    implicit none
    class(bdf_time_scheme_t), intent(inout) :: this
    real(kind=rp), intent(inout), dimension(10) :: dt
    real(kind=rp), dimension(4) :: coeffs_old

    associate(n => this%n, coeffs => this%coeffs, coeffs_d => this%coeffs_d)
      ! will be used to check whether the coefficients changed
      coeffs_old = coeffs
      
      ! Increment the order of the scheme if below time_order
      n = n + 1
      n = min(n, this%time_order)

      call rzero(coeffs, 4)
      
      ! Note, these are true coeffs, multiplied by dt(1)
      select case (n)
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
      end select
      
      if (c_associated(coeffs_d)) then
         if (maxval(abs(coeffs - coeffs_old)) .gt. 1e-10_rp) then
            call device_memcpy(coeffs, coeffs_d, 4, HOST_TO_DEVICE)
         end if
      end if
    end associate
    
  end subroutine bdf_time_scheme_set_coeffs

  !> Contructor for ext_bdf_scheme_t
  subroutine ext_bdf_scheme_init(this, torder)
     class(ext_bdf_scheme_t) :: this
     integer :: torder !< Desired order of the scheme: 1, 2 or 3.
     call time_scheme_init(this%ext, torder)
     call time_scheme_init(this%bdf, torder)
  end subroutine ext_bdf_scheme_init

  !> Destructor for ext_bdf_scheme_t
  subroutine ext_bdf_scheme_free(this)
     class(ext_bdf_scheme_t) :: this
     call time_scheme_free(this%ext)
     call time_scheme_free(this%bdf)
  end subroutine ext_bdf_scheme_free
  
end module ext_bdf_scheme
