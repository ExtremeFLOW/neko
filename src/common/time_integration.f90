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
module time_integration
  use neko_config
  use num_types
  use math
  use utils
  use device
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Class storing time-integration coefficients for the explicit extrapolation
  !! and Backward-differencing schemes. 
  type, public :: time_integration_t
     real(kind=rp), dimension(10) :: ext
     real(kind=rp), dimension(10) :: bdf
     integer :: nab = 0
     integer :: nbd = 0
     integer :: time_order  !< Default is 3
     type(c_ptr) :: ext_d = C_NULL_PTR !< dev. ptr for coefficients
     type(c_ptr) :: bdf_d = C_NULL_PTR !< dev. ptr for coefficients
   contains
     procedure, pass(this) :: set_bd => time_integration_set_bdf
     procedure, pass(this) :: set_abbd => time_integration_set_ext
     procedure, pass(this) :: set_time_order => time_integration_set_time_order
     final :: time_integration_free
  end type time_integration_t


contains

  subroutine time_integration_free(this)
    type(time_integration_t), intent(inout) :: this

    if (c_associated(this%ext_d)) then
       call device_free(this%ext_d)
    end if
       
    if (c_associated(this%bdf_d)) then
       call device_free(this%bdf_d)
    end if
    
  end subroutine time_integration_free

  subroutine time_integration_set_time_order(this,torder)
    integer, intent(in) :: torder
    class(time_integration_t), intent(inout) :: this
    if(torder .le. 3 .and. torder .gt. 0) then
       this%time_order = torder
    else
       this%time_order = 3
       call neko_warning('Invalid time order, defaulting to 3')
    end if

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_map(this%ext, this%ext_d, 10)
       call device_map(this%bdf, this%bdf_d, 10)
    end if

  end subroutine time_integration_set_time_order

  !>Compute backward-differentiation coefficients of order NBD
  subroutine time_integration_set_bdf(this, dtbd)
    class(time_integration_t), intent(inout) :: this
    real(kind=rp), intent(inout), dimension(10) :: dtbd
    real(kind=rp), dimension(10,10) :: bdmat
    real(kind=rp), dimension(10) :: bdrhs
    real(kind=rp), dimension(10) :: bd_old
    real(kind=rp) :: bdf
    integer, parameter :: ldim = 10
    integer, dimension(10) :: ir, ic
    integer :: ibd, nsys, i

    associate(nbd => this%nbd, bd => this%bdf, bdf_d => this%bdf_d)
      bd_old = bd
      nbd = nbd + 1
      nbd = min(nbd, this%time_order)
      call rzero(bd, 10)
      if (nbd .eq. 1) then
         bd(1) = 1.0_rp
         bdf = 1.0_rp
      else if (nbd .ge. 2) then
         nsys = nbd + 1
         call bdsys(bdmat, bdrhs, dtbd, nbd, ldim)
         call lu(bdmat, nsys, ldim, ir, ic)
         call solve(bdrhs, bdmat, 1, nsys, ldim, ir, ic)
         do i = 1, nbd
            bd(i) = bdrhs(i)
         end do
         bdf = bdrhs(nbd + 1)
      endif
      
      !Normalize
      do ibd = nbd, 1, -1
         bd(ibd + 1) = bd(ibd)
      end do
      bd(1) = 1.0_rp
      do ibd= 1, nbd + 1
         bd(ibd) = bd(ibd)/bdf
      end do

      if (c_associated(bdf_d)) then
         if (maxval(abs(bd - bd_old)) .gt. 1e-10_rp) then
            call device_memcpy(bd, bdf_d, 10, HOST_TO_DEVICE)
         end if
      end if
    end associate
    
  end subroutine time_integration_set_bdf

  !>
  !! Compute Adams-Bashforth coefficients (order NAB, less or equal to 3)
  !!   
  !! NBD .EQ. 1
  !! Standard Adams-Bashforth coefficients 
  !!
  !! NBD .GT. 1
  !! Modified Adams-Bashforth coefficients to be used in con-
  !! junction with Backward Differentiation schemes (order NBD)
  !!
  subroutine time_integration_set_ext(this, dtlag)
    class(time_integration_t), intent(inout)  :: this
    real(kind=rp), intent(inout), dimension(10) :: dtlag
    real(kind=rp) :: dt0, dt1, dt2, dts, dta, dtb, dtc, dtd, dte
    real(kind=rp), dimension(10) :: ab_old
    associate(nab => this%nab, nbd => this%nbd, ext => this%ext, ext_d => this%ext_d)
      ab_old = ext
      nab = nab + 1
      nab = min(nab, this%time_order)
    
      dt0 = dtlag(1)
      dt1 = dtlag(2)
      dt2 = dtlag(3)
      call rzero(ext, 10)
      
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
    
  end subroutine time_integration_set_ext

  !
  ! CLEAN UP THIS MESS BELOW, USE LAPACK OR SIMILAR
  !


  
  subroutine bdsys (a, b, dt, nbd, ldim)
    integer :: ldim, j, n, k, i, nsys, nbd
    real(kind=rp) ::  A(ldim,9),B(9),DT(9)
    real(kind=rp) :: SUMDT

    CALL RZERO (A, ldim**2)
    N = NBD + 1
    DO J = 1, NBD
       A(1,J) = 1.0_rp
    end DO
    A(1,NBD+1) = 0.0_rp
    B(1) = 1.0_rp
    DO J = 1, NBD
       SUMDT = 0.0_rp
       DO  K = 1, J
          SUMDT = SUMDT + DT(K)
       end DO
       A(2,J) = SUMDT
    end DO
    A(2,NBD+1) = -DT(1)
    B(2) = 0.0_rp
    DO I = 3, NBD + 1
       DO J = 1, NBD
          SUMDT = 0.0_rp
          DO K = 1, J
             SUMDT = SUMDT + DT(K)
          end DO
          A(I,J) = SUMDT**(I-1)
       end DO
       A(I,NBD+1) = 0.0_rp
       B(I) = 0.0_rp
    end DO
      
    end subroutine bdsys


    SUBROUTINE LU(A, N, ldim, IR, IC)
      integer :: n, ldim, IR(10), IC(10)
      real(kind=rp) :: A(ldim,10), xmax, ymax, B, Y, C
      integer :: i, j, k, l, m, icm, irl, k1

      DO I = 1, N
         IR(I) = I
         IC(I) = I
      end DO
      K = 1
      L = K
      M = K
      XMAX = ABS(A(K,K))
      DO I = K, N
         DO J = K, N
            Y = ABS(A(I,J))
            IF(XMAX .GE. Y) GOTO 100
            XMAX = Y
            L = I
            M = J
100      END DO
      END DO
      DO K = 1, N
         IRL = IR(L)
         IR(L) = IR(K)
         IR(K) = IRL
         ICM = IC(M)
         IC(M) = IC(K)
         IC(K) = ICM
         IF(L .EQ. K) GOTO 300
         DO J = 1, N
            B = A(K,J)
            A(K,J) = A(L,J)
            A(L,J) = B
         END DO
300      IF(M .EQ. K) GOTO 500
         DO I = 1, N
            B = A(I,K)
            A(I,K) = A(I,M)
            A(I,M) = B
         END DO
500      C = 1.0_rp / A(K,K)
         A(K,K) = C
         IF(K .EQ. N) GOTO 1000
         K1 = K + 1
         XMAX = ABS(A(K1,K1))
         L = K1
         M = K1
         DO I = K1, N
            A(I,K) = C * A(I,K)
         END DO
         DO I = K1, N
            B = A(I,K)
            DO J = K1, N
               A(I,J) = A(I,J) - B * A(K,J)
               Y = ABS(A(I,J))
               IF(XMAX .GE. Y) GOTO 800
               XMAX = Y
               L = I
               M = J
800         END DO
         END DO
1000  END DO
    end subroutine lu
   
    SUBROUTINE SOLVE(F, A, K, N, ldim, IR, IC)
      integer :: IR(10),IC(10), N, N1, k, kk, i, j, ldim, ICM, URL, K1, ICI
      integer :: I1, IRI,IRL, IT
      real(kind=rp) ::  A(ldim,10), F(ldim,10), G(2000), B, Y

      N1 = N + 1
      DO KK = 1, K
         DO I = 1, N
            IRI = IR(I)
            G(I) = F(IRI,KK)
         END DO
         DO I = 2, N
            I1 = I - 1
            B = G(I)
            DO J = 1, I1
               B = B - A(I,J) * G(J)
            END DO
            G(I) = B
         END DO
         DO IT = 1, N
            I = N1 - IT
            I1 = I + 1
            B = G(I)
            IF(I .EQ. N) GOTO 701
            DO J = I1, N
               B = B - A(I,J) * G(J)
            END DO
701         G(I) = B * A(I,I)
         END DO
         DO I = 1, N
            ICI = IC(I)
            F(ICI,KK) = G(I)
         END DO
      END DO
    END SUBROUTINE SOLVE

  

  
end module time_integration
