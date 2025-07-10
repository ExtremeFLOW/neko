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
!==============================================================================
!
!>     LIBRARY ROUTINES FOR SPECTRAL METHODS
!!
!!     March 1989
!!
!!     For questions, comments or suggestions, please contact:
!!
!!     Einar Malvin Ronquist
!!     Room 3-243
!!     Department of Mechanical Engineering
!!     Massachusetts Institute of Technology
!!     77 Massachusetts Avenue
!!     Cambridge, MA 0299
!!     U.S.A.
!
!------------------------------------------------------------------------------
!
!     ABBRIVIATIONS:
!
!     M   - Set of mesh points
!     Z   - Set of collocation/quadrature points
!     W   - Set of quadrature weights
!     H   - Lagrangian interpolant
!     D   - Derivative operator
!     I   - Interpolation operator
!     GL  - Gauss Legendre
!     GLL - Gauss-Lobatto Legendre
!     GJ  - Gauss Jacobi
!     GLJ - Gauss-Lobatto Jacobi
!
!
!     MAIN ROUTINES:
!
!     Points and weights:
!
!     ZWGL      Compute Gauss Legendre points and weights
!     ZWGLL     Compute Gauss-Lobatto Legendre points and weights
!     ZWGJ      Compute Gauss Jacobi points and weights (general)
!     ZWGLJ     Compute Gauss-Lobatto Jacobi points and weights (general)
!
!     Lagrangian interpolants:
!
!     HGL       Compute Gauss Legendre Lagrangian interpolant
!     HGLL      Compute Gauss-Lobatto Legendre Lagrangian interpolant
!     HGJ       Compute Gauss Jacobi Lagrangian interpolant (general)
!     HGLJ      Compute Gauss-Lobatto Jacobi Lagrangian interpolant (general)
!
!     Derivative operators:
!
!     DGLL      Compute Gauss-Lobatto Legendre derivative matrix
!     DGLLGL    Compute derivative matrix for a staggered mesh (GLL->GL)
!     DGJ       Compute Gauss Jacobi derivative matrix (general)
!     DGLJ      Compute Gauss-Lobatto Jacobi derivative matrix (general)
!     DGLJGJ    Compute derivative matrix for a staggered mesh (GLJ->GJ) (general)
!
!     Interpolation operators:
!
!     IGLM      Compute interpolation operator GL  -> M
!     IGLLM     Compute interpolation operator GLL -> M
!     IGJM      Compute interpolation operator GJ  -> M  (general)
!     IGLJM     Compute interpolation operator GLJ -> M  (general)
!
!     Other:
!
!     PNLEG     Compute Legendre polynomial of degree N
!     legendre_poly Compute Legendre polynomial of degree 0-N
!     PNDLEG    Compute derivative of Legendre polynomial of degree N
!
!     Comments:
!
!     Note that many of the above routines exist in both single and
!     double precision. If the name of the single precision routine is
!     SUB, the double precision version is called SUBD. In most cases
!     all the "low-level" arithmetic is done in double precision, even
!     for the single precsion versions.
!
!     Useful references:
!
! [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
!     Providence, Rhode Island, 1939.
! [2] Abramowitz & Stegun: Handbook of Mathematical Functions,
!     Dover, New York, 1972.
! [3] Canuto, Hussaini, Quarteroni & Zang: Spectral Methods in Fluid
!     Dynamics, Springer-Verlag, 1988.
!
!
!==============================================================================
module speclib
  use num_types, only : rp, xp
  use utils, only : neko_error
  use, intrinsic :: iso_fortran_env, only : stderr => error_unit
  implicit none

contains
  !> Generate `NP` Gauss Legendre points `Z` and weights `W`
  !! associated with Jacobi polynomial \f$ P(N)(\alpha=0, \beta=0) \f$.
  !! The polynomial degree `N = NP-1`.
  !! @param Z Quadrature points.
  !! @param W Quadrature weights.
  !! @param NP Number of quadrature points.
  subroutine ZWGL(Z, W, NP)
    integer, intent(in) :: NP
    real(kind=rp), intent(inout) :: Z(NP), W(NP)
    real(kind=rp) ALPHA, BETA
    ALPHA = 0.0_rp
    BETA = 0.0_rp
    call ZWGJ(Z, W, NP, ALPHA, BETA)
  end subroutine ZWGL

  subroutine ZWGLL(Z, W, NP)
!--------------------------------------------------------------------
!
!     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
!     The polynomial degree N=NP-1.
!     Z and W are in single precision, but all the arithmetic
!     operations are done in double precision.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NP
    real(kind=rp), intent(inout) :: Z(NP), W(NP)
    real(kind=rp) ALPHA, BETA
    ALPHA = 0.0_rp
    BETA = 0.0_rp
    call ZWGLJ(Z, W, NP, ALPHA, BETA)
  end subroutine ZWGLL

  subroutine ZWGJ(Z, W, NP, ALPHA, BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NP
    real(kind=rp), intent(inout) :: Z(NP), W(NP)
    real(kind=rp), intent(in) :: ALPHA, BETA

    integer, parameter :: NMAX = 84
    integer, parameter :: NZD = NMAX

    real(kind=xp) ZD(NZD), WD(NZD), ALPHAD, BETAD
    integer :: I, NPMAX

    NPMAX = NZD
    if (NP .gt. NPMAX) then
       write (stderr, *) 'Too large polynomial degree in ZWGJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here NP=', NP
       call neko_error
    end if
    ALPHAD = real(ALPHA, kind=xp)
    BETAD = real(BETA, kind=xp)
    call ZWGJD(ZD, WD, NP, ALPHAD, BETAD)
    do I = 1, NP
       Z(I) = real(ZD(I), kind=rp)
       W(I) = real(WD(I), kind=rp)
    end do
  end subroutine ZWGJ

  subroutine ZWGJD(Z, W, NP, ALPHA, BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!
!--------------------------------------------------------------------

    integer, intent(in) :: NP
    real(kind=xp), intent(inout) :: Z(NP), W(NP)
    real(kind=xp), intent(in) :: ALPHA, BETA

    real(kind=xp) :: DN, ONE, TWO, APB
    real(kind=xp) :: FAC1, FAC2, FAC3, FNORM
    real(kind=xp) :: RCOEF, P, PD, PM1, PDM1, PM2, PDM2
    real(kind=xp) :: DNP1, DNP2
    integer :: N, NP1, NP2, I

    N = NP - 1
    DN = real(N, kind=xp)
    ONE = 1.0_xp
    TWO = 2.0_xp
    APB = ALPHA + BETA

    if (NP .le. 0) then
       write (stderr, *) 'ZWGJD: Minimum number of Gauss points is 1', np
       call neko_error
    end if
    if ((ALPHA .le. -ONE) .or. (BETA .le. -ONE)) then
       write (stderr, *) 'ZWGJD: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    if (NP .eq. 1) then
       Z(1) = (BETA-ALPHA)/(APB+TWO)
       W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO) &
            * TWO**(APB+ONE)
       return
    end if

    call JACG (Z, NP, ALPHA, BETA)

    NP1 = N+1
    NP2 = N+2
    DNP1 = real(NP1, kind=xp)
    DNP2 = real(NP2, kind=xp)
    FAC1 = DNP1+ALPHA+BETA+ONE
    FAC2 = FAC1+DNP1
    FAC3 = FAC2+ONE
    FNORM = PNORMJ(NP1, ALPHA, BETA)
    RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
    do I = 1, NP
       call JACOBF (P, PD, PM1, PDM1, PM2, PDM2, NP2, ALPHA, BETA, Z(I))
       W(I) = -RCOEF/(P*PDM1)
    end do
  end subroutine ZWGJD

  subroutine ZWGLJ(Z, W, NP, ALPHA, BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.
!
!--------------------------------------------------------------------

    integer, intent(in) :: NP
    real(kind=rp), intent(inout) :: Z(NP), W(NP)
    real(kind=rp), intent(in) :: ALPHA, BETA

    integer, parameter :: NMAX = 84
    integer, parameter :: NZD = NMAX

    real(kind=xp) ZD(NZD), WD(NZD), ALPHAD, BETAD
    integer :: I, NPMAX

    NPMAX = NZD
    if (NP .gt. NPMAX) then
       write (stderr, *) 'Too large polynomial degree in ZWGLJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here NP=', NP
       call neko_error
    end if
    ALPHAD = real(ALPHA, kind=xp)
    BETAD = real(BETA, kind=xp)
    call ZWGLJD(ZD, WD, NP, ALPHAD, BETAD)
    do I = 1, NP
       Z(I) = real(ZD(I), kind=rp)
       W(I) = real(WD(I), kind=rp)
    end do
  end subroutine ZWGLJ

  subroutine ZWGLJD(Z, W, NP, ALPHA, BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!
!--------------------------------------------------------------------

    integer, intent(in) :: NP
    real(kind=xp), intent(inout) :: Z(NP), W(NP)
    real(kind=xp), intent(in) :: ALPHA, BETA

    real(kind=xp) :: ONE, TWO, ALPG, BETG
    real(kind=xp) :: P, PD, PM1, PDM1, PM2, PDM2
    integer :: N, NM1, I

    N = NP-1
    NM1 = N-1
    ONE = 1.0_xp
    TWO = 2.0_xp

    if (NP .le. 1) then
       write (stderr, *) 'ZWGLJD: Minimum number of Gauss-Lobatto points is 2'
       write (stderr, *) 'ZWGLJD: alpha, beta:', alpha, beta, np
       call neko_error
    end if
    if ((ALPHA .le. -ONE) .or. (BETA .le. -ONE)) then
       write (stderr, *) 'ZWGLJD: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    if (NM1 .gt. 0) then
       ALPG = ALPHA+ONE
       BETG = BETA+ONE
       call ZWGJD(Z(2), W(2), NM1, ALPG, BETG)
    end if
    Z(1) = -ONE
    Z(NP) = ONE
    do I = 2, NP-1
       W(I) = W(I)/(ONE-Z(I)**2)
    end do
    call JACOBF (P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(1))
    W(1) = ENDW1(N, ALPHA, BETA)/(TWO*PD)
    call JACOBF (P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(NP))
    W(NP) = ENDW2(N, ALPHA, BETA)/(TWO*PD)

  end subroutine ZWGLJD

  real(kind=xp) function ENDW1(N, ALPHA, BETA)

    real(kind=xp), intent(in) :: ALPHA, BETA
    integer, intent(in) :: N

    real(kind=xp) :: ZERO, ONE, TWO, THREE, FOUR, APB
    real(kind=xp) :: F1, F2, F3, FINT1, FINT2
    real(kind=xp) :: A1, A2, A3, DI, ABN, ABNN
    integer :: I

    ZERO = 0.0_xp
    ONE = 1.0_xp
    TWO = 2.0_xp
    THREE = 3.0_xp
    FOUR = 4.0_xp
    APB = ALPHA+BETA
    if (N .eq. 0) then
       ENDW1 = ZERO
       return
    end if
    F1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
    F1 = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
    if (N .eq. 1) then
       ENDW1 = F1
       return
    end if
    FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
    FINT1 = FINT1*TWO**(APB+TWO)
    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
    FINT2 = FINT2*TWO**(APB+THREE)
    F2 = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2) &
         * (APB+THREE)/FOUR
    if (N .eq. 2) then
       ENDW1 = F2
       return
    end if
    do I = 3, N
       DI = real(I-1, kind=xp)
       ABN = ALPHA+BETA+DI
       ABNN = ABN+DI
       A1 = -(TWO*(DI + ALPHA)*(DI + BETA)) / (ABN*ABNN*(ABNN + ONE))
       A2 = (TWO*(ALPHA - BETA)) / (ABNN*(ABNN + TWO))
       A3 = (TWO*(ABN + ONE)) / ((ABNN + TWO)*(ABNN + ONE))
       F3 = -(A2*F2 + A1*F1) / A3
       F1 = F2
       F2 = F3
    end do
    ENDW1 = F3
  end function ENDW1

  real(kind=xp) function ENDW2(N, ALPHA, BETA)

    real(kind=xp), intent(in) :: ALPHA, BETA
    integer, intent(in) :: N

    real(kind=xp) :: ZERO, ONE, TWO, THREE, FOUR, APB
    real(kind=xp) :: F1, F2, F3, FINT1, FINT2
    real(kind=xp) :: A1, A2, A3, DI, ABN, ABNN
    integer :: I

    ZERO = 0.0_xp
    ONE = 1.0_xp
    TWO = 2.0_xp
    THREE = 3.0_xp
    FOUR = 4.0_xp
    APB = ALPHA+BETA
    if (N .eq. 0) then
       ENDW2 = ZERO
       return
    end if
    F1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
    F1 = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
    if (N .eq. 1) then
       ENDW2 = F1
       return
    end if
    FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
    FINT1 = FINT1*TWO**(APB+TWO)
    FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
    FINT2 = FINT2*TWO**(APB+THREE)
    F2 = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2) &
         * (APB+THREE)/FOUR
    if (N .eq. 2) then
       ENDW2 = F2
       return
    end if
    do I = 3, N
       DI = ((I-1))
       ABN = ALPHA+BETA+DI
       ABNN = ABN+DI
       A1 = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
       A2 = (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
       A3 = (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
       F3 = -(A2*F2+A1*F1)/A3
       F1 = F2
       F2 = F3
    end do
    ENDW2 = F3
  end function ENDW2

  real(kind=xp) function GAMMAF (X)
    real(kind=xp), intent(in) :: X
    real(kind=xp) :: ZERO, HALF, ONE, TWO, FOUR, PI

    ZERO = 0.0_xp
    HALF = 0.5_xp
    ONE = 1.0_xp
    TWO = 2.0_xp
    FOUR = 4.0_xp
    PI = FOUR*atan(ONE)
    GAMMAF = ONE
    if (X .eq. -HALF) GAMMAF = -TWO*sqrt(PI)
    if (X .eq. HALF) GAMMAF = sqrt(PI)
    if (X .eq. ONE ) GAMMAF = ONE
    if (X .eq. TWO ) GAMMAF = ONE
    if (X .eq. 1.5_xp) GAMMAF = sqrt(PI) / 2.0_xp
    if (X .eq. 2.5_xp) GAMMAF = 1.5_xp * sqrt(PI) / 2.0_xp
    if (X .eq. 3.5_xp) GAMMAF = 0.5_xp * (2.5_xp * (1.5_xp * sqrt(PI)))
    if (X .eq. 3.0_xp) GAMMAF = 2.0_xp
    if (X .eq. 4.0_xp) GAMMAF = 6.0_xp
    if (X .eq. 5.0_xp) GAMMAF = 24.0_xp
    if (X .eq. 6.0_xp) GAMMAF = 120.0_xp
  end function GAMMAF

  real(kind=xp) function PNORMJ (N, ALPHA, BETA)
    real(kind=xp), intent(in) :: ALPHA, BETA
    integer, intent(in) :: N

    real(kind=xp) :: ONE, TWO, DN, DINDX
    real(kind=xp) :: CONST, PROD, FRAC
    integer :: I

    ONE = 1.0_xp
    TWO = 2.0_xp
    DN = real(N, kind=xp)
    CONST = ALPHA + BETA + ONE
    if (N .le. 1) then
       PROD = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
       PROD = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
       PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
       return
    end if
    PROD = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
    PROD = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
    PROD = PROD*(ONE+ALPHA)*(TWO+ALPHA)
    PROD = PROD*(ONE+BETA)*(TWO+BETA)
    do I = 3, N
       DINDX = real(I, kind=xp)
       FRAC = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
       PROD = PROD*FRAC
    end do
    PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
  end function PNORMJ

  subroutine JACG (XJAC, NP, ALPHA, BETA)
!--------------------------------------------------------------------
!
!     Compute NP Gauss points XJAC, which are the zeros of the
!     Jacobi polynomial J(NP) with parameters ALPHA and BETA.
!     ALPHA and BETA determines the specific type of Gauss points.
!     Examples:
!     ALPHA = BETA =  0.0  ->  Legendre points
!     ALPHA = BETA = -0.5  ->  Chebyshev points
!
!--------------------------------------------------------------------


    integer, intent(in) :: NP
    real(kind=xp), intent(inout) :: XJAC(NP)
    real(kind=xp), intent(in) :: ALPHA, BETA

    integer, parameter :: KSTOP = 10
    real(kind=rp), parameter :: EPS = 1.0e-12_rp

    real(kind=xp) :: one, DTH, X, X1, X2, XLAST, DELX, XMIN
    real(kind=xp) :: P, PD, PM1, PDM1, PM2, PDM2
    real(kind=xp) :: RECSUM, SWAP
    integer :: I, J, K, N, JM, JMIN

    N = NP-1
    one = 1.0_xp
    DTH = 4.0_xp*atan(one)/(2.0_xp*real(N, kind=xp) + 2.0_xp)
    do J = 1, NP
       if (J .eq. 1) then
          X = cos((2.0_xp*(real(J, kind=xp)-1.0_xp)+1.0_xp)*DTH)
       else
          X1 = cos((2.0_xp*(real(J, kind=xp)-1.0_xp)+1.0_xp)*DTH)
          X2 = XLAST
          X = (X1+X2)/2.0_xp
       end if
       do K = 1, KSTOP
          call JACOBF (P, PD, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, X)
          RECSUM = 0.0_xp
          JM = J-1
          do I = 1, JM
             RECSUM = RECSUM+1.0_xp/(X-XJAC(NP-I+1))
          end do
          DELX = -P/(PD-RECSUM*P)
          X = X+DELX
          if (abs(DELX) .lt. EPS) exit
       end do

       XJAC(NP-J+1) = X
       XLAST = X
    end do
    do I = 1, NP
       XMIN = 2.
       do J = I, NP
          if (XJAC(J) .lt. XMIN) then
             XMIN = XJAC(J)
             JMIN = J
          end if
       end do
       if (JMIN .ne. I) then
          SWAP = XJAC(I)
          XJAC(I) = XJAC(JMIN)
          XJAC(JMIN) = SWAP
       end if
    end do
  end subroutine JACG

  subroutine JACOBF (POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2, N, ALP, BET, X)
!--------------------------------------------------------------------
!
!     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
!     of degree N at X.
!
!--------------------------------------------------------------------

    real(kind=xp), intent(inout) :: POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2
    real(kind=xp), intent(in) :: ALP, BET, X
    integer, intent(in) :: N

    real(kind=xp) :: APB, POLYL, PDERL, POLYN, PDERN
    real(kind=xp) :: PSAVE, PDSAVE
    real(kind=xp) :: A1, A2, A3, A4, B3
    real(kind=xp) :: DK
    integer :: K

    APB = ALP+BET
    POLY = 1.0_xp
    PDER = 0.0_xp
    if (N .eq. 0) return
    POLYL = POLY
    PDERL = PDER
    POLY = (ALP-BET+(APB+2.0_xp)*X)/2.0_xp
    PDER = (APB+2.0_xp)/2.0_xp
    if (N .eq. 1) return
    do K = 2, N
       DK = real(K, kind=xp)
       A1 = 2.0_xp*DK*(DK+APB)*(2.0_xp*DK+APB-2.0_xp)
       A2 = (2.0_xp*DK+APB-1.0_xp)*(ALP**2-BET**2)
       B3 = (2.0_xp*DK+APB-2.0_xp)
       A3 = B3*(B3+1.0_xp)*(B3+2.0_xp)
       A4 = 2.0_xp*(DK+ALP-1.0_xp)*(DK+BET-1.0_xp)*(2.0_xp*DK+APB)
       POLYN = ((A2+A3*X)*POLY-A4*POLYL)/A1
       PDERN = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
       PSAVE = POLYL
       PDSAVE = PDERL
       POLYL = POLY
       POLY = POLYN
       PDERL = PDER
       PDER = PDERN
    end do
    POLYM1 = POLYL
    PDERM1 = PDERL
    POLYM2 = PSAVE
    PDERM2 = PDSAVE
  end subroutine JACOBF

  real(kind=xp) function HGJ (II, Z, ZGJ, NP, ALPHA, BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGJ through
!     the NP Gauss Jacobi points ZGJ at the point Z.
!     Single precision version.
!
!---------------------------------------------------------------------
    integer, intent(in) :: NP, II
    real(kind=xp), intent(in) :: Z, ZGJ(NP), ALPHA, BETA

    integer, parameter :: NMAX = 84
    integer, parameter :: NZD = NMAX

    real(kind=xp) ZD, ZGJD(NZD)
    integer :: I, NPMAX

    NPMAX = NZD
    if (NP .gt. NPMAX) then
       write (stderr, *) 'Too large polynomial degree in HGJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here NP=', NP
       call neko_error
    end if
    ZD = Z
    do I = 1, NP
       ZGJD(I) = ZGJ(I)
    end do
    HGJ = HGJD (II, ZD, ZGJD, NP, ALPHA, BETA)
  end function HGJ

  real(kind=xp) function HGJD (II, Z, ZGJ, NP, ALPHA, BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGJD through
!     the NZ Gauss-Lobatto Jacobi points ZGJ at the point Z.
!     Double precision version.
!
!---------------------------------------------------------------------

    integer, intent(in) :: NP, II
    real(kind=xp), intent(in) :: Z, ZGJ(NP), ALPHA, BETA

    real(kind=xp) :: EPS, ONE, ZI, DZ
    real(kind=xp) :: PZ, PDZ, PZI, PDZI, PM1, PDM1, PM2, PDM2

    EPS = 1.0e-5_xp
    ONE = 1.0_xp
    ZI = ZGJ(II)
    DZ = Z-ZI
    if (abs(DZ) .lt. EPS) then
       HGJD = ONE
       return
    end if
    call JACOBF (PZI, PDZI, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, ZI)
    call JACOBF (PZ, PDZ, PM1, PDM1, PM2, PDM2, NP, ALPHA, BETA, Z)
    HGJD = PZ/(PDZI*(Z-ZI))
  end function HGJD

  real(kind=xp) function HGLJ (II, Z, ZGLJ, NP, ALPHA, BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGLJ through
!     the NZ Gauss-Lobatto Jacobi points ZGLJ at the point Z.
!     Single precision version.
!
!---------------------------------------------------------------------
    integer, parameter :: NMAX = 84
    integer, parameter :: NZD = NMAX

    integer, intent(in) :: NP, II
    real(kind=xp), intent(in) :: Z, ZGLJ(NP), ALPHA, BETA
    real(kind=xp) ZD, ZGLJD(NZD)
    integer :: I, NPMAX

    NPMAX = NZD
    if (NP .gt. NPMAX) then
       write (stderr, *) 'Too large polynomial degree in HGLJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here NP=', NP
       call neko_error
    end if
    ZD = Z
    do I = 1, NP
       ZGLJD(I) = ZGLJ(I)
    end do
    HGLJ = HGLJD (II, ZD, ZGLJD, NP, ALPHA, BETA)
  end function HGLJ

  real(kind=xp) function HGLJD (I, Z, ZGLJ, NP, ALPHA, BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGLJD through
!     the NZ Gauss-Lobatto Jacobi points ZJACL at the point Z.
!     Double precision version.
!
!---------------------------------------------------------------------

    integer, intent(in) :: NP, I
    real(kind=xp), intent(in) :: Z, ZGLJ(NP), ALPHA, BETA

    real(kind=xp) :: EPS, ONE, ZI, DZ, DN
    real(kind=xp) :: P, PD, PI, PDI, PM1, PDM1, PM2, PDM2
    real(kind=xp) :: EIGVAL, CONST
    integer :: N

    EPS = 1.0e-5_xp
    ONE = 1.0_xp
    ZI = ZGLJ(I)
    DZ = Z-ZI
    if (abs(DZ) .lt. EPS) then
       HGLJD = ONE
       return
    end if
    N = NP-1
    DN = real(N, kind=xp)
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
    call JACOBF (PI, PDI, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, ZI)
    CONST = EIGVAL*PI+ALPHA*(ONE+ZI)*PDI-BETA*(ONE-ZI)*PDI
    call JACOBF (P, PD, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z)
    HGLJD = (ONE-Z**2)*PD/(CONST*(Z-ZI))
  end function HGLJD

  subroutine DGJ (D, DT, Z, NZ, NZD, ALPHA, BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!
!-----------------------------------------------------------------
    integer, intent(in) :: NZ, NZD
    real(kind=xp), intent(inout) :: D(NZD, NZD), DT(NZD, NZD)
    real(kind=xp), intent(in) :: Z(NZ), ALPHA, BETA

    integer, parameter :: NMAX = 84
    integer, parameter :: NZDD = NMAX

    real(kind=xp) :: DD(NZDD, NZDD), DTD(NZDD, NZDD), ZD(NZDD)
    integer :: I, J

    if (NZ .le. 0) then
       write (stderr, *) 'DGJ: Minimum number of Gauss points is 1'
       call neko_error
    end if
    if (NZ .gt. NMAX) then
       write (stderr, *) 'Too large polynomial degree in DGJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here Nz=', Nz
       call neko_error
    end if
    if ((ALPHA .le. -1.0_xp) .or. (BETA .le. -1.0_xp)) then
       write (stderr, *) 'DGJ: Alpha and Beta must be greater than -1'
       call neko_error
    end if
    do I = 1, NZ
       ZD(I) = Z(I)
    end do
    call DGJD (DD, DTD, ZD, NZ, NZDD, ALPHA, BETA)
    do I = 1, NZ
       do J = 1, NZ
          D(I, J) = DD(I, J)
          DT(I, J) = DTD(I, J)
       end do
    end do
  end subroutine DGJ

  subroutine DGJD (D, DT, Z, NZ, NZD, ALPHA, BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!
!-----------------------------------------------------------------

    integer, intent(in) :: NZ, NZD
    real(kind=xp), intent(inout) :: D(NZD, NZD), DT(NZD, NZD)
    real(kind=xp), intent(in) :: Z(NZ), ALPHA, BETA

    real(kind=xp) :: DN, ONE, TWO
    real(kind=xp) :: PDI, PDJ, PI, PJ, PM1, PDM1, PM2, PDM2
    integer :: I, J, N

    N = NZ-1
    DN = real(N, kind=xp)
    ONE = 1.0_xp
    TWO = 2.0_xp

    if (NZ .le. 1) then
       write (stderr, *) 'DGJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
    end if
    if ((ALPHA .le. -ONE) .or. (BETA .le. -ONE)) then
       write (stderr, *) 'DGJD: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    do I = 1, NZ
       do J = 1, NZ
          call JACOBF (PI, PDI, PM1, PDM1, PM2, PDM2, NZ, ALPHA, BETA, Z(I))
          call JACOBF (PJ, PDJ, PM1, PDM1, PM2, PDM2, NZ, ALPHA, BETA, Z(J))
          if (I .ne. J) D(I, J) = PDI/(PDJ*(Z(I)-Z(J)))
          if (I .eq. J) D(I, J) = ((ALPHA+BETA+TWO)*Z(I)+ALPHA-BETA)/ &
               (TWO*(ONE-Z(I)**2))
          DT(J, I) = D(I, J)
       end do
    end do
  end subroutine DGJD

  subroutine DGLJ (D, DT, Z, NZ, NZD, ALPHA, BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!
!-----------------------------------------------------------------
    integer, parameter :: NMAX = 84
    integer, parameter :: NZDD = NMAX
    integer, intent(in) :: NZ, NZD
    real(kind=xp), intent(inout) :: D(NZD, NZD), DT(NZD, NZD)
    real(kind=xp), intent(in) :: Z(NZ), ALPHA, BETA

    real(kind=xp) :: DD(NZDD, NZDD), DTD(NZDD, NZDD), ZD(NZDD)
    integer :: I, J

    if (NZ .le. 1) then
       write (stderr, *) 'DGLJ: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
    end if
    if (NZ .gt. NMAX) then
       write (stderr, *) 'Too large polynomial degree in DGLJ'
       write (stderr, *) 'Maximum polynomial degree is', NMAX
       write (stderr, *) 'Here NZ=', NZ
       call neko_error
    end if
    if ((ALPHA .le. -1.0_xp) .or. (BETA .le. -1.0_xp)) then
       write (stderr, *) 'DGLJ: Alpha and Beta must be greater than -1'
       call neko_error
    end if
    do I = 1, NZ
       ZD(I) = Z(I)
    end do
    call DGLJD (DD, DTD, ZD, NZ, NZDD, ALPHA, BETA)
    do I = 1, NZ
       do J = 1, NZ
          D(I, J) = DD(I, J)
          DT(I, J) = DTD(I, J)
       end do
    end do
  end subroutine DGLJ

  subroutine DGLJD (D, DT, Z, NZ, NZD, ALPHA, BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!
!-----------------------------------------------------------------

    integer, intent(in) :: NZ, NZD
    real(kind=xp), intent(inout) :: D(NZD, NZD), DT(NZD, NZD)
    real(kind=xp), intent(in) :: Z(NZ), ALPHA, BETA

    real(kind=xp) :: DN, ONE, TWO, EIGVAL
    real(kind=xp) :: PDI, PDJ, PI, PJ, PM1, PDM1, PM2, PDM2
    real(kind=xp) :: CI, CJ
    integer :: I, J, N

    N = NZ-1
    DN = real(N, kind=xp)
    ONE = 1.0_xp
    TWO = 2.0_xp
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

    if (NZ .le. 1) then
       write (stderr, *) 'DGLJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
    end if
    if ((ALPHA .le. -ONE) .or. (BETA .le. -ONE)) then
       write (stderr, *) 'DGLJD: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    do I = 1, NZ
       do J = 1, NZ
          call JACOBF (PI, PDI, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(I))
          call JACOBF (PJ, PDJ, PM1, PDM1, PM2, PDM2, N, ALPHA, BETA, Z(J))
          CI = EIGVAL*PI-(BETA*(ONE-Z(I))-ALPHA*(ONE+Z(I)))*PDI
          CJ = EIGVAL*PJ-(BETA*(ONE-Z(J))-ALPHA*(ONE+Z(J)))*PDJ
          if (I .ne. J) D(I, J) = CI/(CJ*(Z(I)-Z(J)))
          if ((I .eq. J) .and. (I .ne. 1) .and. (I .ne. NZ)) &
               D(I, J) = (ALPHA*(ONE+Z(I))-BETA*(ONE-Z(I)))/ &
               (TWO*(ONE-Z(I)**2))
          if ((I .eq. J) .and. (I .eq. 1)) &
               D(I, J) = (EIGVAL+ALPHA)/(TWO*(BETA+TWO))
          if ((I .eq. J) .and. (I .eq. NZ)) &
               D(I, J) = -(EIGVAL+BETA)/(TWO*(ALPHA+TWO))
          DT(J, I) = D(I, J)
       end do
    end do
  end subroutine DGLJD

  subroutine DGLL (D, DT, Z, NZ, NZD)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Legendre points Z.
!     Note: D and DT are square matrices.
!
!-----------------------------------------------------------------

    integer, intent(in) :: NZ, NZD
    real(kind=rp), intent(inout) :: D(NZD, NZD), DT(NZD, NZD)
    real(kind=rp), intent(in) :: Z(NZ)

    integer, parameter :: NMAX = 84

    real(kind=xp) :: d0, FN
    integer :: I, J, N

    N = NZ-1
    if (NZ .gt. NMAX) then
       write (stderr, *) 'Subroutine DGLL'
       write (stderr, *) 'Maximum polynomial degree =', NMAX
       write (stderr, *) 'Polynomial degree         =', NZ
    end if
    if (NZ .eq. 1) then
       D(1, 1) = 0.0_rp
       return
    end if
    FN = real(N, kind=xp)
    d0 = FN*(FN + 1.0_xp)/4.0_xp
    do I = 1, NZ
       do J = 1, NZ
          D(I, J) = 0.0_rp
          if (I .ne. J) D(I, J) = PNLEG(real(Z(I), xp), N)/ &
               (PNLEG(real(Z(J), xp), N)*(Z(I)-Z(J)))
          if ((I .eq. J) .and. (I .eq. 1)) D(I, J) = -d0
          if ((I .eq. J) .and. (I .eq. NZ)) D(I, J) = d0
          DT(J, I) = D(I, J)
       end do
    end do
  end subroutine DGLL

  real(kind=xp) function HGLL (I, Z, ZGLL, NZ)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant L through
!     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
!
!---------------------------------------------------------------------

    integer, intent(in) :: I, NZ
    real(kind=xp), intent(in) :: ZGLL(NZ), Z

    real(kind=xp) :: EPS, DZ
    real(kind=xp) :: ALFAN
    integer :: N

    EPS = 1.0e-5_xp
    DZ = Z - ZGLL(I)
    if (abs(DZ) .lt. EPS) then
       HGLL = 1.0_xp
       return
    end if
    N = NZ - 1
    ALFAN = real(N, kind=xp)*(real(N, kind=xp) + 1.0_xp)
    HGLL = -(1.0_xp - Z*Z)*PNDLEG(Z, N) / &
         (ALFAN*PNLEG(ZGLL(I), N)*(Z - ZGLL(I)))
  end function HGLL

  real(kind=xp) function HGL (I, Z, ZGL, NZ)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGL through
!     the NZ Gauss Legendre points ZGL at the point Z.
!
!---------------------------------------------------------------------
    integer, intent(in) :: I, NZ
    real(kind=xp), intent(in) :: ZGL(NZ), Z
    real(kind=xp) :: EPS, DZ

    integer :: N

    EPS = 1.0e-5_xp
    DZ = Z - ZGL(I)
    if (abs(DZ) .lt. EPS) then
       HGL = 1.0_xp
       return
    end if
    N = NZ-1
    HGL = PNLEG(Z, NZ)/(PNDLEG(ZGL(I), NZ)*(Z-ZGL(I)))
  end function HGL

  real(kind=xp) function PNLEG (Z, N)
!---------------------------------------------------------------------
!
!     Compute the value of the Nth order Legendre polynomial at Z.
!     (Simpler than JACOBF)
!     Based on the recursion formula for the Legendre polynomials.
!
!---------------------------------------------------------------------
!
!     This next statement is to overcome the underflow bug in the i860.
!     It can be removed at a later date.  11 Aug 1990   pff.
!

    real(kind=xp), intent(in) :: Z
    integer, intent(in) :: N

    real(kind=xp) :: P1, P2, P3, FK
    integer :: K

    P1 = 1.0_xp
    if (N .eq. 0) then
       PNLEG = P1
       return
    end if
    P2 = Z
    P3 = P2
    do K = 1, N-1
       FK = real(K, kind=xp)
       P3 = ((2.0_xp*FK + 1.0_xp)*Z*P2 - FK*P1)/(FK + 1.0_xp)
       P1 = P2
       P2 = P3
    end do
    PNLEG = P3
  end function PNLEG

  !> Evaluate Legendre polynomials of degrees 0-N at point x
  !! and store in array L
  subroutine legendre_poly(L, x, N)
    integer, intent(in) :: N
    real(kind=rp), intent(inout) :: L(0:N)
    real(kind=rp), intent(in) :: x

    real(kind=rp) :: DJ
    integer :: j

    L(0) = 1.0_rp
    if (N .eq. 0) return
    L(1) = x

    do j = 1, N-1
       DJ = real(j, kind=rp)
       L(j+1) = ((2.0_rp*DJ + 1.0_rp)*x*L(j) - DJ*L(j-1)) / (DJ + 1.0_rp)
    end do
  end subroutine legendre_poly

  real(kind=xp) function PNDLEG (Z, N)
!----------------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     (Simpler than JACOBF)
!     Based on the recursion formula for the Legendre polynomials.
!
!----------------------------------------------------------------------

    real(kind=xp), intent(in) :: Z
    integer, intent(in) :: N

    real(kind=xp) :: P1, P2, P3, P1D, P2D, P3D, FK
    integer :: K

    if (N .eq. 0) then
       PNDLEG = 0.0_xp
       return
    end if

    P1 = 1.0_xp
    P2 = Z
    P1D = 0.0_xp
    P2D = 1.0_xp
    P3D = 1.0_xp
    do K = 1, N-1
       FK = real(K, kind=xp)
       P3 = ((2.0_xp*FK + 1.0_xp)*Z*P2 - FK*P1) / (FK + 1.0_xp)
       P3D = ((2.0_xp*FK + 1.0_xp)*P2 + (2.0_xp*FK + 1.0_xp)*Z*P2D - FK*P1D) / &
            (FK + 1.0_xp)
       P1 = P2
       P2 = P3
       P1D = P2D
       P2D = P3D
    end do
    PNDLEG = P3D
  end function PNDLEG

  subroutine DGLLGL(D, DT, ZM1, ZM2, IM12, NZM1, NZM2, ND1, ND2)
!-----------------------------------------------------------------------
!
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Legendre mesh (M1), and evaluate its
!     derivative on a Guass Legendre mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLLGL).
!     Note: D and DT are rectangular matrices.
!
!-----------------------------------------------------------------------

    integer, intent(in) :: NZM1, NZM2, ND1, ND2
    real(kind=xp), intent(inout) :: D(ND2, ND1), DT(ND1, ND2)
    real(kind=xp), intent(in) :: ZM1(ND1), ZM2(ND2), IM12(ND2, ND1)

    real(kind=xp) EPS, ZP, ZQ
    integer :: IP, JQ, NM1

    if (NZM1 .eq. 1) then
       D (1, 1) = 0.0_xp
       DT(1, 1) = 0.0_xp
       return
    end if
    EPS = 1.0E-6_xp
    NM1 = NZM1-1
    do IP = 1, NZM2
       do JQ = 1, NZM1
          ZP = ZM2(IP)
          ZQ = ZM1(JQ)
          if ((abs(ZP) .lt. EPS) .and. (abs(ZQ) .lt. EPS)) then
             D(IP, JQ) = 0.0_xp
          else
             D(IP, JQ) = (PNLEG(ZP, NM1) / PNLEG(ZQ, NM1) - IM12(IP, JQ)) / &
                  (ZP - ZQ)
          end if
          DT(JQ, IP) = D(IP, JQ)
       end do
    end do
  end subroutine DGLLGL

  subroutine DGLJGJ(D, DT, ZGL, ZG, IGLG, NPGL, NPG, ND1, ND2, ALPHA, BETA)
!-----------------------------------------------------------------------
!
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
!     derivative on a Guass Jacobi mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLJGJ).
!     Note: D and DT are rectangular matrices.
!     Single precision version.
!
!-----------------------------------------------------------------------
    integer, intent(in) :: NPGL, NPG, ND1, ND2
    real(kind=xp), intent(inout) :: D(ND2, ND1), DT(ND1, ND2)
    real(kind=xp), intent(in) :: ZGL(ND1), ZG(ND2), IGLG(ND2, ND1), ALPHA, BETA

    integer, parameter :: NMAX = 84
    integer, parameter :: NDD = NMAX

    real(kind=xp) DD(NDD, NDD), DTD(NDD, NDD)
    real(kind=xp) ZGD(NDD), ZGLD(NDD), IGLGD(NDD, NDD)
    integer :: I, J

    if (NPGL .le. 1) then
       write(6, *) 'DGLJGJ: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
    end if
    if (NPGL .gt. NMAX) then
       write(6, *) 'Polynomial degree too high in DGLJGJ'
       write(6, *) 'Maximum polynomial degree is', NMAX
       write(6, *) 'Here NPGL=', NPGL
       call neko_error
    end if
    if ((ALPHA .le. -1.0_xp) .or. (BETA .le. -1.0_xp)) then
       write(6, *) 'DGLJGJ: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    do I = 1, NPG
       ZGD(I) = ZG(I)
       do J = 1, NPGL
          IGLGD(I, J) = IGLG(I, J)
       end do
    end do
    do I = 1, NPGL
       ZGLD(I) = ZGL(I)
    end do
    call DGLJGJD(DD, DTD, ZGLD, ZGD, IGLGD, NPGL, NPG, NDD, NDD, ALPHA, BETA)
    do I = 1, NPG
       do J = 1, NPGL
          D(I, J) = DD(I, J)
          DT(J, I) = DTD(J, I)
       end do
    end do
  end subroutine DGLJGJ

  subroutine DGLJGJD(D, DT, ZGL, ZG, IGLG, NPGL, NPG, ND1, ND2, ALPHA, BETA)
!-----------------------------------------------------------------------
!
!     Compute the (one-dimensional) derivative matrix D and its
!     transpose DT associated with taking the derivative of a variable
!     expanded on a Gauss-Lobatto Jacobi mesh (M1), and evaluate its
!     derivative on a Guass Jacobi mesh (M2).
!     Need the one-dimensional interpolation operator IM12
!     (see subroutine IGLJGJ).
!     Note: D and DT are rectangular matrices.
!     Double precision version.
!
!-----------------------------------------------------------------------

    integer, intent(in) :: NPGL, NPG, ND1, ND2
    real(kind=xp), intent(inout) :: D(ND2, ND1), DT(ND1, ND2)
    real(kind=xp), intent(in) :: ZGL(ND1), ZG(ND2), IGLG(ND2, ND1), ALPHA, BETA

    real(kind=xp) :: EPS, ONE, TWO, EIGVAL, DN
    real(kind=xp) :: PDI, PDJ, PI, PJ, PM1, PDM1, PM2, PDM2
    real(kind=xp) :: DZ, FACI, FACJ, CONST
    integer :: I, J, NGL

    if (NPGL .le. 1) then
       write(6, *) 'DGLJGJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
    end if
    if ((ALPHA .le. -1.) .or. (BETA .le. -1.)) then
       write(6, *) 'DGLJGJD: Alpha and Beta must be greater than -1'
       call neko_error
    end if

    EPS = 1.0e-6_xp
    ONE = 1.0_xp
    TWO = 2.0_xp
    NGL = NPGL-1
    DN = real(NGL, kind=xp)
    EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

    do I = 1, NPG
       do J = 1, NPGL
          DZ = abs(ZG(I)-ZGL(J))
          if (DZ .lt. EPS) then
             D(I, J) = (ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I)))/ &
                  (TWO*(ONE-ZG(I)**2))
          else
             call JACOBF (PI, PDI, PM1, PDM1, PM2, PDM2, NGL, ALPHA, BETA, ZG(I))
             call JACOBF (PJ, PDJ, PM1, PDM1, PM2, PDM2, NGL, ALPHA, BETA, ZGL(J))
             FACI = ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I))
             FACJ = ALPHA*(ONE+ZGL(J))-BETA*(ONE-ZGL(J))
             CONST = EIGVAL*PJ+FACJ*PDJ
             D(I, J) = ((EIGVAL*PI+FACI*PDI)*(ZG(I)-ZGL(J)) &
                  -(ONE-ZG(I)**2)*PDI)/(CONST*(ZG(I)-ZGL(J))**2)
          end if
          DT(J, I) = D(I, J)
       end do
    end do
  end subroutine DGLJGJD

  subroutine IGLM (I12, IT12, Z1, Z2, NZ1, NZ2, ND1, ND2)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss Legendre points.
!     Z2 : NZ2 points on mesh M.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NZ1, NZ2, ND1, ND2
    real(kind=xp), intent(inout) :: I12(ND2, ND1), IT12(ND1, ND2)
    real(kind=xp), intent(in) :: Z1(ND1), Z2(ND2)
    real(kind=xp) :: ZI
    integer :: I, J

    if (NZ1 .eq. 1) then
       I12 (1, 1) = 1.0_xp
       IT12(1, 1) = 1.0_xp
       return
    end if
    do I = 1, NZ2
       ZI = Z2(I)
       do J = 1, NZ1
          I12 (I, J) = HGL(J, ZI, Z1, NZ1)
          IT12(J, I) = I12(I, J)
       end do
    end do
  end subroutine IGLM

  subroutine IGLLM (I12, IT12, Z1, Z2, NZ1, NZ2, ND1, ND2)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss-Lobatto Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss-Lobatto Legendre points.
!     Z2 : NZ2 points on mesh M.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NZ1, NZ2, ND1, ND2
    real(kind=xp), intent(inout) :: I12(ND2, ND1), IT12(ND1, ND2)
    real(kind=xp), intent(in) :: Z1(ND1), Z2(ND2)
    real(kind=xp) :: ZI
    integer :: I, J

    if (NZ1 .eq. 1) then
       I12 (1, 1) = 1.0_xp
       IT12(1, 1) = 1.0_xp
       return
    end if
    do I = 1, NZ2
       ZI = Z2(I)
       do J = 1, NZ1
          I12 (I, J) = HGLL(J, ZI, Z1, NZ1)
          IT12(J, I) = I12(I, J)
       end do
    end do
  end subroutine IGLLM

  subroutine IGJM (I12, IT12, Z1, Z2, NZ1, NZ2, ND1, ND2, ALPHA, BETA)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss Jacobi mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss Jacobi points.
!     Z2 : NZ2 points on mesh M.
!     Single precision version.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NZ1, NZ2, ND1, ND2
    real(kind=xp), intent(inout) :: I12(ND2, ND1), IT12(ND1, ND2)
    real(kind=xp), intent(in) :: Z1(ND1), Z2(ND2), ALPHA, BETA
    real(kind=xp) :: ZI
    integer :: I, J

    if (NZ1 .eq. 1) then
       I12 (1, 1) = 1.0_xp
       IT12(1, 1) = 1.0_xp
       return
    end if
    do I = 1, NZ2
       ZI = Z2(I)
       do J = 1, NZ1
          I12 (I, J) = HGJ(J, ZI, Z1, NZ1, ALPHA, BETA)
          IT12(J, I) = I12(I, J)
       end do
    end do
  end subroutine IGJM

  subroutine IGLJM (I12, IT12, Z1, Z2, NZ1, NZ2, ND1, ND2, ALPHA, BETA)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss-Lobatto Jacobi mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss-Lobatto Jacobi points.
!     Z2 : NZ2 points on mesh M.
!     Single precision version.
!
!--------------------------------------------------------------------
    integer, intent(in) :: NZ1, NZ2, ND1, ND2
    real(kind=xp), intent(inout) :: I12(ND2, ND1), IT12(ND1, ND2)
    real(kind=xp), intent(in) :: Z1(ND1), Z2(ND2), ALPHA, BETA
    real(kind=xp) :: ZI
    integer :: I, J

    if (NZ1 .eq. 1) then
       I12 (1, 1) = 1.0_xp
       IT12(1, 1) = 1.0_xp
       return
    end if
    do I = 1, NZ2
       ZI = Z2(I)
       do J = 1, NZ1
          I12 (I, J) = HGLJ(J, ZI, Z1, NZ1, ALPHA, BETA)
          IT12(J, I) = I12(I, J)
       end do
    end do
  end subroutine IGLJM
end module speclib
