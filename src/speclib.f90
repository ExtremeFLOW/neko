!==============================================================================
!
!     LIBRARY ROUTINES FOR SPECTRAL METHODS
!
!     March 1989
!
!     For questions, comments or suggestions, please contact:
!
!     Einar Malvin Ronquist
!     Room 3-243
!     Department of Mechanical Engineering
!     Massachusetts Institute of Technology
!     77 Massachusetts Avenue
!     Cambridge, MA 0299
!     U.S.A.
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
  use num_types
  use utils


contains
!
!--------------------------------------------------------------------
      SUBROUTINE ZWGL (Z,W,NP)
!--------------------------------------------------------------------
!
!     Generate NP Gauss Legendre points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
!     The polynomial degree N=NP-1.
!     Z and W are in single precision, but all the arithmetic
!     operations are done in double precision.
!
!--------------------------------------------------------------------
      REAL(KIND=RP) Z(1),W(1), ALPHA, BETA
      ALPHA = 0.
      BETA  = 0.
      CALL ZWGJ (Z,W,NP,ALPHA,BETA)
      RETURN
      END SUBROUTINE

      SUBROUTINE ZWGLL (Z,W,NP)
!--------------------------------------------------------------------
!
!     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha=0,beta=0).
!     The polynomial degree N=NP-1.
!     Z and W are in single precision, but all the arithmetic
!     operations are done in double precision.
!
!--------------------------------------------------------------------
      REAL(KIND=RP) Z(1),W(1), ALPHA, BETA
      ALPHA = 0.
      BETA  = 0.
      CALL ZWGLJ (Z,W,NP,ALPHA,BETA)
      RETURN
      END SUBROUTINE

      SUBROUTINE ZWGJ (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.
!
!--------------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZD = NMAX)
      REAL(KIND=RP)  ZD(NZD),WD(NZD),ALPHAD,BETAD
      REAL(KIND=RP) Z(1),W(1),ALPHA,BETA

      NPMAX = NZD
      IF (NP.GT.NPMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in ZWGJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NP=',NP
         call neko_error
      ENDIF
      ALPHAD = ALPHA
      BETAD  = BETA
      CALL ZWGJD (ZD,WD,NP,ALPHAD,BETAD)
      DO 100 I=1,NP
         Z(I) = ZD(I)
         W(I) = WD(I)
 100  CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE ZWGJD (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!
!--------------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  Z(1),W(1),ALPHA,BETA

      N     = NP-1
      DN    = ((N))
      ONE   = 1.
      TWO   = 2.
      APB   = ALPHA+BETA

      IF (NP.LE.0) THEN
         WRITE (6,*) 'ZWGJD: Minimum number of Gauss points is 1',np
         call neko_error
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'ZWGJD: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      IF (NP.EQ.1) THEN
         Z(1) = (BETA-ALPHA)/(APB+TWO)
         W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO) &
              * TWO**(APB+ONE)
         RETURN
      ENDIF

      CALL JACG (Z,NP,ALPHA,BETA)

      NP1   = N+1
      NP2   = N+2
      DNP1  = ((NP1))
      DNP2  = ((NP2))
      FAC1  = DNP1+ALPHA+BETA+ONE
      FAC2  = FAC1+DNP1
      FAC3  = FAC2+ONE
      FNORM = PNORMJ(NP1,ALPHA,BETA)
      RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
      DO 100 I=1,NP
         CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z(I))
         W(I) = -RCOEF/(P*PDM1)
 100  CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE ZWGLJ (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Single precision version.
!
!--------------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZD = NMAX)
      REAL(KIND=RP)  ZD(NZD),WD(NZD),ALPHAD,BETAD
      REAL(KIND=RP) Z(1),W(1),ALPHA,BETA

      NPMAX = NZD
      IF (NP.GT.NPMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in ZWGLJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NP=',NP
         call neko_error
      ENDIF
      ALPHAD = ALPHA
      BETAD  = BETA
      CALL ZWGLJD (ZD,WD,NP,ALPHAD,BETAD)
      DO 100 I=1,NP
         Z(I) = ZD(I)
         W(I) = WD(I)
 100  CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE ZWGLJD (Z,W,NP,ALPHA,BETA)
!--------------------------------------------------------------------
!
!     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
!     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
!     The polynomial degree N=NP-1.
!     Double precision version.
!
!--------------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  Z(NP),W(NP),ALPHA,BETA

      N     = NP-1
      NM1   = N-1
      ONE   = 1.
      TWO   = 2.

      IF (NP.LE.1) THEN
       WRITE (6,*) 'ZWGLJD: Minimum number of Gauss-Lobatto points is 2'
       WRITE (6,*) 'ZWGLJD: alpha,beta:',alpha,beta,np
       call neko_error
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'ZWGLJD: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      IF (NM1.GT.0) THEN
         ALPG  = ALPHA+ONE
         BETG  = BETA+ONE
         CALL ZWGJD (Z(2),W(2),NM1,ALPG,BETG)
      ENDIF
      Z(1)  = -ONE
      Z(NP) =  ONE
      DO 100  I=2,NP-1
         W(I) = W(I)/(ONE-Z(I)**2)
 100  CONTINUE
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(1))
      W(1)  = ENDW1 (N,ALPHA,BETA)/(TWO*PD)
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(NP))
      W(NP) = ENDW2 (N,ALPHA,BETA)/(TWO*PD)

!      RETURN
      END SUBROUTINE

      REAL(KIND=RP)  FUNCTION ENDW1 (N,ALPHA,BETA)
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW1 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW1 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2) &
           * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW1 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW1  = F3
      RETURN
      END FUNCTION

      REAL(KIND=RP)  FUNCTION ENDW2 (N,ALPHA,BETA)
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  ALPHA,BETA
      ZERO  = 0.
      ONE   = 1.
      TWO   = 2.
      THREE = 3.
      FOUR  = 4.
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW2 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW2 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2) &
           * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW2 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = ((I-1))
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW2  = F3
      RETURN
      END FUNCTION

      REAL(KIND=RP)  FUNCTION GAMMAF (X)
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  X
      ZERO = 0.0
      HALF = 0.5
      ONE  = 1.0
      TWO  = 2.0
      FOUR = 4.0
      PI   = FOUR*ATAN(ONE)
      GAMMAF = ONE
      IF (X.EQ.-HALF) GAMMAF = -TWO*SQRT(PI)
      IF (X.EQ. HALF) GAMMAF =  SQRT(PI)
      IF (X.EQ. ONE ) GAMMAF =  ONE
      IF (X.EQ. TWO ) GAMMAF =  ONE
      IF (X.EQ. 1.5  ) GAMMAF =  SQRT(PI)/2.
      IF (X.EQ. 2.5) GAMMAF =  1.5*SQRT(PI)/2.
      IF (X.EQ. 3.5) GAMMAF =  0.5*(2.5*(1.5*SQRT(PI)))
      IF (X.EQ. 3. ) GAMMAF =  2.
      IF (X.EQ. 4. ) GAMMAF = 6.
      IF (X.EQ. 5. ) GAMMAF = 24.
      IF (X.EQ. 6. ) GAMMAF = 120.
      RETURN
      END FUNCTION

      REAL(KIND=RP)  FUNCTION PNORMJ (N,ALPHA,BETA)
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  ALPHA,BETA
      ONE   = 1.
      TWO   = 2.
      DN    = ((N))
      CONST = ALPHA+BETA+ONE
      IF (N.LE.1) THEN
         PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
         PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
         PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
         RETURN
      ENDIF
      PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
      PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
      PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
      PROD  = PROD*(ONE+BETA)*(TWO+BETA)
      DO 100 I=3,N
         DINDX = ((I))
         FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
         PROD  = PROD*FRAC
 100  CONTINUE
      PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
      RETURN
      END FUNCTION

      SUBROUTINE JACG (XJAC,NP,ALPHA,BETA)
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
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  XJAC(1)
      DATA KSTOP /10/
      DATA EPS/1.0e-12/
      N   = NP-1
      one = 1.
      DTH = 4.*ATAN(one)/(2.*((N))+2.)
      DO 40 J=1,NP
         IF (J.EQ.1) THEN
            X = COS((2.*(((J))-1.)+1.)*DTH)
         ELSE
            X1 = COS((2.*(((J))-1.)+1.)*DTH)
            X2 = XLAST
            X  = (X1+X2)/2.
         ENDIF
         DO 30 K=1,KSTOP
            CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X)
            RECSUM = 0.
            JM = J-1
            DO 29 I=1,JM
               RECSUM = RECSUM+1./(X-XJAC(NP-I+1))
 29         CONTINUE
            DELX = -P/(PD-RECSUM*P)
            X    = X+DELX
            IF (ABS(DELX) .LT. EPS) GOTO 31
 30      CONTINUE
 31      CONTINUE
         XJAC(NP-J+1) = X
         XLAST        = X
 40   CONTINUE
      DO 200 I=1,NP
         XMIN = 2.
         DO 100 J=I,NP
            IF (XJAC(J).LT.XMIN) THEN
               XMIN = XJAC(J)
               JMIN = J
            ENDIF
 100     CONTINUE
         IF (JMIN.NE.I) THEN
            SWAP = XJAC(I)
            XJAC(I) = XJAC(JMIN)
            XJAC(JMIN) = SWAP
         ENDIF
 200  CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE JACOBF (POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,N,ALP,BET,X)
!--------------------------------------------------------------------
!
!     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
!     of degree N at X.
!
!--------------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      APB  = ALP+BET
      POLY = 1.
      PDER = 0.
      IF (N .EQ. 0) RETURN
      POLYL = POLY
      PDERL = PDER
      POLY  = (ALP-BET+(APB+2.)*X)/2.
      PDER  = (APB+2.)/2.
      IF (N .EQ. 1) RETURN
      DO 20 K=2,N
         DK = ((K))
         A1 = 2.*DK*(DK+APB)*(2.*DK+APB-2.)
         A2 = (2.*DK+APB-1.)*(ALP**2-BET**2)
         B3 = (2.*DK+APB-2.)
         A3 = B3*(B3+1.)*(B3+2.)
         A4 = 2.*(DK+ALP-1.)*(DK+BET-1.)*(2.*DK+APB)
         POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
         PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
         PSAVE  = POLYL
         PDSAVE = PDERL
         POLYL  = POLY
         POLY   = POLYN
         PDERL  = PDER
         PDER   = PDERN
 20   CONTINUE
      POLYM1 = POLYL
      PDERM1 = PDERL
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      RETURN
      END SUBROUTINE

      REAL(KIND=RP) FUNCTION HGJ (II,Z,ZGJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGJ through
!     the NP Gauss Jacobi points ZGJ at the point Z.
!     Single precision version.
!
!---------------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZD = NMAX)
      REAL(KIND=RP)  ZD,ZGJD(NZD),ALPHAD,BETAD
      REAL(KIND=RP)  Z,ZGJ(1),ALPHA,BETA
      NPMAX = NZD
      IF (NP.GT.NPMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in HGJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NP=',NP
         call neko_error
      ENDIF
      ZD = Z
      DO 100 I=1,NP
         ZGJD(I) = ZGJ(I)
 100  CONTINUE
      ALPHAD = ALPHA
      BETAD  = BETA
      HGJ    = HGJD (II,ZD,ZGJD,NP,ALPHAD,BETAD)
      RETURN
      END FUNCTION

      REAL(KIND=RP)  FUNCTION HGJD (II,Z,ZGJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGJD through
!     the NZ Gauss-Lobatto Jacobi points ZGJ at the point Z.
!     Double precision version.
!
!---------------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  Z,ZGJ(1),ALPHA,BETA
      EPS = 1.e-5
      ONE = 1.
      ZI  = ZGJ(II)
      DZ  = Z-ZI
      IF (ABS(DZ).LT.EPS) THEN
         HGJD = ONE
         RETURN
      ENDIF
      CALL JACOBF (PZI,PDZI,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,ZI)
      CALL JACOBF (PZ,PDZ,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,Z)
      HGJD  = PZ/(PDZI*(Z-ZI))
      RETURN
      END FUNCTION

      REAL(KIND=RP) FUNCTION HGLJ (II,Z,ZGLJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGLJ through
!     the NZ Gauss-Lobatto Jacobi points ZGLJ at the point Z.
!     Single precision version.
!
!---------------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZD = NMAX)
      REAL(KIND=RP)  ZD,ZGLJD(NZD),ALPHAD,BETAD
      REAL(KIND=RP)  Z,ZGLJ(1),ALPHA,BETA
      NPMAX = NZD
      IF (NP.GT.NPMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in HGLJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NP=',NP
         call neko_error
      ENDIF
      ZD = Z
      DO 100 I=1,NP
         ZGLJD(I) = ZGLJ(I)
 100  CONTINUE
      ALPHAD = ALPHA
      BETAD  = BETA
      HGLJ   = HGLJD (II,ZD,ZGLJD,NP,ALPHAD,BETAD)
      RETURN
      END FUNCTION

      REAL(KIND=RP)  FUNCTION HGLJD (I,Z,ZGLJ,NP,ALPHA,BETA)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGLJD through
!     the NZ Gauss-Lobatto Jacobi points ZJACL at the point Z.
!     Double precision version.
!
!---------------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  Z,ZGLJ(1),ALPHA,BETA
      EPS = 1.e-5
      ONE = 1.
      ZI  = ZGLJ(I)
      DZ  = Z-ZI
      IF (ABS(DZ).LT.EPS) THEN
         HGLJD = ONE
         RETURN
      ENDIF
      N      = NP-1
      DN     = ((N))
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)
      CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,ZI)
      CONST  = EIGVAL*PI+ALPHA*(ONE+ZI)*PDI-BETA*(ONE-ZI)*PDI
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z)
      HGLJD  = (ONE-Z**2)*PD/(CONST*(Z-ZI))
      RETURN
      END FUNCTION

      SUBROUTINE DGJ (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!
!-----------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZDD = NMAX)
      REAL(KIND=RP)  DD(NZDD,NZDD),DTD(NZDD,NZDD),ZD(NZDD),ALPHAD,BETAD
      REAL(KIND=RP) D(NZD,NZD),DT(NZD,NZD),Z(1),ALPHA,BETA

      IF (NZ.LE.0) THEN
         WRITE (6,*) 'DGJ: Minimum number of Gauss points is 1'
         call neko_error
      ENDIF
      IF (NZ .GT. NMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in DGJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here Nz=',Nz
         call neko_error
      ENDIF
      IF ((ALPHA.LE.-1.).OR.(BETA.LE.-1.)) THEN
         WRITE (6,*) 'DGJ: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF
      ALPHAD = ALPHA
      BETAD  = BETA
      DO 100 I=1,NZ
         ZD(I) = Z(I)
 100  CONTINUE
      CALL DGJD (DD,DTD,ZD,NZ,NZDD,ALPHAD,BETAD)
      DO I=1,NZ
      DO J=1,NZ
         D(I,J)  = DD(I,J)
         DT(I,J) = DTD(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGJD (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!
!-----------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  D(NZD,NZD),DT(NZD,NZD),Z(1),ALPHA,BETA
      N    = NZ-1
      DN   = ((N))
      ONE  = 1.
      TWO  = 2.

      IF (NZ.LE.1) THEN
       WRITE (6,*) 'DGJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'DGJD: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      DO I=1,NZ
      DO J=1,NZ
         CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z(I))
         CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,NZ,ALPHA,BETA,Z(J))
         IF (I.NE.J) D(I,J) = PDI/(PDJ*(Z(I)-Z(J)))
         IF (I.EQ.J) D(I,J) = ((ALPHA+BETA+TWO)*Z(I)+ALPHA-BETA)/ &
              (TWO*(ONE-Z(I)**2))
         DT(J,I) = D(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGLJ (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Single precision version.
!
!-----------------------------------------------------------------
      PARAMETER (NMAX=84)
      PARAMETER (NZDD = NMAX)
      REAL(KIND=RP)  DD(NZDD,NZDD),DTD(NZDD,NZDD),ZD(NZDD),ALPHAD,BETAD
      REAL(KIND=RP) D(NZD,NZD),DT(NZD,NZD),Z(1),ALPHA,BETA

      IF (NZ.LE.1) THEN
       WRITE (6,*) 'DGLJ: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
      ENDIF
      IF (NZ .GT. NMAX) THEN
         WRITE (6,*) 'Too large polynomial degree in DGLJ'
         WRITE (6,*) 'Maximum polynomial degree is',NMAX
         WRITE (6,*) 'Here NZ=',NZ
         call neko_error
      ENDIF
      IF ((ALPHA.LE.-1.).OR.(BETA.LE.-1.)) THEN
         WRITE (6,*) 'DGLJ: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF
      ALPHAD = ALPHA
      BETAD  = BETA
      DO 100 I=1,NZ
         ZD(I) = Z(I)
 100  CONTINUE
      CALL DGLJD (DD,DTD,ZD,NZ,NZDD,ALPHAD,BETAD)
      DO I=1,NZ
      DO J=1,NZ
         D(I,J)  = DD(I,J)
         DT(I,J) = DTD(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGLJD (D,DT,Z,NZ,NZD,ALPHA,BETA)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Jacobi points Z.
!     Note: D and DT are square matrices.
!     Double precision version.
!
!-----------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  D(NZD,NZD),DT(NZD,NZD),Z(1),ALPHA,BETA
      N    = NZ-1
      DN   = ((N))
      ONE  = 1.
      TWO  = 2.
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

      IF (NZ.LE.1) THEN
       WRITE (6,*) 'DGLJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'DGLJD: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      DO I=1,NZ
      DO J=1,NZ
         CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(I))
         CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(J))
         CI = EIGVAL*PI-(BETA*(ONE-Z(I))-ALPHA*(ONE+Z(I)))*PDI
         CJ = EIGVAL*PJ-(BETA*(ONE-Z(J))-ALPHA*(ONE+Z(J)))*PDJ
         IF (I.NE.J) D(I,J) = CI/(CJ*(Z(I)-Z(J)))
         IF ((I.EQ.J).AND.(I.NE.1).AND.(I.NE.NZ)) &
              D(I,J) = (ALPHA*(ONE+Z(I))-BETA*(ONE-Z(I)))/ &
              (TWO*(ONE-Z(I)**2))
         IF ((I.EQ.J).AND.(I.EQ.1)) &
              D(I,J) =  (EIGVAL+ALPHA)/(TWO*(BETA+TWO)) 
         IF ((I.EQ.J).AND.(I.EQ.NZ)) &
              D(I,J) = -(EIGVAL+BETA)/(TWO*(ALPHA+TWO)) 
         DT(J,I) = D(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGLL (D,DT,Z,NZ,NZD)
!-----------------------------------------------------------------
!
!     Compute the derivative matrix D and its transpose DT
!     associated with the Nth order Lagrangian interpolants
!     through the NZ Gauss-Lobatto Legendre points Z.
!     Note: D and DT are square matrices.
!
!-----------------------------------------------------------------
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      PARAMETER (NMAX=84)
      REAL(KIND=RP) D(NZD,NZD),DT(NZD,NZD),Z(1)
      N  = NZ-1
      IF (NZ .GT. NMAX) THEN
         WRITE (6,*) 'Subroutine DGLL'
         WRITE (6,*) 'Maximum polynomial degree =',NMAX
         WRITE (6,*) 'Polynomial degree         =',NZ
      ENDIF
      IF (NZ .EQ. 1) THEN
         D(1,1) = 0.
         RETURN
      ENDIF
      FN = (N)
      d0 = FN*(FN+1.)/4.
      DO I=1,NZ
      DO J=1,NZ
         D(I,J) = 0.
         IF  (I.NE.J) D(I,J) = PNLEG(Z(I),N)/ &
                             (PNLEG(Z(J),N)*(Z(I)-Z(J)))
         IF ((I.EQ.J).AND.(I.EQ.1))  D(I,J) = -d0
         IF ((I.EQ.J).AND.(I.EQ.NZ)) D(I,J) =  d0
         DT(J,I) = D(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      REAL(KIND=RP) FUNCTION HGLL (I,Z,ZGLL,NZ)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant L through
!     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
!
!---------------------------------------------------------------------
      REAL(KIND=RP) ZGLL(1), EPS, DZ, Z
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      EPS = 1.E-5
      DZ = Z - ZGLL(I)
      IF (ABS(DZ) .LT. EPS) THEN
         HGLL = 1.
         RETURN
      ENDIF
      N = NZ - 1
      ALFAN = (N)*((N)+1.)
      HGLL = - (1.-Z*Z)*PNDLEG(Z,N)/ &
           (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))
      RETURN
      END FUNCTION

      REAL(KIND=RP) FUNCTION HGL (I,Z,ZGL,NZ)
!---------------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant HGL through
!     the NZ Gauss Legendre points ZGL at the point Z.
!
!---------------------------------------------------------------------
      REAL(KIND=RP) ZGL(1), Z, EPS, DZ
      EPS = 1.E-5
      DZ = Z - ZGL(I)
      IF (ABS(DZ) .LT. EPS) THEN
         HGL = 1.
         RETURN
      ENDIF
      N = NZ-1
      HGL = PNLEG(Z,NZ)/(PNDLEG(ZGL(I),NZ)*(Z-ZGL(I)))
      RETURN
      END FUNCTION

      REAL(KIND=RP) FUNCTION PNLEG (Z,N)
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
      REAL(KIND=RP) Z, P1, P2, P3
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      IF(ABS(Z) .LT. 1.0E-25) Z = 0.0


      P1   = 1.
      IF (N.EQ.0) THEN
         PNLEG = P1
         RETURN
      ENDIF
      P2   = Z
      P3   = P2
      DO 10 K = 1, N-1
         FK  = (K)
         P3  = ((2.*FK+1.)*Z*P2 - FK*P1)/(FK+1.)
         P1  = P2
         P2  = P3
 10   CONTINUE
      PNLEG = P3
      if (n.eq.0) pnleg = 1.
      RETURN
      END FUNCTION

      REAL(KIND=RP) FUNCTION PNDLEG (Z,N)
!----------------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     (Simpler than JACOBF)
!     Based on the recursion formula for the Legendre polynomials.
!
!----------------------------------------------------------------------
      REAL(KIND=RP) P1, P2, P1D, P2D, P3D, Z
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      P1   = 1.
      P2   = Z
      P1D  = 0.
      P2D  = 1.
      P3D  = 1.
      DO 10 K = 1, N-1
         FK  = (K)
         P3  = ((2.*FK+1.)*Z*P2 - FK*P1)/(FK+1.)
         P3D = ((2.*FK+1.)*P2 + (2.*FK+1.)*Z*P2D - FK*P1D)/(FK+1.)
         P1  = P2
         P2  = P3
         P1D = P2D
         P2D = P3D
 10   CONTINUE
      PNDLEG = P3D
      IF (N.eq.0) pndleg = 0.
      RETURN
      END FUNCTION

      SUBROUTINE DGLLGL (D,DT,ZM1,ZM2,IM12,NZM1,NZM2,ND1,ND2)
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
      REAL(KIND=RP) D(ND2,ND1), DT(ND1,ND2), ZM1(ND1), ZM2(ND2), IM12(ND2,ND1)
      REAL(KIND=RP) EPS, ZP, ZQ  
      IF (NZM1.EQ.1) THEN
        D (1,1) = 0.
        DT(1,1) = 0.
        RETURN
      ENDIF
      EPS = 1.E-6
      NM1 = NZM1-1
      DO IP = 1, NZM2
         DO JQ = 1, NZM1
            ZP = ZM2(IP)
            ZQ = ZM1(JQ)
            IF ((ABS(ZP) .LT. EPS).AND.(ABS(ZQ) .LT. EPS)) THEN
                D(IP,JQ) = 0.
            ELSE
                D(IP,JQ) = (PNLEG(ZP,NM1)/PNLEG(ZQ,NM1) &
                     -IM12(IP,JQ))/(ZP-ZQ)
            ENDIF
            DT(JQ,IP) = D(IP,JQ)
         END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGLJGJ (D,DT,ZGL,ZG,IGLG,NPGL,NPG,ND1,ND2,ALPHA,BETA)
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
      REAL(KIND=RP) D(ND2,ND1), DT(ND1,ND2), ZGL(ND1), ZG(ND2), IGLG(ND2,ND1)
      PARAMETER (NMAX=84)
      PARAMETER (NDD = NMAX)
      REAL(KIND=RP)  DD(NDD,NDD), DTD(NDD,NDD)
      REAL(KIND=RP)  ZGD(NDD), ZGLD(NDD), IGLGD(NDD,NDD)
      REAL(KIND=RP)  ALPHAD, BETAD

      IF (NPGL.LE.1) THEN
       WRITE(6,*) 'DGLJGJ: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
      ENDIF
      IF (NPGL.GT.NMAX) THEN
         WRITE(6,*) 'Polynomial degree too high in DGLJGJ'
         WRITE(6,*) 'Maximum polynomial degree is',NMAX
         WRITE(6,*) 'Here NPGL=',NPGL
         call neko_error
      ENDIF
      IF ((ALPHA.LE.-1.).OR.(BETA.LE.-1.)) THEN
         WRITE(6,*) 'DGLJGJ: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      ALPHAD = ALPHA
      BETAD  = BETA
      DO I=1,NPG
         ZGD(I) = ZG(I)
         DO J=1,NPGL
            IGLGD(I,J) = IGLG(I,J)
         END DO
      END DO
      DO 200 I=1,NPGL
         ZGLD(I) = ZGL(I)
 200  CONTINUE
      CALL DGLJGJD (DD,DTD,ZGLD,ZGD,IGLGD,NPGL,NPG,NDD,NDD,ALPHAD,BETAD)
      DO I=1,NPG
      DO J=1,NPGL
         D(I,J)  = DD(I,J)
         DT(J,I) = DTD(J,I)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE DGLJGJD (D,DT,ZGL,ZG,IGLG,NPGL,NPG,ND1,ND2,ALPHA,BETA)
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
      IMPLICIT REAL(KIND=RP)  (A-H,O-Z)
      REAL(KIND=RP)  D(ND2,ND1), DT(ND1,ND2), ZGL(ND1), ZG(ND2)
      REAL(KIND=RP)  IGLG(ND2,ND1), ALPHA, BETA

      IF (NPGL.LE.1) THEN
       WRITE(6,*) 'DGLJGJD: Minimum number of Gauss-Lobatto points is 2'
       call neko_error
      ENDIF
      IF ((ALPHA.LE.-1.).OR.(BETA.LE.-1.)) THEN
         WRITE(6,*) 'DGLJGJD: Alpha and Beta must be greater than -1'
         call neko_error
      ENDIF

      EPS    = 1.e-6
      ONE    = 1.
      TWO    = 2.
      NGL    = NPGL-1
      DN     = ((NGL))
      EIGVAL = -DN*(DN+ALPHA+BETA+ONE)

      DO I=1,NPG
      DO J=1,NPGL
         DZ = ABS(ZG(I)-ZGL(J))
         IF (DZ.LT.EPS) THEN
            D(I,J) = (ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I)))/ &
                 (TWO*(ONE-ZG(I)**2))
         ELSE
            CALL JACOBF (PI,PDI,PM1,PDM1,PM2,PDM2,NGL,ALPHA,BETA,ZG(I))
            CALL JACOBF (PJ,PDJ,PM1,PDM1,PM2,PDM2,NGL,ALPHA,BETA,ZGL(J))
            FACI   = ALPHA*(ONE+ZG(I))-BETA*(ONE-ZG(I))
            FACJ   = ALPHA*(ONE+ZGL(J))-BETA*(ONE-ZGL(J))
            CONST  = EIGVAL*PJ+FACJ*PDJ
            D(I,J) = ((EIGVAL*PI+FACI*PDI)*(ZG(I)-ZGL(J)) &
                 -(ONE-ZG(I)**2)*PDI)/(CONST*(ZG(I)-ZGL(J))**2)
         ENDIF
         DT(J,I) = D(I,J)
      END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE IGLM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss Legendre points.
!     Z2 : NZ2 points on mesh M.
!
!--------------------------------------------------------------------
      REAL(KIND=RP) I12(ND2,ND1),IT12(ND1,ND2),Z1(ND1),Z2(ND2), ZI
      IF (NZ1 .EQ. 1) THEN
         I12 (1,1) = 1.
         IT12(1,1) = 1.
         RETURN
      ENDIF
      DO I=1,NZ2
         ZI = Z2(I)
         DO J=1,NZ1
            I12 (I,J) = HGL(J,ZI,Z1,NZ1)
            IT12(J,I) = I12(I,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE IGLLM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2)
!----------------------------------------------------------------------
!
!     Compute the one-dimensional interpolation operator (matrix) I12
!     ands its transpose IT12 for interpolating a variable from a
!     Gauss-Lobatto Legendre mesh (1) to a another mesh M (2).
!     Z1 : NZ1 Gauss-Lobatto Legendre points.
!     Z2 : NZ2 points on mesh M.
!
!--------------------------------------------------------------------
      REAL(KIND=RP) I12(ND2,ND1),IT12(ND1,ND2),Z1(ND1),Z2(ND2),ZI
      IF (NZ1 .EQ. 1) THEN
         I12 (1,1) = 1.
         IT12(1,1) = 1.
         RETURN
      ENDIF
      DO I=1,NZ2
         ZI = Z2(I)
         DO J=1,NZ1
            I12 (I,J) = HGLL(J,ZI,Z1,NZ1)
            IT12(J,I) = I12(I,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE IGJM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2,ALPHA,BETA)
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
      REAL(KIND=RP) I12(ND2,ND1),IT12(ND1,ND2),Z1(ND1),Z2(ND2),ZI,ALPHA,BETA
      IF (NZ1 .EQ. 1) THEN
         I12 (1,1) = 1.
         IT12(1,1) = 1.
         RETURN
      ENDIF
      DO I=1,NZ2
         ZI = Z2(I)
         DO J=1,NZ1
            I12 (I,J) = HGJ(J,ZI,Z1,NZ1,ALPHA,BETA)
            IT12(J,I) = I12(I,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE

      SUBROUTINE IGLJM (I12,IT12,Z1,Z2,NZ1,NZ2,ND1,ND2,ALPHA,BETA)
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
      REAL(KIND=RP) I12(ND2,ND1),IT12(ND1,ND2),Z1(ND1),Z2(ND2),ZI,ALPHA,BETA
      IF (NZ1 .EQ. 1) THEN
         I12 (1,1) = 1.
         IT12(1,1) = 1.
         RETURN
      ENDIF
      DO I=1,NZ2
         ZI = Z2(I)
         DO J=1,NZ1
            I12 (I,J) = HGLJ(J,ZI,Z1,NZ1,ALPHA,BETA)
            IT12(J,I) = I12(I,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE
end module
