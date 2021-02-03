!> Adams-Bashforth coefs for Backward Differentiation schemes
module abbdf
  use num_types
  use math
  implicit none

  !> AB-BDF coefficeints
  type abbdf_t
     real(kind=dp), dimension(10) :: ab
     real(kind=dp), dimension(10) :: bd
     integer :: nab = 0
     integer :: nbd = 0
   contains
     procedure, pass(this) :: set_bd => abbdf_set_bd
     procedure, pass(this) :: set_abbd => abbdf_set_abbd
  end type abbdf_t


contains

  !>Compute backward-differentiation coefficients of order NBD
  subroutine abbdf_set_bd(this, dtbd)
    class(abbdf_t), intent(inout) :: this
    real(kind=dp), intent(inout), dimension(10) :: dtbd
    real(kind=dp), dimension(10,10) :: bdmat
    real(kind=dp), dimension(10) :: bdrhs
    real(kind=dp) :: bdf
    integer, parameter :: ldim = 10
    integer, dimension(10) :: ir, ic
    integer :: ibd, nsys, i

    associate(nbd => this%nbd, bd => this%bd)
    
      nbd = nbd + 1
      nbd = min(nbd, 3)

      if (nbd .eq. 1) then
         bd(1) = 1d0
         bdf = 1d0
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
      bd(1) = 1d0
      do ibd= 1, nbd + 1
         bd(ibd) = bd(ibd)/bdf
      end do
    end associate
  end subroutine abbdf_set_bd

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
  subroutine abbdf_set_abbd(this, dtlag)
    class(abbdf_t), intent(inout)  :: this
    real(kind=dp), intent(inout), dimension(10) :: dtlag
    real(kind=dp) :: dt0, dt1, dt2, dts, dta, dtb, dtc, dtd, dte

    associate(nab => this%nab, nbd => this%nbd, ab => this%ab)

      nab = nab + 1
      nab = min(nab, 3)
    
      dt0 = dtlag(1)
      dt1 = dtlag(2)
      dt2 = dtlag(3)
      
      if (nab .eq. 1) then
         ab(1) = 1.0
      else if (nab .eq. 2) then
         dta =  dt0/dt1
         if (nbd .eq. 1) then
            ab(2) = -0.5 * dta
            ab(1) =  1d0 - ab(2)
         else if (nbd .eq. 2) then
            ab(2) = -dta
            ab(1) =  1d0 - ab(2)
         endif
      else if (nab .eq. 3) then
         dts =  dt1 + dt2
         dta =  dt0 / dt1
         dtb =  dt1 / dt2
         dtc =  dt0 / dt2
         dtd =  dts / dt1
         dte =  dt0 / dts
         if (nbd .eq. 1) then
            ab(3) =  dte*( 0.5*dtb + dtc/3. )
            ab(2) = -0.5*dta - ab(3)*dtd
            ab(1) =  1.0 - ab(2) - ab(3)
         elseif (nbd .eq. 2) then
            ab(3) =  2./3.*dtc*(1./dtd + dte)
            ab(2) = -dta - ab(3)*dtd
            ab(1) =  1.0 - ab(2) - ab(3)
         elseif (nbd .eq. 3) then
            ab(3) =  dte * (dtb + dtc)
            ab(2) = -dta * (1.0 + dtb + dtc)
            ab(1) =  1.0 - ab(2) - ab(3)
         endif
      endif
    end associate
    
  end subroutine abbdf_set_abbd

  !
  ! CLEAN UP THIS MESS BELOW, USE LAPACK OR SIMILAR
  !


  
  subroutine bdsys (a,b,dt,nbd,ldim)
    real(kind=dp) ::  A(ldim,9),B(9),DT(9)
    integer :: ldim, j, n, k, i, nsys, nbd
    real(kind=dp) :: SUMDT
    CALL RZERO (A,ldim**2)
    N = NBD+1
    DO J=1,NBD
       A(1,J) = 1d0
    end DO
    A(1,NBD+1) = 0d0
    B(1) = 1.
    DO J=1,NBD
       SUMDT = 0d0
       DO  K=1,J
          SUMDT = SUMDT+DT(K)
       end DO
       A(2,J) = SUMDT
    end DO
    A(2,NBD+1) = -DT(1)
    B(2) = 0d0
    DO I=3,NBD+1
       DO J=1,NBD
          SUMDT = 0d0
          DO K=1,J
             SUMDT = SUMDT+DT(K)
          end DO
          A(I,J) = SUMDT**(I-1)
       end DO
       A(I,NBD+1) = 0d0
       B(I) = 0d0
    end DO
      
    end subroutine bdsys


    SUBROUTINE LU(A,N,ldim,IR,IC)
      integer :: n, ldim, IR(10), IC(10)
      real(kind=dp) :: A(ldim,10), xmax, ymax, B, Y, C
      integer :: i, j, k, l, m, icm, irl, k1
      DO I=1,N
         IR(I)=I
         IC(I)=I
      end DO
      K=1
      L=K
      M=K
      XMAX=ABS(A(K,K))
      DO 100 I=K,N
         DO 100 J=K,N
            Y=ABS(A(I,J))
            IF(XMAX.GE.Y) GOTO 100
            XMAX=Y
            L=I
            M=J
100   CONTINUE
      DO 1000 K=1,N
         IRL=IR(L)
         IR(L)=IR(K)
         IR(K)=IRL
         ICM=IC(M)
         IC(M)=IC(K)
         IC(K)=ICM
         IF(L.EQ.K) GOTO 300
         DO 200 J=1,N
            B=A(K,J)
            A(K,J)=A(L,J)
            A(L,J)=B
200      CONTINUE
300      IF(M.EQ.K) GOTO 500
         DO 400 I=1,N
            B=A(I,K)
            A(I,K)=A(I,M)
            A(I,M)=B
400      CONTINUE
500      C=1d0/A(K,K)
         A(K,K)=C
         IF(K.EQ.N) GOTO 1000
         K1=K+1
         XMAX=ABS(A(K1,K1))
         L=K1
         M=K1
         DO 600 I=K1,N
            A(I,K)=C*A(I,K)
600      CONTINUE
         DO 800 I=K1,N
            B=A(I,K)
            DO 800 J=K1,N
               A(I,J)=A(I,J)-B*A(K,J)
               Y=ABS(A(I,J))
               IF(XMAX.GE.Y) GOTO 800
               XMAX=Y
               L=I
               M=J
800     CONTINUE
1000  CONTINUE
    end subroutine lu
   
    SUBROUTINE SOLVE(F,A,K,N,ldim,IR,IC)
      real(kind=dp) ::  A(ldim,10),F(ldim,10), G(2000), B, Y
      integer :: IR(10),IC(10), N, N1, k, kk, i, j, ldim, ICM, URL, K1, ICI
      integer :: I1, IRI,IRL, IT
        

!      IF (N.GT.2000) THEN
!         write(6,*) 'Abort IN Subrtouine SOLVE: N>2000, N=',N
!         call exitt
!      ENDIF
      N1=N+1
      DO 1000 KK=1,K
         DO 100 I=1,N
            IRI=IR(I)
            G(I)=F(IRI,KK)
100      CONTINUE
         DO 400 I=2,N
            I1=I-1
            B=G(I)
            DO 300 J=1,I1
               B=B-A(I,J)*G(J)
300         CONTINUE
            G(I)=B
400      CONTINUE
         DO 700 IT=1,N
            I=N1-IT
            I1=I+1
            B=G(I)
            IF(I.EQ.N) GOTO 701
            DO 600 J=I1,N
               B=B-A(I,J)*G(J)
600         CONTINUE
701         G(I)=B*A(I,I)
700      CONTINUE
         DO 900 I=1,N
            ICI=IC(I)
            F(ICI,KK)=G(I)
900      CONTINUE
1000  CONTINUE
      RETURN
      END

  

  
end module abbdf
