!> Defines a domain as a subset of facets in a mesh
module curve
  use stack
  use utils
  use structs
  use point
  use hex
  implicit none

  ! Maybe should be moved somewhere else
  type :: curve_t
     type(struct_curve_t), allocatable :: curve_el(:)
     type(stack_curve_t), private :: scratch
     integer :: size = 0
     logical, private :: finalized = .false.
   contains
     procedure, pass(z) :: init => curve_element_init
     procedure, pass(z) :: free => curve_element_free
     procedure, pass(z) :: finalize => curve_element_finalize
     procedure, pass(z) :: add_element => curve_element_add
!     procedure, pass(z) :: apply_xyz => curve_apply_xyz
  end type curve_t
  
contains

  !> Initialize a curved domain
  subroutine curve_element_init(z, size)
    class(curve_t), intent(inout) :: z
    integer, optional :: size

    call curve_element_free(z)

    if (present(size)) then
       call z%scratch%init(size)
    else
       call z%scratch%init()
    end if
    
  end subroutine curve_element_init

  !> Deallocate a domain
  subroutine curve_element_free(z)
    class(curve_t), intent(inout) :: z
    if (allocated(z%curve_el)) then
       deallocate(z%curve_el)
    end if

    z%finalized = .false.
    z%size = 0

    call z%scratch%free()
    
  end subroutine curve_element_free

  !> Finalize a domain list
  !! @details Create a static list of (facet,el) tuples
  subroutine curve_element_finalize(z)
    class(curve_t), intent(inout) :: z
    type(struct_curve_t), pointer :: tp(:)
    integer :: i
    
    if (.not. z%finalized) then

       allocate(z%curve_el(z%scratch%size()))
       
       tp => z%scratch%array()
       do i = 1, z%scratch%size()
          z%curve_el(i) = tp(i)
       end do

       z%size = z%scratch%size()

       call z%scratch%clear()

       z%finalized = .true.
       
    end if
    
  end subroutine curve_element_finalize

  !> Add a (facet, el) tuple to an unfinalized domain
  subroutine curve_element_add(z, el_idx, curve_data, curve_type )
    class(curve_t), intent(inout) :: z
    real(kind=dp), dimension(5,8), intent(inout) :: curve_data
    integer, dimension(8), intent(inout) :: curve_type
    integer, intent(in) :: el_idx
    type(struct_curve_t) :: c_el

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if
    c_el%curve_data = curve_data
    c_el%curve_type = curve_type
    c_el%el_idx = el_idx
    call z%scratch%push(c_el)
  end subroutine curve_element_add

!  !> Add a (facet, el) tuple to an unfinalized domain
!  subroutine curve_apply_xyz(this, x, y, z, Xh, elements, nel)
!    class(curve_t), intent(inout) :: this
!    integer, intent(in) :: lx, nel
!    class(element_t)
!    integer :: el_idx
!    real(kind=dp), dimension(Xh%lxyz, nel) :: x, y, z
!    do i = 1, this%size
!       el_idx = this%curve_el(i)%el_idx
!       do j = 1, 12
!          if (this%curve_el(i)%curve_type(j) .eq. 1) then
!             !call sphere_surface(j, this%curve_el(i)%curve_data(:,j),x(1,el_idx), y(1,el_idx), z(1, el_idx)) 
!          else if (this%curve_el(i)%curve_type(j) .eq. 2) then
!             !call generate_surface(j, this%curve_el(i)%curve_data(:,j), x(1,1,1,el_idx), y(1,1,1,el_idx), z(1,1,1, el_idx)) 
!          end if
!       end do
!       do j = 1, 12
!          if (this%curve_el(i)%curve_type(j) .eq. 3) then
!             call arc_surface(j, this%curve_el(i)%curve_data(:,j),x(1,1,1,el_idx), y(1,1,1,el_idx), z(1,1,1, el_idx)) 
!          end if
!       end do
!
!    end do
!  end subroutine curve_apply_xyz


!  subroutine sphere_surface(face, curve_data, x, y, z, Xh) 
!!  Major overhaul.
!!  Martin Karp 1 March 2021 15:30:15
!!
!!  5 Aug 1988 19:29:52 
!!
!!  Program to generate spherical shell elements for NEKTON
!!  input.  Paul F. Fischer
!  
!    NXM1 = NX-1
!    NYM1 = NY-1
!    NXY  = NX*NZ
!    NXY3 = 3*NX*NZ
!    XCTR   = curve_data(1,face)
!    YCTR   = curve_data(2,face)
!    ZCTR   = curve_data(3,face)
!    RADIUS = curve_data(4,face)
!    IFACE  = EFACE1(IFCE)
!
!    Generate (normalized) corner vectors XCV(1,i,j):
!
!    CALL CRN3D(XCV,XC(1,IE),YC(1,IE),ZC(1,IE),CURVE(1,IFCE,IE),IFACE)
!
!    Generate edge vectors on the sphere RR=1.0, 
!    for (r,s) = (-1,*),(1,*),(*,-1),(*,1)      
!
!    CALL EDG3D(XYSRF,XCV(1,1,1),XCV(1,1,2), 1, 1, 1,NY,NX,NY)
!    CALL EDG3D(XYSRF,XCV(1,2,1),XCV(1,2,2),NX,NX, 1,NY,NX,NY)
!    CALL EDG3D(XYSRF,XCV(1,1,1),XCV(1,2,1), 1,NX, 1, 1,NX,NY)
!    CALL EDG3D(XYSRF,XCV(1,1,2),XCV(1,2,2), 1,NX,NY,NY,NX,NY)
!
!    ivtx = cface(1,ifce)
!    ivto = cface(2,ifce)
!    vout(1) = xc(ivtx,ie)-xc(ivto,ie)
!    vout(2) = yc(ivtx,ie)-yc(ivto,ie)
!    vout(3) = zc(ivtx,ie)-zc(ivto,ie)
!
!    vsph(1) = xc(ivtx,ie)-xctr
!    vsph(2) = yc(ivtx,ie)-yctr
!    vsph(3) = zc(ivtx,ie)-zctr
!    ifconcv = .true.
!    sign    = DOT(vsph,vout,3)
!    if (sign.gt.0) ifconcv = .false.
!    write(6,*) 'THIS IS SIGN:',sign
!
!    DO 200 J=2,NYM1
!       CALL CROSS(VN1,XYSRF(1,1,J),XYSRF(1,NX,J))
!       DO 200 I=2,NXM1
!          CALL CROSS(VN2,XYSRF(1,I,1),XYSRF(1,I,NY))
!          if (ifconcv) then
!          IF (IFACE.EQ.1.OR.IFACE.EQ.4.OR.IFACE.EQ.5) THEN
!             CALL CROSS(XYSRF(1,I,J),VN2,VN1)
!          ELSE
!             CALL CROSS(XYSRF(1,I,J),VN1,VN2)
!          ENDIF
!20  CONTINUE
!
!    Normalize all vectors to the unit sphere.
!
!    DO 300 I=1,NXY
!       CALL NORM3D(XYSRF(1,I,1))
!30  CONTINUE
!
!    Scale by actual radius
!
!    CALL CMULT(XYSRF,RADIUS,NXY3)
!
!    Add back the sphere center offset
!
!    DO 400 I=1,NXY
!       XYSRF(1,I,1)=XYSRF(1,I,1)+XCTR
!       XYSRF(2,I,1)=XYSRF(2,I,1)+YCTR
!       XYSRF(3,I,1)=XYSRF(3,I,1)+ZCTR
!40  CONTINUE
!
!
!    Transpose data, if necessary
!
!    IF (IFACE.EQ.1.OR.IFACE.EQ.4.OR.IFACE.EQ.5) THEN
!       DO 500 J=1  ,NY
!       DO 500 I=J+1,NX
!          TMP=XYSRF(1,I,J)
!          XYSRF(1,I,J)=XYSRF(1,J,I)
!          XYSRF(1,J,I)=TMP
!          TMP=XYSRF(2,I,J)
!          XYSRF(2,I,J)=XYSRF(2,J,I)
!          XYSRF(2,J,I)=TMP
!          TMP=XYSRF(3,I,J)
!          XYSRF(3,I,J)=XYSRF(3,J,I)
!          XYSRF(3,J,I)=TMP
!50     CONTINUE
!    ENDIF
!
!    Compute surface deflection and perturbation due to face IFACE
!
!    CALL DSSET(NX,NY,NZ)
!    JS1    = SKPDAT(1,IFACE)
!    JF1    = SKPDAT(2,IFACE)
!    JSKIP1 = SKPDAT(3,IFACE)
!    JS2    = SKPDAT(4,IFACE)
!    JF2    = SKPDAT(5,IFACE)
!    JSKIP2 = SKPDAT(6,IFACE)
!
!    IOPP(1) = NX-1
!    IOPP(2) = NX*(NY-1)
!    IOPP(3) = NX*NY*(NZ-1)
!    NXX(1)  = NX
!    NXX(2)  = NY
!    NXX(3)  = NZ
!    IDIR    = 2*MOD(IFACE,2) - 1
!    IFC2    = (IFACE+1)/2
!    DELT    = 0.0
!    I=0
!    DO 700 J2=JS2,JF2,JSKIP2
!    DO 700 J1=JS1,JF1,JSKIP1
!       I=I+1
!       JOPP = J1 + IOPP(IFC2)*IDIR
!       X2(1) = XML(J1,J2,1,IE)
!       X2(2) = YML(J1,J2,1,IE)
!       X2(3) = ZML(J1,J2,1,IE)
!
!       DX(1) = XYSRF(1,I,1)-X2(1)
!       DX(2) = XYSRF(2,I,1)-X2(2)
!       DX(3) = XYSRF(3,I,1)-X2(3)
!
!       NXS = NXX(IFC2)
!       JOFF = (J1-JOPP)/(NXS-1)
!       DO 600 IX = 2,NXS
!          J = JOPP + JOFF*(IX-1)
!          ZETA = 0.5*(ZGML(IX,IFC2) + 1.0)
!          XML(J,J2,1,IE) = XML(J,J2,1,IE)+DX(1)*ZETA
!          YML(J,J2,1,IE) = YML(J,J2,1,IE)+DX(2)*ZETA
!          ZML(J,J2,1,IE) = ZML(J,J2,1,IE)+DX(3)*ZETA
!60     CONTINUE
!70  CONTINUE
!
!    return
!  end


end module curve
