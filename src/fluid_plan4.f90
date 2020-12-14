!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4
  use fluid_method
  use ax_helm
  implicit none

  type, extends(fluid_scheme_t) :: fluid_plan4_t
     type(field_t) :: u_e
     type(field_t) :: v_e
     type(field_t) :: w_e

     real(kind=dp), allocatable :: p_res(:)
     real(kind=dp), allocatable :: u_res(:,:,:,:)
     real(kind=dp), allocatable :: v_res(:,:,:,:)
     real(kind=dp), allocatable :: w_res(:,:,:,:)
     real(kind=dp), allocatable :: p_old(:,:,:,:)
     real(kind=dp), allocatable :: u_old(:,:,:,:)
     real(kind=dp), allocatable :: v_old(:,:,:,:)
     real(kind=dp), allocatable :: w_old(:,:,:,:)

     type(field_t) :: dp
     type(field_t) :: du
     type(field_t) :: dv
     type(field_t) :: dw

     type(field_t) :: wa1
     type(field_t) :: wa2
     type(field_t) :: wa3

     type(field_t) :: ta1
     type(field_t) :: ta2
     type(field_t) :: ta3
     type(field_t) :: ta4
     
     !> @todo move this to a scratch space
     type(field_t) :: work1
     type(field_t) :: work2

     type(ax_helm_t) :: Ax 

     real(kind=dp) :: tpres
     integer :: ncalls = 0
     integer :: niter = 1000


     !> Time variables
     real(kind=dp) :: ab(10), bd(10),dt_old(10)
!     real(kind=dp) :: dt = 1e-3
     real(kind=dp) :: t_old
     integer :: nab, nbd
   contains
     procedure, pass(this) :: init => fluid_plan4_init
     procedure, pass(this) :: free => fluid_plan4_free
     procedure, pass(this) :: step => fluid_plan4_step
  end type fluid_plan4_t

contains

  subroutine fluid_plan4_init(this, msh, lx, param, vel, prs)    
    class(fluid_plan4_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: param        
    character(len=80), intent(inout) :: vel
    character(len=80), intent(inout) :: prs

    call this%free()
    
    !> Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, solver_vel=vel, solver_prs=prs)
    
    !> Initialize variables specific to this plan
    allocate(this%p_res(this%dm_Xh%n_dofs))
    allocate(this%u_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%p_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%u_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
            
    call field_init(this%u_e, this%dm_Xh, 'u_e')
    call field_init(this%v_e, this%dm_Xh, 'v_e')
    call field_init(this%w_e, this%dm_Xh, 'w_e')
    
    call field_init(this%wa1, this%dm_Xh, 'wa1')
    call field_init(this%wa2, this%dm_Xh, 'wa2')
    call field_init(this%wa3, this%dm_Xh, 'wa3')

    call field_init(this%ta1, this%dm_Xh, 'ta1')
    call field_init(this%ta2, this%dm_Xh, 'ta2')
    call field_init(this%ta3, this%dm_Xh, 'ta3')
    call field_init(this%ta4, this%dm_Xh, 'ta3')
    
    call field_init(this%du, this%dm_Xh, 'du')
    call field_init(this%dv, this%dm_Xh, 'dv')
    call field_init(this%dw, this%dm_Xh, 'dw')
    call field_init(this%dp, this%dm_Xh, 'dp')

    call field_init(this%work1, this%dm_Xh, 'work1')
    call field_init(this%work2, this%dm_Xh, 'work2')

  end subroutine fluid_plan4_init

  subroutine fluid_plan4_free(this)
    class(fluid_plan4_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call field_free(this%u_e)
    call field_free(this%v_e)
    call field_free(this%w_e)
    
    call field_free(this%wa1)
    call field_free(this%wa2)
    call field_free(this%wa3)

    call field_free(this%ta1)
    call field_free(this%ta2)
    call field_free(this%ta3)
    call field_free(this%ta4)

    call field_free(this%du)
    call field_free(this%dv)
    call field_free(this%dw)
    call field_free(this%dp)

    call field_free(this%work1)
    call field_free(this%work2)

    if (allocated(this%p_res)) then
       deallocate(this%p_res)
    end if
    
    if (allocated(this%u_res)) then
       deallocate(this%u_res)
    end if
    
    if (allocated(this%v_res)) then
       deallocate(this%v_res)
    end if
    
    if (allocated(this%w_res)) then
       deallocate(this%w_res)
    end if
    
    if (allocated(this%u_old)) then
       deallocate(this%u_old)
    end if
       
    if (allocated(this%v_old)) then
       deallocate(this%v_old)
    end if
    
    if (allocated(this%w_old)) then
       deallocate(this%w_old)
    end if
    
    if (allocated(this%p_old)) then
       deallocate(this%p_old)
    end if

  end subroutine fluid_plan4_free
  
  subroutine fluid_plan4_step(this, t, tstep)
    class(fluid_plan4_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: t
    integer, intent(inout) :: tstep

    integer :: n, iter
    n = this%dm_Xh%n_dofs

    if (this%ncalls .eq. 0) then
       this%tpres=0.0
    end if

    
    this%ncalls = this%ncalls + 1
    call settime(t, this%params%dt, this%t_old, this%dt_old,&
                 this%ab, this%bd, this%nab, this%nbd, tstep)
    
    ! compute explicit contributions bfx,bfy,bfz 
    ! how should we handle this?A
    !It seems like mane of the operators are in navier1.f
    !Mybae time for a navier.f90? Or operators.f90
    call this%f_Xh%eval()
    call plan4_sumab(this%u_e%x,this%u%x,this%u_old,n,this%ab,this%nab)
    call plan4_sumab(this%v_e%x,this%v%x,this%v_old,n,this%ab,this%nab)

    if (this%dm_Xh%msh%gdim .eq. 3) then
       call plan4_sumab(this%w_e%x,this%w%x,this%w_old,n,this%ab,this%nab)
    end if
    
    !shopuld we have this or not?
    !if(iflomach) call opcolv(bfx,bfy,bfz,vtrans)

    ! add user defined divergence to qtl 
    ! Where do we save usrdiv? Is this necessary?
    ! Seems to be added for stabiity reasons in some cases.
    ! call add2 (qtl,usrdiv,ntot1)
    !lagvel, we keep several old velocity?
    !call lagvel
    call copy(this%p_old, this%p%x, n)
    call copy(this%u_old, this%u%x, n)
    call copy(this%v_old, this%v%x, n)
    call copy(this%w_old, this%w%x, n)

    ! mask Dirichlet boundaries (velocity)
    call this%bc_apply_vel()
    call this%bc_apply_prs()
    ! compute pressure
    call plan4_pres_setup(this%c_Xh%h1, this%params%rho, &
                          this%dm_Xh%n_dofs, this%c_Xh%ifh2)    
    call plan4_pres_residual(this%p,this%p_res, this%u_e, this%v_e, this%w_e, &
                             this%ta1, this%ta2, this%ta3, &
                             this%wa1, this%wa2, this%wa3, &
                             this%work1, this%work2, this%f_Xh, &
                             this%c_Xh, this%gs_Xh, this%bclst_prs, &
                             this%Ax, &
                             this%params%Re, this%params%rho)

    !Sets tolerances
    !call ctolspl  (tolspl,respr)
    !!OBSERVE we do not solve anything 
    !!bclist is input to the krylov solver, when bcs are inplace uncomment all the solve
    !statement!

    iter = this%ksp_prs%solve(this%Ax, this%dp, this%p_res, n, &
         this%c_Xh, this%bclst_prs, this%gs_Xh, this%niter)
    call add2(this%p%x,this%dp%x,n)
    !call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)
    !We only need to update h2 once I think then use the flag to switch on/off
    call plan4_vel_setup(this%c_Xh%h1, this%c_Xh%h2, &
         this%params%Re, this%params%rho, this%bd(1), &
         this%params%dt, this%dm_Xh%n_dofs, this%c_Xh%ifh2)
    
    call plan4_vel_residual(this%Ax, this%u, this%v, this%w, &
                            this%u_res, this%v_res, this%w_res, &
                            this%p, this%ta1%x, this%ta2%x, this%ta3%x, this%ta4%x, &
                            this%f_Xh, this%c_Xh, &
                            this%msh, this%Xh, this%dm_Xh%n_dofs)
    
    iter = this%ksp_vel%solve(this%Ax, this%du, this%u_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)
    iter = this%ksp_vel%solve(this%Ax, this%dv, this%v_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)
    iter = this%ksp_vel%solve(this%Ax, this%dw, this%w_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)

    call opadd2cm(this%u%x,this%v%x,this%w%x,this%du%x,this%dv%x,this%dw%x,1d0,n,this%msh%gdim)

  end subroutine fluid_plan4_step
  
  subroutine plan4_pres_setup(h1, rho, n, ifh2)
    integer, intent(inout) :: n
    real(kind=dp), intent(inout) :: h1(n)
    real(kind=dp), intent(inout) :: rho
    logical, intent(inout) :: ifh2
    call cmult(h1, 1d0 /rho, n)
    ifh2 = .false.
  end subroutine plan4_pres_setup

  subroutine plan4_vel_setup(h1, h2, Re, rho, bd, dt, n, ifh2)
    integer, intent(inout) :: n
    real(kind=dp), intent(inout) :: h1(n)
    real(kind=dp), intent(inout) :: h2(n)
    real(kind=dp), intent(inout) :: Re
    real(kind=dp), intent(inout) :: rho
    real(kind=dp), intent(inout) :: bd
    real(kind=dp), intent(inout) :: dt
    logical, intent(inout) :: ifh2
    real(kind=dp) :: dtbd    
    dtbd = rho * (bd / dt)
    h1 = (1d0 / Re)
    call cmult(h2, dtbd, n)
    ifh2 = .true.
  end subroutine plan4_vel_setup

  subroutine plan4_vel_residual(Ax, u, v, w, u_res, v_res, w_res, &
                                p, ta1, ta2, ta3, ta4, f_Xh, c_Xh, msh, Xh, n)
    type(ax_helm_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    real(kind=dp), intent(inout) :: u_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: v_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: w_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: ta1(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: ta2(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: ta3(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: ta4(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    integer, intent(inout) :: n

    real(kind=dp) :: scl

    call Ax%compute(u_res, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res, v%x, c_Xh, msh, Xh)
    if (msh%gdim .eq. 3) then
       call Ax%compute(w_res, w%x, c_Xh, msh, Xh)
    end if

    call opchsign(u_res, v_res, w_res, msh%gdim, n)

    scl = -1./3.

    !call col3    (this%ta4,this%vdiff,qtl,ntot)
    call rzero(ta4,c_xh%dof%n_dofs)
    call add2s1  (ta4,p%x,scl,c_Xh%dof%n_dofs)
    call opgrad  (ta1,ta2,ta3,ta4,c_Xh)
    !if(IFAXIS) then
    !  !   CALL COL2 (TA2, OMASK,NTOT)
    !  !   CALL COL2 (TA3, OMASK,NTOT)
    !  !endif
    !call opsub2  (u_res,v_res,w_res,ta1,ta2,ta3)
    call opadd2cm(u_res, v_res, w_res, ta1, ta2, ta3, -1d0, n, msh%gdim)
    call opadd2cm(u_res, v_res, w_res, f_Xh%u, f_Xh%v, f_Xh%w, 1d0, n, msh%gdim)


  end subroutine plan4_vel_residual

  subroutine plan4_pres_residual(p,p_res,u_e, v_e, w_e, ta1, ta2, ta3, &
       wa1, wa2, wa3, work1, work2, f_Xh, c_xh, gs_Xh, bclst_prs, Ax, Re, rho)
    type(field_t), intent(inout) :: p
    real(kind=dp), intent(inout) :: p_res(p%dof%n_dofs)
    type(field_t), intent(inout) :: u_e
    type(field_t), intent(inout) :: v_e
    type(field_t), intent(inout) :: w_e
    type(field_t), intent(inout) :: ta1
    type(field_t), intent(inout) :: ta2
    type(field_t), intent(inout) :: ta3
    type(field_t), intent(inout) :: wa1
    type(field_t), intent(inout) :: wa2
    type(field_t), intent(inout) :: wa3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(bc_list_t), intent(inout) :: bclst_prs
    class(Ax_t), intent(inout) :: Ax
    real(kind=dp), intent(inout) :: Re
    real(kind=dp), intent(inout) :: rho
    real(kind=dp) :: scl
    integer :: i
    integer :: n, gdim

    n = c_Xh%dof%size()
    gdim = c_Xh%msh%gdim
    
     call plan4_op_curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
     call plan4_op_curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)
     call opcolv(wa1%x, wa2%x, wa3%x, c_Xh%B, gdim, n)
     scl = -4d0 / 3d0
     call opadd2cm (wa1%x, wa2%x, wa3%x, ta1%x, ta2%x, ta3%x, scl, n, gdim)

     work1%x = (1d0 / Re) / rho
     call opcolv(wa1%x, wa2%x, wa3%x, work1%x, gdim, n)

     !BOUNDARY CONDITION, DIRICHLET PRESSURE!
     !call bcdirpr (pr)

     
     call Ax%compute(p_res,p%x,c_Xh,p%msh,p%Xh)
     call chsign  (p_res,n)

     do i=1,n
        ta1%x(i,1,1,1) = f_Xh%u(i,1,1,1) / rho - wa1%x(i,1,1,1)
        ta2%x(i,1,1,1) = f_Xh%v(i,1,1,1) / rho - wa2%x(i,1,1,1)
        ta3%x(i,1,1,1) = f_Xh%w(i,1,1,1) / rho - wa3%x(i,1,1,1)
     enddo
     
     print *, sum(ta1%x)
     !Need to consider cyclic bcs here...
     call gs_op(gs_Xh, ta1, GS_OP_ADD) 
     call gs_op(gs_Xh, ta2, GS_OP_ADD) 
     call gs_op(gs_Xh, ta3, GS_OP_ADD) 

     do i=1,n
        ta1%x(i,1,1,1) = ta1%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
        ta2%x(i,1,1,1) = ta2%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
        ta3%x(i,1,1,1) = ta3%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
     enddo

     if (gdim .eq. 3) then
         call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
         call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
         call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
         do i=1,n
            p_res(i) = p_res(i)+wa1%x(i,1,1,1)+wa2%x(i,1,1,1)+wa3%x(i,1,1,1)
         enddo
      else
         call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
         call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)

         do i=1,n
            p_res(i) = p_res(i)+wa1%x(i,1,1,1)+wa2%x(i,1,1,1)
         enddo
      endif

    ! APPLY boundary conditions! 
     ! call bclist_apply(this%p_res,this%bclist)

!      DO 100 IEL=1,NELV
!         DO 300 IFC=1,NFACES
!            CALL RZERO  (W1(1,IEL),NXYZ1)
!            CALL RZERO  (W2(1,IEL),NXYZ1)
!            IF (ldim.EQ.3)
!     $      CALL RZERO  (W3(1,IEL),NXYZ1)
!            CB = CBC(IFC,IEL,IFIELD)
!            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v'.or.
!     $         cb.eq.'MV '.or.cb.eq.'mv '.or.cb.eq.'shl') then
!               CALL FACCL3
!     $         (W1(1,IEL),VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
!               CALL FACCL3
!     $         (W2(1,IEL),VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
!               IF (ldim.EQ.3)
!     $          CALL FACCL3
!     $         (W3(1,IEL),VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
!            ELSE IF (CB(1:3).EQ.'SYM') THEN
!               CALL FACCL3
!     $         (W1(1,IEL),TA1(1,IEL),UNX(1,1,IFC,IEL),IFC)
!               CALL FACCL3
!     $         (W2(1,IEL),TA2(1,IEL),UNY(1,1,IFC,IEL),IFC)
!               IF (ldim.EQ.3)
!     $          CALL FACCL3
!     $         (W3(1,IEL),TA3(1,IEL),UNZ(1,1,IFC,IEL),IFC)
!            ENDIF
!            CALL ADD2   (W1(1,IEL),W2(1,IEL),NXYZ1)
!            IF (ldim.EQ.3)
!     $      CALL ADD2   (W1(1,IEL),W3(1,IEL),NXYZ1)
!            CALL FACCL2 (W1(1,IEL),AREA(1,1,IFC,IEL),IFC)
!            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v'.or.
!     $         cb.eq.'MV '.or.cb.eq.'mv ') then
!              CALL CMULT(W1(1,IEL),dtbd,NXYZ1)
!            endif
!            CALL SUB2 (RESPR(1,IEL),W1(1,IEL),NXYZ1)
!  300    CONTINUE
!  100 CONTINUE

!     Assure that the residual is orthogonal to (1,1,...,1)T 
!     (only if all Dirichlet b.c.)
!     REALLY NOT sure when we shoudl do this. Results for poisson was not ecourgaging
!      CALL ORTHO (RESPR)


  end subroutine plan4_pres_residual

  subroutine plan4_op_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(inout)  :: c_Xh 
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    endif
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
           call sub3(w2%x, work1%x, work2%x, n)
        else
           call rzero (work1%x, n)
           call dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
           call sub3(w2%x, work1%x, work2%x, n)
        endif
        !     this%work1=dv/dx ; this%work2=du/dy
        call dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
        call dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
        call sub3(w3%x, work1%x, work2%x, n)
        !!    BC dependent, how should we handle this? they do Avg at bndry
!     if (ifavg .and. .not. ifcyclic) then
!
!         ifielt = ifield
!         ifield = 1
!       
!         call opcolv  (w1,w2,w3,bm1)
!         call opdssum (w1,w2,w3)
!         call opcolv  (w1,w2,w3,binvm1)
!
!         ifield = ifielt
!
!      endif
  end subroutine plan4_op_curl
  
  !> Sum up AB/BDF contributions 
  subroutine plan4_sumab(v_e,v,v_old,n,ab,nab)
    integer, intent(in) :: n, nab
    real(kind=dp), dimension(n), intent(inout) :: v_e,v, v_old
    real(kind=dp), dimension(3), intent(in) :: ab
    real(kind=dp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    call add3s2(v_e,v,v_old,ab0,ab1,n)
    !should we have this right now we only save one old velocity
    !if(nab .eq. 3) call add2s2(v_e,v_old(1,2),ab2,n)
  end subroutine plan4_sumab
!>     Store old time steps and compute new time step, time and timef.
!!     Set time-dependent coefficients in time-stepping schemes.
!!     @note this really should be placed somewhere else. Right now dt, etc is hardcoded. 
  subroutine settime(t, dt, t_old, dt_old, ab, bd, nab, nbd, tstep)
    real(kind=dp), intent(inout) :: t
    real(kind=dp), intent(inout) :: dt
    real(kind=dp), intent(inout) :: t_old
    real(kind=dp), dimension(10), intent(inout) :: dt_old
    real(kind=dp), dimension(10), intent(inout) :: ab
    real(kind=dp), dimension(10), intent(inout) :: bd
    integer, intent(inout) :: nab
    integer, intent(inout) :: nbd
    integer, intent(in) :: tstep
    integer :: i


!    Set time step.
      do i=10,2,-1
         dt_old(i) = dt_old(i-1)
      end do
      ! SHOULD be implemented!
      !call setdt
      dt_old(1) = dt
      if (tstep .eq. 1) then
         dt_old(2) = dt
      end if
     
      !IF (ISTEP.EQ.1 .and. irst.le.0) DTLAG(2) = DT

      !Set time.

      t_old = t
      t = t + dt

      !Set coefficients in AB/BD-schemes.

      !call SETORDBD
      !hardcoded for now
      nbd = 1
      call rzero (bd, 10)
      call setbd(bd, dt_old, nbd)
      ! This can also be varied, hardcoded for now
      nab = min(tstep, 3)
      !dont know what irst is really
      !IF (ISTEP.lt.NAB.and.irst.le.0) NAB = ISTEP
      call rzero(ab, 10)
      call setabbd(ab, dt_old, nab,nbd)

  end subroutine settime
!>
!!     Compute Adams-Bashforth coefficients (order NAB, less or equal to 3)
!!     
!!     NBD .EQ. 1
!!     Standard Adams-Bashforth coefficients 
!!
!!     NBD .GT. 1
!!     Modified Adams-Bashforth coefficients to be used in con-
!!     junction with Backward Differentiation schemes (order NBD)
!!
  subroutine setabbd (ab,dtlag,nab,nbd)
    INTEGER, intent(in) :: nab, nbd
    REAL(KIND=DP), intent(inout), dimension(NAB) :: AB(NAB),DTLAG(NAB)
    REAL(KIND=DP) :: DT0, DT1, DT2, DTS, DTA, DTB, DTC, DTD, DTE
    DT0 = DTLAG(1)
    DT1 = DTLAG(2)
    DT2 = DTLAG(3)
    IF ( NAB.EQ.1 ) THEN
        AB(1) = 1.0
    ELSEIF ( NAB.EQ.2 ) THEN
        DTA =  DT0/DT1
        IF ( NBD.EQ.1 ) THEN
        AB(2) = -0.5*DTA
        AB(1) =  1.0 - AB(2)
        ELSEIF ( NBD.EQ.2 ) THEN
        AB(2) = -DTA
        AB(1) =  1.0 - AB(2)
        ENDIF
    ELSEIF ( NAB.EQ.3 ) THEN
        DTS =  DT1 + DT2
        DTA =  DT0 / DT1
        DTB =  DT1 / DT2
        DTC =  DT0 / DT2
        DTD =  DTS / DT1
        DTE =  DT0 / DTS
        IF ( NBD.EQ.1 ) THEN
        AB(3) =  DTE*( 0.5*DTB + DTC/3. )
        AB(2) = -0.5*DTA - AB(3)*DTD
        AB(1) =  1.0 - AB(2) - AB(3)
        ELSEIF ( NBD.EQ.2 ) THEN
        AB(3) =  2./3.*DTC*(1./DTD + DTE)
        AB(2) = -DTA - AB(3)*DTD
        AB(1) =  1.0 - AB(2) - AB(3)
        ELSEIF ( NBD.EQ.3 ) THEN
        AB(3) =  DTE*(DTB + DTC)
        AB(2) = -DTA*(1.0 + DTB + DTC)
        AB(1) =  1.0 - AB(2) - AB(3)
        ENDIF
    ENDIF
  end subroutine setabbd

!>Compute bacward-differentiation coefficients of order NBD
! @note this need to be fixed as well...
  subroutine setbd (bd,dtbd,nbd)
    REAL(kind=dp), intent(inout), dimension(10) :: bd, dtbd
    INTEGER, intent(in) :: nbd
    REAL(kind=dp) :: BDMAT(10,10),BDRHS(10), BDF
    INTEGER :: IR(10),IC(10)
    INTEGER :: IBD, ldim = 10
    !  IF (NBD.EQ.1) THEN
         BD(1) = 1d0
         BDF   = 1d0
    !  ELSEIF (NBD.GE.2) THEN
    !     NSYS = NBD+1
    !     CALL BDSYS (BDMAT,BDRHS,DTBD,NBD,ldim)
    !     CALL LU    (BDMAT,NSYS,ldim,IR,IC)
    !     CALL SOLVE (BDRHS,BDMAT,1,NSYS,ldim,IR,IC)
    !     DO 30 I=1,NBD
    !        BD(I) = BDRHS(I)
    !     end do
    !     BDF = BDRHS(NBD+1)
    !  ENDIF
    !Normalize
      DO IBD=NBD,1,-1
         BD(IBD+1) = BD(IBD)
      end do
      BD(1) = 1d0
      DO IBD=1,NBD+1
         BD(IBD) = BD(IBD)/BDF
      end do
  end subroutine setbd
end module fluid_plan4
