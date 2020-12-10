!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4
  use fluid_method
  use ax_helm
  implicit none

  type, extends(fluid_scheme_t) :: fluid_plan4_t
     !>@todo Add plan4 related data, ax, precon ect
     real(kind=dp), allocatable :: u_e(:,:,:,:)
     real(kind=dp), allocatable :: v_e(:,:,:,:)
     real(kind=dp), allocatable :: w_e(:,:,:,:)
     real(kind=dp), allocatable :: p_res(:,:)
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

     real(kind=dp), allocatable :: bfx(:,:,:,:)
     real(kind=dp), allocatable :: bfy(:,:,:,:)
     real(kind=dp), allocatable :: bfz(:,:,:,:)

     real(kind=dp), allocatable :: w1(:)
     real(kind=dp), allocatable :: wa1(:)
     real(kind=dp), allocatable :: wa2(:)
     real(kind=dp), allocatable :: wa3(:)
     real(kind=dp), allocatable :: ta1(:,:)
     real(kind=dp), allocatable :: ta2(:,:)
     real(kind=dp), allocatable :: ta3(:,:)
     real(kind=dp), allocatable :: ta4(:,:)


     !> @todo move this to a scratch space
     real(kind=dp), allocatable :: work1(:,:,:,:)
     real(kind=dp), allocatable :: work2(:,:,:,:)
     
     real(kind=dp), allocatable :: vtrans(:,:,:,:) !< Inverted Mass matrix/volume matrix
     real(kind=dp), allocatable :: vdiff(:,:,:,:) !< Inverted Mass matrix/volume matrix
 
     type(ax_helm_t) :: ax 

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
    
    !> Initialize variables sepcific to this plan
    allocate(this%u_e(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_e(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_e(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%p_res(this%Xh%lxyz,this%msh%nelv))
    allocate(this%u_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%p_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%u_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_old(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%bfx(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%bfy(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%bfz(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%w1 (this%dm_Xh%n_dofs))
    allocate(this%wa1(this%dm_Xh%n_dofs))
    allocate(this%wa2(this%dm_Xh%n_dofs))
    allocate(this%wa3(this%dm_Xh%n_dofs))
    allocate(this%ta1(this%Xh%lxyz,this%msh%nelv))
    allocate(this%ta2(this%Xh%lxyz,this%msh%nelv))
    allocate(this%ta3(this%Xh%lxyz,this%msh%nelv))
    allocate(this%ta4(this%Xh%lxyz,this%msh%nelv))

    allocate(this%work1(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%work2(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%vtrans(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%vdiff(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))

    call field_init(this%du, this%dm_Xh, 'du')
    call field_init(this%dv, this%dm_Xh, 'dv')
    call field_init(this%dw, this%dm_Xh, 'dw')
    call field_init(this%dp, this%dm_Xh, 'dp')

  end subroutine fluid_plan4_init

  subroutine fluid_plan4_free(this)
    class(fluid_plan4_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()
    
    call field_free(this%du)
    call field_free(this%dv)
    call field_free(this%dw)
    call field_free(this%dp)

    if (allocated(this%u_e)) then
       deallocate(this%u_e)
    end if
    
    if (allocated(this%v_e)) then
       deallocate(this%v_e)
    end if
    
    if (allocated(this%w_e)) then
       deallocate(this%w_e)
    end if
    
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
    
    if (allocated(this%p_res)) then
       deallocate(this%p_res)
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
    
    if (allocated(this%bfx)) then
       deallocate(this%bfx)
    end if
    
    if (allocated(this%bfy)) then
       deallocate(this%bfy)
    end if

    if (allocated(this%bfz)) then
       deallocate(this%bfz)
    end if
    
    if (allocated(this%w1)) then
       deallocate(this%w1)
    end if
    
    if (allocated(this%wa1)) then
       deallocate(this%wa1)
    end if
    
    if (allocated(this%wa2)) then
       deallocate(this%wa2)
    end if
    
    if (allocated(this%wa3)) then
       deallocate(this%wa3)
    end if
    
    if (allocated(this%ta1)) then
       deallocate(this%ta1)
    end if
    
    if (allocated(this%ta2)) then
       deallocate(this%ta2)
    end if
    
    if (allocated(this%ta3)) then
       deallocate(this%ta3)
    end if
    
    if (allocated(this%ta4)) then
       deallocate(this%ta4)
    end if
    
    if (allocated(this%work1)) then
       deallocate(this%work1)
    end if
    
    if (allocated(this%work2)) then
       deallocate(this%work2)
    end if
    
    if (allocated(this%vtrans)) then
       deallocate(this%vtrans)
    end if
    
    if (allocated(this%vdiff)) then
       deallocate(this%vdiff)
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
    call plan4_sumab(this%u_e,this%u%x,this%u_old,n,this%ab,this%nab)
    call plan4_sumab(this%v_e,this%v%x,this%v_old,n,this%ab,this%nab)

    if (this%dm_Xh%msh%gdim .eq. 3) then
       call plan4_sumab(this%w_e,this%w%x,this%w_old,n,this%ab,this%nab)
    end if
    
    !shopuld we have this or not?
    !if(iflomach) call opcolv(bfx,bfy,bfz,vtrans)

    ! add user defined divergence to qtl 
    ! Where do we save usrdiv? Is this necessary?
    ! Seems to be added for stabiity reasons in some cases.
    ! call add2 (qtl,usrdiv,ntot1)
    !lagvel, we keep several old velocity?
    !call lagvel
    call copy(this%p_old,this%p%x,n)
    call copy(this%u_old,this%u%x,n)
    call copy(this%v_old,this%v%x,n)
    call copy(this%w_old,this%w%x,n)

    ! mask Dirichlet boundaries
    ! Add our new boundry apply
    !call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask) 

    ! compute pressure
    call plan4_pres_setup(this)
    call plan4_pres_residual(this)
    !Sets tolerances
    !call ctolspl  (tolspl,respr)
    !!OBSERVE we do not solve anything 
    !!bclist is input to the krylov solver, when bcs are inplace uncomment all the solve
    !statement!
    !iter = this%ksp_prs%solve(this%Ax,this%dp, this%p_res, n, this%c_Xh, this%bclst, this%gs_Xh, this%niter)
    call add2(this%p%x,this%dp%x,n)
    call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)
    !We only need to update h2 once I think then use the flag to switch on/off
    call plan4_vel_setup(this) 
    call plan4_vel_residual(this)
    
    iter = this%ksp_vel%solve(this%Ax,this%du, this%u_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)
    iter = this%ksp_vel%solve(this%Ax,this%dv, this%v_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)
    iter = this%ksp_vel%solve(this%Ax,this%dw, this%w_res, n, &
         this%c_Xh, this%bclst_vel, this%gs_Xh, this%niter)

    call opadd2cm(this%u%x,this%v%x,this%w%x,this%du%x,this%dv%x,this%dw%x,1d0,n,this%msh%gdim)

  end subroutine fluid_plan4_step
  
  subroutine plan4_pres_setup(this)
    type(fluid_plan4_t) :: this
    call invers2(this%c_Xh%h1,this%vtrans,this%dm_Xh%n_dofs)
    this%c_Xh%ifh2 = .false.

  end subroutine plan4_pres_setup

  subroutine plan4_vel_setup(this)
    type(fluid_plan4_t) :: this
    real(kind=dp) :: dtbd
    dtbd = this%bd(1)/this%params%dt
    call copy(this%c_Xh%h1,this%vdiff(1,1,1,1),this%dm_Xh%n_dofs)
    call cmult2(this%c_Xh%h2,this%vtrans(1,1,1,1),dtbd,this%dm_Xh%n_dofs)
    this%c_Xh%ifh2 = .true.
  end subroutine plan4_vel_setup

  subroutine plan4_vel_residual(this) 
    type(fluid_plan4_t) :: this
    real(kind=dp) :: scl

    call this%Ax%compute(this%u_res, this%u%x, this%c_Xh, this%msh, this%Xh)
    call this%Ax%compute(this%v_res, this%v%x, this%c_Xh, this%msh, this%Xh)
    if (this%msh%gdim .eq. 3) then
       call this%Ax%compute(this%w_res, this%w%x, this%c_Xh, this%msh,t his%Xh)
    end if
    call opchsign(this%u_res, this%v_res, this%w_res, &
         this%msh%gdim, this%dm_Xh%n_dofs)

    !scl = -1./3.

    !call col3    (this%ta4,this%vdiff,qtl,ntot)
    !call rzero(this%ta4,this%dm_Xh%n_dofs)
    !call add2s1  (this%ta4,this%p,scl,this%dm_Xh%n_dofs)
    !call opgrad  (this%ta1,this%ta2,this%ta3,this%ta4,this%c_Xh)
    !if(IFAXIS) then
    !  !   CALL COL2 (TA2, OMASK,NTOT)
    !  !   CALL COL2 (TA3, OMASK,NTOT)
    !  !endif
    !call opsub2  (this%u_res,this%v_res,this%w_res,this%ta1,this%ta2,this%ta3)
    call opadd2cm(this%u_res, this%v_res, this%w_res, this%bfx, &
         this%bfy, this%bfz,1d0, this%dm_Xh%n_dofs, this%msh%gdim)


  end subroutine plan4_vel_residual

  subroutine plan4_pres_residual(this)
     type(fluid_plan4_t), intent(inout), target :: this
     type(coef_t), pointer :: coef
     real(kind=dp) :: scl
     integer :: i, n
      
     coef => this%c_Xh
     n  = this%dm_Xh%n_dofs
     call plan4_op_curl (this%ta1,this%ta2,this%ta3,this%u_e,this%u_e,this%w_e, this)
     call plan4_op_curl  (this%wa1,this%wa2,this%wa3,this%ta1,this%ta2,this%ta3, this)
     call opcolv   (this%wa1,this%wa2,this%wa3,coef%B,this%msh%gdim, n)
     scl = -4./3. 
     call opadd2cm (this%wa1,this%wa2,this%wa3,this%ta1,this%ta2,this%ta3,scl,n,this%msh%gdim)
     call invcol3  (this%w1,this%vdiff,this%vtrans,n)
     call opcolv   (this%wa1,this%wa2,this%wa3,this%w1, this%msh%gdim, n)

     !BOUNDARY CONDITION, DIRICHLET PRESSURE!
     !call bcdirpr (pr)

     call this%Ax%compute(this%p_res,this%p%x,this%c_Xh,this%msh,this%Xh)
     call chsign  (this%p_res,n)

     do i=1,n
        this%ta1(i,1) = this%bfx(i,1,1,1)/this%vtrans(i,1,1,1)-this%wa1(i)
        this%ta2(i,1) = this%bfy(i,1,1,1)/this%vtrans(i,1,1,1)-this%wa2(i)
        this%ta3(i,1) = this%bfz(i,1,1,1)/this%vtrans(i,1,1,1)-this%wa3(i)
     enddo
     
     !Need to consider cyclic bcs here...
     call gs_op_vector(this%gs_Xh, this%ta1, n, GS_OP_ADD) 
     call gs_op_vector(this%gs_Xh, this%ta2, n, GS_OP_ADD) 
     call gs_op_vector(this%gs_Xh, this%ta3, n, GS_OP_ADD) 

     do i=1,n
        this%ta1(i,1) = this%ta1(i,1)*this%c_Xh%Binv(i,1,1,1)
        this%ta2(i,1) = this%ta2(i,1)*this%c_Xh%Binv(i,1,1,1)
        this%ta3(i,1) = this%ta3(i,1)*this%c_Xh%Binv(i,1,1,1)
     enddo

     if (this%msh%gdim .eq. 3) then
         call cdtp    (this%wa1,this%ta1,coef%drdx,coef%dsdx,coef%dtdx,coef)
         call cdtp    (this%wa2,this%ta2,coef%drdy,coef%dsdy,coef%dtdy,coef)
         call cdtp    (this%wa3,this%ta3,coef%drdz,coef%dsdz,coef%dtdz,coef)
         do i=1,n
            this%p_res(i,1) = this%p_res(i,1)+this%wa1(i)+this%wa2(i)+this%wa3(i)
         enddo
      else
         call cdtp    (this%wa1,this%ta1,coef%drdx,coef%dsdx,coef%dtdx,coef)
         call cdtp    (this%wa2,this%ta2,coef%drdy,coef%dsdy,coef%dtdy,coef)

         do i=1,n
            this%p_res(i,1) = this%p_res(i,1)+this%wa1(i)+this%wa2(i)
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

  subroutine plan4_op_curl(w1,w2,w3,u1,u2,u3,this)
    type(fluid_plan4_t), intent(inout), target :: this
    real(kind=dp), dimension(this%dm_Xh%n_dofs), intent(inout) :: w1,w2,w3,u1,u2,u3
    type(coef_t), pointer :: coef
    integer :: lxyz, n
    coef => this%c_Xh
    lxyz  = this%Xh%lx*this%Xh%ly*this%Xh%lz
    n = this%dm_Xh%n_dofs
!     this%work1=dw/dy ; this%work2=dv/dz
call dudxyz(this%work1,u3,coef%drdy,coef%dsdy,coef%dtdy,coef)
        if (this%msh%gdim .eq. 3) then
           call dudxyz(this%work2,u2,coef%drdz,coef%dsdz,coef%dtdz,this%c_Xh)
           call sub3(w1,this%work1,this%work2,n)
        else
           call copy(w1,this%work1,n)
        endif
!     this%work1=du/dz ; this%work2=dw/dx
        if (this%msh%gdim .eq. 3) then
           call dudxyz(this%work1,u1,coef%drdz,coef%dsdz,coef%dtdz,this%c_Xh)
           call dudxyz(this%work2,u3,coef%drdx,coef%dsdx,coef%dtdx,this%c_Xh)
           call sub3(w2,this%work1,this%work2,n)
        else
           call rzero (this%work1,n)
           call dudxyz(this%work2,u3,coef%drdx,coef%dsdx,coef%dtdx,this%c_Xh)
           call sub3(w2,this%work1,this%work2,n)
        endif
!     this%work1=dv/dx ; this%work2=du/dy
        call dudxyz(this%work1,u2,coef%drdx,coef%dsdx,coef%dtdx,this%c_Xh)
        call dudxyz(this%work2,u1,coef%drdy,coef%dsdy,coef%dtdy,this%c_Xh)
        call sub3(w3,this%work1,this%work2,n)
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
