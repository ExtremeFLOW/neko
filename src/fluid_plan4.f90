!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4
  use fluid_method
  use facet_normal
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
     real(kind=dp), allocatable :: ulag(:,:,:,:,:)
     real(kind=dp), allocatable :: vlag(:,:,:,:,:)
     real(kind=dp), allocatable :: wlag(:,:,:,:,:)

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

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs

     !> Time variables
     real(kind=dp) :: ab(10), bd(10),dt_old(10)
     real(kind=dp), allocatable :: abx1(:,:,:,:), aby1(:,:,:,:), abz1(:,:,:,:)
     real(kind=dp), allocatable :: abx2(:,:,:,:), aby2(:,:,:,:), abz2(:,:,:,:)
     real(kind=dp) :: t_old
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
    
    ! Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, solver_vel=vel, solver_prs=prs)
    
    ! Initialize variables specific to this plan
    allocate(this%p_res(this%dm_Xh%n_dofs))
    allocate(this%u_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%v_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%w_res(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%ulag(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv,2))
    allocate(this%vlag(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv,2))
    allocate(this%wlag(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv,2))
    
    
    allocate(this%abx1(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%aby1(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%abz1(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    
    allocate(this%abx2(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%aby2(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
    allocate(this%abz2(this%Xh%lx,this%Xh%ly,this%Xh%lz,this%msh%nelv))
            
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

    ! Initialize velocity surface terms in pressure rhs
    call this%bc_prs_surface%init(this%dm_Xh)
    call this%bc_prs_surface%mark_zone(msh%inlet)
    call this%bc_prs_surface%finalize()
    call this%bc_prs_surface%set_coef(this%c_Xh)
    
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
    
    if (allocated(this%ulag)) then
       deallocate(this%ulag)
    end if
       
    if (allocated(this%vlag)) then
       deallocate(this%vlag)
    end if
    
    if (allocated(this%wlag)) then
       deallocate(this%wlag)
    end if
    
    if (allocated(this%abx1)) then
       deallocate(this%abx1)
    end if
    if (allocated(this%aby1)) then
       deallocate(this%aby1)
    end if
    if (allocated(this%abz1)) then
       deallocate(this%abz1)
    end if
    if (allocated(this%abx2)) then
       deallocate(this%abx2)
    end if
    if (allocated(this%aby2)) then
       deallocate(this%aby2)
    end if
    if (allocated(this%abz2)) then
       deallocate(this%abz2)
    end if
  end subroutine fluid_plan4_free
  
  subroutine fluid_plan4_step(this, t, tstep)
    class(fluid_plan4_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: t
    integer, intent(inout) :: tstep
    integer tt
    integer :: n, iter, i, nab, nbd, niter
    n = this%dm_Xh%n_dofs
    tt = tstep
    niter = 1000

    call settime(t, this%params%dt, this%t_old, this%dt_old,&
                 this%ab, this%bd, nab, nbd, tt)
    
    call plan4_sumab(this%u_e%x,this%u%x,this%ulag,n,this%ab,nab)
    call plan4_sumab(this%v_e%x,this%v%x,this%vlag,n,this%ab,nab)
    if (this%dm_Xh%msh%gdim .eq. 3) then
       call plan4_sumab(this%w_e%x,this%w%x,this%wlag,n,this%ab,nab)
    end if

    call this%f_Xh%eval()
    
    call makeabf(this%ta1, this%ta2, this%ta3,&
                 this%abx1, this%aby1, this%abz1,&
                 this%abx2, this%aby2, this%abz2, &
                 this%f_Xh%u, this%f_Xh%v, this%f_Xh%w,&
                 this%params%rho, this%ab, n, this%msh%gdim)
    call makebdf(this%ta1, this%ta2, this%ta3,&
                 this%wa1%x, this%wa2%x, this%wa3%x,&
                 this%c_Xh%h2, this%ulag, this%vlag, this%wlag, &
                 this%f_Xh%u, this%f_Xh%v, this%f_Xh%w,&
                 this%u, this%v, this%w,&
                 this%c_Xh%B, this%params%rho, this%params%dt, &
                 this%bd, nbd, n, this%msh%gdim)

    
    do i = 3-1,2,-1
       call copy(this%ulag(1,1,1,1,i), this%ulag(1,1,1,1,i-1), n)
       call copy(this%vlag(1,1,1,1,i), this%vlag(1,1,1,1,i-1), n)
       call copy(this%wlag(1,1,1,1,i), this%wlag(1,1,1,1,i-1), n)
    end do
    
    call copy(this%ulag, this%u%x, n)
    call copy(this%vlag, this%v%x, n)
    call copy(this%wlag, this%w%x, n)
    

    ! mask Dirichlet boundaries (velocity)
    call this%bc_apply_vel()

    ! compute pressure
    call this%bc_apply_prs()
    call plan4_pres_setup(this%c_Xh%h1, this%c_Xh%h2, this%params%rho, &
                          this%dm_Xh%n_dofs, this%c_Xh%ifh2)    
    call plan4_pres_residual(this%p,this%p_res, this%u, this%v, this%w, &
                             this%u_e, this%v_e, this%w_e, &
                             this%ta1, this%ta2, this%ta3, &
                             this%wa1, this%wa2, this%wa3, &
                             this%work1, this%work2, this%f_Xh, &
                             this%c_Xh, this%gs_Xh, this%bc_prs_surface, &
                             this%Ax,this%bd(1), this%params%dt, &
                             this%params%Re, this%params%rho)

    !Sets tolerances
    !call ctolspl  (tolspl,respr)
    call gs_op_vector(this%gs_Xh, this%p_res, n, GS_OP_ADD) 
    call bc_list_apply_scalar(this%bclst_prs, this%p_res, this%p%dof%n_dofs)
    select type(pcp => this%pc_prs)
    type is(jacobi_t)
       call jacobi_set_d(pcp,this%c_Xh, this%dm_Xh, this%gs_Xh)
    end select
    write(*,*) "PRES"
    iter = this%ksp_prs%solve(this%Ax, this%dp, this%p_res, n, &
         this%c_Xh, this%bclst_prs, this%gs_Xh, niter)    
    call add2(this%p%x,this%dp%x,n)
!    call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)
    
    !We only need to update h2 once I think then use the flag to switch on/off
    call plan4_vel_setup(this%c_Xh%h1, this%c_Xh%h2, &
         this%params%Re, this%params%rho, this%bd(1), &
         this%params%dt, this%dm_Xh%n_dofs, this%c_Xh%ifh2)
    
    call plan4_vel_residual(this%Ax, this%u, this%v, this%w, &
                            this%u_res, this%v_res, this%w_res, &
                            this%p, this%ta1%x, this%ta2%x, this%ta3%x, this%ta4%x, &
                            this%f_Xh, this%c_Xh, &
                            this%msh, this%Xh, this%dm_Xh%n_dofs)

    call gs_op_vector(this%gs_Xh, this%u_res, n, GS_OP_ADD) 
    call gs_op_vector(this%gs_Xh, this%v_res, n, GS_OP_ADD) 
    call gs_op_vector(this%gs_Xh, this%w_res, n, GS_OP_ADD) 
    call bc_list_apply_vector(this%bclst_res,&
         this%u_res, this%v_res, this%w_res, this%dm_Xh%n_dofs)
    select type(pcp =>this%pc_vel)
    type is(jacobi_t)
       call jacobi_set_d(pcp,this%c_Xh, this%dm_Xh, this%gs_Xh)
    end select

    write(*,*) 'U'
    iter = this%ksp_vel%solve(this%Ax, this%du, this%u_res, n, &
         this%c_Xh, this%bclst_res, this%gs_Xh, niter)
    write(*,*) 'V'
    iter = this%ksp_vel%solve(this%Ax, this%dv, this%v_res, n, &
         this%c_Xh, this%bclst_res, this%gs_Xh, niter)
    write(*,*) 'W'
    iter = this%ksp_vel%solve(this%Ax, this%dw, this%w_res, n, &
         this%c_Xh, this%bclst_res, this%gs_Xh, niter)

    call opadd2cm(this%u%x,this%v%x,this%w%x,this%du%x,this%dv%x,this%dw%x,1d0,n,this%msh%gdim)

  end subroutine fluid_plan4_step
  
  subroutine plan4_pres_setup(h1, h2, rho, n, ifh2)
    integer, intent(inout) :: n
    real(kind=dp), intent(inout) :: h1(n)
    real(kind=dp), intent(inout) :: h2(n)
    real(kind=dp), intent(inout) :: rho
    logical, intent(inout) :: ifh2
    call rone(h1, n)
    call cmult(h1, 1d0 /rho, n)
    call rzero(h2, n)
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
    h2 = dtbd
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

    scl = -1d0/3d0

    call rzero(ta4,c_xh%dof%n_dofs)
    call add2s1  (ta4,p%x,scl,c_Xh%dof%n_dofs)
    call opgrad  (ta1,ta2,ta3,ta4,c_Xh)

    call opadd2cm(u_res, v_res, w_res, ta1, ta2, ta3, -1d0, n, msh%gdim)

    call opadd2cm(u_res, v_res, w_res, f_Xh%u, f_Xh%v, f_Xh%w, 1d0, n, msh%gdim)

  end subroutine plan4_vel_residual

  subroutine plan4_pres_residual(p, p_res, u, v, w, u_e, v_e, w_e, &
       ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, f_Xh, c_xh, gs_Xh, &
       bc_prs_surface, Ax, bd, dt, Re, rho)
    type(field_t), intent(inout) :: p
    real(kind=dp), intent(inout) :: p_res(p%dof%n_dofs)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
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
    type(facet_normal_t), intent(inout) :: bc_prs_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=dp), intent(inout) :: bd
    real(kind=dp), intent(inout) :: dt
    real(kind=dp), intent(inout) :: Re
    real(kind=dp), intent(inout) :: rho
    real(kind=dp) :: scl, dtbd
    integer :: i, idx(4)
    integer :: n, gdim, glb_n,m, k

    n = c_Xh%dof%size()
    gdim = c_Xh%msh%gdim
    glb_n = n / p%msh%nelv * p%msh%glb_nelv
    
    call plan4_op_curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call plan4_op_curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)
    call opcolv(wa1%x, wa2%x, wa3%x, c_Xh%B, gdim, n)
    scl = -4d0 / 3d0
    ta1 = 0d0
    ta2 = 0d0
    ta3 = 0d0
    call opadd2cm (wa1%x, wa2%x, wa3%x, ta1%x, ta2%x, ta3%x, scl, n, gdim)

    work1%x = (1d0 / Re) / rho
    call opcolv(wa1%x, wa2%x, wa3%x, work1%x, gdim, n)

    call Ax%compute(p_res,p%x,c_Xh,p%msh,p%Xh)
    call chsign  (p_res,n)

    do i=1,n
       ta1%x(i,1,1,1) = f_Xh%u(i,1,1,1) / rho - wa1%x(i,1,1,1)
       ta2%x(i,1,1,1) = f_Xh%v(i,1,1,1) / rho - wa2%x(i,1,1,1)
       ta3%x(i,1,1,1) = f_Xh%w(i,1,1,1) / rho - wa3%x(i,1,1,1)
    enddo
     
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


    !
    ! Surface velocity terms
    !
    dtbd = bd/dt
    call rzero(ta1%x, n)
    call rzero(ta2%x, n)
    call rzero(ta3%x, n)
    call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)
    call add2(ta1%x, ta2%x, n)
    call add2(ta1%x, ta3%x, n)    
    call cmult(ta1%x, dtbd, n)
    call sub2(p_res, ta1%x, n)
        
!    call ortho(p_res,n,glb_n) ! Orthogonalize wrt null space, if present

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
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine plan4_op_curl
  
  !> Sum up AB/BDF contributions 
  subroutine plan4_sumab(v,vv,vvlag,n,ab,nab)
    integer, intent(in) :: n, nab
    real(kind=dp), dimension(n), intent(inout) :: v,vv
    real(kind=dp), dimension(n,2), intent(inout) :: vvlag
    real(kind=dp), dimension(3), intent(in) :: ab
    real(kind=dp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    call add3s2(v,vv,vvlag(1,1),ab0,ab1,n)
    !should we have this right now we only save one old velocity
    if(nab .eq. 3) call add2s2(v,vvlag(1,2),ab2,n)
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
!      if ((tstep .eq. 0) .or. (tstep .eq. 1)) nbd = 1
!      if ((tstep .gt. 1) .and. (tstep .le. 2)) nbd = tstep
!      if (tstep.gt.2) nbd = 3
      !      call rzero (bd, 10)
      nbd = min(tstep, 3)
      call setbd(bd, dt_old, nbd)
      ! This can also be varied, hardcoded for now
      nab = min(tstep, 3)
      !dont know what irst is really
      !IF (ISTEP.lt.NAB.and.irst.le.0) NAB = ISTEP
!      if(tstep .lt. 3) NAB = tstep
!      if(tstep .ge. 3) NAB = 3
 !     call rzero(ab, 10)
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
    REAL(KIND=DP), intent(inout), dimension(NAB) :: AB(10),DTLAG(10)
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
    integer :: I, NSYS
    IF (NBD.EQ.1) THEN
         BD(1) = 1d0
         BDF   = 1d0
      ELSEIF (NBD.GE.2) THEN
         NSYS = NBD+1
         CALL BDSYS (BDMAT,BDRHS,DTBD,NBD,ldim)
         CALL LU    (BDMAT,NSYS,ldim,IR,IC)
         CALL SOLVE (BDRHS,BDMAT,1,NSYS,ldim,IR,IC)
         DO I=1,NBD
            BD(I) = BDRHS(I)
         end do
         BDF = BDRHS(NBD+1)
      ENDIF
    !Normalize
      DO IBD=NBD,1,-1
         BD(IBD+1) = BD(IBD)
      end do
      BD(1) = 1d0
      DO IBD=1,NBD+1
         BD(IBD) = BD(IBD)/BDF
      end do
    end subroutine setbd


    subroutine bdsys (a,b,dt,nbd,ldim)
      real(kind=dp) ::  A(ldim,9),B(9),DT(9)
      integer :: ldim, j, n, k, i, nsys, nbd
      real(kind=dp) :: SUMDT
      CALL RZERO (A,ldim**2)
      N = NBD+1
      DO  J=1,NBD
         A(1,J) = 1.
      end DO
      A(1,NBD+1) = 0.
      B(1) = 1.
      DO J=1,NBD
         SUMDT = 0.
         DO  K=1,J
            SUMDT = SUMDT+DT(K)
         end DO
         A(2,J) = SUMDT
      end DO
      A(2,NBD+1) = -DT(1)
      B(2) = 0.
      DO I=3,NBD+1
         DO  J=1,NBD
            SUMDT = 0.
            DO  K=1,J
               SUMDT = SUMDT+DT(K)
            end DO
            A(I,J) = SUMDT**(I-1)
         end DO
         A(I,NBD+1) = 0.
         B(I) = 0.
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
100     CONTINUE
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
200     CONTINUE
300     IF(M.EQ.K) GOTO 500
      DO 400 I=1,N
      B=A(I,K)
      A(I,K)=A(I,M)
       A(I,M)=B
400    CONTINUE
500     C=1./A(K,K)
      A(K,K)=C
      IF(K.EQ.N) GOTO 1000
      K1=K+1
      XMAX=ABS(A(K1,K1))
      L=K1
      M=K1
      DO 600 I=K1,N
       A(I,K)=C*A(I,K)
600     CONTINUE
      DO 800 I=K1,N
      B=A(I,K)
      DO 800 J=K1,N
      A(I,J)=A(I,J)-B*A(K,J)
      Y=ABS(A(I,J))
      IF(XMAX.GE.Y) GOTO 800
      XMAX=Y
      L=I
      M=J
800    CONTINUE
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
100     CONTINUE
      DO 400 I=2,N
      I1=I-1
      B=G(I)
      DO 300 J=1,I1
        B=B-A(I,J)*G(J)
300     CONTINUE
        G(I)=B
400     CONTINUE
      DO 700 IT=1,N
      I=N1-IT
      I1=I+1
      B=G(I)
      IF(I.EQ.N) GOTO 701
      DO 600 J=I1,N
        B=B-A(I,J)*G(J)
600     CONTINUE
701     G(I)=B*A(I,I)
700     CONTINUE
      DO 900 I=1,N
      ICI=IC(I)
        F(ICI,KK)=G(I)
900     CONTINUE
1000    CONTINUE
      RETURN
      END
    

 subroutine makebdf(ta1, ta2, ta3, tb1, tb2, tb3, h2, ulag, vlag, wlag, &
                    bfx, bfy, bfz, u, v, w, B, rho, dt,bd, nbd, n, gdim)
!
!     Add contributions to F from lagged BD terms.
!
    integer, intent(inout) :: n, nbd, gdim
    type(field_t) :: ta1, ta2, ta3, u, v, w
    real(kind=dp), intent(inout) :: TB1(n)
    real(kind=dp), intent(inout) :: TB2(n)
    real(kind=dp), intent(inout) :: TB3(n)
    real(kind=dp), intent(inout) :: BFX(n)
    real(kind=dp), intent(inout) :: BFY(n)
    real(kind=dp), intent(inout) :: BFZ(n)
    real(kind=dp), intent(inout) :: H2 (n), B(n)
    real(kind=dp), intent(inout) :: ulag(n,nbd), vlag(n,nbd), wlag(n,nbd)
    real(kind=dp), intent(inout) :: dt, rho, bd(10)
    real(kind=dp) :: const
    integer :: ilag
    CONST = rho /DT
    call rone(h2, n)
    call cmult(h2,const,n)
    CALL OPCOLV3c (TB1,TB2,TB3,u%x,v%x,w%x,B,bd(2),n,gdim)
    DO ILAG=2,NBD
       CALL OPCOLV3c(TA1%x,TA2%x,TA3%x,ulag (1,ILAG-1),&
                                      vlag (1,ILAG-1),&
                                      wlag (1,ILAG-1),&
                                      B,bd(ilag+1),n,gdim)
       CALL OPADD2cm  (TB1,TB2,TB3,TA1%x,TA2%x,TA3%x, 1d0, n, gdim)
   ENDDO
   CALL OPADD2col (BFX,BFY,BFZ,TB1,TB2,TB3,h2, n, gdim)
  END subroutine makebdf
 subroutine makeabf(ta1, ta2, ta3, abx1, aby1, abz1, abx2, aby2, abz2, &
                    bfx, bfy, bfz, rho, ab, n, gdim)
!-----------------------------------------------------------------------
!
!     Sum up contributions to kth order extrapolation scheme.
!
!-----------------------------------------------------------------------
    type(field_t), intent(inout) :: ta1, ta2, ta3
    real(kind=dp), intent(inout) :: rho, ab(10)
     integer, intent(in) :: n, gdim
    real(kind=dp), intent(inout) :: BFX(n)
    real(kind=dp), intent(inout) :: BFY(n)
    real(kind=dp), intent(inout) :: BFZ(n)
    real(kind=dp), intent(inout) :: abx1(n), aby1(n), abz1(n)
    real(kind=dp), intent(inout) :: abx2(n), aby2(n), abz2(n)
    real(kind=dp) :: ab0, ab1, ab2
      AB0 = AB(1)
      AB1 = AB(2)
      AB2 = AB(3)
      CALL ADD3S2 (ta1%x,ABX1,ABX2,AB1,AB2,n)
      CALL ADD3S2 (TA2%x,ABY1,ABY2,AB1,AB2,n)
      CALL COPY   (ABX2,ABX1,n)
      CALL COPY   (ABY2,ABY1,N)
      CALL COPY   (ABX1,BFX,N)
      CALL COPY   (ABY1,BFY,N)
      CALL ADD2S1 (BFX,ta1%x,AB0,N)
      CALL ADD2S1 (BFY,TA2%x,AB0,N)
      CALL Cmult   (BFX,rho,N)          ! multiply by density
      CALL Cmult   (BFY,rho,N)
      IF (gdim.EQ.3) THEN
         CALL ADD3S2 (TA3%x,ABZ1,ABZ2,AB1,AB2,N)
         CALL COPY   (ABZ2,ABZ1,N)
         CALL COPY   (ABZ1,BFZ,N)
         CALL ADD2S1 (BFZ,TA3%x,AB0,N)
         CALL cmult   (BFZ,rho,N)
      ENDIF
      END subroutine makeabf
end module fluid_plan4
