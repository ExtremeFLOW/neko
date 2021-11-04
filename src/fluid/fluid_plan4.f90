!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4  
  use ax_helm_fctry
  use fluid_method
  use facet_normal
  use neko_config
  use abbdf
  use projection
  use logger
  use advection
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan  
  implicit none
  private

  type, public, extends(fluid_scheme_t) :: fluid_plan4_t
     type(field_t) :: u_e
     type(field_t) :: v_e
     type(field_t) :: w_e

     real(kind=rp), allocatable :: p_res(:)
     real(kind=rp), allocatable :: u_res(:,:,:,:)
     real(kind=rp), allocatable :: v_res(:,:,:,:)
     real(kind=rp), allocatable :: w_res(:,:,:,:)
     real(kind=rp), allocatable :: ulag(:,:,:,:,:)
     real(kind=rp), allocatable :: vlag(:,:,:,:,:)
     real(kind=rp), allocatable :: wlag(:,:,:,:,:)

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
     
     !> @todo move this to a scratch space
     type(field_t) :: work1
     type(field_t) :: work2

     class(ax_t), allocatable :: Ax
     
     type(projection_t) :: proj

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_residual   !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_residual  

     class(advection_t), allocatable :: adv 

     ! Time variables
     real(kind=rp), allocatable :: abx1(:,:,:,:), aby1(:,:,:,:), abz1(:,:,:,:)
     real(kind=rp), allocatable :: abx2(:,:,:,:), aby2(:,:,:,:), abz2(:,:,:,:)

     ! Vol_flow
     
     integer :: flow_dir !< these two should be moved to params
     logical :: avflow 
     real(kind=rp) :: flow_rate 
     real(kind=rp) :: dtlag = 0d0
     real(kind=rp) :: bdlag = 0d0!< Really quite pointless since we do not vary the timestep
     type(field_t) :: u_vol, v_vol, w_vol, p_vol
     real(kind=rp) :: domain_length, base_flow
   contains
     procedure, pass(this) :: init => fluid_plan4_init
     procedure, pass(this) :: free => fluid_plan4_free
     procedure, pass(this) :: step => fluid_plan4_step
  end type fluid_plan4_t

contains

  subroutine fluid_plan4_init(this, msh, lx, param)    
    class(fluid_plan4_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: param     

    call this%free()
    
    ! Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, .true., .true.)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax)
    call advection_factory(this%adv)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      allocate(this%p_res(dm_Xh%n_dofs))
      allocate(this%u_res(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%v_res(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%w_res(Xh_lx, Xh_ly, Xh_lz, nelv))
      
      allocate(this%ulag(Xh_lx, Xh_ly, Xh_lz, nelv,2))
      allocate(this%vlag(Xh_lx, Xh_ly, Xh_lz, nelv,2))
      allocate(this%wlag(Xh_lx, Xh_ly, Xh_lz, nelv,2))
          
      allocate(this%abx1(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%aby1(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%abz1(Xh_lx, Xh_ly, Xh_lz, nelv))
    
      allocate(this%abx2(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%aby2(Xh_lx, Xh_ly, Xh_lz, nelv))
      allocate(this%abz2(Xh_lx, Xh_ly, Xh_lz, nelv))
      
      call rzero(this%abx1, dm_Xh%n_dofs)
      call rzero(this%aby1, dm_Xh%n_dofs)
      call rzero(this%abz1, dm_Xh%n_dofs)
      
      call rzero(this%abx2, dm_Xh%n_dofs)
      call rzero(this%aby2, dm_Xh%n_dofs)
      call rzero(this%abz2, dm_Xh%n_dofs)
            
      call field_init(this%u_e, dm_Xh, 'u_e')
      call field_init(this%v_e, dm_Xh, 'v_e')
      call field_init(this%w_e, dm_Xh, 'w_e')
    
      call field_init(this%wa1, dm_Xh, 'wa1')
      call field_init(this%wa2, dm_Xh, 'wa2')
      call field_init(this%wa3, dm_Xh, 'wa3')

      call field_init(this%ta1, dm_Xh, 'ta1')
      call field_init(this%ta2, dm_Xh, 'ta2')
      call field_init(this%ta3, dm_Xh, 'ta3')
    
      call field_init(this%du, dm_Xh, 'du')
      call field_init(this%dv, dm_Xh, 'dv')
      call field_init(this%dw, dm_Xh, 'dw')
      call field_init(this%dp, dm_Xh, 'dp')

      call field_init(this%work1, dm_Xh, 'work1')
      call field_init(this%work2, dm_Xh, 'work2')

    end associate
    
    ! Initialize velocity surface terms in pressure rhs
    call this%bc_prs_surface%init(this%dm_Xh)
    call this%bc_prs_surface%mark_zone(msh%inlet)
    call this%bc_prs_surface%finalize()
    call this%bc_prs_surface%set_coef(this%c_Xh)

    ! Initialize boundary condition for velocity residual
    call this%bc_vel_residual%init(this%dm_Xh)
    call this%bc_vel_residual%mark_zone(msh%inlet)
    call this%bc_vel_residual%mark_zone(msh%wall)
    call this%bc_vel_residual%finalize()
    call this%bc_vel_residual%set_g(0.0_rp)
    call bc_list_init(this%bclst_vel_residual)
    call bc_list_add(this%bclst_vel_residual, this%bc_vel_residual)

    !Intialize projection space thingy
    call this%proj%init(this%dm_Xh%n_dofs, param%proj_dim)

    !Initialize vol_flow (if there is a forced voume flow)
    this%flow_dir = param%vol_flow_dir
    this%avflow = param%avflow
    this%flow_rate = param%flow_rate
    
    call field_init(this%u_vol, this%dm_Xh, 'u_vol')
    call field_init(this%v_vol, this%dm_Xh, 'v_vol')
    call field_init(this%w_vol, this%dm_Xh, 'w_vol')
    call field_init(this%p_vol, this%dm_Xh, 'p_vol')

    ! Add lagged term to checkpoint
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)    

  end subroutine fluid_plan4_init

  subroutine fluid_plan4_free(this)
    class(fluid_plan4_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call this%bc_prs_surface%free()  
    call bc_list_free(this%bclst_vel_residual)
    call this%proj%free()
   
    call field_free(this%u_e)
    call field_free(this%v_e)
    call field_free(this%w_e)
    
    call field_free(this%wa1)
    call field_free(this%wa2)
    call field_free(this%wa3)

    call field_free(this%ta1)
    call field_free(this%ta2)
    call field_free(this%ta3)

    call field_free(this%du)
    call field_free(this%dv)
    call field_free(this%dw)
    call field_free(this%dp)
    
    call field_free(this%u_vol)
    call field_free(this%v_vol)
    call field_free(this%w_vol)
    call field_free(this%p_vol)

    call field_free(this%work1)
    call field_free(this%work2)

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
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
  
  subroutine fluid_plan4_step(this, t, tstep, ab_bdf)
    class(fluid_plan4_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    type(abbdf_t), intent(inout) :: ab_bdf
    integer, intent(inout) :: tstep
    integer :: n, i, niter
    type(ksp_monitor_t) :: ksp_results(4)
    real(kind=rp), parameter :: one = 1.0
    n = this%dm_Xh%n_dofs
    niter = 1000

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         ta1 => this%ta1, ta2 => this%ta2, ta3 => this%ta3, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         params => this%params, msh => this%msh)

      call fluid_plan4_sumab(u_e%x, u%x,this%ulag,n,ab_bdf%ab,ab_bdf%nab)
      call fluid_plan4_sumab(v_e%x, v%x,this%vlag,n,ab_bdf%ab,ab_bdf%nab)
      if (msh%gdim .eq. 3) then
         call fluid_plan4_sumab(w_e%x, w%x,this%wlag,n,ab_bdf%ab,ab_bdf%nab)
      end if

      call f_Xh%eval()
      call opcolv(f_Xh%u, f_Xh%v, f_Xh%w, c_Xh%B, msh%gdim, n)
      call this%adv%apply(ta1, ta2, ta3, &
                 this%u, this%v, this%w, &
                 f_Xh%u, f_Xh%v, f_Xh%w, &
                 Xh, this%c_Xh, msh%nelv, &
                 dm_Xh%n_dofs, msh%gdim)
   
      call makeabf(ta1, ta2, ta3,&
                  this%abx1, this%aby1, this%abz1,&
                  this%abx2, this%aby2, this%abz2, &
                  f_Xh%u, f_Xh%v, f_Xh%w,&
                  params%rho, ab_bdf%ab, n, msh%gdim)
      call makebdf(ta1, ta2, ta3,&
                   this%wa1%x, this%wa2%x, this%wa3%x,&
                   c_Xh%h2, this%ulag, this%vlag, this%wlag, &
                   f_Xh%u, f_Xh%v, f_Xh%w, u, v, w,&
                   c_Xh%B, params%rho, params%dt, &
                   ab_bdf%bd, ab_bdf%nbd, n, msh%gdim)

    
      do i = 3-1,2,-1
         call copy(this%ulag(1,1,1,1,i), this%ulag(1,1,1,1,i-1), n)
         call copy(this%vlag(1,1,1,1,i), this%vlag(1,1,1,1,i-1), n)
         call copy(this%wlag(1,1,1,1,i), this%wlag(1,1,1,1,i-1), n)
      end do
    
      call copy(this%ulag, u%x, n)
      call copy(this%vlag, v%x, n)
      call copy(this%wlag, w%x, n)
      

      ! mask Dirichlet boundaries (velocity)
      call this%bc_apply_vel()
      
      ! compute pressure
      call this%bc_apply_prs()
      call fluid_plan4_pres_setup(c_Xh%h1, c_Xh%h2, params%rho, &
                                  dm_Xh%n_dofs, c_Xh%ifh2)    
      call fluid_plan4_pres_residual(p, p_res, u, v, w, &
                                     u_e, v_e, w_e, &
                                     ta1, ta2, ta3, &
                                     this%wa1, this%wa2, this%wa3, &
                                     this%work1, this%work2, f_Xh, &
                                     c_Xh, gs_Xh, this%bc_prs_surface, &
                                     Ax, ab_bdf%bd(1), params%dt, &
                                     params%Re, params%rho)

      !Sets tolerances
      !call ctolspl  (tolspl,respr)
      call gs_op_vector(gs_Xh, p_res, n, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_prs, p_res, p%dof%n_dofs)

      if( tstep .gt. 5) call this%proj%project_on(p_res, Ax, c_Xh, &
                                this%bclst_prs, gs_Xh, n)
      call this%pc_prs%update()
      ksp_results(1) = this%ksp_prs%solve(Ax, dp, p_res, n, c_Xh, &
                                this%bclst_prs, gs_Xh, niter)    
      if( tstep .gt. 5) call this%proj%project_back(dp%x, Ax, c_Xh, &
                                  this%bclst_prs, gs_Xh, n)
      call add2(p%x, dp%x,n)
      !    call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)
    
      !We only need to update h2 once I think then use the flag to switch on/off
      call fluid_plan4_vel_setup(c_Xh%h1, c_Xh%h2, &
                                 params%Re, params%rho, ab_bdf%bd(1), &
                                 params%dt, dm_Xh%n_dofs, c_Xh%ifh2)
    
      call fluid_plan4_vel_residual(Ax, u, v, w, &
                                    u_res, v_res, w_res, &
                                    p, ta1%x, ta2%x, ta3%x, &
                                    f_Xh, c_Xh, msh, Xh, dm_Xh%n_dofs)

      call gs_op_vector(gs_Xh, u_res, n, GS_OP_ADD) 
      call gs_op_vector(gs_Xh, v_res, n, GS_OP_ADD) 
      call gs_op_vector(gs_Xh, w_res, n, GS_OP_ADD) 

      call bc_list_apply_vector(this%bclst_vel_residual,&
                                u_res, v_res, w_res, dm_Xh%n_dofs)
      call this%pc_vel%update()

      ksp_results(2) = this%ksp_vel%solve(Ax, du, u_res, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      ksp_results(3) = this%ksp_vel%solve(Ax, dv, v_res, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      ksp_results(4) = this%ksp_vel%solve(Ax, dw, w_res, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      
      call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, one, n, msh%gdim)
     
      if (this%flow_dir .ne. 0) then
         call plan4_vol_flow(this, ab_bdf, niter)
      end if
      
      call fluid_step_info(tstep, t, params%dt, ksp_results)
      
    end associate
  end subroutine fluid_plan4_step
  
  subroutine fluid_plan4_pres_setup(h1, h2, rho, n, ifh2)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: h1(n)
    real(kind=rp), intent(inout) :: h2(n)
    real(kind=rp), intent(in) :: rho
    real(kind=rp), parameter :: one = 1.0
    logical, intent(inout) :: ifh2
    call rone(h1, n)
    call cmult(h1, one /rho, n)
    call rzero(h2, n)
    ifh2 = .false.
  end subroutine fluid_plan4_pres_setup

  subroutine fluid_plan4_vel_setup(h1, h2, Re, rho, bd, dt, n, ifh2)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: h1(n)
    real(kind=rp), intent(inout) :: h2(n)
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    logical, intent(inout) :: ifh2
    real(kind=rp), parameter :: one = 1.0
    real(kind=rp) :: dtbd    
    dtbd = rho * (bd / dt)
    h1 = (one / Re)
    h2 = dtbd
    ifh2 = .true.
  end subroutine fluid_plan4_vel_setup

  subroutine fluid_plan4_vel_residual(Ax, u, v, w, u_res, v_res, w_res, &
       p, ta1, ta2, ta3, f_Xh, c_Xh, msh, Xh, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh    
    type(field_t), intent(inout) :: p, u, v, w
    real(kind=rp), intent(inout) :: u_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: v_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: w_res(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: ta1(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: ta2(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: ta3(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    integer, intent(in) :: n
    
    call Ax%compute(u_res, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res, v%x, c_Xh, msh, Xh)
    if (msh%gdim .eq. 3) then
       call Ax%compute(w_res, w%x, c_Xh, msh, Xh)
    end if

    call opchsign(u_res, v_res, w_res, msh%gdim, n)

    call opgrad(ta1, ta2, ta3, p%x, c_Xh)

    call opadd2cm(u_res, v_res, w_res, ta1, ta2, ta3, -1.0_rp, n, msh%gdim)

    call opadd2cm(u_res, v_res, w_res, &
                  f_Xh%u, f_Xh%v, f_Xh%w, 1.0_rp, n, msh%gdim)

  end subroutine fluid_plan4_vel_residual

  subroutine fluid_plan4_pres_residual(p, p_res, u, v, w, u_e, v_e, w_e, &
       ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, f_Xh, c_xh, gs_Xh, &
       bc_prs_surface, Ax, bd, dt, Re, rho)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(inout) :: wa1, wa2, wa3
    type(field_t), intent(inout) :: work1, work2
    real(kind=rp), intent(inout) :: p_res(p%dof%n_dofs)
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(inout) :: bc_prs_surface
    class(Ax_t), intent(inout) :: Ax
    real(kind=rp), intent(inout) :: bd
    real(kind=rp), intent(in) :: dt
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp) :: dtbd
    integer :: n, gdim
    integer :: i        

    n = c_Xh%dof%size()
    gdim = c_Xh%msh%gdim
    
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)
    call opcolv(wa1%x, wa2%x, wa3%x, c_Xh%B, gdim, n)

    work1 = (1.0_rp / Re) / rho
    call opcolv(wa1%x, wa2%x, wa3%x, work1%x, gdim, n)

    call Ax%compute(p_res,p%x,c_Xh,p%msh,p%Xh)
    call chsign(p_res, n)

    do i = 1, n
       ta1%x(i,1,1,1) = f_Xh%u(i,1,1,1) / rho - wa1%x(i,1,1,1)
       ta2%x(i,1,1,1) = f_Xh%v(i,1,1,1) / rho - wa2%x(i,1,1,1)
       ta3%x(i,1,1,1) = f_Xh%w(i,1,1,1) / rho - wa3%x(i,1,1,1)
    enddo
     
     !Need to consider cyclic bcs here...
    call gs_op(gs_Xh, ta1, GS_OP_ADD) 
    call gs_op(gs_Xh, ta2, GS_OP_ADD) 
    call gs_op(gs_Xh, ta3, GS_OP_ADD) 

    do i = 1, n
       ta1%x(i,1,1,1) = ta1%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta2%x(i,1,1,1) = ta2%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
       ta3%x(i,1,1,1) = ta3%x(i,1,1,1) * c_Xh%Binv(i,1,1,1)
    enddo

    if (gdim .eq. 3) then
       call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
       call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       do i = 1, n
          p_res(i) = p_res(i) + wa1%x(i,1,1,1) + wa2%x(i,1,1,1) + wa3%x(i,1,1,1)
       enddo
    else
       call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)

       do i = 1, n
          p_res(i) = p_res(i)+wa1%x(i,1,1,1)+wa2%x(i,1,1,1)
       enddo
    endif


    !
    ! Surface velocity terms
    !
    dtbd = bd / dt
    call rzero(ta1%x, n)
    call rzero(ta2%x, n)
    call rzero(ta3%x, n)
    call bc_prs_surface%apply_surfvec(ta1%x, ta2%x, ta3%x, u%x, v%x, w%x, n)
    call add2(ta1%x, ta2%x, n)
    call add2(ta1%x, ta3%x, n)    
    call cmult(ta1%x, dtbd, n)
    call sub2(p_res, ta1%x, n)

        
!    call ortho(p_res,n,glb_n) ! Orthogonalize wrt null space, if present

  end subroutine fluid_plan4_pres_residual

  !> Sum up AB/BDF contributions 
  subroutine fluid_plan4_sumab(v,vv,vvlag,n,ab,nab)
    integer, intent(in) :: n, nab
    real(kind=rp), dimension(n), intent(inout) :: v, vv
    real(kind=rp), dimension(n,2), intent(inout) :: vvlag
    real(kind=rp), dimension(3), intent(in) :: ab
    real(kind=rp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    call add3s2(v,vv,vvlag(1,1),ab0,ab1,n)
    if(nab .eq. 3) call add2s2(v,vvlag(1,2),ab2,n)
  end subroutine fluid_plan4_sumab
  
  !> Add contributions to F from lagged BD terms.
  subroutine makebdf(ta1, ta2, ta3, tb1, tb2, tb3, h2, ulag, vlag, wlag, &
                     bfx, bfy, bfz, u, v, w, B, rho, dt, bd, nbd, n, gdim)
    integer, intent(in) :: n, nbd, gdim
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(in) :: u, v, w
    real(kind=rp), intent(inout) :: tb1(n), tb2(n), tb3(n)
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(inout) :: h2(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(inout) :: ulag(n,nbd), vlag(n,nbd), wlag(n,nbd)
    real(kind=rp), intent(in) :: dt, rho, bd(10)
    real(kind=rp) :: const
    integer :: ilag
    
    const = rho / dt
    call rone(h2, n)
    call cmult(h2,const,n)
    call opcolv3c(tb1, tb2, tb3, u%x, v%x, w%x, B, bd(2), n, gdim)
    do ilag = 2, nbd
       call opcolv3c(ta1%x, ta2%x, ta3%x, &
                     ulag(1,ilag-1), vlag(1,ilag-1), wlag(1,ilag-1), &
                     B, bd(ilag+1), n, gdim)
       call opadd2cm(tb1, tb2, tb3, ta1%x, ta2%x, ta3%x, 1.0_rp, n, gdim)
    end do
    call opadd2col(bfx, bfy, bfz, tb1, tb2, tb3, h2, n, gdim)
  end subroutine makebdf

  !> Sum up contributions to kth order extrapolation scheme.
  subroutine makeabf(ta1, ta2, ta3, abx1, aby1, abz1, abx2, aby2, abz2, &
                     bfx, bfy, bfz, rho, ab, n, gdim)
    type(field_t), intent(inout) :: ta1, ta2, ta3
    real(kind=rp), intent(inout) :: rho, ab(10)
     integer, intent(in) :: n, gdim
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(inout) :: abx1(n), aby1(n), abz1(n)
    real(kind=rp), intent(inout) :: abx2(n), aby2(n), abz2(n)
    real(kind=rp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)
    call add3s2(ta1%x, abx1, abx2, ab1, ab2, n)
    call add3s2(ta2%x, aby1, aby2, ab1, ab2, n)
    call copy(abx2, abx1, n)
    call copy(aby2, aby1, n)
    call copy(abx1, bfx, n)
    call copy(aby1, bfy, n)
    call add2s1(bfx, ta1%x, ab0, n)
    call add2s1(bfy, ta2%x, ab0, n)
    call cmult(bfx, rho, n)          ! multiply by density
    call cmult(bfy, rho, n)
    if (gdim.eq.3) then
       call add3s2(ta3%x, abz1, abz2, ab1, ab2, n)
       call copy(abz2, abz1, n)
       call copy(abz1, bfz, n)
       call add2s1(bfz, ta3%x, ab0, n)
       call cmult(bfz, rho, n)
    end if
  end subroutine makeabf
  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual 
  subroutine fluid_step_info(step, t, dt, ksp_results)
    type(ksp_monitor_t), intent(in) :: ksp_results(4)
    integer, intent(in) :: step
    real(kind=rp), intent(in) :: t, dt
    character(len=LOG_SIZE) :: log_buf
    integer :: i


    call neko_log%message('Pressure')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(1)%iter, &
         ksp_results(1)%res_start, ksp_results(1)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('X-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(2)%iter, &
         ksp_results(2)%res_start, ksp_results(2)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('Y-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(3)%iter, &
         ksp_results(3)%res_start, ksp_results(3)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('Z-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ', &
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(4)%iter, &
         ksp_results(4)%res_start, ksp_results(4)%res_final
    call neko_log%message(log_buf)

    ! Check for divergence
    do i = 1, 4
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_log%error("Fluid solver diverged")
          stop
       end if
    end do
    
  end subroutine fluid_step_info

  subroutine plan4_compute_vol_flow(this, ab_bdf, niter)

!     Compute pressure and velocity using fractional step method.
!     (Tombo splitting scheme).

    type(fluid_plan4_t), intent(inout) :: this
    type(abbdf_t), intent(inout) :: ab_bdf
    integer, intent(in) :: niter
    integer :: n
    real(kind=rp) :: xlmin, xlmax
    real(kind=rp) :: ylmin, ylmax
    real(kind=rp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_result


    associate(c => this%c_Xh, p_vol => this%p_vol, p_res => this%p_res, &
         u_res => this%u_res, v_res => this%v_res, w_res=>this%w_res, &
         msh => this%msh)
      
      n = this%dm_Xh%size()
      xlmin = glmin(this%dm_Xh%x,n)
      xlmax = glmax(this%dm_Xh%x,n)
      ylmin = glmin(this%dm_Xh%y,n)          !  for Y!
      ylmax = glmax(this%dm_Xh%y,n)
      zlmin = glmin(this%dm_Xh%z,n)          !  for Z!
      zlmax = glmax(this%dm_Xh%z,n)
      if (this%flow_dir.eq.1) this%domain_length = xlmax - xlmin
      if (this%flow_dir.eq.2) this%domain_length = ylmax - ylmin
      if (this%flow_dir.eq.3) this%domain_length = zlmax - zlmin
      
      call fluid_plan4_pres_setup(c%h1, c%h2, this%params%rho, n, c%ifh2)
      !   Compute pressure 

      if (this%flow_dir .eq. 1) then
         call cdtp(p_res, c%h1, c%drdx, c%dsdx, c%dtdx, c)
      end if
      
      if (this%flow_dir .eq. 2) then
         call cdtp(p_res, c%h1, c%drdy, c%dsdy, c%dtdy, c)
      end if
    
      if (this%flow_dir .eq. 3) then
         call cdtp(p_res, c%h1, c%drdz, c%dsdz, c%dtdz, c)
      end if
    
      !call ortho    (respr)

      call gs_op_vector(this%gs_Xh, p_res, n, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_prs, p_res, n)
      call this%pc_prs%update()
      ksp_result = this%ksp_prs%solve(this%Ax, p_vol, p_res, n, c, &
                                      this%bclst_prs, this%gs_Xh, niter)    
      
      !   Compute velocity
      
      call opgrad(u_res, v_res, w_res, p_vol%x, c)
      call opchsign(u_res, v_res, w_res, msh%gdim, n)
      call copy(this%ta1%x, c%B, n)
      call copy(this%ta2%x, c%B, n)
      call copy(this%ta3%x, c%B, n)
      call bc_list_apply_vector(this%bclst_vel,&
                                this%ta1%x, this%ta2%x, this%ta3%x,n)

      if (this%flow_dir.eq.1) then
         call add2(u_res, this%ta1%x,n) ! add forcing
      else if (this%flow_dir.eq.2) then
         call add2(v_res, this%ta2%x,n)
      else if (this%flow_dir.eq.3) then
         call add2(w_res, this%ta3%x,n)
      end if
      

      call fluid_plan4_vel_setup(c%h1, c%h2, &
                                 this%params%Re, this%params%rho,&
                                 ab_bdf%bd(1), &
                                 this%params%dt, n, c%ifh2)
      call gs_op_vector(this%gs_Xh, u_res, n, GS_OP_ADD) 
      call gs_op_vector(this%gs_Xh, v_res, n, GS_OP_ADD) 
      call gs_op_vector(this%gs_Xh, w_res, n, GS_OP_ADD) 
      
      call bc_list_apply_vector(this%bclst_vel,&
                                u_res, v_res, w_res, this%dm_Xh%n_dofs)
      call this%pc_vel%update()

      ksp_result = this%ksp_vel%solve(this%Ax, this%u_vol, u_res, n, &
           c, this%bclst_vel_residual, this%gs_Xh, niter)
      ksp_result = this%ksp_vel%solve(this%Ax, this%v_vol, v_res, n, &
           c, this%bclst_vel_residual, this%gs_Xh, niter)
      ksp_result = this%ksp_vel%solve(this%Ax, this%w_vol, w_res, n, &
           c, this%bclst_vel_residual, this%gs_Xh, niter)
      
      if (this%flow_dir.eq.1) then
         this%base_flow = glsc2(this%u_vol%x, c%B, n) / this%domain_length
      end if
      
      if (this%flow_dir.eq.2) then
         this%base_flow = glsc2(this%v_vol%x, c%B, n) / this%domain_length
      end if
      
      if (this%flow_dir.eq.3) then
         this%base_flow = glsc2(this%w_vol%x, c%B, n) / this%domain_length
      end if    

    end associate
    
  end subroutine  plan4_compute_vol_flow

  subroutine plan4_vol_flow(this, ab_bdf, niter)
!     Adust flow volume at end of time step to keep flow rate fixed by
!     adding an appropriate multiple of the linear solution to the Stokes
!     problem arising from a unit forcing in the X-direction.  This assumes
!     that the flow rate in the X-direction is to be fixed (as opposed to Y-
!     or Z-) *and* that the periodic boundary conditions in the X-direction
!     occur at the extreme left and right ends of the mesh.
!
!     pff 6/28/98
      type(fluid_plan4_t), intent(inout) :: this
      type(abbdf_t), intent(inout) :: ab_bdf
      integer, intent(in) :: niter
      real(kind=rp) :: ifcomp, flow_rate, xsec
      real(kind=rp) :: current_flow, delta_flow, base_flow, scale
      integer :: n, ierr

      n = this%dm_Xh%n_dofs

!     If either dt or the backwards difference coefficient change,
!     then recompute base flow solution corresponding to unit forcing:

      ifcomp = 0.0_rp

      if (this%params%dt .ne. this%dtlag .or. ab_bdf%bd(1) .ne. this%bdlag) then
         ifcomp = 1.0_rp
      end if

      this%dtlag = this%params%dt
      this%bdlag = ab_bdf%bd(1)

      call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      if (ifcomp .gt. 0d0) then
         call plan4_compute_vol_flow(this, ab_bdf, niter)
      end if

      if (this%flow_dir .eq. 1) then
         current_flow=glsc2(this%u%x,this%c_Xh%B,n)/this%domain_length  ! for X
      else if (this%flow_dir .eq. 2) then
         current_flow=glsc2(this%v%x,this%c_Xh%B,n)/this%domain_length  ! for Y
      else if (this%flow_dir .eq. 3) then
         current_flow=glsc2(this%w%x,this%c_Xh%B,n)/this%domain_length  ! for Z
      end if
      
      if (this%avflow) then
         xsec = this%c_Xh%volume / this%domain_length
         flow_rate = this%flow_rate*xsec
      endif  

      delta_flow = flow_rate-current_flow

!     Note, this scale factor corresponds to FFX, provided FFX has
!     not also been specified in userf.   If ffx is also specified
!     in userf then the true FFX is given by ffx_userf + scale.

      scale = delta_flow/this%base_flow

      call add2s2(this%u%x,this%u_vol%x,scale,n)
      call add2s2(this%v%x,this%v_vol%x,scale,n)
      call add2s2(this%w%x,this%w_vol%x,scale,n)
      call add2s2(this%p%x,this%p_vol%x,scale,n)
  end subroutine plan4_vol_flow

end module fluid_plan4
