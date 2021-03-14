!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
module fluid_plan4  
  use fluid_method
  use facet_normal
  use neko_config
  use ax_helm_sx
  use ax_helm
  use abbdf
  use projection
  use log
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

     class(ax_t), allocatable :: Ax
     
     type(projection_t) :: proj

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_residual   !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_residual  

     !> Time variables
     real(kind=dp), allocatable :: abx1(:,:,:,:), aby1(:,:,:,:), abz1(:,:,:,:)
     real(kind=dp), allocatable :: abx2(:,:,:,:), aby2(:,:,:,:), abz2(:,:,:,:)

     !> all the shit for vol_flow
     
     integer :: flow_dir !> these two should be moved to params
     logical :: avflow 
     real(kind=dp) :: flow_rate 
     real(kind=dp) :: dtlag = 0d0
     real(kind=dp) :: bdlag = 0d0!> Really quite pointless since we do not vary the timestep
     type(field_t) :: u_vol, v_vol, w_vol, p_vol
     real(kind=dp) :: domain_length, base_flow

     integer :: niter = 1000

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

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(ax_helm_sx_t::this%Ax)
    else
       allocate(ax_helm_t::this%Ax)
    end if
    
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

    ! Initialize boundary condition for velocity residual
    call this%bc_vel_residual%init(this%dm_Xh)
    call this%bc_vel_residual%mark_zone(msh%inlet)
    call this%bc_vel_residual%mark_zone(msh%wall)
    call this%bc_vel_residual%finalize()
    call this%bc_vel_residual%set_g(0d0)
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

  end subroutine fluid_plan4_init

  subroutine fluid_plan4_free(this)
    class(fluid_plan4_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call this%bc_prs_surface%free()  
    call bc_list_free(this%bclst_vel_residual)
    call projection_free(this%proj)
   
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
    real(kind=dp), intent(inout) :: t
    type(abbdf_t), intent(inout) :: ab_bdf
    integer, intent(inout) :: tstep
    integer tt
    integer :: n, iter, i, niter
    type(ksp_monitor_t) :: ksp_results(4)
    n = this%dm_Xh%n_dofs
    niter = 1000

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         ta1 => this%ta1, ta2 => this%ta2, ta3 => this%ta3, ta4 => this%ta4, &
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
      call advab(ta1, ta2, ta3, &
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
                                    p, ta1%x, ta2%x, ta3%x, ta4%x, &
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
      
      call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x,1d0,n,msh%gdim)
     
      if (this%flow_dir .ne. 0) call plan4_vol_flow(this, ab_bdf)
      call fluid_step_info(tstep, t, params%dt, ksp_results)
    end associate
  end subroutine fluid_plan4_step
  
  subroutine fluid_plan4_pres_setup(h1, h2, rho, n, ifh2)
    integer, intent(inout) :: n
    real(kind=dp), intent(inout) :: h1(n)
    real(kind=dp), intent(inout) :: h2(n)
    real(kind=dp), intent(inout) :: rho
    logical, intent(inout) :: ifh2
    call rone(h1, n)
    call cmult(h1, 1d0 /rho, n)
    call rzero(h2, n)
    ifh2 = .false.
  end subroutine fluid_plan4_pres_setup

  subroutine fluid_plan4_vel_setup(h1, h2, Re, rho, bd, dt, n, ifh2)
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
  end subroutine fluid_plan4_vel_setup

  subroutine fluid_plan4_vel_residual(Ax, u, v, w, u_res, v_res, w_res, &
       p, ta1, ta2, ta3, ta4, f_Xh, c_Xh, msh, Xh, n)
    class(ax_t), intent(in) :: Ax
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

  end subroutine fluid_plan4_vel_residual

  subroutine fluid_plan4_pres_residual(p, p_res, u, v, w, u_e, v_e, w_e, &
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
    
    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)
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

  end subroutine fluid_plan4_pres_residual

  !> Sum up AB/BDF contributions 
  subroutine fluid_plan4_sumab(v,vv,vvlag,n,ab,nab)
    integer, intent(in) :: n, nab
    real(kind=dp), dimension(n), intent(inout) :: v,vv
    real(kind=dp), dimension(n,2), intent(inout) :: vvlag
    real(kind=dp), dimension(3), intent(in) :: ab
    real(kind=dp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    call add3s2(v,vv,vvlag(1,1),ab0,ab1,n)
    if(nab .eq. 3) call add2s2(v,vvlag(1,2),ab2,n)
  end subroutine fluid_plan4_sumab
  
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

  subroutine advab(ta1, ta2, ta3, vx, vy, vz, bfx, bfy, bfz, Xh, coef, nelv, n, gdim)
!---------------------------------------------------------------
!
!     Eulerian scheme, add convection term to forcing function 
!     at current time step.
!
!---------------------------------------------------------------
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: ta1, ta2, ta3, vx, vy, vz
    integer, intent(inout) :: nelv, n, gdim
    real(kind=dp), intent(inout), dimension(n) :: bfx, bfy, bfz
    call rzero  (ta1%x,n)
    call rzero  (ta2%x,n)
    call rzero  (ta3%x,n)
    CALL CONV1  (TA1%x,vx%x, vx%x, vy%x, vz%x, Xh, coef, nelv, gdim)
    CALL CONV1  (TA2%x,vy%x, vx%x, vy%x, vz%x, Xh, coef, nelv, gdim)
    CALL SUBCOL3 (BFX,coef%B,TA1%x,N)
    CALL SUBCOL3 (BFY,coef%B,TA2%x,N)
    IF (gdim.EQ.2) THEN
       CALL RZERO (TA3%x,N)
    ELSE
       CALL CONV1  (TA3%x,vz%x, vx%X, vy%x, vz%x, Xh, coef, nelv, gdim)
       CALL SUBCOL3 (BFZ,coef%B,TA3%x,N)
    ENDIF
  END subroutine advab
  !> Prints for prs, velx, vely, velz the following:
  !> Number of iterations, start residual, end residual 
  subroutine fluid_step_info(step, t, dt, ksp_results)
    type(ksp_monitor_t), intent(in) :: ksp_results(4)
    integer, intent(in) :: step
    real(kind=dp), intent(in) :: t, dt
    character(len=LOG_SIZE) :: log_buf


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

  end subroutine fluid_step_info

  subroutine plan4_compute_vol_flow(this, ab_bdf)

!     Compute pressure and velocity using fractional step method.
!     (Tombo splitting scheme).

    type(fluid_plan4_t), intent(inout) :: this
    type(abbdf_t), intent(inout) :: ab_bdf
    integer :: n
    real(kind=dp) :: xlmin, xlmax
    real(kind=dp) :: ylmin, ylmax
    real(kind=dp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_result


    associate(c => this%c_Xh) 
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

    if (this%flow_dir .eq. 1) call cdtp(this%p_res,c%h1,c%drdx,c%dsdx,c%dtdx,c)
    if (this%flow_dir .eq. 2) call cdtp(this%p_res,c%h1,c%drdy,c%dsdy,c%dtdy,c)
    if (this%flow_dir .eq. 3) call cdtp(this%p_res,c%h1,c%drdz,c%dsdz,c%dtdz,c)

    !call ortho    (respr)

    call gs_op_vector(this%gs_Xh, this%p_res, n, GS_OP_ADD) 
    call bc_list_apply_scalar(this%bclst_prs, this%p_res, n)
    call this%pc_prs%update()
    ksp_result = this%ksp_prs%solve(this%Ax, this%p_vol, this%p_res, n, c, &
                              this%bclst_prs, this%gs_Xh, this%niter)    

!   Compute velocity

    call opgrad(this%u_res,this%v_res,this%w_res,this%p_vol%x,c)
    call opchsign(this%u_res,this%v_res,this%w_res, this%msh%gdim, n)
    call copy(this%ta1%x, c%B, n)
    call copy(this%ta2%x, c%B, n)
    call copy(this%ta3%x, c%B, n)
    call bc_list_apply_vector(this%bclst_vel,&
                              this%ta1%x, this%ta2%x, this%ta3%x,n)

    if (this%flow_dir.eq.1) call add2(this%u_res,this%ta1%x,n) ! add forcing
    if (this%flow_dir.eq.2) call add2(this%v_res,this%ta2%x,n)
    if (this%flow_dir.eq.3) call add2(this%w_res,this%ta3%x,n)


    call fluid_plan4_vel_setup(c%h1, c%h2, &
                               this%params%Re, this%params%rho,&
                               ab_bdf%bd(1), &
                               this%params%dt, n, c%ifh2)
    call gs_op_vector(this%gs_Xh, this%u_res, n, GS_OP_ADD) 
    call gs_op_vector(this%gs_Xh, this%v_res, n, GS_OP_ADD) 
    call gs_op_vector(this%gs_Xh, this%w_res, n, GS_OP_ADD) 

    call bc_list_apply_vector(this%bclst_vel,&
                              this%u_res, this%v_res, this%w_res, this%dm_Xh%n_dofs)
    call this%pc_vel%update()

    ksp_result = this%ksp_vel%solve(this%Ax, this%u_vol, this%u_res, n, &
         c, this%bclst_vel_residual, this%gs_Xh, this%niter)
    ksp_result = this%ksp_vel%solve(this%Ax, this%v_vol, this%v_res, n, &
         c, this%bclst_vel_residual, this%gs_Xh, this%niter)
    ksp_result = this%ksp_vel%solve(this%Ax, this%w_vol, this%w_res, n, &
         c, this%bclst_vel_residual, this%gs_Xh, this%niter)
      
    if (this%flow_dir.eq.1) this%base_flow = glsc2(this%u_vol%x,c%B,n)/this%domain_length
    if (this%flow_dir.eq.2) this%base_flow = glsc2(this%v_vol%x,c%B,n)/this%domain_length
    if (this%flow_dir.eq.3) this%base_flow = glsc2(this%w_vol%x,c%B,n)/this%domain_length
  end associate
  end subroutine  plan4_compute_vol_flow

  subroutine plan4_vol_flow(this, ab_bdf)
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
      real(kind=dp) :: ifcomp(1), flow_rate, xsec
      real(kind=dp) :: current_flow, delta_flow, base_flow, scale
      integer :: n

      n = this%dm_Xh%n_dofs

!     If either dt or the backwards difference coefficient change,
!     then recompute base flow solution corresponding to unit forcing:

      ifcomp(1) = 0d0
      if (this%params%dt .ne. this%dtlag .or. ab_bdf%bd(1) .ne. this%bdlag) ifcomp= (/1d0/)
      this%dtlag = this%params%dt
      this%bdlag = ab_bdf%bd(1)

      if (glsum(ifcomp,1) .gt. 0d0) call plan4_compute_vol_flow(this, ab_bdf)

      if (this%flow_dir .eq. 1) current_flow=glsc2(this%u%x,this%c_Xh%B,n)/this%domain_length  ! for X
      if (this%flow_dir .eq. 2) current_flow=glsc2(this%v%x,this%c_Xh%B,n)/this%domain_length  ! for Y
      if (this%flow_dir .eq. 3) current_flow=glsc2(this%w%x,this%c_Xh%B,n)/this%domain_length  ! for Z

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
