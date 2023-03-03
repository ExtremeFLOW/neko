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
!> Classic Nek5000 PN/PN formulation for fluids
!! Splitting scheme A.G. Tomboulides et al.
!! Journal of Sci.Comp.,Vol. 12, No. 2, 1998
!! @deprecated The plan4 solver has been retired,
!! bugs will not be fixed, and new features are not added.
module fluid_plan4
  use ax_helm_fctry
  use field_series    
  use fluid_method
  use facet_normal
  use neko_config
  use fluid_aux    
  use ext_bdf_scheme
  use projection
  use logger
  use advection
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan  
  implicit none
  private

  type, public, extends(fluid_scheme_t) :: fluid_plan4_t
     type(field_t) :: u_e, v_e, w_e

     type(field_t) :: p_res, u_res, v_res, w_res

     type(field_series_t) :: ulag, vlag, wlag

     type(field_t) :: dp, du, dv, dw

     type(field_t) :: wa1, wa2, wa3
     type(field_t) :: ta1, ta2, ta3
     
     !> @todo move this to a scratch space
     type(field_t) :: work1, work2

     class(ax_t), allocatable :: Ax
     
     type(projection_t) :: proj

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_residual   !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_residual  

     class(advection_t), allocatable :: adv 

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2

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
    class(fluid_plan4_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), target, intent(inout) :: param     
    character(len=13), parameter :: scheme = 'plan4 (Pn/Pn)'

    call this%free()
    
    ! Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, .true., .true., scheme)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call field_init(this%p_res, dm_Xh, "p_res")
      call field_init(this%u_res, dm_Xh, "u_res")
      call field_init(this%v_res, dm_Xh, "v_res")
      call field_init(this%w_res, dm_Xh, "w_res")            
      call field_init(this%abx1, dm_Xh, "abx1")
      call field_init(this%aby1, dm_Xh, "aby1")
      call field_init(this%abz1, dm_Xh, "abz1")

      call field_init(this%abx2, dm_Xh, "abx2")
      call field_init(this%aby2, dm_Xh, "aby2")
      call field_init(this%abz2, dm_Xh, "abz2")
                  
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

      call this%ulag%init(this%u, 2)
      call this%vlag%init(this%v, 2)
      call this%wlag%init(this%w, 2)
      
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
    call this%proj%init(this%dm_Xh%size(), param%proj_prs_dim)

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
    call advection_factory(this%adv, this%c_Xh, param%dealias, param%lxd)

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

    call field_free(this%p_res)        
    call field_free(this%u_res)
    call field_free(this%v_res)
    call field_free(this%w_res)
    
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

    call field_free(this%abx1)
    call field_free(this%aby1)
    call field_free(this%abz1)

    call field_free(this%abx2)
    call field_free(this%aby2)
    call field_free(this%abz2)
    
    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if
    
    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()
    
  end subroutine fluid_plan4_free
  
  subroutine fluid_plan4_step(this, t, tstep, ext_bdf)
    class(fluid_plan4_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
    integer, intent(inout) :: tstep
    integer :: n, niter
    type(ksp_monitor_t) :: ksp_results(4)
    real(kind=rp), parameter :: one = 1.0
    n = this%dm_Xh%size()
    niter = 1000

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         ta1 => this%ta1, ta2 => this%ta2, ta3 => this%ta3, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         params => this%params, msh => this%msh)

      call fluid_plan4_sumab(u_e%x, u%x, ulag ,n, ext_bdf%ext, ext_bdf%nab)
      call fluid_plan4_sumab(v_e%x, v%x, vlag ,n, ext_bdf%ext, ext_bdf%nab)
      if (msh%gdim .eq. 3) then
         call fluid_plan4_sumab(w_e%x, w%x, wlag,n, ext_bdf%ext, ext_bdf%nab)
      end if

      call f_Xh%eval(t)
      call opcolv(f_Xh%u, f_Xh%v, f_Xh%w, c_Xh%B, msh%gdim, n)
      call this%adv%apply(this%u, this%v, this%w, &
                 f_Xh%u, f_Xh%v, f_Xh%w, &
                 Xh, this%c_Xh, dm_Xh%size())
   
      call makeabf(ta1, ta2, ta3,&
                  this%abx1, this%aby1, this%abz1,&
                  this%abx2, this%aby2, this%abz2, &
                  f_Xh%u, f_Xh%v, f_Xh%w,&
                  params%rho, ext_bdf%ext, n, msh%gdim)
      call makebdf(ta1, ta2, ta3,&
                   this%wa1, this%wa2, this%wa3,&
                   c_Xh%h2, ulag, vlag, wlag, &
                   f_Xh%u, f_Xh%v, f_Xh%w, u, v, w,&
                   c_Xh%B, params%rho, params%dt, &
                   ext_bdf%bdf, ext_bdf%nbd, n, msh%gdim)

      call ulag%update()
      call vlag%update()
      call wlag%update()

      ! mask Dirichlet boundaries (velocity)
      call this%bc_apply_vel()
      
      ! compute pressure
      call this%bc_apply_prs()
      call fluid_plan4_pres_setup(c_Xh%h1, c_Xh%h2, params%rho, &
                                  dm_Xh%size(), c_Xh%ifh2)    
      call fluid_plan4_pres_residual(p, p_res%x, u, v, w, &
                                     u_e, v_e, w_e, &
                                     ta1, ta2, ta3, &
                                     this%wa1, this%wa2, this%wa3, &
                                     this%work1, this%work2, f_Xh, &
                                     c_Xh, gs_Xh, this%bc_prs_surface, &
                                     Ax, ext_bdf%bdf(1), params%dt, &
                                     params%Re, params%rho)

      !Sets tolerances
      !call ctolspl  (tolspl,respr)
      call gs_op(gs_Xh, p_res, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_prs, p_res%x, p%dof%size())

      if( tstep .gt. 5) call this%proj%project_on(p_res%x, c_Xh, n)
      call this%pc_prs%update()
      ksp_results(1) = this%ksp_prs%solve(Ax, dp, p_res%x, n, c_Xh, &
                                this%bclst_prs, gs_Xh, niter)    
      if( tstep .gt. 5) call this%proj%project_back(dp%x, Ax, c_Xh, &
                                  this%bclst_prs, gs_Xh, n)
      call add2(p%x, dp%x,n)
      !    call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)
    
      !We only need to update h2 once I think then use the flag to switch on/off
      call fluid_plan4_vel_setup(c_Xh%h1, c_Xh%h2, &
                                 params%Re, params%rho, ext_bdf%bdf(1), &
                                 params%dt, dm_Xh%size(), c_Xh%ifh2)
    
      call fluid_plan4_vel_residual(Ax, u, v, w, &
                                    u_res, v_res, w_res, &
                                    p, ta1, ta2, ta3, &
                                    f_Xh, c_Xh, msh, Xh, dm_Xh%size())

      call gs_op(gs_Xh, u_res, GS_OP_ADD) 
      call gs_op(gs_Xh, v_res, GS_OP_ADD) 
      call gs_op(gs_Xh, w_res, GS_OP_ADD) 

      call bc_list_apply_vector(this%bclst_vel_residual,&
                                u_res%x, v_res%x, w_res%x, dm_Xh%size())
      call this%pc_vel%update()

      ksp_results(2) = this%ksp_vel%solve(Ax, du, u_res%x, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      ksp_results(3) = this%ksp_vel%solve(Ax, dv, v_res%x, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      ksp_results(4) = this%ksp_vel%solve(Ax, dw, w_res%x, n, &
           c_Xh, this%bclst_vel_residual, gs_Xh, niter)
      
      call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, one, n, msh%gdim)
     
      if (this%flow_dir .ne. 0) then
         call plan4_vol_flow(this, ext_bdf, niter)
      end if
      
      call fluid_step_info(tstep, t, params%dt, ksp_results)
      
    end associate
  end subroutine fluid_plan4_step
  
  subroutine fluid_plan4_pres_setup(h1, h2, rho, n, ifh2)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: h1(n)
    real(kind=rp), intent(inout) :: h2(n)
    real(kind=rp), intent(in) :: rho
    real(kind=rp), parameter :: one = 1.0_rp
    logical, intent(inout) :: ifh2
    call cfill(h1, one /rho, n)
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
    real(kind=rp), parameter :: one = 1.0_rp
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
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(source_t), intent(inout) :: f_Xh
    type(coef_t), intent(inout) :: c_Xh
    integer, intent(in) :: n
    
    call Ax%compute(u_res%x, u%x, c_Xh, msh, Xh)
    call Ax%compute(v_res%x, v%x, c_Xh, msh, Xh)
    if (msh%gdim .eq. 3) then
       call Ax%compute(w_res%x, w%x, c_Xh, msh, Xh)
    end if

    call opchsign(u_res%x, v_res%x, w_res%x, msh%gdim, n)

    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)

    call opadd2cm(u_res%x, v_res%x, w_res%x, &
         ta1%x, ta2%x, ta3%x, -1.0_rp, n, msh%gdim)

    call opadd2cm(u_res%x, v_res%x, w_res%x, &
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
    real(kind=rp), intent(inout) :: p_res(p%dof%size())
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
          p_res(i) = p_res(i) + wa1%x(i,1,1,1) + wa2%x(i,1,1,1)
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
    type(field_series_t), intent(inout) :: vvlag
    real(kind=rp), dimension(3), intent(in) :: ab
    real(kind=rp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    call add3s2(v,vv,vvlag%lf(1)%x,ab0,ab1,n)
    if(nab .eq. 3) call add2s2(v,vvlag%lf(2)%x,ab2,n)
  end subroutine fluid_plan4_sumab
  
  !> Add contributions to F from lagged BD terms.
  subroutine makebdf(ta1, ta2, ta3, tb1, tb2, tb3, h2, ulag, vlag, wlag, &
                     bfx, bfy, bfz, u, v, w, B, rho, dt, bd, nbd, n, gdim)
    integer, intent(in) :: n, nbd, gdim
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: tb1, tb2, tb3
    type(field_series_t), intent(in) :: ulag, vlag, wlag        
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(inout) :: h2(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, rho, bd(10)
    real(kind=rp) :: const
    integer :: ilag
    
    const = rho / dt
    call cfill(h2,const,n)
    call opcolv3c(tb1%x, tb2%x, tb3%x, u%x, v%x, w%x, B, bd(2), n, gdim)
    do ilag = 2, nbd
       call opcolv3c(ta1%x, ta2%x, ta3%x, &
                     ulag%lf(ilag-1)%x, vlag%lf(ilag-1)%x, wlag%lf(ilag-1)%x, &
                     B, bd(ilag+1), n, gdim)
       call opadd2cm(tb1%x, tb2%x, tb3%x, ta1%x, ta2%x, ta3%x, 1.0_rp, n, gdim)
    end do
    call opadd2col(bfx, bfy, bfz, tb1%x, tb2%x, tb3%x, h2, n, gdim)
  end subroutine makebdf

  !> Sum up contributions to kth order extrapolation scheme.
  subroutine makeabf(ta1, ta2, ta3, abx1, aby1, abz1, abx2, aby2, abz2, &
                     bfx, bfy, bfz, rho, ab, n, gdim)
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(inout) :: abx1, aby1, abz1
    type(field_t), intent(inout) :: abx2, aby2, abz2
    real(kind=rp), intent(inout) :: rho, ab(10)
    integer, intent(in) :: n, gdim
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)
    call add3s2(ta1%x, abx1%x, abx2%x, ab1, ab2, n)
    call add3s2(ta2%x, aby1%x, aby2%x, ab1, ab2, n)
    call copy(abx2%x, abx1%x, n)
    call copy(aby2%x, aby1%x, n)
    call copy(abx1%x, bfx, n)
    call copy(aby1%x, bfy, n)
    call add2s1(bfx, ta1%x, ab0, n)
    call add2s1(bfy, ta2%x, ab0, n)
    call cmult(bfx, rho, n)          ! multiply by density
    call cmult(bfy, rho, n)
    if (gdim.eq.3) then
       call add3s2(ta3%x, abz1%x, abz2%x, ab1, ab2, n)
       call copy(abz2%x, abz1%x, n)
       call copy(abz1%x, bfz, n)
       call add2s1(bfz, ta3%x, ab0, n)
       call cmult(bfz, rho, n)
    end if
  end subroutine makeabf
  
  subroutine plan4_compute_vol_flow(this, ext_bdf, niter)

!     Compute pressure and velocity using fractional step method.
!     (Tombo splitting scheme).

    type(fluid_plan4_t), intent(inout) :: this
    type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
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
         call cdtp(p_res%x, c%h1, c%drdx, c%dsdx, c%dtdx, c)
      end if
      
      if (this%flow_dir .eq. 2) then
         call cdtp(p_res%x, c%h1, c%drdy, c%dsdy, c%dtdy, c)
      end if
    
      if (this%flow_dir .eq. 3) then
         call cdtp(p_res%x, c%h1, c%drdz, c%dsdz, c%dtdz, c)
      end if
    
      !call ortho    (respr)

      call gs_op(this%gs_Xh, p_res, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_prs, p_res%x, n)
      call this%pc_prs%update()
      ksp_result = this%ksp_prs%solve(this%Ax, p_vol, p_res%x, n, c, &
                                      this%bclst_prs, this%gs_Xh, niter)    
      
      !   Compute velocity
      
      call opgrad(u_res%x, v_res%x, w_res%x, p_vol%x, c)
      call opchsign(u_res%x, v_res%x, w_res%x, msh%gdim, n)
      call copy(this%ta1%x, c%B, n)
      call copy(this%ta2%x, c%B, n)
      call copy(this%ta3%x, c%B, n)
      call bc_list_apply_vector(this%bclst_vel,&
                                this%ta1%x, this%ta2%x, this%ta3%x,n)

      if (this%flow_dir.eq.1) then
         call add2(u_res%x, this%ta1%x,n) ! add forcing
      else if (this%flow_dir.eq.2) then
         call add2(v_res%x, this%ta2%x,n)
      else if (this%flow_dir.eq.3) then
         call add2(w_res%x, this%ta3%x,n)
      end if
      

      call fluid_plan4_vel_setup(c%h1, c%h2, &
                                 this%params%Re, this%params%rho,&
                                 ext_bdf%bdf(1), &
                                 this%params%dt, n, c%ifh2)
      call gs_op(this%gs_Xh, u_res, GS_OP_ADD) 
      call gs_op(this%gs_Xh, v_res, GS_OP_ADD) 
      call gs_op(this%gs_Xh, w_res, GS_OP_ADD) 
      
      call bc_list_apply_vector(this%bclst_vel,&
                                u_res%x, v_res%x, w_res%x, this%dm_Xh%size())
      call this%pc_vel%update()

      ksp_result = this%ksp_vel%solve(this%Ax, this%u_vol, u_res%x, n, &
           c, this%bclst_vel_residual, this%gs_Xh, niter)
      ksp_result = this%ksp_vel%solve(this%Ax, this%v_vol, v_res%x, n, &
           c, this%bclst_vel_residual, this%gs_Xh, niter)
      ksp_result = this%ksp_vel%solve(this%Ax, this%w_vol, w_res%x, n, &
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

  subroutine plan4_vol_flow(this, ext_bdf, niter)
!     Adust flow volume at end of time step to keep flow rate fixed by
!     adding an appropriate multiple of the linear solution to the Stokes
!     problem arising from a unit forcing in the X-direction.  This assumes
!     that the flow rate in the X-direction is to be fixed (as opposed to Y-
!     or Z-) *and* that the periodic boundary conditions in the X-direction
!     occur at the extreme left and right ends of the mesh.
!
!     pff 6/28/98
      type(fluid_plan4_t), intent(inout) :: this
      type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
      integer, intent(in) :: niter
      real(kind=rp) :: ifcomp, flow_rate, xsec
      real(kind=rp) :: current_flow, delta_flow, base_flow, scale
      integer :: n, ierr

      n = this%dm_Xh%size()

!     If either dt or the backwards difference coefficient change,
!     then recompute base flow solution corresponding to unit forcing:

      ifcomp = 0.0_rp

      if (this%params%dt .ne. this%dtlag .or. ext_bdf%bdf(1) .ne. this%bdlag) then
         ifcomp = 1.0_rp
      end if

      this%dtlag = this%params%dt
      this%bdlag = ext_bdf%bdf(1)

      call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      if (ifcomp .gt. 0d0) then
         call plan4_compute_vol_flow(this, ext_bdf, niter)
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
