!> Accelerator version of the PN/PN formulation for fluids
module device_fluid_plan4
  use ax_helm_fctry
  use fluid_plan4
  use fluid_method
  use field_series
  use facet_normal
  use fluid_aux
  use ext_bdf_scheme
  use projection
  use logger
  use advection
  use device_math
  use device_mathops
  implicit none
  private

  type, public, extends(fluid_plan4_t) :: device_fluid_plan4_t
   contains
     procedure, pass(this) :: step => device_fluid_plan4_step
  end type device_fluid_plan4_t

contains

  subroutine device_fluid_plan4_step(this, t, tstep, ext_bdf)
    class(device_fluid_plan4_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
    integer, intent(inout) :: tstep
    integer :: n, niter
    type(ksp_monitor_t) :: ksp_results(4)
    real(kind=rp), parameter :: one = 1.0
    n = this%dm_Xh%size()
    niter = 3000

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         ta1 => this%ta1, ta2 => this%ta2, ta3 => this%ta3, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         params => this%params, msh => this%msh)

      call device_fluid_plan4_sumab(u_e%x, u%x, ulag ,n, ext_bdf%ext, ext_bdf%nab)
      call device_fluid_plan4_sumab(v_e%x, v%x, vlag ,n, ext_bdf%ext, ext_bdf%nab)
      if (msh%gdim .eq. 3) then
         call device_fluid_plan4_sumab(w_e%x, w%x, wlag,n, ext_bdf%ext, ext_bdf%nab)
      end if

      call f_Xh%eval(t)

      call device_opcolv(f_Xh%u_d, f_Xh%v_d, f_Xh%w_d, c_Xh%B_d, msh%gdim, n)


      call this%adv%apply(this%u, this%v, this%w, &
                 f_Xh%u, f_Xh%v, f_Xh%w, &
                 Xh, this%c_Xh, dm_Xh%size())

      call device_makeabf(ta1, ta2, ta3,&
                          this%abx1, this%aby1, this%abz1,&
                          this%abx2, this%aby2, this%abz2, &
                          f_Xh%u_d, f_Xh%v_d, f_Xh%w_d,&
                          params%rho, ext_bdf%ext, n, msh%gdim)
      call device_makebdf(ta1, ta2, ta3,&
                          this%wa1, this%wa2, this%wa3,&
                          c_Xh%h2_d, ulag, vlag, wlag, &
                          f_Xh%u_d, f_Xh%v_d, f_Xh%w_d, u, v, w,&
                          c_Xh%B_d, params%rho, params%dt, &
                          ext_bdf%bdf, ext_bdf%nbd, n, msh%gdim)

      call ulag%update()
      call vlag%update()
      call wlag%update()


      ! mask Dirichlet boundaries (velocity)
      call this%bc_apply_vel()

      ! compute pressure
      call this%bc_apply_prs()
      call device_fluid_plan4_pres_setup(c_Xh%h1_d, c_Xh%h2_d, params%rho, &
                                         dm_Xh%size(), c_Xh%ifh2)    
      call device_fluid_plan4_pres_residual(p, p_res, u, v, w, &
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
      call device_add2(p%x_d, dp%x_d, n)

      !    call ortho(this%p%x,n,this%Xh%lxyz*this%msh%glb_nelv)

      !We only need to update h2 once I think then use the flag to switch on/off
      call device_fluid_plan4_vel_setup(c_Xh%h1_d, c_Xh%h2_d, &
                                        params%Re, params%rho, ext_bdf%bdf(1), &
                                        params%dt, dm_Xh%size(), c_Xh%ifh2)
    
      call device_fluid_plan4_vel_residual(Ax, u, v, w, &
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

      call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
                           du%x_d, dv%x_d, dw%x_d, one, n, msh%gdim)

     
      call fluid_step_info(tstep, t, params%dt, ksp_results)
     
    end associate
  end subroutine device_fluid_plan4_step
  
  subroutine device_fluid_plan4_pres_setup(h1_d, h2_d, rho, n, ifh2)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: h1_d
    type(c_ptr), intent(inout) :: h2_d
    real(kind=rp), intent(in) :: rho
    real(kind=rp), parameter :: one = 1.0_rp
    logical, intent(inout) :: ifh2
    call device_cfill(h1_d, one/rho, n)
    call device_rzero(h2_d, n)
    ifh2 = .false.
  end subroutine device_fluid_plan4_pres_setup

  subroutine device_fluid_plan4_vel_setup(h1_d, h2_d, Re, rho, bd, dt, n, ifh2)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: h1_d
    type(c_ptr), intent(inout) :: h2_d
    real(kind=rp), intent(in) :: Re
    real(kind=rp), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    logical, intent(inout) :: ifh2
    real(kind=rp), parameter :: one = 1.0_rp
    real(kind=rp) :: dtbd
    dtbd = rho * (bd / dt)
    call device_cfill(h1_d, (one / Re), n)
    call device_cfill(h2_d, dtbd, n)
    ifh2 = .true.
  end subroutine device_fluid_plan4_vel_setup

  subroutine device_fluid_plan4_vel_residual(Ax, u, v, w, u_res, v_res, w_res, &
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

    call device_opchsign(u_res%x_d, v_res%x_d, w_res%x_d, msh%gdim, n)
    call opgrad(ta1%x, ta2%x, ta3%x, p%x, c_Xh)
    
    call device_opadd2cm(u_res%x_d, v_res%x_d, w_res%x_d, &
                         ta1%x_d, ta2%x_d, ta3%x_d, -1.0_rp, n, msh%gdim)

    call device_opadd2cm(u_res%x_d, v_res%x_d, w_res%x_d, &
                         f_Xh%u_d, f_Xh%v_d, f_Xh%w_d, 1.0_rp, n, msh%gdim)

  end subroutine device_fluid_plan4_vel_residual

  subroutine device_fluid_plan4_pres_residual(p, p_res, u, v, w, u_e, v_e, w_e, &
       ta1, ta2, ta3, wa1, wa2, wa3, work1, work2, f_Xh, c_xh, gs_Xh, &
       bc_prs_surface, Ax, bd, dt, Re, rho)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(inout) :: wa1, wa2, wa3
    type(field_t), intent(inout) :: p_res, work1, work2
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
 
    n = c_Xh%dof%size()
    gdim = c_Xh%msh%gdim

    call curl(ta1, ta2, ta3, u_e, v_e, w_e, work1, work2, c_Xh)
    call curl(wa1, wa2, wa3, ta1, ta2, ta3, work1, work2, c_Xh)

    call device_opcolv(wa1%x_d, wa2%x_d, wa3%x_d, c_Xh%B_d, gdim, n)
    
    call device_cfill(work1%x_d, (1.0_rp / Re) / rho, n)
       
    call device_opcolv(wa1%x_d, wa2%x_d, wa3%x_d, work1%x_d, gdim, n)

    call Ax%compute(p_res%x,p%x,c_Xh,p%msh,p%Xh)

    call device_cmult(p_res%x_d, -1.0_rp, n)

    call device_add3s2(ta1%x_d, f_Xh%u_d, wa1%x_d, 1.0_rp / rho, -1.0_rp, n)
    call device_add3s2(ta2%x_d, f_Xh%v_d, wa2%x_d, 1.0_rp / rho, -1.0_rp, n)
    call device_add3s2(ta3%x_d, f_Xh%w_d, wa3%x_d, 1.0_rp / rho, -1.0_rp, n)
                 
     !Need to consider cyclic bcs here...
    call gs_op(gs_Xh, ta1, GS_OP_ADD) 
    call gs_op(gs_Xh, ta2, GS_OP_ADD) 
    call gs_op(gs_Xh, ta3, GS_OP_ADD) 

    call device_opcolv(ta1%x_d, ta2%x_d, ta3%x_D, c_Xh%Binv_d, gdim, n)

    if (gdim .eq. 3) then
       call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
       call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

       call device_add2(wa1%x_d,wa2%x_d,n)
       call device_add2(wa1%x_d,wa3%x_d,n)
       call device_add2(p_res%x_d,wa1%x_d,n)
    else
       call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)

       call device_add2(wa1%x_d,wa2%x_d,n)
       call device_add2(p_res%x_d,wa1%x_d,n)
    endif



    !
    ! Surface velocity terms
    !
    dtbd = bd / dt
    call device_rzero(ta1%x_d, n)
    call device_rzero(ta2%x_d, n)
    call device_rzero(ta3%x_d, n)
    call bc_prs_surface%apply_surfvec_dev(ta1%x_d, ta2%x_d, ta3%x_d, u%x_d, v%x_d, w%x_d)
    call device_add2(ta1%x_d, ta2%x_d, n)
    call device_add2(ta1%x_d, ta3%x_d, n)    
    call device_cmult(ta1%x_d, dtbd, n)
    call device_sub2(p_res%x_d, ta1%x_d, n)

        
!    call ortho(p_res,n,glb_n) ! Orthogonalize wrt null space, if present

  end subroutine device_fluid_plan4_pres_residual

  !> Sum up AB/BDF contributions 
  subroutine device_fluid_plan4_sumab(v,vv,vvlag,n,ab,nab)
    integer, intent(in) :: n, nab
    real(kind=rp), dimension(n), intent(inout) :: v, vv
    type(field_series_t), intent(inout) :: vvlag
    real(kind=rp), dimension(3), intent(in) :: ab
    real(kind=rp) :: ab0, ab1, ab2
    type(c_ptr) :: v_d, vv_d

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)

    v_d = device_get_ptr(v)
    vv_d = device_get_ptr(vv)
       
    call device_add3s2(v_d, vv_d, vvlag%lf(1)%x_d, ab0, ab1, n)
    if(nab .eq. 3) call device_add2s2(v_d, vvlag%lf(2)%x_d, ab2, n)

  end subroutine device_fluid_plan4_sumab
  
  !> Add contributions to F from lagged BD terms.
  subroutine device_makebdf(ta1, ta2, ta3, tb1, tb2, tb3, &
                            h2_d, ulag, vlag, wlag, bfx_d, bfy_d, bfz_d, &
                            u, v, w, b_d, rho, dt, bd, nbd, n, gdim)
    integer, intent(in) :: n, nbd, gdim
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: tb1, tb2, tb3
    type(field_series_t), intent(in) :: ulag, vlag, wlag
    type(c_ptr) :: bfx_d, bfy_d, bfz_d
    type(c_ptr) :: h2_d
    type(c_ptr) :: B_d
    real(kind=rp), intent(in) :: dt, rho, bd(10)
    real(kind=rp) :: const
    integer :: ilag
    
    const = rho / dt
    call device_cfill(h2_d, const, n)
    call device_opcolv3c(tb1%x_d, tb2%x_d, tb3%x_d, &
                         u%x_d, v%x_d, w%x_d, B_d, bd(2), n, gdim)
    do ilag = 2, nbd
       call device_opcolv3c(ta1%x_d, ta2%x_d, ta3%x_d, &
                            ulag%lf(ilag-1)%x_d, vlag%lf(ilag-1)%x_d, &
                            wlag%lf(ilag-1)%x_d, B_d, bd(ilag+1), n, gdim)
       call device_opadd2cm(tb1%x_d, tb2%x_d, tb3%x_d, &
                            ta1%x_d, ta2%x_D, ta3%x_d, 1.0_rp, n, gdim)
    end do
    call device_opadd2col(bfx_d, bfy_d, bfz_d, &
                          tb1%x_d, tb2%x_d, tb3%x_d, h2_d, n, gdim)
  end subroutine device_makebdf

  !> Sum up contributions to kth order extrapolation scheme.
  subroutine device_makeabf(ta1, ta2, ta3, abx1, aby1, abz1, abx2, aby2, abz2, &
                     bfx_d, bfy_d, bfz_d, rho, ab, n, gdim)
    type(field_t), intent(inout) :: ta1, ta2, ta3
    type(field_t), intent(inout) :: abx1, aby1, abz1
    type(field_t), intent(inout) :: abx2, aby2, abz2
    real(kind=rp), intent(inout) :: rho, ab(10)
    integer, intent(in) :: n, gdim
    type(c_ptr) :: bfx_d, bfy_d, bfz_d
    real(kind=rp) :: ab0, ab1, ab2

    ab0 = ab(1)
    ab1 = ab(2)
    ab2 = ab(3)
    call device_add3s2(ta1%x_d, abx1%x_d, abx2%x_d, ab1, ab2, n)
    call device_add3s2(ta2%x_d, aby1%x_d, aby2%x_d, ab1, ab2, n)
    call device_copy(abx2%x_d, abx1%x_d, n)
    call device_copy(aby2%x_d, aby1%x_d, n)
    call device_copy(abx1%x_d, bfx_d, n)
    call device_copy(aby1%x_d, bfy_d, n)
    call device_add2s1(bfx_d, ta1%x_d, ab0, n)
    call device_add2s1(bfy_d, ta2%x_d, ab0, n)
    call device_cmult(bfx_d, rho, n)          ! multiply by density
    call device_cmult(bfy_d, rho, n)
    if (gdim.eq.3) then
       call device_add3s2(ta3%x_d, abz1%x_d, abz2%x_d, ab1, ab2, n)
       call device_copy(abz2%x_d, abz1%x_d, n)
       call device_copy(abz1%x_d, bfz_d, n)
       call device_add2s1(bfz_d, ta3%x_d, ab0, n)
       call device_cmult(bfz_d, rho, n)
    end if
  end subroutine device_makeabf
  
end module device_fluid_plan4
