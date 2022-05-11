module fluid_pnpn
  use pnpn_res_fctry
  use ax_helm_fctry
  use fluid_abbdf_fctry
  use fluid_method
  use field_series  
  use facet_normal
  use device_math
  use device_mathops
  use fluid_aux    
  use abbdf
  use projection
  use logger
  use advection
  implicit none
  private

  
  type, public, extends(fluid_scheme_t) :: fluid_pnpn_t
     type(field_t) :: u_e, v_e, w_e

     type(field_t) :: p_res, u_res, v_res, w_res

     type(field_series_t) :: ulag, vlag, wlag

     type(field_t) :: dp, du, dv, dw

     type(field_t) :: wa1, wa2, wa3
     type(field_t) :: ta1, ta2, ta3
     
     !> @todo move this to a scratch space
     type(field_t) :: work1, work2

     class(ax_t), allocatable :: Ax
     
     type(projection_t) :: proj_prs
     type(projection_t) :: proj_u
     type(projection_t) :: proj_v
     type(projection_t) :: proj_w

     type(facet_normal_t) :: bc_prs_surface !< Surface term in pressure rhs
     type(facet_normal_t) :: bc_sym_surface !< Surface term in pressure rhs
     type(dirichlet_t) :: bc_vel_residual   !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_du   !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_dv   !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_dw   !< Dirichlet condition vel. res.
     type(dirichlet_t) :: bc_dp   !< Dirichlet condition vel. res.
     type(non_normal_t) :: bc_vel_residual_non_normal   !< Dirichlet condition vel. res.
     type(bc_list_t) :: bclst_vel_residual  
     type(bc_list_t) :: bclst_du
     type(bc_list_t) :: bclst_dv
     type(bc_list_t) :: bclst_dw
     type(bc_list_t) :: bclst_dp  

     class(advection_t), allocatable :: adv 

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2
     
     !> Pressure residual
     class(pnpn_prs_res_t), allocatable :: prs_res

     !> Velocity residual
     class(pnpn_vel_res_t), allocatable :: vel_res

     !> Summation of AB/BDF contributions
     class(fluid_sumab_t), allocatable :: sumab

     !> Contributions to kth order extrapolation scheme
     class(fluid_makeabf_t), allocatable :: makeabf

     !> Contributions to F from lagged BD terms
     class(fluid_makebdf_t), allocatable :: makebdf
     
   contains
     procedure, pass(this) :: init => fluid_pnpn_init
     procedure, pass(this) :: free => fluid_pnpn_free
     procedure, pass(this) :: step => fluid_pnpn_step
  end type fluid_pnpn_t

contains
  
  subroutine fluid_pnpn_init(this, msh, lx, param)    
    class(fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), target, intent(inout) :: param    

    call this%free()
    
    ! Setup velocity and pressure fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, .true., .true.)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax)

    ! Setup backend dependent prs residual routines
    call pnpn_prs_res_factory(this%prs_res)

    ! Setup backend dependent vel residual routines
    call pnpn_vel_res_factory(this%vel_res)

    ! Setup backend dependent summation of AB/BDF
    call fluid_sumab_fctry(this%sumab)

    ! Setup backend dependent summation of extrapolation scheme
    call fluid_makeabf_fctry(this%makeabf)

    ! Setup backend depenent contributions to F from lagged BD terms
    call fluid_makebdf_fctry(this%makebdf)
    
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
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
                        'v', this%params%bc_labels)
    call this%bc_prs_surface%mark_zones_from_list(msh%labeled_zones,&
        'o+dong', this%params%bc_labels)
    call this%bc_prs_surface%finalize()
    call this%bc_prs_surface%set_coef(this%c_Xh)
    ! Initialize symmetry surface terms in pressure rhs
    call this%bc_sym_surface%init(this%dm_Xh)
    call this%bc_sym_surface%mark_zone(msh%sympln)
    call this%bc_sym_surface%mark_zones_from_list(msh%labeled_zones,&
                        'sym', this%params%bc_labels)
    call this%bc_sym_surface%finalize()
    call this%bc_sym_surface%set_coef(this%c_Xh)
    ! Initialize dirichlet bcs for velocity residual
    call this%bc_vel_residual_non_normal%init(this%dm_Xh)
    call this%bc_vel_residual_non_normal%mark_zone(msh%outlet_normal)
    call this%bc_vel_residual_non_normal%mark_zones_from_list(msh%labeled_zones,&
                        'on', this%params%bc_labels)
    call this%bc_vel_residual_non_normal%mark_zones_from_list(msh%labeled_zones,&
                        'on+dong', this%params%bc_labels)
    call this%bc_vel_residual_non_normal%finalize()
    call this%bc_vel_residual_non_normal%init_msk(this%c_Xh)    

    call this%bc_dp%init(this%dm_Xh)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones,&
                        'o', this%params%bc_labels)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones,&
                        'on', this%params%bc_labels)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones,&
                        'on+dong', this%params%bc_labels)
    call this%bc_dp%mark_zones_from_list(msh%labeled_zones,&
                        'o+dong', this%params%bc_labels)
    call this%bc_dp%finalize()
    call this%bc_dp%set_g(0.0_rp)
    call bc_list_init(this%bclst_dp)
    call bc_list_add(this%bclst_dp, this%bc_dp)

    call this%bc_vel_residual%init(this%dm_Xh)
    call this%bc_vel_residual%mark_zone(msh%inlet)
    call this%bc_vel_residual%mark_zone(msh%wall)
    call this%bc_vel_residual%mark_zones_from_list(msh%labeled_zones,&
                        'v', this%params%bc_labels)
    call this%bc_vel_residual%mark_zones_from_list(msh%labeled_zones,&
                        'w', this%params%bc_labels)
    call this%bc_vel_residual%finalize()
    call this%bc_vel_residual%set_g(0.0_rp)
    call bc_list_init(this%bclst_vel_residual)
    call bc_list_add(this%bclst_vel_residual, this%bc_vel_residual)
    call bc_list_add(this%bclst_vel_residual, this%bc_vel_residual_non_normal)
    call bc_list_add(this%bclst_vel_residual, this%bc_sym)

    !Initialize bcs for u, v, w velocity components
    call bc_list_init(this%bclst_du)
    call bc_list_add(this%bclst_du, this%bc_vel_residual)
    call this%bc_du%init(this%dm_Xh)
    if (this%bc_vel_residual_non_normal%xaxis_msk(0) .gt. 0) then
       call this%bc_du%mark_facets(this%bc_vel_residual_non_normal%marked_facet)
    end if
    if (this%bc_sym%xaxis_msk(0) .gt. 0) then
       call this%bc_du%mark_facets(this%bc_sym%marked_facet)
    end if
    call this%bc_du%finalize()
    call this%bc_du%set_g(0.0_rp)
    call bc_list_add(this%bclst_du, this%bc_du)

    call bc_list_init(this%bclst_dv)
    call bc_list_add(this%bclst_dv, this%bc_vel_residual)
    call this%bc_dv%init(this%dm_Xh)
    if (this%bc_vel_residual_non_normal%yaxis_msk(0) .gt. 0) then
       call this%bc_dv%mark_facets(this%bc_vel_residual_non_normal%marked_facet)
    end if
    if (this%bc_sym%yaxis_msk(0) .gt. 0) then
       call this%bc_dv%mark_facets(this%bc_sym%marked_facet)
    end if
    call this%bc_dv%finalize()
    call this%bc_dv%set_g(0.0_rp)
    call bc_list_add(this%bclst_dv, this%bc_dv)

    call bc_list_init(this%bclst_dw)
    call bc_list_add(this%bclst_dw, this%bc_vel_residual)
    call this%bc_dw%init(this%dm_Xh)
    if (this%bc_vel_residual_non_normal%zaxis_msk(0) .gt. 0) then
       call this%bc_dw%mark_facets(this%bc_vel_residual_non_normal%marked_facet)
    end if
    if (this%bc_sym%zaxis_msk(0) .gt. 0) then
       call this%bc_dw%mark_facets(this%bc_sym%marked_facet)
    end if
    call this%bc_dw%finalize()
    call this%bc_dw%set_g(0.0_rp)
    call bc_list_add(this%bclst_dw, this%bc_dw)

    !Intialize projection space thingy
    call this%proj_prs%init(this%dm_Xh%n_dofs, param%proj_prs_dim)
    if (param%proj_vel_dim .gt. 0) then
       call this%proj_u%init(this%dm_Xh%n_dofs, param%proj_vel_dim)
       call this%proj_v%init(this%dm_Xh%n_dofs, param%proj_vel_dim)
       call this%proj_w%init(this%dm_Xh%n_dofs, param%proj_vel_dim)
    end if

    ! Add lagged term to checkpoint
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)    
    call advection_factory(this%adv, this%c_Xh, param%dealias, param%lxd)

  end subroutine fluid_pnpn_init

  subroutine fluid_pnpn_free(this)
    class(fluid_pnpn_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call this%bc_prs_surface%free() 
    call this%bc_sym_surface%free()  
    call bc_list_free(this%bclst_vel_residual)
    call bc_list_free(this%bclst_dp)
    call this%proj_prs%free()
    call this%proj_u%free()
    call this%proj_v%free()
    call this%proj_w%free()
   
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

    if (allocated(this%prs_res)) then
       deallocate(this%prs_res)
    end if
    
    if (allocated(this%vel_res)) then
       deallocate(this%vel_res)
    end if

    if (allocated(this%sumab)) then
       deallocate(this%sumab)
    end if

    if (allocated(this%makeabf)) then
       deallocate(this%makeabf)
    end if

    if (allocated(this%makebdf)) then
       deallocate(this%makebdf)
    end if
    
    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()
    
  end subroutine fluid_pnpn_free

  subroutine fluid_pnpn_step(this, t, tstep, ab_bdf)
    class(fluid_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    type(abbdf_t), intent(inout) :: ab_bdf
    integer, intent(inout) :: tstep
    integer :: n, niter
    type(ksp_monitor_t) :: ksp_results(4)
    n = this%dm_Xh%n_dofs
    niter = 1000

    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         ta1 => this%ta1, ta2 => this%ta2, ta3 => this%ta3, &
         wa1 => this%wa1, wa2 => this%wa2, wa3 => this%wa3, &
         u_res =>this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         params => this%params, msh => this%msh, prs_res => this%prs_res, &
         vel_res => this%vel_res, sumab => this%sumab, &
         makeabf => this%makeabf, makebdf => this%makebdf)
         
      

      call sumab%compute(u_e, v_e, w_e, u, v, w, &
                         ulag, vlag, wlag, ab_bdf%ab, ab_bdf%nab)
     
      call f_Xh%eval()

      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_opcolv(f_Xh%u_d, f_Xh%v_d, f_Xh%w_d, c_Xh%B_d, msh%gdim, n)
      else
         call opcolv(f_Xh%u, f_Xh%v, f_Xh%w, c_Xh%B, msh%gdim, n)
      end if

      call this%adv%apply(this%u, this%v, this%w, &
                          f_Xh%u, f_Xh%v, f_Xh%w, &
                          Xh, this%c_Xh, dm_Xh%n_dofs)
   
      call makeabf%compute(ta1, ta2, ta3,&
                           this%abx1, this%aby1, this%abz1,&
                           this%abx2, this%aby2, this%abz2, &
                           f_Xh%u, f_Xh%v, f_Xh%w,&
                           params%rho, ab_bdf%ab, n)
      
      call makebdf%compute(ta1, ta2, ta3, this%wa1, this%wa2, this%wa3,&
                           ulag, vlag, wlag, f_Xh%u, f_Xh%v, f_Xh%w, &
                           u, v, w, c_Xh%B, params%rho, params%dt, &
                           ab_bdf%bd, ab_bdf%nbd, n)

      call ulag%update()
      call vlag%update()
      call wlag%update()
      !> We assume that no change of boundary conditions 
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply dirichlet
      call this%bc_apply_vel()
      call this%bc_apply_prs()

      ! compute pressure
      call prs_res%compute(p, p_res, u, v, w, u_e, v_e, w_e, &
                           ta1, ta2, ta3, wa1, wa2, wa3, &
                           this%work1, this%work2, f_Xh, &
                           c_Xh, gs_Xh, this%bc_prs_surface, &
                           this%bc_sym_surface, Ax, ab_bdf%bd(1), &
                           params%dt, params%Re, params%rho)

      call gs_op(gs_Xh, p_res, GS_OP_ADD) 
      call bc_list_apply_scalar(this%bclst_dp, p_res%x, p%dof%n_dofs)

      if( tstep .gt. 5 .and. params%proj_prs_dim .gt. 0) call this%proj_prs%project_on(p_res%x, c_Xh, n)
      call this%pc_prs%update()
      ksp_results(1) = this%ksp_prs%solve(Ax, dp, p_res%x, n, c_Xh, &
                                this%bclst_dp, gs_Xh, niter)    
      if( tstep .gt. 5 .and. params%proj_prs_dim .gt. 0) call this%proj_prs%project_back(dp%x, Ax, c_Xh, &
                                  this%bclst_dp, gs_Xh, n)

      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_add2(p%x_d, dp%x_d,n)
      else
         call add2(p%x, dp%x,n)
      end if
      

      ! compute velocity
      call vel_res%compute(Ax, u, v, w, &
                           u_res, v_res, w_res, &
                           p, ta1, ta2, ta3, &
                           f_Xh, c_Xh, msh, Xh, &
                           params%Re, params%rho, ab_bdf%bd(1), &
                           params%dt, dm_Xh%n_dofs)
      
      call gs_op(gs_Xh, u_res, GS_OP_ADD) 
      call gs_op(gs_Xh, v_res, GS_OP_ADD) 
      call gs_op(gs_Xh, w_res, GS_OP_ADD) 

      call bc_list_apply_vector(this%bclst_vel_residual,&
                                u_res%x, v_res%x, w_res%x, dm_Xh%n_dofs)
      
      if (tstep .gt. 5 .and. params%proj_vel_dim .gt. 0) then 
         call this%proj_u%project_on(u_res%x, c_Xh, n)
         call this%proj_v%project_on(v_res%x, c_Xh, n)
         call this%proj_w%project_on(w_res%x, c_Xh, n)
      end if

      call this%pc_vel%update()

      ksp_results(2) = this%ksp_vel%solve(Ax, du, u_res%x, n, &
           c_Xh, this%bclst_du, gs_Xh, niter)
      ksp_results(3) = this%ksp_vel%solve(Ax, dv, v_res%x, n, &
           c_Xh, this%bclst_dv, gs_Xh, niter)
      ksp_results(4) = this%ksp_vel%solve(Ax, dw, w_res%x, n, &
           c_Xh, this%bclst_dw, gs_Xh, niter)

      if (tstep .gt. 5 .and. params%proj_vel_dim .gt. 0) then
         call this%proj_u%project_back(du%x, Ax, c_Xh, &
                                  this%bclst_du, gs_Xh, n)
         call this%proj_v%project_back(dv%x, Ax, c_Xh, &
                                  this%bclst_dv, gs_Xh, n)
         call this%proj_w%project_back(dw%x, Ax, c_Xh, &
                                  this%bclst_dw, gs_Xh, n)
      end if
      
      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if
      
      call fluid_step_info(tstep, t, params%dt, ksp_results)
      
    end associate
  end subroutine fluid_pnpn_step

  
end module fluid_pnpn
