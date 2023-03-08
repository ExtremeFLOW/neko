! Copyright (c) 2022-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Modular version of the Classic Nek5000 Pn/Pn formulation for scalars
module scalar_pnpn
  use scalar_residual_fctry
  use ax_helm_fctry
  use rhs_maker_fctry
  use scalar
  use field_series  
  use facet_normal
  use device_math
  use device_mathops
  use scalar_aux    
  use ext_bdf_scheme
  use projection
  use logger
  use advection
  use profiler
  implicit none
  private


  type, public, extends(scalar_scheme_t) :: scalar_pnpn_t
     
     type(field_t) :: s_res

     type(field_series_t) :: slag

     type(field_t) :: ds

     type(field_t) :: wa1
     type(field_t) :: ta1

     
     class(ax_t), allocatable :: Ax

     type(projection_t) :: proj_s

     type(dirichlet_t) :: bc_res   !< Dirichlet condition for scala 
     type(bc_list_t) :: bclst_ds

     class(advection_t), allocatable :: adv 

     ! Time variables
     type(field_t) :: abx1
     type(field_t) :: abx2

     !> Residual
     class(scalar_residual_t), allocatable :: res

     !> Contributions to kth order extrapolation scheme
     class(rhs_maker_ext_t), allocatable :: makeext

     !> Contributions to F from lagged BD terms
     class(rhs_maker_bdf_t), allocatable :: makebdf

   contains
     procedure, pass(this) :: init => scalar_pnpn_init
     procedure, pass(this) :: free => scalar_pnpn_free
     procedure, pass(this) :: step => scalar_pnpn_step
  end type scalar_pnpn_t

contains

  subroutine scalar_pnpn_init(this, msh, coef, gs, param)    
    class(scalar_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    type(coef_t), target, intent(inout) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(param_t), target, intent(inout) :: param
    integer :: i
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'

    call this%free()

    ! Setup fields on the space \f$ Xh \f$
    call this%scheme_init(msh, coef, gs, param, scheme)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax)

    ! Setup backend dependent scalar residual routines
    call scalar_residual_factory(this%res)

    ! Setup backend dependent summation of extrapolation scheme
    call rhs_maker_ext_fctry(this%makeext)

    ! Setup backend depenent contributions to F from lagged BD terms
    call rhs_maker_bdf_fctry(this%makebdf)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call field_init(this%s_res, dm_Xh, "s_res")

      call field_init(this%abx1, dm_Xh, "abx1")

      call field_init(this%abx2, dm_Xh, "abx2")

      call field_init(this%wa1, dm_Xh, 'wa1')

      call field_init(this%ta1, dm_Xh, 'ta1')

      call field_init(this%ds, dm_Xh, 'ds')

      call this%slag%init(this%s, 2)

    end associate

    ! Initialize dirichlet bcs for scalar residual
    ! todo: look that this works
    call this%bc_res%init(this%dm_Xh)
    do i = 1, this%n_dir_bcs
       call this%bc_res%mark_facets(this%dir_bcs(i)%marked_facet)
    end do
    call this%bc_res%finalize()
    call this%bc_res%set_g(0.0_rp)
    call bc_list_init(this%bclst_ds)
    call bc_list_add(this%bclst_ds, this%bc_res)

    ! @todo not param stuff again, using velocity stuff
    ! Intialize projection space thingy
    if (param%proj_vel_dim .gt. 0) then
       call this%proj_s%init(this%dm_Xh%size(), param%proj_vel_dim)
    end if

    ! Add lagged term to checkpoint
    ! @todo Init chkp object, note, adding 3 slags
    ! call this%chkp%add_lag(this%slag, this%slag, this%slag)    
    
    ! @todo add dealiasing here, now hardcoded to false
    call advection_factory(this%adv, this%c_Xh, .false., param%lxd)

  end subroutine scalar_pnpn_init

  subroutine scalar_pnpn_free(this)
    class(scalar_pnpn_t), intent(inout) :: this

    !Deallocate scalar field
    call this%scheme_free()

    call bc_list_free(this%bclst_ds)
    call this%proj_s%free()

    call field_free(this%s_res)        

    call field_free(this%wa1)

    call field_free(this%ta1)

    call field_free(this%ds)

    call field_free(this%abx1)
    call field_free(this%abx2)

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    if (allocated(this%res)) then
       deallocate(this%res)
    end if

    if (allocated(this%makeext)) then
       deallocate(this%makeext)
    end if

    if (allocated(this%makebdf)) then
       deallocate(this%makebdf)
    end if


    call this%slag%free()

  end subroutine scalar_pnpn_free

  subroutine scalar_pnpn_step(this, t, tstep, ext_bdf)
    class(scalar_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
    integer, intent(inout) :: tstep
    integer :: n, niter
    type(ksp_monitor_t) :: ksp_results(1)
    n = this%dm_Xh%size()
    niter = 1000
    
    call profiler_start_region('Scalar')
    associate(u => this%u, v => this%v, w => this%w, s => this%s, &
         ds => this%ds, &
         ta1 => this%ta1, &
         wa1 => this%wa1, &
         s_res =>this%s_res, &
         Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         slag => this%slag, &
         params => this%params, msh => this%msh, res => this%res, &
         makeext => this%makeext, makebdf => this%makebdf)

      ! evaluate the source term and scale with the mass matrix
      call f_Xh%eval(t)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(f_Xh%s_d, c_Xh%B_d, n)
      else
         call col2(f_Xh%s, c_Xh%B, n)
      end if

      call this%adv%apply_scalar(u, v, w, s, f_Xh%s, &
                                 Xh, this%c_Xh, dm_Xh%size())

      call makeext%compute_scalar(ta1, this%abx1, this%abx2, f_Xh%s, &
           params%rho, ext_bdf%ext, n)

      call makebdf%compute_scalar(ta1, wa1, slag, f_Xh%s, s, c_Xh%B, &
           params%rho, params%dt, ext_bdf%bdf, ext_bdf%nbd, n)

      call slag%update()
      !> We assume that no change of boundary conditions 
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply dirichlet
      call this%bc_apply()

      ! compute scalar residual
      call profiler_start_region('Scalar residual')
      call res%compute(Ax, s,  s_res, f_Xh, c_Xh, msh, Xh, params%Pr, &
          params%Re, params%rho, ext_bdf%bdf(1), params%dt, &
          dm_Xh%size())

      call gs_op(gs_Xh, s_res, GS_OP_ADD) 

      call bc_list_apply_scalar(this%bclst_ds,&
           s_res%x, dm_Xh%size())
      call profiler_end_region

      if (tstep .gt. 5 .and. params%proj_vel_dim .gt. 0) then 
         call this%proj_s%project_on(s_res%x, c_Xh, n)
      end if

      call this%pc%update()
      call profiler_start_region('Scalar solve')
      ksp_results(1) = this%ksp%solve(Ax, ds, s_res%x, n, &
           c_Xh, this%bclst_ds, gs_Xh, niter)
      call profiler_end_region

      if (tstep .gt. 5 .and. params%proj_vel_dim .gt. 0) then
         call this%proj_s%project_back(ds%x, Ax, c_Xh, &
              this%bclst_ds, gs_Xh, n)
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(s%x_d, ds%x_d, 1.0_rp, n)
      else
         call add2s2(s%x, ds%x, 1.0_rp, n)
      end if

      call scalar_step_info(tstep, t, params%dt, ksp_results)

    end associate
    call profiler_end_region
  end subroutine scalar_pnpn_step


end module scalar_pnpn
