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
  use scalar_residual_fctry, only : scalar_residual_factory
  use ax_helm_fctry, only: ax_helm_factory
  use rhs_maker_fctry
  use scalar_scheme, only : scalar_scheme_t
  use dirichlet, only : dirichlet_t
  use field, only : field_t
  use bc, only : bc_list_t, bc_list_init, bc_list_free, bc_list_apply_scalar, &
                 bc_list_add
  use mesh, only : mesh_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_residual, only :scalar_residual_t
  use ax_product, only : ax_t
  use field_series  
  use facet_normal
  use device_math
  use device_mathops
  use scalar_aux    
  use time_scheme_controller
  use projection
  use math
  use logger
  use advection
  use profiler
  use json_utils, only: json_get, json_get_or_default
  use json_module, only : json_file
  use user_intf, only : user_t
  use material_properties, only : material_properties_t
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
     !> Constructor.
     procedure, pass(this) :: init => scalar_pnpn_init
     !> Destructor.
     procedure, pass(this) :: free => scalar_pnpn_free
     procedure, pass(this) :: step => scalar_pnpn_step
  end type scalar_pnpn_t

contains

  !> Constructor.
  !! @param msh The mesh.
  !! @param coef The coefficients.
  !! @param gs The gather-scatter.
  !! @param params The case parameter file in json.
  !! @param user Type with user-defined procedures.
  subroutine scalar_pnpn_init(this, msh, coef, gs, params, user, &
                              material_properties)
    class(scalar_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    type(coef_t), target, intent(inout) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(material_properties_t), intent(inout) :: material_properties
    integer :: i
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    ! Variables for retrieving json parameters
    logical :: found, logical_val
    integer :: integer_val

    call this%free()

    ! Initiliaze base type.
    call this%scheme_init(msh, coef, gs, params, scheme, user, &
                          material_properties)

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

      call this%s_res%init(dm_Xh, "s_res")

      call this%abx1%init(dm_Xh, "abx1")

      call this%abx2%init(dm_Xh, "abx2")

      call this%wa1%init(dm_Xh, 'wa1')

      call this%ta1%init(dm_Xh, 'ta1')

      call this%ds%init(dm_Xh, 'ds')

      call this%slag%init(this%s, 2)

    end associate

    ! Initialize dirichlet bcs for scalar residual
    ! todo: look that this works
    call this%bc_res%init(this%dm_Xh)
    do i = 1, this%n_dir_bcs
       call this%bc_res%mark_facets(this%dir_bcs(i)%marked_facet)
    end do

    ! Check for user bcs
    if (this%user_bc%msk(0) .gt. 0) then
       call this%bc_res%mark_facets(this%user_bc%marked_facet)
    end if
    call this%bc_res%finalize()
    call this%bc_res%set_g(0.0_rp)
    call bc_list_init(this%bclst_ds)
    call bc_list_add(this%bclst_ds, this%bc_res)

    ! @todo not param stuff again, using velocity stuff
    ! Intialize projection space thingy
    if (this%projection_dim .gt. 0) then
       call this%proj_s%init(this%dm_Xh%size(), this%projection_dim)
    end if

    ! Add lagged term to checkpoint
    ! @todo Init chkp object, note, adding 3 slags
    ! call this%chkp%add_lag(this%slag, this%slag, this%slag)    
    
    ! Uses sthe same parameter as the fluid to set dealiasing
    call json_get(params, 'case.numerics.dealias', logical_val)
    call params%get('case.numerics.dealiased_polynomial_order', integer_val, &
                    found)
    if (.not. found) then
       call json_get(params, 'case.numerics.polynomial_order', integer_val)
       integer_val =  3.0_rp / 2.0_rp * (integer_val + 1) - 1
    end if
    call advection_factory(this%adv, this%c_Xh, logical_val, integer_val + 1)

  end subroutine scalar_pnpn_init

  subroutine scalar_pnpn_free(this)
    class(scalar_pnpn_t), intent(inout) :: this

    !Deallocate scalar field
    call this%scheme_free()

    call bc_list_free(this%bclst_ds)
    call this%proj_s%free()

    call this%s_res%free()

    call this%wa1%free()

    call this%ta1%free()

    call this%ds%free()

    call this%abx1%free()
    call this%abx2%free()

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

  subroutine scalar_pnpn_step(this, t, tstep, dt, ext_bdf)
    class(scalar_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    ! Number of degrees of freedom
    integer :: n
    ! Linear solver results monitor
    type(ksp_monitor_t) :: ksp_results(1)

    n = this%dm_Xh%size()
    
    call profiler_start_region('Scalar')
    associate(u => this%u, v => this%v, w => this%w, s => this%s, &
         cp => this%cp, lambda => this%lambda, rho => this%rho, &
         ds => this%ds, &
         ta1 => this%ta1, &
         wa1 => this%wa1, &
         s_res =>this%s_res, &
         Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         slag => this%slag, &
         projection_dim => this%projection_dim, &
         ksp_maxiter => this%ksp_maxiter, &
         msh => this%msh, res => this%res, &
         makeext => this%makeext, makebdf => this%makebdf)

      ! Evaluate the source term and scale with the mass matrix.
      call f_Xh%eval(t)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(f_Xh%s_d, c_Xh%B_d, n)
      else
         call col2(f_Xh%s, c_Xh%B, n)
      end if

      ! Add the advection operators to the right-hans-side.
      call this%adv%compute_scalar(u, v, w, s, f_Xh%s, &
                                   Xh, this%c_Xh, dm_Xh%size())

      call makeext%compute_scalar(ta1, this%abx1, this%abx2, f_Xh%s, &
           rho, ext_bdf%advection_coeffs, n)

      call makebdf%compute_scalar(ta1, wa1, slag, f_Xh%s, s, c_Xh%B, &
           rho, dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)

      call slag%update()
      !> We assume that no change of boundary conditions 
      !! occurs between elements. I.e. we do not apply gsop here like in Nek5000
      !> Apply dirichlet
      call this%bc_apply()

      ! Compute scalar residual.
      call profiler_start_region('Scalar residual')
      call res%compute(Ax, s,  s_res, f_Xh, c_Xh, msh, Xh, lambda, rho * cp, &
          ext_bdf%diffusion_coeffs(1), dt, &
          dm_Xh%size())

      call gs_Xh%op(s_res, GS_OP_ADD) 

      call bc_list_apply_scalar(this%bclst_ds,&
           s_res%x, dm_Xh%size())
      call profiler_end_region

      if (tstep .gt. 5 .and. projection_dim .gt. 0) then 
         call this%proj_s%project_on(s_res%x, c_Xh, n)
      end if

      call this%pc%update()
      call profiler_start_region('Scalar solve')
      ksp_results(1) = this%ksp%solve(Ax, ds, s_res%x, n, &
           c_Xh, this%bclst_ds, gs_Xh, ksp_maxiter)
      call profiler_end_region

      if (tstep .gt. 5 .and. projection_dim .gt. 0) then 
         call this%proj_s%project_back(ds%x, Ax, c_Xh, &
              this%bclst_ds, gs_Xh, n)
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(s%x_d, ds%x_d, 1.0_rp, n)
      else
         call add2s2(s%x, ds%x, 1.0_rp, n)
      end if

      call scalar_step_info(tstep, t, dt, ksp_results)

    end associate
    call profiler_end_region
  end subroutine scalar_pnpn_step


end module scalar_pnpn
