! Copyright (c) 2022, The Neko Authors
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

! todo: module name
module scalar
  use gather_scatter
  use mean_sqr_flow    
  use neko_config
  use parameters
  use checkpoint
  use mean_flow
  use num_types
  use source
  use source_scalar
  use field
  use space
  use dofmap
  use krylov
  use coefs
  use wall
  use inflow
  use usr_inflow
  use blasius
  use dirichlet
  use dong_outflow
  use symmetry
  use non_normal
  use krylov_fctry
  use precon_fctry
  use bc
  use mesh
  use math
  use ext_bdf_scheme
  use mathops
  use operators
  use logger
  use field_registry

  type, abstract :: scalar_scheme_t
     type(field_t), pointer :: u         !< x-component of Velocity
     type(field_t), pointer :: v         !< y-component of Velocity
     type(field_t), pointer :: w         !< z-component of Velocity
     type(field_t), pointer :: s         !< the scalar
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$
     type(source_scalar_t) :: f_Xh     !< Source term associated with \f$ X_h \f$
     class(ksp_t), allocatable  :: ksp         !< Krylov solver
     class(pc_t), allocatable :: pc            !< Preconditioner
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     type(dirichlet_t) :: bc_inflow           !< Dirichlet condition for scalar
     type(bc_list_t) :: bclst                  !< List of boundary conditions
     type(param_t), pointer :: params          !< Parameters          
     type(mesh_t), pointer :: msh => null()    !< Mesh
     type(chkp_t) :: chkp                      !< Checkpoint
   contains
     procedure, pass(this) :: scalar_scheme_init
     procedure, pass(this) :: scheme_free => scalar_scheme_free
     procedure, pass(this) :: validate => scalar_scheme_validate
     procedure, pass(this) :: bc_apply => scalar_scheme_bc_apply
     procedure, pass(this) :: set_source => scalar_scheme_set_source
     procedure(scalar_method_init), pass(this), deferred :: init
     procedure(scalar_method_free), pass(this), deferred :: free
     procedure(scalar_method_step), pass(this), deferred :: step
     generic :: scheme_init => scalar_scheme_init
  end type scalar_scheme_t

  !> Abstract interface to initialize a scalar formulation
  abstract interface
     subroutine scalar_method_init(this, msh, lx, param)
       import scalar_scheme_t
       import param_t
       import mesh_t
       class(scalar_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(inout) :: msh       
       integer, intent(inout) :: lx
       type(param_t), target, intent(inout) :: param              
     end subroutine scalar_method_init
  end interface

  !> Abstract interface to dealocate a scalar formulation
  abstract interface
     subroutine scalar_method_free(this)
       import scalar_scheme_t
       class(scalar_scheme_t), intent(inout) :: this
     end subroutine scalar_method_free
  end interface
  
  !> Abstract interface to compute a time-step
  abstract interface
     subroutine scalar_method_step(this, t, tstep, ext_bdf)
       import scalar_scheme_t
       import ext_bdf_scheme_t
       import rp
       class(scalar_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
     end subroutine scalar_method_step
  end interface

contains
    
  !> Initialize common data for the current scheme
  subroutine scalar_scheme_init_common(this, msh, lx, params, scheme)
    class(scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    character(len=*), intent(in) :: scheme
    type(param_t), target, intent(inout) :: params
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    
    call neko_log%section('Scalar')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'lx         : ', lx
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'lx         : ', lx
    else
       write(log_buf, '(A, I3)') 'lx         : ', lx
    end if

    ! can reuse the same as the velocity
    if (msh%gdim .eq. 2) then
       call space_init(this%Xh, GLL, lx, lx)
    else
       call space_init(this%Xh, GLL, lx, lx, lx)
    end if

    this%dm_Xh = dofmap_t(msh, this%Xh)

    this%params => params

    this%msh => msh

    call gs_init(this%gs_Xh, this%dm_Xh)

    call coef_init(this%c_Xh, this%gs_Xh)

    call source_scalar_init(this%f_Xh, this%dm_Xh)

    !
    ! Setup scalar boundary conditions
    !
    call bc_list_init(this%bclst)

    !if (trim(params%scalar_inflow) .eq. "default") then
    !   allocate(inflow_t::this%bc_inflow)
    !else
    !   call neko_error('Invalid Inflow condition')
    !end if
    
    call this%bc_inflow%init(this%dm_Xh)
    ! todo: add scalar zones to mesh
    call this%bc_inflow%mark_zone(msh%inlet)
    ! todo: add scalar labeled zones to mesh
    call this%bc_inflow%mark_zones_from_list(msh%labeled_zones,&
                        't', this%params%bc_labels)
    call this%bc_inflow%finalize()
    !call this%bc_inflow%set_inflow(params%uinf)
    call bc_list_add(this%bclst, this%bc_inflow)

    call bc_list_add(this%bclst, this%bc_wall)
  end subroutine scalar_scheme_init_common

  !> Initialize all velocity related components of the current scheme
  subroutine scalar_scheme_init(this, msh, lx, params, kspv_init, scheme)
    class(scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), target, intent(inout) :: params
    logical :: kspv_init
    character(len=*), intent(in) :: scheme

    integer :: idx(4)
    real :: dx,dy,dz

    call scalar_scheme_init_common(this, msh, lx, params, scheme)
    
    call neko_field_registry%add_field(this%dm_Xh, 's')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')
    this%s => neko_field_registry%get_field('s')

    ! todo: note, adhoc init
    ! ADHOC INITAL VALUE SECTION FOR THE SCALAR

    do i = 1, this%dm_Xh%size()
      idx = nonlinear_index(i, this%Xh%lx,this%Xh%lx,this%Xh%lx)
      dx = this%dm_Xh%x(idx(1), idx(2), idx(3), idx(4))
      dy = this%dm_Xh%y(idx(1), idx(2), idx(3), idx(4))
      dz = this%dm_Xh%z(idx(1), idx(2), idx(3), idx(4))
      this%s%x(idx(1), idx(2), idx(3), idx(4)) = &
         exp(-(dx**2 + dy**2 + dz**2))
    end do

      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_memcpy(this%s%x, this%s%x_d, this%s%dof%size(), &
                            HOST_TO_DEVICE)
      end if

    if (kspv_init) then
       ! todo parameter file ksp tol should be added
       call scalar_scheme_solver_factory(this%ksp, this%dm_Xh%size(), &
            params%ksp_vel, params%abstol_vel)
       call scalar_scheme_precon_factory(this%pc, this%ksp, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst, params%pc_vel)
    end if

    call neko_log%end_section()
  end subroutine scalar_scheme_init


  !> Deallocate a scalar formulation
  subroutine scalar_scheme_free(this)
    class(scalar_scheme_t), intent(inout) :: this

!    if (allocated(this%bc_inflow)) then
!       call this%bc_inflow%free()
!    end if

    write(*,*) "this Xh"
    call space_free(this%Xh)    

    if (allocated(this%ksp)) then
       call krylov_solver_destroy(this%ksp)
       deallocate(this%ksp)
    end if

    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if

    call gs_free(this%gs_Xh)

    call coef_free(this%c_Xh)

    call source_scalar_free(this%f_Xh)

    call bc_list_free(this%bclst)

    nullify(this%params)
    
  end subroutine scalar_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine scalar_scheme_validate(this)
    class(scalar_scheme_t), target, intent(inout) :: this

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%s%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp)) then
       call neko_error('No Krylov solver for velocity defined')
    end if
    
    if (.not. associated(this%f_Xh%eval)) then
       call neko_error('No source term defined')
    end if

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

!    select type(ip => this%bc_inflow)
!    type is(usr_inflow_t)
!       call ip%validate
!    end select

    !
    ! Setup checkpoint structure (if everything is fine)
    !
!    @todo no io for now
!    call this%chkp%init(this%u, this%v, this%w, this%p)

    !
    ! Setup mean flow fields if requested
    !
!    if (this%params%stats_mean_flow) then
!       call this%mean%init(this%u, this%v, this%w, this%p)
!    end if

!    if (this%params%stats_mean_sqr_flow) then
!       call this%mean_sqr%init(this%u, this%v, this%w, this%p)
!    end if

  end subroutine scalar_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! @todo Why can't we call the interface here?
  subroutine scalar_scheme_bc_apply(this)
    class(scalar_scheme_t), intent(inout) :: this
    call bc_list_apply_scalar(this%bclst, this%s%x, this%dm_Xh%size())
  end subroutine scalar_scheme_bc_apply
  
  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine scalar_scheme_solver_factory(ksp, n, solver, abstol)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=20), intent(inout) :: solver
    real(kind=rp) :: abstol

    call krylov_solver_factory(ksp, n, solver, abstol)
    
  end subroutine scalar_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine scalar_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=20) :: pctype
    
    call precon_factory(pc, pctype)
    
    select type(pcp => pc)
    type is(jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is(hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, &
                  bclst, trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)
    
  end subroutine scalar_scheme_precon_factory

  !> Initialize source term
  subroutine scalar_scheme_set_source(this, source_term_type, usr_f, usr_f_vec)
    class(scalar_scheme_t), intent(inout) :: this
    character(len=*) :: source_term_type
    procedure(source_term_pw), optional :: usr_f
    procedure(source_term), optional :: usr_f_vec

    if (trim(source_term_type) .eq. 'noforce') then
       call source_scalar_set_type(this%f_Xh, source_scalar_eval_noforce)
!    else if (trim(source_term_type) .eq. 'user' .and. present(usr_f)) then
!       call source_set_pw_type(this%f_Xh, usr_f)
!    else if (trim(source_term_type) .eq. 'user_vector' .and. present(usr_f_vec)) then
!       call source_set_type(this%f_Xh, usr_f_vec)
!    else if (trim(source_term_type) .eq. '') then
!       if (pe_rank .eq. 0) then
!          call neko_warning('No source term defined, using default (noforce)')
!       end if
!       call source_set_type(this%f_Xh, source_eval_noforce)
    else
       call neko_error('Invalid source term')
    end if

  end subroutine scalar_scheme_set_source

  !> Initialize a user defined inflow condition
!  subroutine scalar_scheme_set_usr_inflow(this, usr_eval)
!    class(scalar_scheme_t), intent(inout) :: this
!    procedure(usr_inflow_eval) :: usr_eval

!    select type(bc_if => this%bc_inflow)
!    type is(usr_inflow_t)
!      call bc_if%set_eval(usr_eval)
!    class default
!      call neko_error("Not a user defined inflow condition")
!    end select
    
!  end subroutine scalar_scheme_set_usr_inflow
     
end module scalar
