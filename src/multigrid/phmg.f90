! Copyright (c) 2024, The Neko Authors
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
!> Hybrid ph-multigrid preconditioner
module phmg
  use num_types, only : rp
  use precon, only : pc_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use field, only : field_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use bc, only : bc_t, bc_list_apply_scalar, bc_list_t, bc_list_add, &
       bc_list_init
  use dirichlet , only : dirichlet_t
  use utils, only : neko_error
  use cheby, only : cheby_t
  use jacobi, only : jacobi_t
  use ax_product, only : ax_t, ax_helm_factory
  use tree_amg_multigrid, only : tamg_solver_t
  use interpolation, only : interpolator_t
  use math, only : copy, col2, add2
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, &
       krylov_solver_factory, krylov_solver_destroy
  implicit none
  private


  type, private :: phmg_lvl_t
     integer :: lvl = -1
     type(space_t) :: Xh
     type(dofmap_t) :: dm_Xh
     type(gs_t) :: gs_h
     type(cheby_t) :: cheby
     type(jacobi_t), allocatable :: jacobi
     type(coef_t) :: coef
     type(bc_list_t) :: bclst
     type(dirichlet_t) :: bc
     type(field_t) :: e, r, w, z
  end type phmg_lvl_t

  type, public :: phmg_hrchy_t
     type(phmg_lvl_t), allocatable :: lvl(:)
  end type phmg_hrchy_t
  
    
  type, public, extends(pc_t) :: phmg_t
     type(tamg_solver_t) :: amg_solver
     integer :: nlvls 
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: w(:)
     type(phmg_hrchy_t) :: phmg_hrchy
     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh
     type(coef_t), pointer :: coef
     type(dofmap_t), pointer :: dof
     type(gs_t), pointer :: gs_h
     type(bc_list_t) :: bclst
     type(dirichlet_t) :: bc
     class(ax_t), allocatable :: ax
     type(cheby_t) :: cheby
     type(interpolator_t), allocatable :: intrp(:)
     type(field_t) :: z
   contains
     procedure, pass(this) :: init => phmg_init
     procedure, pass(this) :: free => phmg_free
     procedure, pass(this) :: solve => phmg_solve
     procedure, pass(this) :: update => phmg_update     
  end type phmg_t

contains

  subroutine phmg_init(this, msh, Xh, coef, dof, gs_h, bclst)
    class(phmg_t), intent(inout), target :: this
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h
    type(bc_list_t), intent(inout), target :: bclst
    integer :: lx_crs, lx_mid, lx_lvls(2)
    integer :: n, i, j

    this%msh => msh
    this%Xh => Xh
    this%coef => coef
    this%dof => dof
    this%gs_h => gs_h

    
    this%nlvls = 2
    lx_crs = 2
    lx_lvls(1) = lx_crs
    if (Xh%lx .lt. 5) then
       lx_mid = max(Xh%lx-1,3)

       if(Xh%lx .le. 2) then
          call neko_error('Polynomial order < 2 not supported for phmg precon')
       end if

    else
       lx_mid = 4
    end if
    lx_lvls(2) = lx_mid
    this%msh => msh

    allocate(this%r(dof%size()))
    allocate(this%w(dof%size()))
    allocate(this%phmg_hrchy%lvl(this%nlvls))
    do i = 1, this%nlvls
       this%phmg_hrchy%lvl(i)%lvl = i
       call this%phmg_hrchy%lvl(i)%Xh%init(GLL, lx_lvls(i), lx_lvls(i), lx_lvls(i))
       call this%phmg_hrchy%lvl(i)%dm_Xh%init(msh, this%phmg_hrchy%lvl(i)%Xh)
       call this%phmg_hrchy%lvl(i)%gs_h%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%coef%init(this%phmg_hrchy%lvl(i)%gs_h)
       call this%phmg_hrchy%lvl(i)%e%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%r%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%w%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%z%init(this%phmg_hrchy%lvl(i)%dm_Xh)


       call this%phmg_hrchy%lvl(i)%cheby%init(this%phmg_hrchy%lvl(i)%dm_Xh%size(), KSP_MAX_ITER)
              
!       call krylov_solver_factory(phmg_hrchy%lvl(i)%, &
!            this%dm_crs%size(), trim(crs_pctype), KSP_MAX_ITER

       this%phmg_hrchy%lvl(i)%coef%ifh2 = coef%ifh2
       call copy(this%phmg_hrchy%lvl(i)%coef%h1, coef%h1, &
            this%phmg_hrchy%lvl(i)%dm_Xh%size())

       call this%phmg_hrchy%lvl(i)%bc%init_base(this%phmg_hrchy%lvl(i)%coef)       
       if (bclst%n .gt. 0 ) then
          do j = 1, bclst%n
             call this%phmg_hrchy%lvl(i)%bc%mark_facets(&
                  bclst%bc(j)%bcp%marked_facet)
          end do
       end if
       call this%phmg_hrchy%lvl(i)%bc%finalize()
       call this%phmg_hrchy%lvl(i)%bc%set_g(0.0_rp)
       call bc_list_init(this%phmg_hrchy%lvl(i)%bclst)
       call bc_list_add(this%phmg_hrchy%lvl(i)%bclst, this%phmg_hrchy%lvl(i)%bc)              
    end do

    ! Create backend specific Ax operator
    call ax_helm_factory(this%ax, full_formulation = .false.)

    call this%cheby%init(dof%size(), KSP_MAX_ITER)

    call this%z%init(dof)


    call this%bc%init_base(this%coef)
    if (bclst%n .gt. 0 ) then
       do j = 1, bclst%n
          call this%bc%mark_facets(bclst%bc(j)%bcp%marked_facet)
       end do
    end if
    call this%bc%finalize()
    call this%bc%set_g(0.0_rp)
    call bc_list_init(this%bclst)
    call bc_list_add(this%bclst, this%bc)              
    
    ! Interpolator Fine + mg levels
    allocate(this%intrp(this%nlvls))
    call this%intrp(2)%init(this%Xh, this%phmg_hrchy%lvl(2)%Xh)
    call this%intrp(1)%init(this%phmg_hrchy%lvl(2)%Xh, this%phmg_hrchy%lvl(1)%Xh)

    call this%amg_solver%init(this%ax, this%phmg_hrchy%lvl(1)%Xh, &
         this%phmg_hrchy%lvl(1)%coef, this%msh, this%phmg_hrchy%lvl(1)%gs_h, 4, &
         this%phmg_hrchy%lvl(1)%bclst, 1)
        
  end subroutine phmg_init

  subroutine phmg_free(this)
    class(phmg_t), intent(inout) :: this

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    
  end subroutine phmg_free

  subroutine phmg_solve(this, z, r, n)
    class(phmg_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(ksp_monitor_t) :: ksp_results
         
    !We should not work with the input
    call copy(this%r, r, n)

    call copy(this%z%x, z, n)

    ksp_results =  this%cheby%solve(this%Ax, this%z, &
         this%r, this%dof%size(), &
         this%coef, this%bclst, &
         this%gs_h, niter = 10)

    call copy(z, this%z%x, n)

    call this%ax%compute(this%w, this%z%x, this%coef, this%msh, this%Xh)

    this%w = this%r - this%w 
    
    ! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
    !Restrict to middle level
    call this%intrp(2)%map(this%phmg_hrchy%lvl(2)%r%x, this%w, this%msh%nelv, &
         this%phmg_hrchy%lvl(2)%Xh)
    call this%phmg_hrchy%lvl(2)%gs_h%op(this%phmg_hrchy%lvl(2)%r%x, &
                                        this%phmg_hrchy%lvl(2)%dm_Xh%size(), &
                                        GS_OP_ADD)
    call col2(this%phmg_hrchy%lvl(2)%r%x, this%phmg_hrchy%lvl(2)%coef%mult, &
              this%phmg_hrchy%lvl(2)%dm_Xh%size())

    call bc_list_apply_scalar(this%phmg_hrchy%lvl(2)%bclst, &
                              this%phmg_hrchy%lvl(2)%r%x, &
                              this%phmg_hrchy%lvl(2)%dm_Xh%size())

    this%phmg_hrchy%lvl(2)%z = 0.0_rp
    ksp_results =  this%phmg_hrchy%lvl(2)%cheby%solve(this%Ax, this%phmg_hrchy%lvl(2)%z, &
         this%phmg_hrchy%lvl(2)%r%x, this%phmg_hrchy%lvl(2)%dm_Xh%size(), &
         this%phmg_hrchy%lvl(2)%coef, this%phmg_hrchy%lvl(2)%bclst, &
         this%phmg_hrchy%lvl(2)%gs_h, niter = 10)

    call this%ax%compute(this%phmg_hrchy%lvl(2)%w%x, this%phmg_hrchy%lvl(2)%z%x, &
         this%phmg_hrchy%lvl(2)%coef, this%msh, this%phmg_hrchy%lvl(2)%Xh)
    
    this%phmg_hrchy%lvl(2)%w%x = this%phmg_hrchy%lvl(2)%r%x - this%phmg_hrchy%lvl(2)%w%x
    
    !restrict residual to crs    
    call this%intrp(1)%map(this%phmg_hrchy%lvl(1)%r%x, this%phmg_hrchy%lvl(2)%w%x, &
                           this%msh%nelv, this%phmg_hrchy%lvl(1)%Xh)
    
    
    ! Crs solve
    call this%phmg_hrchy%lvl(1)%gs_h%op(this%phmg_hrchy%lvl(1)%r%x, &
         this%phmg_hrchy%lvl(1)%dm_Xh%size(), GS_OP_ADD)

    call col2(this%phmg_hrchy%lvl(1)%r%x, this%phmg_hrchy%lvl(1)%coef%mult, &
         this%phmg_hrchy%lvl(1)%dm_Xh%size())
    
    call bc_list_apply_scalar(this%phmg_hrchy%lvl(1)%bclst, &
                              this%phmg_hrchy%lvl(1)%r%x, &
                              this%phmg_hrchy%lvl(1)%dm_Xh%size())

    call this%amg_solver%solve(this%phmg_hrchy%lvl(1)%e%x, &
                               this%phmg_hrchy%lvl(1)%r%x, &
                               this%phmg_hrchy%lvl(1)%dm_Xh%size())
    
    call bc_list_apply_scalar(this%phmg_hrchy%lvl(1)%bclst, &
                              this%phmg_hrchy%lvl(1)%e%x,&
                              this%phmg_hrchy%lvl(1)%dm_Xh%size())


    call this%intrp(1)%map(this%phmg_hrchy%lvl(2)%w%x, &
                           this%phmg_hrchy%lvl(1)%e%x, &
                           this%msh%nelv, this%phmg_hrchy%lvl(2)%Xh)

    call add2(this%phmg_hrchy%lvl(2)%e%x, this%phmg_hrchy%lvl(2)%w%x, &
         this%phmg_hrchy%lvl(2)%dm_Xh%size())

    ksp_results =  this%phmg_hrchy%lvl(2)%cheby%solve(this%Ax, this%phmg_hrchy%lvl(2)%e, &
         this%phmg_hrchy%lvl(2)%r%x, this%phmg_hrchy%lvl(2)%dm_Xh%size(), &
         this%phmg_hrchy%lvl(2)%coef, this%phmg_hrchy%lvl(2)%bclst, &
         this%phmg_hrchy%lvl(2)%gs_h, niter = 10)
                 
    call this%intrp(2)%map(this%w, this%phmg_hrchy%lvl(2)%e%x, &
                           this%msh%nelv, this%Xh)
    call add2(z, this%w, this%dof%size())
    call this%gs_h%op(z, this%dof%size(), GS_OP_ADD)
    call col2(z, this%coef%mult, this%dof%size())

    call copy(this%z%x, z, n)

    ksp_results =  this%cheby%solve(this%Ax, this%z, &
         this%r, this%dof%size(), &
         this%coef, this%bclst, &
         this%gs_h, niter = 10)

    call copy(z, this%z%x, n)


  end subroutine phmg_solve

  subroutine phmg_update(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_update
  
end module phmg
