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
  use bc, only : bc_t
  use bc_list, only : bc_list_t
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
     type(space_t), pointer :: Xh
     type(dofmap_t), pointer :: dm_Xh
     type(gs_t), pointer :: gs_h
     type(cheby_t) :: cheby
     type(jacobi_t) :: jacobi
     type(coef_t), pointer :: coef
     type(bc_list_t) :: bclst
     type(dirichlet_t) :: bc
     type(field_t) :: r, w, z
  end type phmg_lvl_t

  type, public :: phmg_hrchy_t
     type(phmg_lvl_t), allocatable :: lvl(:)
  end type phmg_hrchy_t
  
    
  type, public, extends(pc_t) :: phmg_t
     type(tamg_solver_t) :: amg_solver
     integer :: nlvls 
     type(phmg_hrchy_t) :: phmg_hrchy
     class(ax_t), allocatable :: ax
     type(interpolator_t), allocatable :: intrp(:)
     type(mesh_t), pointer :: msh
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
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof
    type(gs_t), intent(inout), target :: gs_h
    type(bc_list_t), intent(inout), target :: bclst
    integer :: lx_crs, lx_mid
    integer, allocatable :: lx_lvls(:)
    integer :: n, i, j
    class(bc_t), pointer :: bc_j
    
    this%msh => msh

    this%nlvls = Xh%lx - 1

    allocate(lx_lvls(0:this%nlvls - 1))
    lx_lvls(0) = Xh%lx
    do i = 1, this%nlvls -1
       lx_lvls(i) = Xh%lx - i
    end do

    allocate(this%phmg_hrchy%lvl(0:this%nlvls - 1))

    this%phmg_hrchy%lvl(0)%lvl = 0
    this%phmg_hrchy%lvl(0)%Xh => Xh
    this%phmg_hrchy%lvl(0)%coef => coef
    this%phmg_hrchy%lvl(0)%dm_Xh => dof
    this%phmg_hrchy%lvl(0)%gs_h => gs_h
    
    do i = 1, this%nlvls - 1
       allocate(this%phmg_hrchy%lvl(i)%Xh)
       allocate(this%phmg_hrchy%lvl(i)%dm_Xh)
       allocate(this%phmg_hrchy%lvl(i)%gs_h)
       allocate(this%phmg_hrchy%lvl(i)%coef)

       this%phmg_hrchy%lvl(i)%lvl = i
       call this%phmg_hrchy%lvl(i)%Xh%init(GLL, lx_lvls(i), lx_lvls(i), lx_lvls(i))
       call this%phmg_hrchy%lvl(i)%dm_Xh%init(msh, this%phmg_hrchy%lvl(i)%Xh)
       call this%phmg_hrchy%lvl(i)%gs_h%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%coef%init(this%phmg_hrchy%lvl(i)%gs_h)
    end do

    do i = 0, this%nlvls - 1
       call this%phmg_hrchy%lvl(i)%r%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%w%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%z%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       
       call this%phmg_hrchy%lvl(i)%cheby%init(this%phmg_hrchy%lvl(i)%dm_Xh%size(), KSP_MAX_ITER)

       call this%phmg_hrchy%lvl(i)%jacobi%init(this%phmg_hrchy%lvl(i)%coef, &
                                               this%phmg_hrchy%lvl(i)%dm_Xh, &
                                               this%phmg_hrchy%lvl(i)%gs_h)
              
       this%phmg_hrchy%lvl(i)%coef%ifh2 = coef%ifh2
       call copy(this%phmg_hrchy%lvl(i)%coef%h1, coef%h1, &
            this%phmg_hrchy%lvl(i)%dm_Xh%size())

       call this%phmg_hrchy%lvl(i)%bc%init_base(this%phmg_hrchy%lvl(i)%coef)
       if (bclst%size() .gt. 0 ) then
          do j = 1, bclst%size()
             bc_j => bclst%get(j)
             call this%phmg_hrchy%lvl(i)%bc%mark_facets(&
                  bc_j%marked_facet)
          end do
       end if
       call this%phmg_hrchy%lvl(i)%bc%finalize()
       call this%phmg_hrchy%lvl(i)%bc%set_g(0.0_rp)
       call this%phmg_hrchy%lvl(i)%bclst%init()
       call this%phmg_hrchy%lvl(i)%bclst%append(this%phmg_hrchy%lvl(i)%bc)
    end do

    ! Create backend specific Ax operator
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Interpolator Fine + mg levels
    allocate(this%intrp(this%nlvls - 1))
    do i = 1, this%nlvls -1     
       call this%intrp(i)%init(this%phmg_hrchy%lvl(i-1)%Xh, this%phmg_hrchy%lvl(i)%Xh)
    end do
    
    call this%amg_solver%init(this%ax, this%phmg_hrchy%lvl(this%nlvls -1)%Xh, &
                              this%phmg_hrchy%lvl(this%nlvls -1)%coef, this%msh, &
                              this%phmg_hrchy%lvl(this%nlvls-1)%gs_h, 4, &
                              this%phmg_hrchy%lvl(this%nlvls -1)%bclst, 1)
        
  end subroutine phmg_init

  subroutine phmg_free(this)
    class(phmg_t), intent(inout) :: this

    
  end subroutine phmg_free

  subroutine phmg_solve(this, z, r, n)
    class(phmg_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(ksp_monitor_t) :: ksp_results


    associate( mglvl => this%phmg_hrchy%lvl)
      !We should not work with the input
      call copy(mglvl(0)%r%x, r, n)
      
      mglvl(0)%z%x = 0.0_rp
      mglvl(0)%w%x = 0.0_rp

      call phmg_mg_cycle(mglvl(0)%z, mglvl(0)%r, mglvl(0)%w, 0, this%nlvls -1, &
           mglvl, this%intrp, this%msh, this%Ax, this%amg_solver)

      call copy(z, mglvl(0)%z%x, n)
      
    end associate

  end subroutine phmg_solve

  subroutine phmg_update(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_update


  recursive subroutine phmg_mg_cycle(z, r, w, lvl, clvl, &
                                     mg, intrp, msh, Ax, amg_solver)
    type(ksp_monitor_t) :: ksp_results
    integer :: lvl, clvl
    type(phmg_lvl_t) :: mg(0:clvl)
    type(interpolator_t) :: intrp(1:clvl)
    type(tamg_solver_t), intent(inout) :: amg_solver
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t) :: z, r, w
    integer :: i


    ksp_results =  mg(lvl)%cheby%solve(Ax, z, &
                                       r%x, mg(lvl)%dm_Xh%size(), &
                                       mg(lvl)%coef, mg(lvl)%bclst, &
                                       mg(lvl)%gs_h, niter = 15)

    call Ax%compute(w%x, z%x, mg(lvl)%coef, msh, mg(lvl)%Xh)
      
    w%x = r%x - w%x

    call intrp(lvl+1)%map(mg(lvl+1)%r%x, w%x, msh%nelv, mg(lvl+1)%Xh)

    call mg(lvl+1)%gs_h%op(mg(lvl+1)%r%x, mg(lvl+1)%dm_Xh%size(), GS_OP_ADD)
    
    call col2(mg(lvl+1)%r%x, mg(lvl+1)%coef%mult, mg(lvl+1)%dm_Xh%size())


    call mg(lvl+1)%bclst%apply_scalar( &
                              mg(lvl+1)%r%x, &
                              mg(lvl+1)%dm_Xh%size())
      
    mg(lvl+1)%z%x = 0.0_rp      
    if (lvl+1 .eq. clvl) then
       
       call amg_solver%solve(mg(lvl+1)%z%x, &
                             mg(lvl+1)%r%x, &
                             mg(lvl+1)%dm_Xh%size())
      
       call mg(lvl+1)%bclst%apply_scalar( &
                                 mg(lvl+1)%z%x,&
                                 mg(lvl+1)%dm_Xh%size())
    else
       call phmg_mg_cycle(mg(lvl+1)%z, mg(lvl+1)%r, mg(lvl+1)%w, lvl+1, &
            clvl, mg, intrp, msh, Ax, amg_solver)
    end if

    call intrp(lvl+1)%map(w%x, mg(lvl+1)%z%x, msh%nelv, mg(lvl)%Xh)

    call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD)
    
    call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
    
    call mg(lvl)%bclst%apply_scalar(w%x, mg(lvl)%dm_Xh%size())
       
    z%x = z%x + w%x
    
    ksp_results =  mg(lvl)%cheby%solve(Ax, z, &
                                       r%x, mg(lvl)%dm_Xh%size(), &
                                       mg(lvl)%coef, mg(lvl)%bclst, &
                                       mg(lvl)%gs_h, niter = 15)      

    

  end subroutine phmg_mg_cycle
  
end module phmg
