!> Krylov preconditioner
module hsmg
  use math
  use utils
  use precon
  use ax_product
  use gather_scatter
  use fast3d
  use bc
  use cg
  use dirichlet
  use fdm
  use schwarz
  use ax_helm
  use gmres
  use tensor
  use jacobi
  implicit none

  !Struct to arrange our multigridlevels
  type, private :: multigrid_t
     type(dofmap_t), pointer :: dof
     type(gs_t), pointer  :: gs_h
     type(space_t), pointer :: Xh
     type(coef_t), pointer :: coef
     type(bc_list_t), pointer :: bclst
     type(schwarz_t), pointer :: schwarz
     type(field_t), pointer :: e
  end type multigrid_t
 
  type, public, extends(pc_t) :: hsmg_t
     type(mesh_t), pointer :: msh
     integer :: nlvls !< Number of levels in the multigrid
     type(multigrid_t), allocatable :: grids(:) !< array for multigrids
     type(gs_t) :: gs_crs, gs_mg !< gather scatter for lower levels
     type(space_t) :: Xh_crs, Xh_mg !< spaces for lower levels
     type(dofmap_t) :: dm_crs, dm_mg 
     type(coef_t) :: c_crs, c_mg 
     type(dirichlet_t) :: bc_crs, bc_mg
     type(bc_list_t) :: bclst_crs, bclst_mg
     type(schwarz_t) :: schwarz, schwarz_mg, schwarz_crs !< Schwarz decompostions
     type(field_t) :: e, e_mg, e_crs !< Solve fields
     type(cg_t) :: crs_solver !< Solver for course problem
     integer :: niter = 10 !< Number of iter of crs sovlve
     type(jacobi_t) :: pc_crs !< Some basic precon for crs
     type(ax_helm_t) :: ax !< Matrix for crs solve
     real(kind=dp), allocatable :: jh(:,:) !< Interpolator crs -> fine
     real(kind=dp), allocatable :: jht(:,:)!< Interpolator crs -> fine, transpose
     real(kind=dp), allocatable :: r(:)!< Residual work array
     real(kind=dp), allocatable :: jhfc(:,:)!< Interpolator fine -> crs
     real(kind=dp), allocatable :: jhfct(:,:) !< Interpolator fine -> crs, transpose
     real(kind=dp), allocatable :: w(:) !< work array
  contains
     procedure, pass(this) :: init => hsmg_init
     procedure, pass(this) :: free => hsmg_free
     procedure, pass(this) :: solve => hsmg_solve
     procedure, pass(this) :: update => hsmg_set_h
  end type hsmg_t
  
  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
contains
  !> @note I do not think we actually use the same grids as they do in the original!
  subroutine hsmg_init(this, msh, Xh, coef, dof, gs_h, bclst)
    class(hsmg_t), intent(inout) :: this
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh
    type(coef_t), intent(inout), target :: coef
    type(dofmap_t), intent(inout), target :: dof
    type(gs_t), intent(inout), target :: gs_h 
    type(bc_list_t), intent(inout), target :: bclst
    integer :: lx, n
    
    call this%free()
    if(Xh%lx .lt. 5) then
       call neko_error('Insufficient number of GLL points for hsmg precon. Minimum degree 4 and 5 GLL points required.')
    end if
    this%nlvls = 3 
    this%msh => msh
    allocate(this%grids(this%nlvls))
    allocate(this%jh(Xh%lxy,this%nlvls))
    allocate(this%jht(Xh%lxy,this%nlvls))
    allocate(this%jhfc(Xh%lxy,this%nlvls))
    allocate(this%jhfct(Xh%lxy,this%nlvls))
    allocate(this%w(dof%n_dofs))
    allocate(this%r(dof%n_dofs))

    ! Compute all elements as if they are deformed
    call mesh_all_deformed(msh)

    n = dof%n_dofs
    call field_init(this%e, dof,'work array')
    
    call space_init(this%Xh_crs, GLL, 2, 2, 2)
    this%dm_crs = dofmap_t(msh, this%Xh_crs) 
    call gs_init(this%gs_crs, this%dm_crs)
    call field_init(this%e_crs, this%dm_crs,'work crs')
    call coef_init(this%c_crs, this%gs_crs)
    
    call this%pc_crs%init(this%c_crs, this%dm_crs, this%gs_crs)
    call this%crs_solver%init(this%dm_crs%n_dofs, M= this%pc_crs)

    call this%bc_crs%init(this%dm_crs)
    call this%bc_crs%mark_zone(msh%outlet)
    call this%bc_crs%finalize()
    call this%bc_crs%set_g(0d0)
    call bc_list_init(this%bclst_crs)
    call bc_list_add(this%bclst_crs, this%bc_crs)

    call space_init(this%Xh_mg, GLL, 4, 4, 4)
    this%dm_mg = dofmap_t(msh, this%Xh_mg) 
    call gs_init(this%gs_mg, this%dm_mg)
    call field_init(this%e_mg, this%dm_mg,'work midl')
    call coef_init(this%c_mg, this%gs_mg)
    
    call this%bc_mg%init(this%dm_mg)
    call this%bc_mg%mark_zone(msh%outlet)
    call this%bc_mg%finalize()
    call this%bc_mg%set_g(0d0)
    call bc_list_init(this%bclst_mg)
    call bc_list_add(this%bclst_mg, this%bc_mg)

    call this%schwarz%init(Xh, dof, gs_h, bclst, msh)
    call this%schwarz_mg%init(this%Xh_mg, this%dm_mg, this%gs_mg,&
                              this%bclst_mg, msh)


    call hsmg_fill_grid(dof, gs_h, Xh, coef, bclst, this%schwarz, &
                        this%e, this%grids, 3) 
    call hsmg_fill_grid(this%dm_mg, this%gs_mg, this%Xh_mg, this%c_mg, &
                        this%bclst_mg, this%schwarz_mg, this%e_mg, &
                        this%grids, 2) 
    call hsmg_fill_grid(this%dm_crs, this%gs_crs, this%Xh_crs, &
                        this%c_crs, this%bclst_crs, this%schwarz_crs, &
                        this%e_crs, this%grids, 1) 

    ! Interpolation between levels
    ! Interpolation operators
    call hsmg_setup_intp(this%grids, this%jh, this%jht, &
                     this%jhfc, this%jhfct, this%nlvls) 
                 
    call hsmg_set_h(this)

  end subroutine hsmg_init
  
  subroutine hsmg_set_h(this)
   class(hsmg_t), intent(inout) :: this
    integer :: i
   !Yeah I dont really know what to do here. For incompressible flow not much happens
!    do i = this%nlvls,2,-1
!       call hsmg_intp_fc(this%grids(i-1),this%grids(i), this%jhfc(1,i-1),this%jhfct(1,i-1))
!    end do
    this%grids(1)%coef%ifh2 = .false.
    call copy(this%grids(1)%coef%h1,this%grids(3)%coef%h1,this%grids(1)%dof%n_dofs)
  end subroutine hsmg_set_h

  subroutine hsmg_intp_fc(grid_c,grid_f,jhfc, jhfct) ! l is coarse level
    type(multigrid_t), intent(inout) :: grid_c, grid_f
    real(kind=dp), intent(inout), dimension(grid_c%Xh%lx*grid_f%Xh%lx) :: jhfc, jhfct
    integer :: nc, nf
    nc = grid_c%Xh%lx
    nf = grid_f%Xh%lx
    call tnsr3d(grid_c%coef%h1,nc,grid_f%coef%h1,nf,jhfc,jhfct, jhfct, grid_c%dof%msh%nelv)
    call tnsr3d(grid_c%coef%h2,nc,grid_f%coef%h2,nf,jhfc,jhfct, jhfct, grid_c%dof%msh%nelv)
  end subroutine hsmg_intp_fc

  subroutine hsmg_fill_grid(dof, gs_h, Xh, coef, bclst, schwarz, e, grids, l) 
    type(dofmap_t), target, intent(in):: dof
    type(gs_t), target, intent(in) :: gs_h
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(bc_list_t), target, intent(in) :: bclst
    type(schwarz_t), target, intent(in) :: schwarz
    type(field_t), target, intent(in) :: e
    integer, intent(in) :: l
    type(multigrid_t), intent(inout), dimension(l) :: grids


    grids(l)%dof => dof
    grids(l)%gs_h => gs_h
    grids(l)%Xh => Xh
    grids(l)%coef => coef
    grids(l)%bclst => bclst
    grids(l)%schwarz => schwarz
    grids(l)%e => e

  end subroutine hsmg_fill_grid


  subroutine hsmg_setup_intp(grids, jh, jht, jhfc, jhfct, lvls)
    integer, intent(in) :: lvls
    type(multigrid_t) :: grids(lvls)
    real(kind=dp), intent(inout), dimension(grids(lvls)%Xh%lxy,lvls) :: jh, jht
    real(kind=dp), intent(inout), dimension(grids(lvls)%Xh%lxy,lvls) :: jhfc, jhfct
    integer l,nf,nc

    do l=1,lvls-1

       nf=grids(l+1)%Xh%lx
       nc=grids(l)%Xh%lx

       !Standard multigrid coarse-to-fine interpolation
       call hsmg_setup_intpm(jh(1,l),grids(l+1)%Xh%zg,grids(l)%Xh%zg,nf,nc)
       call trsp(jht(1,l),nc,jh(1,l),nf)

       !Fine-to-coarse interpolation for variable-coefficient operators
       call hsmg_setup_intpm(jhfc(1,l),grids(l)%Xh%zg,grids(l+1)%Xh%zg,nc,nf)
       call trsp(jhfct(1,l),nf,jhfc(1,l),nc)

    enddo
  end subroutine hsmg_setup_intp

  subroutine hsmg_setup_intpm(jh,zf,zc,nf,nc)
    integer, intent(in) :: nf,nc
    real(kind=dp), intent(inout) :: jh(nf,nc),zf(nf),zc(nc)
    real(kind=dp) ::  w(2*(nf+nc)+4)
    integer :: i, j
    do i=1,nf
       call fd_weights_full(zf(i),zc,nc-1,1,w)
       do j=1,nc
          jh(i,j)=w(j)
       enddo
    enddo
  end subroutine hsmg_setup_intpm
  subroutine hsmg_free(this)
    class(hsmg_t), intent(inout) :: this
    if (allocated(this%grids)) deallocate(this%grids)
    if (allocated(this%jh)) deallocate(this%jh)
    if (allocated(this%jht)) deallocate(this%jht)
    if (allocated(this%jhfc)) deallocate(this%jhfc)
    if (allocated(this%jhfct)) deallocate(this%jhfct)
    if (allocated(this%w)) deallocate(this%w)
    if (allocated(this%r)) deallocate(this%r)
    call this%schwarz%free()
    call this%schwarz_mg%free()
    call coef_free(this%c_crs)
    call coef_free(this%c_mg)
    call field_free(this%e)
    call field_free(this%e_mg)
    call field_free(this%e_crs)
    call this%pc_crs%free()
    call gs_free(this%gs_crs)
    call gs_free(this%gs_mg)
    call this%crs_solver%free()
    


  end subroutine hsmg_free

  !> The h1mg preconditioner from Nek5000.
  subroutine hsmg_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(hsmg_t), intent(inout) :: this
    real(kind=dp), dimension(n), intent(inout) :: z
    real(kind=dp), dimension(n), intent(inout) :: r
    integer :: i
    type(ksp_monitor_t) :: crs_info
    
    !We should not work with the input 
    call copy(this%r, r, n)

    !OVERLAPPING Schwarz exchange and solve
    call this%grids(3)%schwarz%compute(z,this%r)      
    ! DOWNWARD Leg of V-cycle, we are pretty hardcoded here but w/e
    !In original code they do not do col2 but only on faces
    call col2(this%r, this%grids(3)%coef%mult, &
                 this%grids(3)%dof%n_dofs)
    !Restrict to middle level
    call tnsr1_3d(this%r,this%grids(2)%Xh%lx,this%grids(3)%Xh%lx,&
                  this%jht(1,2),this%jh(1,2), this%jh(1,2), &
                  this%msh%nelv)
    call gs_op_vector(this%grids(2)%gs_h, this%r, &
                      this%grids(2)%dof%n_dofs, GS_OP_ADD)
    !OVERLAPPING Schwarz exchange and solve
    call this%grids(2)%schwarz%compute(this%grids(2)%e%x,this%r)  
                                                     !         T
    call col2(this%r, this%grids(2)%coef%mult, this%grids(2)%dof%n_dofs)
    !restrict residual to crs
    call tnsr1_3d(this%r,this%grids(1)%Xh%lx,this%grids(2)%Xh%lx,&
                  this%jht(1,1),this%jh(1,1), this%jh(1,1), &
                  this%msh%nelv)
    !Crs solve
    call gs_op_vector(this%grids(1)%gs_h, this%r, &
                      this%grids(1)%dof%n_dofs, GS_OP_ADD)
    call bc_list_apply_scalar(this%grids(1)%bclst, this%r, &
                              this%grids(1)%dof%n_dofs)
    crs_info = this%crs_solver%solve(this%Ax, this%grids(1)%e, this%r, &
                                 this%grids(1)%dof%n_dofs, &
                                 this%grids(1)%coef, &
                                 this%grids(1)%bclst, &
                                 this%grids(1)%gs_h, this%niter)    
    call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x,&
                              this%grids(1)%dof%n_dofs)
    ! UNWIND.  No smoothing.    
    !to middle level
    call tnsr3d(this%w,this%grids(2)%Xh%lx,this%grids(1)%e%x, &
                this%grids(1)%Xh%lx,this%jh(1,1), &
                this%jht(1,1), this%jht(1,1), this%msh%nelv)
    call add2(this%grids(2)%e%x,this%w,this%grids(2)%dof%n_dofs)
    !to original grid
    call tnsr3d(this%w,this%grids(3)%Xh%lx,this%grids(2)%e%x, &
                this%grids(2)%Xh%lx,this%jh(1,2), &
                this%jht(1,2), this%jht(1,2), this%msh%nelv)
    call add2(z, this%w, this%grids(3)%dof%n_dofs)
    call gs_op_vector(this%grids(3)%gs_h, z, &
                      this%grids(3)%dof%n_dofs, GS_OP_ADD)
    call col2(z, this%grids(3)%coef%mult, this%grids(3)%dof%n_dofs)
  end subroutine hsmg_solve
end module hsmg
