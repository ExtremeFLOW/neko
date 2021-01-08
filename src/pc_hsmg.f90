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
  implicit none
  
  type, public, extends(pc_t) :: hsmg_t
     type(mesh_t), pointer :: msh
     integer :: nlvls 
     type(multigrid_t), allocatable :: grids(:)
     type(gs_t) :: gs_crs, gs_mg
     type(space_t) :: Xh_crs, Xh_mg
     type(dofmap_t) :: dm_crs, dm_mg
     type(coef_t) :: c_crs, c_mg
     type(dirichlet_t) :: bc_crs, bc_mg
     type(bc_list_t) :: bclst_crs, bclst_mg
     type(schwarz_t) :: schwarz, schwarz_mg, schwarz_crs
     type(field_t) :: e, e_mg, e_crs
     type(cg_t) :: crs_solver
     type(ax_helm_t) :: ax
     real(kind=dp), allocatable :: jh(:,:)
     real(kind=dp), allocatable :: jht(:,:)
     real(kind=dp), allocatable :: jhfc(:,:)
     real(kind=dp), allocatable :: jhfct(:,:) 
     real(kind=dp), allocatable :: w(:)
     integer :: niter = 50
  contains
     procedure, pass(this) :: init => hsmg_init
     procedure, pass(this) :: free => hsmg_free
     procedure, pass(this) :: solve => hsmg_solve
  end type hsmg_t
  !This is very preliminary
  type, private :: multigrid_t
    type(dofmap_t), pointer :: dof
    type(gs_t), pointer  :: gs_h
    type(space_t), pointer :: Xh
    type(coef_t), pointer :: coef
    type(bc_list_t), pointer :: bclst
    type(schwarz_t), pointer :: schwarz
    type(field_t), pointer :: e
  end type multigrid_t 
  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
contains
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
    this%nlvls = 3 
    this%msh => msh
    allocate(this%grids(this%nlvls))

    allocate(this%jh(Xh%lxy,this%nlvls))
    allocate(this%jht(Xh%lxy,this%nlvls))
    allocate(this%jhfc(Xh%lxy,this%nlvls))
    allocate(this%jhfct(Xh%lxy,this%nlvls))
    allocate(this%w(dof%n_dofs))

    ! Compute all elements as if they are deformed
    call mesh_all_deformed(msh)

    n = dof%n_dofs
    !call rone (coef%h1,n)
    !call rzero(coef%h2,n)
    !call rzero(coef%h2inv,n)
    call field_init(this%e, dof,'work array')
    
    call space_init(this%Xh_crs, GLL, 2, 2, 2)
    this%dm_crs = dofmap_t(msh, this%Xh_crs) 
    call gs_init(this%gs_crs, this%dm_crs)
    call field_init(this%e_crs, this%dm_crs,'work crs')
    call coef_init(this%c_crs, this%gs_crs)

    call this%crs_solver%init(this%dm_crs%n_dofs)

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
    call this%schwarz_crs%init(this%Xh_crs, this%dm_crs, this%gs_crs,&
                              this%bclst_crs, msh)


    call fill_grid(dof, gs_h, Xh, coef, bclst, this%schwarz, this%e, this%grids, 3) 
    call fill_grid(this%dm_mg, this%gs_mg, this%Xh_mg, this%c_mg, this%bclst_mg, this%schwarz_mg, this%e_mg, this%grids, 2) 
    call fill_grid(this%dm_crs, this%gs_crs, this%Xh_crs, this%c_crs, this%bclst_crs, this%schwarz_crs,this%e_crs, this%grids, 1) 

    ! INterpoaltion between levels
    call hsmg_setup_intp(this%grids, this%jh, this%jht, &
                     this%jhfc, this%jhfct, this%nlvls)  ! Interpolation operators
                 
    !  setup fast diagonolization method
    !  call h1mg_setup_fdm    ! set up fast diagonalization method
    !  call h1mg_setup_schwarz_wt(.false.)
    !  call hsmg_setup_solve  ! set up the solver
    
    !  l=mg_h1_lmax
    !fucking bloat probably
    call hsmg_set_h(this)
    !  call mg_set_h1  (p_h1 ,l)
    !  call mg_set_h2  (p_h2 ,l)
    !  call mg_set_gb  (p_g,p_b,l)
    !  call mg_set_msk (p_msk,l)

  end subroutine hsmg_init
  
  subroutine hsmg_set_h(this)
    type(hsmg_t), intent(inout) :: this
    integer :: i

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
    call hsmg_tnsr3d(grid_c%coef%h1,nc,grid_f%coef%h1,nf,jhfc,jhfct, jhfct, grid_c%dof%msh%nelv)
    call hsmg_tnsr3d(grid_c%coef%h2,nc,grid_f%coef%h2,nf,jhfc,jhfct, jhfct, grid_c%dof%msh%nelv)
  end subroutine hsmg_intp_fc

  subroutine fill_grid(dof, gs_h, Xh, coef, bclst, schwarz, e, grids, l) 
    type(dofmap_t), target, intent(in):: dof
    type(gs_t), target, intent(in) :: gs_h
    type(space_t), target, intent(in) :: Xh
    type(coef_t), target, intent(in) :: coef
    type(bc_list_t), target, intent(in) :: bclst
    type(schwarz_t), target, intent(in) :: schwarz
    type(field_t), target, intent(in) :: e
    type(multigrid_t), intent(inout), dimension(l) :: grids
    integer, intent(in) :: l

    grids(l)%dof => dof
    grids(l)%gs_h => gs_h
    grids(l)%Xh => Xh
    grids(l)%coef => coef
    grids(l)%bclst => bclst
    grids(l)%schwarz => schwarz
    grids(l)%e => e

  end subroutine fill_grid


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
       call transpose(jht(1,l),nc,jh(1,l),nf)

       !Fine-to-coarse interpolation for variable-coefficient operators
       call hsmg_setup_intpm(jhfc(1,l),grids(l)%Xh%zg,grids(l+1)%Xh%zg,nc,nf)
       call transpose(jhfct(1,l),nf,jhfc(1,l),nc)

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
  end subroutine hsmg_free

  !> The hsmg preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine hsmg_solve(this, z, r, n)
    integer, intent(inout) :: n
    class(hsmg_t), intent(inout) :: this
    real(kind=dp), dimension(n), intent(inout) :: z
    real(kind=dp), dimension(n), intent(inout) :: r
    integer :: i, l, iter


      call this%grids(3)%schwarz%compute(z,r)                ! z := sigma W M       rhs

      !prolly unnecessary
      !call copy(r,rhs,n)                              ! r  := rhs

      do l = this%nlvls-1,2,-1                        ! DOWNWARD Leg of V-cycle
                                                      !          T
      call col2(r, this%grids(l+1)%coef%mult, this%grids(l+1)%dof%n_dofs)
      !if (.not. this%msh%gdim .eq. 3) then
      !   call hsmg_tnsr1_2d(v,nv,nu,A,At)
      !else
         call hsmg_tnsr1_3d(r,this%grids(l)%Xh%lx,this%grids(l+1)%Xh%lx,&
                            this%jht(1,l),this%jh(1,l), this%jh(1,l), &
                            this%msh%nelv)
      !endif
    
         !call hsmg_rstr(this,r)                   ! r   :=  J r
                                                      !  l         l+1
         call gs_op_vector(this%grids(l)%gs_h, r, this%grids(l)%dof%n_dofs, GS_OP_ADD)
!        OVERLAPPING Schwarz exchange and solve:
         call this%grids(l)%schwarz%compute(this%grids(l)%e%x,r)                ! z := sigma W M       rhs
                                                      !  l              l
      enddo
                                                     !         T
     call col2(r, this%grids(2)%coef%mult, this%grids(2)%dof%n_dofs)
     call hsmg_tnsr1_3d(r,this%grids(1)%Xh%lx,this%grids(2)%Xh%lx,&
                        this%jht(1,1),this%jh(1,1), this%jh(1,1), &
                        this%msh%nelv)
                                                      !  l         l+1
     call gs_op_vector(this%grids(1)%gs_h, r, this%grids(1)%dof%n_dofs, GS_OP_ADD)
     call bc_list_apply_scalar(this%grids(1)%bclst, r, this%grids(1)%dof%n_dofs)
     !call hsmg_coarse_solve ( e(is) , r )            ! e  := A   r
     iter = this%crs_solver%solve(this%Ax, this%grids(1)%e, r, &
                                  this%grids(1)%dof%n_dofs, &
                                  this%grids(1)%coef, this%grids(1)%bclst, &
                                  this%grids(1)%gs_h, this%niter)    
     call bc_list_apply_scalar(this%grids(1)%bclst, this%grids(1)%e%x, this%grids(1)%dof%n_dofs)
     ! call h1mg_mask(e(is)),mg_imask(p_msk),nel)       !  1     1   1
     do l = 2,this%nlvls-1                           ! UNWIND.  No smoothing.
        !call hsmg_intp (this%w,this%grids%e(l))                 ! w   :=  J e
        call hsmg_tnsr3d(this%w,this%grids(l)%Xh%lx,this%grids(l-1)%e%x, &
                         this%grids(l-1)%Xh%lx,this%jh(1,l-1), &
                         this%jht(1,l-1), this%jht(1,l-1), this%msh%nelv)
        call add2(this%grids(l)%e%x,this%w,this%grids(l)%dof%n_dofs)
     enddo

     l = this%nlvls 
    !call hsmg_intp(w,e(im),l-1)                     ! w   :=  J e
     call hsmg_tnsr3d(this%w,this%grids(l)%Xh%lx,this%grids(l-1)%e%x, &
                      this%grids(l-1)%Xh%lx,this%jh(1,l-1), &
                      this%jht(1,l-1), this%jht(1,l-1), this%msh%nelv)
     call add2(z, this%w, this%grids(l)%dof%n_dofs)
     call gs_op_vector(this%grids(l)%gs_h, z, this%grids(l)%dof%n_dofs, GS_OP_ADD)
     call col2(z, this%grids(l)%coef%mult, this%grids(l)%dof%n_dofs)
      !call dsavg(z) ! Emergency hack --- to ensure continuous z!

  end subroutine hsmg_solve

  subroutine hsmg_tnsr3d(v,nv,u,nu,A,Bt,Ct, nelv)
    integer, intent(inout) :: nv,nu, nelv
    real(kind=dp), intent(inout) :: v(nv*nv*nv,nelv),u(nu*nu*nu,nelv),A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=dp) :: work(0:(nu+nv)**3),work2(0:nu*nv*nv)
    integer :: ie, i, nunu, nvnu, nvnv
    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    do ie=1,nelv
       call mxm(A,nv,u(1,ie),nu,work,nunu)
       do i=0,nu-1
          call mxm(work(nvnu*i),nv,Bt,nu,work2(nv*nv*i),nv)
       enddo
       call mxm(work2,nvnv,Ct,nu,v(1,ie),nv)
    enddo
  end subroutine hsmg_tnsr3d

  subroutine hsmg_tnsr1_3d(v,nv,nu,A,Bt,Ct, nelv) ! v = [C (x) B (x) A] u
    integer, intent(in) :: nv,nu, nelv
    real(kind=dp), intent(inout) :: v(nv*nv*nv*nelv),A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=dp) :: work(0:(nu+nv)**3),work2(0:(nu+nv)**3)
    integer :: e,e0,ee,es, iu, iv, i, nu3, nv3
    e0=1
    es=1
    ee=nelv

    if (nv.gt.nu) then
       e0=nelv
       es=-1
       ee=1
    endif

    nu3 = nu**3
    nv3 = nv**3

    do e=e0,ee,es
       iu = 1 + (e-1)*nu3
       iv = 1 + (e-1)*nv3
       call mxm(A,nv,v(iu),nu,work,nu*nu)
       do i=0,nu-1
          call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
       enddo
       call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
    enddo

    return
    end subroutine hsmg_tnsr1_3d

end module hsmg
