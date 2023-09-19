program poisson
  use neko
  use ax_poisson
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar, iterchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(dirichlet_t) :: dir_bc
  type(bc_list_t) :: bclst
  type(field_t) :: x
  type(ax_poisson_t) :: ax
  type(coef_t) :: coef
  type(cg_t) :: solver
  type(ksp_monitor_t) :: ksp_mon
  integer :: argc, lx, n, n_glb, niter, ierr
  character(len=80) :: suffix
  real(kind=rp), allocatable :: f(:)
  real(kind=rp) :: tol
  tol = -1.0

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./poisson <neko mesh> <N> <Iterations>'
     end if
     stop
  end if
  
  call neko_init 
  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  call get_command_argument(3, iterchar)
  read(lxchar, *) lx
  read(iterchar, *) niter
  
  nmsh_file = file_t(fname)
  call nmsh_file%read(msh)  

  call Xh%init(GLL, lx, lx, lx)

  dm = dofmap_t(msh, Xh)
  call gs_h%init(dm)

  call coef%init(gs_h)
  
  call x%init(dm, "x")

  n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

  call dir_bc%init(dm)
  call dir_bc%set_g(real(0.0d0,rp))
 
  !user specified
  call set_bc(dir_bc, msh)
 
  call dir_bc%finalize()
  call bc_list_init(bclst)
  call bc_list_add(bclst,dir_bc)
  call solver%init(n, abs_tol = tol)

  allocate(f(n))

  !user specified
  call rzero(f,n)
  call set_f(f, coef%mult, dm, n, gs_h)
  call bc_list_apply(bclst,f,n)
  ksp_mon = solver%solve(ax, x, f, n, coef, bclst, gs_h, niter)
  n_glb = Xh%lx * Xh%ly * Xh%lz * msh%glb_nelv
  
  call MPI_Barrier(NEKO_COMM, ierr)

  call set_timer_flop_cnt(0, msh%glb_nelv, x%Xh%lx, niter, n_glb, ksp_mon)
  ksp_mon = solver%solve(ax, x, f, n, coef, bclst, gs_h, niter)
  call set_timer_flop_cnt(1, msh%glb_nelv, x%Xh%lx, niter, n_glb, ksp_mon)
  
  fname = 'out.fld'
  mf =  file_t(fname)
  call mf%write(x)
  deallocate(f)
  call solver%free()
  call dir_bc%free()
  call bc_list_free(bclst)
  call Xh%free()
  call x%free()
  call msh%free() 
  call neko_finalize

end program poisson

subroutine set_timer_flop_cnt(iset, nelt, nx1, niter, n, ksp_mon)
  use comm
  use krylov
  use num_types
  implicit none

  integer :: iset
  integer, intent(inout) :: nelt
  integer, intent(inout) :: nx1
  integer, intent(inout) :: niter
  integer, intent(inout) :: n
  type(ksp_monitor_t), intent(in) :: ksp_mon
  real(kind=dp), save :: time0, time1, mflops, flop_a, flop_cg
  real(kind=dp) :: nxyz, nx


  nx = dble(nx1)
  nxyz = dble(nx1 * nx1 * nx1)
  if (iset .eq. 0) then
     time0 = MPI_Wtime()
  else
     flop_a = (15d0 * nxyz + 12d0 * nx * nxyz) * dble(nelt) * dble(niter)
     flop_cg = dble(niter+1) * 15d0 * dble(n)
     time1 = MPI_Wtime()
     time1 = time1-time0
     if (time1 .gt. 0) mflops = (flop_a+flop_cg)/(1.d6*time1)
     if (pe_rank .eq. 0) then
        write(6,*)
        write(6,1) nelt,pe_size,nx1
        write(6,2) mflops, mflops/pe_size
        write(6,3) flop_a,flop_cg
        write(6,4) ksp_mon%res_start, ksp_mon%res_final
        write(6,5) time1
     endif
1    format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7)
2    format('Tot MFlops = ', 1pe12.4, ', MFlops      = ', e12.4)
3    format('Setup Flop = ', 1pe12.4, ', Solver Flop = ', e12.4)
4    format('Start res  = ', 1pe12.4, ', Final res   = ', e12.4)
5    format('Solve Time = ', e12.4)
  endif

end subroutine set_timer_flop_cnt


