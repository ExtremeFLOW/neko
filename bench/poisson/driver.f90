program poisson
  use neko
  use ax_poisson
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(dirichlet_t) :: dir_bc
  type(bc_list_t) :: bclst
  type(field_t) :: x
  type(ax_poisson_t) :: ax
  type(cg_t) :: solver
  integer :: argc, lx, n, n_glb, niter, ierr, it
  character(len=80) :: suffix
  real(kind=dp), allocatable :: f(:)

  argc = command_argument_count()

  if ((argc .lt. 1) .or. (argc .gt. 1)) then
     write(*,*) 'Usage: ./poisson <N>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, lxchar)
  read(lxchar, *) lx
  
  fname = '512.nmsh'
  nmsh_file = file_t(fname)
  call nmsh_file%read(msh)  
  call mesh_generate_conn(msh)

  call space_init(Xh, GLL, lx, lx, lx)

  dm = dofmap_t(msh, Xh)
  call gs_init(gs_h, dm)
  
  call field_init(x, dm, "x")

  n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

  call dir_bc%init(dm)
  call dir_bc%set_g(0d0)
 
  !user specified
  call set_bc(dir_bc, msh)
 
  call dir_bc%finalize()
  call bc_list_init(bclst,1)
  call bc_list_add(bclst,dir_bc)
  call solver%init(n)

  niter = 100
  allocate(f(n))

  !user specified
  call set_f(f, gs_h%c,dm, n, gs_h)
  
  it = solver%solve(ax,x,f, n, bclst, gs_h, niter)

  n_glb = Xh%lx * Xh%ly * Xh%lz * msh%glb_nelv
  
  call MPI_Barrier(NEKO_COMM, ierr)

  call set_timer_flop_cnt(0, msh%glb_nelv, x%Xh%lx, niter, n_glb)
  it = solver%solve(ax,x,f, n, bclst, gs_h, niter)
  call set_timer_flop_cnt(1, msh%glb_nelv, x%Xh%lx, niter, n_glb)

  fname = 'out.fld'
  mf =  file_t(fname)
  call mf%write(x)
  deallocate(f)
  call solver%free()
  call dir_bc%free()
  call bc_list_free(bclst)
  call space_free(Xh)
  call field_free(x)
  call mesh_free(msh) 
  call neko_finalize

end program poisson

subroutine set_timer_flop_cnt(iset, nelt, nx1, niter, n)
  use comm
  use num_types
  implicit none

  integer :: iset
  integer, intent(inout) :: nelt
  integer, intent(inout) :: nx1
  integer, intent(inout) :: niter
  integer, intent(inout) :: n
  real(kind=dp), save :: time0, time1, mflops, flop_a, flop_cg
  real(kind=dp) :: nxyz, nx

  nx = dble(nx1)
  nxyz = dble(nx1 * nx1 * nx1)
  if (iset .eq. 0) then
     time0 = MPI_Wtime()
  else
     flop_a = (19d0 * nxyz + 12d0 * nx * nxyz) * dble(nelt) * dble(niter)
     flop_cg = dble(niter+1) * 15d0 * dble(n)
     time1 = MPI_Wtime()
     time1 = time1-time0
     if (time1 .gt. 0) mflops = (flop_a+flop_cg)/(1.d6*time1)
     if (pe_rank .eq. 0) then
        write(6,*)
        write(6,1) nelt,pe_size,nx1
        write(6,2) mflops, mflops/pe_size
        write(6,3) flop_a,flop_cg
        write(6,4) time1
     endif
1    format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7)
2    format('Tot MFlops = ', 1pe12.4, ', MFlops      = ', e12.4)
3    format('Setup Flop = ', 1pe12.4, ', Solver Flop = ', e12.4)
4    format('Solve Time = ', e12.4)
  endif

end subroutine set_timer_flop_cnt


