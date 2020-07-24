program nekobone
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(field_t) :: x, w, msk
  integer :: argc, lx, n, n_glb, niter, ierr
  character(len=80) :: suffix
  real(kind=dp), allocatable :: f(:), c(:), r(:), p(:), z(:)
  real(kind=dp), allocatable :: g(:, :, :, :, :)

  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: ./nekobone <neko mesh> <N>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx
  
  nmsh_file = file_t(fname)
  call nmsh_file%read(msh)  
  call mesh_generate_conn(msh)

  call space_init(Xh, 1, GLL, lx, lx, lx)

  dm = dofmap_t(msh, Xh)
  call gs_init(gs_h, dm)
  
  call field_init(x, msh, Xh, "x")
  call field_init(w, msh, Xh, "work")
  call field_init(msk, msh, Xh, "mask")

  msk = 1d0
  call set_mask(msk)

  allocate(g(6, Xh%lx, Xh%ly, Xh%lz, msh%nelv))
  call setup_g(g, Xh%wx, Xh%lx, Xh%ly, Xh%lz, msh%nelv)
  
  niter = 100
  n = Xh%lx * Xh%ly * Xh%lz * msh%nelv
  allocate(f(n), c(n), r(n), p(n), z(n))
  call set_multiplicity(c, n, gs_h)
  call set_f(f, c, n, gs_h)
  
  call cg(x, f, g, c, r, w, p, z, n, msk, niter, gs_h)

  n_glb = Xh%lx * Xh%ly * Xh%lz * msh%glb_nelv
  
  call MPI_Barrier(NEKO_COMM, ierr)

  call set_timer_flop_cnt(0, msh%glb_nelv, x%Xh%lx, niter, n_glb)
  call cg(x, f, g, c, r, w, p, z, n, msk, niter, gs_h)
  call set_timer_flop_cnt(1, msh%glb_nelv, x%Xh%lx, niter, n_glb)
  
  deallocate(f, c, g, r, p, z)
  call space_free(Xh)
  call field_free(x)
  call field_free(msk)
  call mesh_free(msh)
  
  call neko_finalize

end program nekobone

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


