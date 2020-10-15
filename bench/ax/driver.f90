program axbench
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(field_t) :: w
  integer :: argc, lx, m, n, n_glb, niter, ierr
  character(len=80) :: suffix
  real(kind=dp), allocatable :: u(:), ur(:), us(:), ut(:), wk(:)
  real(kind=dp), allocatable :: g(:, :, :, :, :)
  integer :: i
  
  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: ./axbench <neko mesh> <N>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx
  
  nmsh_file = file_t(fname)
  call nmsh_file%read(msh)  

  call space_init(Xh, GLL, lx, lx, lx)
  call field_init(w, msh, Xh, "w")
 
  allocate(g(6, Xh%lx, Xh%ly, Xh%lz, msh%nelv))
  call setup_g(g, Xh%wx, Xh%lx, Xh%ly, Xh%lz, msh%nelv)

  niter = 100
  n = Xh%lx * Xh%ly * Xh%lz
  allocate(ur(n), us(n), ut(n), wk(n))

  m = n * msh%nelv
  allocate(u(m))

  call set_data(u, w%x, m)

  call ax(w, u, g, ur, us, ut, wk, m)

  n_glb = n * msh%glb_nelv
  
  call set_timer_flop_cnt(0, msh%glb_nelv, Xh%lx, niter, n_glb)
  do i = 1, niter
     call ax(w, u, g, ur, us, ut, wk, m)
  end do
  call set_timer_flop_cnt(1, msh%glb_nelv, Xh%lx, niter, n_glb)
  
  deallocate(u,ur,us,ut,wk)
  call space_free(Xh)
  call field_free(w)
  call mesh_free(msh)
  
  call neko_finalize

end program axbench

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
  integer :: ierr  
  real(kind=dp) :: nxyz, nx
  
  nx = dble(nx1)
  nxyz = dble(nx1 * nx1 * nx1)
  call MPI_Barrier(NEKO_COMM, ierr)  
  if (iset .eq. 0) then
     time0 = MPI_Wtime()
  else
     time1 = MPI_Wtime()
     time1 = time1-time0
     flop_a = (19d0 * nxyz + 12d0 * nx * nxyz) * dble(nelt) * dble(niter)
     if (time1 .gt. 0) mflops = (flop_a)/(1.d6*time1)
     if (pe_rank .eq. 0) then
        write(6,*)
        write(6,1) nelt,pe_size,nx1
        write(6,2) mflops, time1
     endif
1    format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7)
2    format('Tot MFlops = ', 1pe12.4, ', Time        = ', e12.4)
  endif

end subroutine set_timer_flop_cnt


