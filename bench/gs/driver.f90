! Gather-scatter benchmark
! Jacob Wahlgren 2022

program gsbench
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  integer :: argc, lx, n, n_glb, ierr, i
  character(len=80) :: suffix
  real(kind=dp), allocatable :: u(:)
  type(c_ptr) :: u_d = C_NULL_PTR
  integer, parameter :: niter = 1000
  real(kind=dp), dimension(niter) :: t
  real(kind=dp) :: mean, error

  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: ./gsbench <neko mesh> <N>'
     stop
  end if

  call neko_init

  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx

  nmsh_file = file_t(fname)
  call nmsh_file%read(msh)
  call mesh_generate_conn(msh)

  call space_init(Xh, GLL, lx, lx, lx)

  dm = dofmap_t(msh, Xh)
  call gs_init(gs_h, dm)

  n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

  allocate(u(n))
  call rzero(u, n)
  call device_map(u, u_d, n)
  call device_memcpy(u, u_d, n, HOST_TO_DEVICE)

  ! warmup
  do i = 1, niter
    call gs_op(gs_h, u, n, GS_OP_ADD)
    call device_sync()
  end do

  n_glb = Xh%lx * Xh%ly * Xh%lz * msh%glb_nelv

  call MPI_Barrier(NEKO_COMM, ierr)

  do i = 1, niter
     t(i) = MPI_Wtime()
     call gs_op(gs_h, u, n, GS_OP_ADD)
     call device_sync()
     t(i) = MPI_Wtime() - t(i)
  end do

  mean = sum(t) / niter
  error = 0
  do i = 1, niter
     error = error + (t(i) - mean)**2
  end do
  error = error / (niter - 1)
  error = sqrt(error)

  if (pe_rank .eq. 0) then
    write(*,*)
    write(6,1) mean, error
  end if
1 format('mean: ', e12.4, ', stddev: ', e12.4)

  stop

  deallocate(u)
  call space_free(Xh)
  call mesh_free(msh)

  call neko_finalize

end program gsbench
