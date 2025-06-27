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

  integer :: op
  integer, parameter :: nops = 3
  integer :: ops(nops) = (/GS_OP_ADD, GS_OP_MIN, GS_OP_MAX/)
  character(len=3) :: op_names(nops) = (/"ADD", "MIN", "MAX"/)

  argc = command_argument_count()

  if ((argc .lt. 2) .or. (argc .gt. 2)) then
     write(*,*) 'Usage: ./gsbench <neko mesh> <N>'
     stop
  end if

  call neko_init

  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx

  call nmsh_file%init(fname)
  call nmsh_file%read(msh)
  call msh%generate_conn()

  call Xh%init(GLL, lx, lx, lx)

  call dm%init(msh, Xh)
  call gs_h%init(dm)

  n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

  allocate(u(n))
  call rzero(u, n)
  call device_map(u, u_d, n)
  call device_memcpy(u, u_d, n, HOST_TO_DEVICE, sync=.false.)

  do op = 1, nops
     ! warmup
     do i = 1, niter
        call gs_h%op(u, n, ops(op))
        call device_sync()
     end do

     call MPI_Barrier(NEKO_COMM, ierr)

     do i = 1, niter
        t(i) = MPI_Wtime()
        call gs_h%op(u, n, ops(op))
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
        write(6,'(A,A,A,e12.4,A,e12.4)') 'GS_OP_', trim(op_names(op)), &
             ' mean: ', mean, ', stddev: ', error
     end if
  end do

  stop

  deallocate(u)
  call Xh%free()
  call msh%free()

  call neko_finalize

end program gsbench
