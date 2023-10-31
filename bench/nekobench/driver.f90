program nekobench
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(coef_t) :: coef
  class(ax_t), allocatable :: ax_helm
  type(field_t) :: f1, f2, f3, f4
  real(kind=rp) :: byte, flop, n_tot, t0, t1, time
  real(kind=rp) :: c1, c2
  integer :: argc, niter, ierr, lx, nelt
  integer :: i
  
  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     write(*,*) 'Usage: ./nekobench <neko mesh> <N> <niter>'
     stop
  end if
  
  call neko_init 
  
  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx
  call get_command_argument(3, lxchar)
  read(lxchar, *) niter
  
  nmsh_file = file_t(fname)
  call nmsh_file%read(msh) 
  
  !Init things
  call Xh%init(GLL, lx, lx, lx)
  dm = dofmap_t(msh,Xh)
  call gs_h%init(dm)
  call f1%init(dm)
  call f2%init(dm)
  call f3%init(dm)
  call f4%init(dm)
  call coef%init(gs_h)
  call ax_helm_factory(ax_helm)
  n_tot = dble(msh%glb_nelv)*dble(niter)*dble(Xh%lxyz) 

  if (pe_rank .eq. 0) then
     write(*,*) 'Nekobench '
     write(*,*) 'lx:', lx
     write(*,*) 'N elements tot:', msh%glb_nelv
     write(*,*) 'N ranks:', pe_size
     write(*,*) 'N iterations:', niter
     write(*,*) 'lx^3*N elements tot*N iter:', n_tot 
  end if
  call MPI_Barrier(NEKO_COMM, ierr)


  f1 = 1.0_rp
  f2 = 1.0_rp
  f3 = 1.0_rp
  
  !Benchmark copy speed f1 = f2
  call device_sync()
  t0 = MPI_Wtime()
  do i = 1, niter
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_copy(f1%x_d,f2%x_d,dm%size())
    else
        call copy(f1%x,f2%x,dm%size())
     end if
  end do
  call device_sync()
  t1 = MPI_Wtime()

  time = t1 - t0
  !Assume double precision
  byte = n_tot*2d0*8d0*1e-9
  if (pe_rank .eq. 0) then
     write(*,*) 'copy Time:',time
     write(*,*) 'copy GByte:',byte
     write(*,*) 'copy GB/s:', byte/time
  end if

  f1 = 1.0_rp
  f2 = 1.0_rp
  f3 = 1.0_rp
  c1 = 1.0_rp
  c2 = 1.0_rp
  !Benchmark vector add f1 = c1*f2 + c2*f3
  call device_sync()
  t0 = MPI_Wtime()
  do i = 1, niter
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_add3s2(f1%x_d,f2%x_d,f3%x_d,c1,c2,dm%size())
    else
        call add3s2(f1%x,f2%x,f3%x,c1,c2,dm%size())
     end if
  end do
  call device_sync()
  t1 = MPI_Wtime()

  time = t1 - t0
  !Assume double precision
  flop = 3d0*n_tot*1e-9
  byte = n_tot*3d0*8d0*1e-9
  if (pe_rank .eq. 0) then
     write(*,*) 'add3s2 Time:',time
     write(*,*) 'add3s2 GFlop:', flop
     write(*,*) 'add3s2 GByte:',byte
     write(*,*) 'add3s2 Gflop/s:', flop/time
     write(*,*) 'add3s2 GB/s:', byte/time
  end if


  f1 = 1.0_rp
  f2 = 1.0_rp
  
  
  !Benchmark Ax, lets not overengineer this...
  ! To not time the autotune...
  call ax_helm%compute(f2%x, f1%x, coef, msh, Xh)       
  call device_sync()
  t0 = MPI_Wtime()
  !h1, h2 are initialized to 1.0 in this case
  do i = 1, niter
     call ax_helm%compute(f2%x, f1%x, coef, msh, Xh)       
  end do
  call device_sync()
  t1 = MPI_Wtime()

  time = t1 - t0
  flop = (19d0 + 12d0*Xh%lx)*n_tot*1e-9
  !Assume double precision
  byte = n_tot*8d0*8d0*1e-9
  if (pe_rank .eq. 0) then
     write(*,*) 'Ax_helm Time:',time
     write(*,*) 'Ax_helm GFlop:', flop
     write(*,*) 'Ax_helm GByte:',byte
     write(*,*) 'Ax_helm Gflop/s:', flop/time
     write(*,*) 'Ax_helm GB/s:', byte/time
  end if

  call neko_finalize  

end program nekobench

