program cache_stream
  use neko
  implicit none

  real(kind=rp), allocatable :: a(:), b(:), c(:)
  type(c_ptr) :: a_d = C_NULL_PTR
  type(c_ptr) :: b_d = C_NULL_PTR
  type(c_ptr) :: c_d = C_NULL_PTR
  integer :: i, j, ierr, niter
  integer :: n, pow, n_max, n_glb
  real(kind=dp) :: time, BW
  logical :: if_device = .false.


  call neko_init 
  if (pe_rank .eq. 0) then
     write(*,*) 'Bandwidth measurements with cache'
     write(*,*) 'Number of PEs', pe_size
  end if

  pow = 27
  n_max = 2**pow

  n_max = n_max/pe_size
  do while (2**(pow-1) > n_max)
     pow = pow-1
  end do
  n_max = 2**pow

  allocate(a(n_max))
  allocate(b(n_max))
  allocate(c(n_max))

  if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
     (NEKO_BCKND_OPENCL .eq. 1)) then

    if_device = .true.
    call device_map(a, a_d, n_max)
    call device_map(b, b_d, n_max)
    call device_map(c, c_d, n_max)
  end if

  call MPI_Barrier(NEKO_COMM, ierr)
  do j = min(10,max(1,pow-10)), pow
     n = 2**j
     niter = 10*2**pow/n
     call MPI_Barrier(NEKO_COMM)
     time = MPI_Wtime()
     do i = 1, niter
        if (if_device) then
           call device_add2s2(a_d,b_d,1.0_rp,n) 
        else 
           call add2s2(a,b,1.0_rp,n) 
        end if
     end do
     call device_sync()
     call MPI_Barrier(NEKO_COMM)
     time = MPI_Wtime() - time
     call MPI_Allreduce(MPI_IN_PLACE, time, 1, &
          MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
     BW = pe_size*3_dp*n*niter/time*8_dp*1e-9_dp
     if (pe_rank .eq. 0) then
        write(6,*) 
        write(6,1) n, niter
        write(6,2) time, BW
1       format('Array length per rank: ',i9, '   Iterations:       ', i10)
2       format('Total time:         ',e12.4, '   Bandwidth GB/s: ', e12.6)
     end if
  end do
  
  deallocate(a)
  deallocate(b)
  deallocate(c)
  call neko_finalize

end program cache_stream
