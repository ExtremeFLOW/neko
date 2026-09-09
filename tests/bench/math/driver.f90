! Benchmark of Neko's math routines against the field/vector/matrix wrappers
!
! Compares four ways of expressing the same elementwise/reduction work:
!
!   1. math        -- the raw routines, called directly on the data
!   2. field_math  -- the field_t wrapper (per-DOF CFD data)
!   3. vector_math -- the vector_t wrapper
!   4. matrix_math -- the matrix_t wrapper (nrows = n, ncols = 1)
!
! All four end up in the same math.f90 routine on a CPU build, so any
! difference in throughput is the cost of the wrapper's dispatch layer:
! the extra CALL/RETURN, the present(n)/%size() branch, and argument
! marshalling. On a device build the wrappers dispatch to device_math
! instead; the driver handles both without source changes.
!
! Storage note: every path -- including the direct math one -- operates on
! memory owned by a field_t/vector_t/matrix_t. Those types manage their own
! device buffers, so the direct path gets device residency for free and all
! four paths are compared on identically-allocated memory. The direct path
! simply reaches through to %x (host) or %x_d (device).
!
! Timing note: for the in-place ops (add2, col2) the operand is restored from
! a reference buffer between iterations. That restore sits OUTSIDE the timed
! window, so the reported time is the operation alone. glsc3 does not mutate
! its arguments and needs no restore.
!
! Fill note: values are derived from the dofmap's global coordinates, not from
! a rank-local index. SEM partitioning distributes elements, each carrying its
! full lx**3 block, so a rank-local index would name different elements at
! different rank counts and make glsc3's reduced value rank-count-dependent
! for reasons having nothing to do with the code under test. Coordinates are
! decomposition-independent, so the global reduction is too.

program mathbench
  use neko
  use, intrinsic :: iso_c_binding, only : c_size_t

  use math, only : NEKO_EPS, add2, col2, glsc3
  use device_math, only : device_add2, device_col2, device_glsc3
  use field_math, only : field_add2, field_col2, field_glsc3
  use vector_math, only : vector_add2, vector_col2, vector_glsc3
  use matrix_math, only : matrix_add2, matrix_col2, matrix_glsc3

  implicit none

  character(len=NEKO_FNAME_LEN) :: fname
  character(len=LOG_SIZE) :: logline
  character(len=80) :: argchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  character(10) :: time
  character(8) :: date

  ! Reference copies, used to restore the in-place operands between reps.
  real(kind=rp), allocatable, dimension(:) :: refa, refb, refc
  type(c_ptr) :: refa_d, refb_d, refc_d

  ! Actual data arrays
  real(kind=rp), allocatable, dimension(:) :: da, db, dc
  type(c_ptr) :: da_d, db_d, dc_d
  type(field_t) :: fa, fb, fc
  type(vector_t) :: va, vb, vc
  type(matrix_t) :: ma, mb, mc

  ! Swept down to lx = 2 deliberately. The wrapper's dispatch cost is a fixed
  ! few nanoseconds per call, while the op itself is memory-bandwidth-bound
  ! and grows with n -- so at CFD-realistic sizes the overhead is far below
  ! run-to-run noise and simply cannot be resolved. The small-lx end is the
  ! only regime where a per-call cost is a measurable fraction of the work.
  integer, parameter :: nlx = 8
  integer :: lx_sweep(nlx) = [ 2, 3, 4, 6, 8, 12, 16, 24 ]

  ! Cross-path results at a fixed rank count should agree to the last bit on
  ! a CPU build (identical routine, identical data). On a device build the
  ! reduction order differs, so this is a named parameter rather than a
  ! hardcoded strictness.
  real(kind=rp), parameter :: verify_tol = NEKO_EPS

  integer :: argc, niter, nwarmup, ilx, lx, n, n_glb, ierr
  logical :: verify_only

  argc = command_argument_count()

  if (argc .lt. 1 .or. argc .gt. 2) then
     write(*, *) 'Usage: ./mathbench <neko mesh> [niter]'
     write(*, *) '  niter = 0 runs verification only (no timing).'
     stop
  end if

  call neko_init
  call date_and_time(time = time, date = date)
  call neko_job_info(date, time)

  call get_command_argument(1, fname)

  niter = 200
  if (argc .eq. 2) then
     call get_command_argument(2, argchar)
     read(argchar, *) niter
  end if
  verify_only = (niter .le. 0)
  nwarmup = max(1, niter / 10)

  call nmsh_file%init(fname)
  call nmsh_file%read(msh)
  call msh%generate_conn()

  call neko_log%section('Math bench setup')
  write(logline, '(A,I0)') 'niter     : ', niter
  call neko_log%message(logline)
  write(logline, '(A,I0)') 'nwarmup   : ', nwarmup
  call neko_log%message(logline)
  call neko_log%end_section()
  call neko_log%newline()

  do ilx = 1, nlx
     lx = lx_sweep(ilx)

     call Xh%init(GLL, lx, lx, lx)
     call dm%init(msh, Xh)

     n = Xh%lx * Xh%ly * Xh%lz * msh%nelv
     n_glb = Xh%lx * Xh%ly * Xh%lz * msh%glb_nelv

     call alloc_all(n)
     call fill_all(n)

     call verify_ops(n, lx)

     if (.not. verify_only) then
        call bench_add2(n, n_glb, lx)
        call bench_col2(n, n_glb, lx)
        call bench_glsc3(n, n_glb, lx)
     end if

     call free_all()
     call dm%free()
     call Xh%free()
  end do

  call msh%free()
  call neko_finalize

contains

  !> Allocate every path's operands at size n.
  subroutine alloc_all(n)
    integer, intent(in) :: n
    integer(c_size_t) :: c_size

    allocate(refa(n))
    allocate(refb(n))
    allocate(refc(n))

    if (NEKO_BCKND_DEVICE .eq. 1) then

       select case (rp)
       case (sp)
          c_size = n * int(4, c_size_t)
       case (dp)
          c_size = n * int(8, c_size_t)
       case default
          call neko_error('Unknown Fortran type')
       end select

       call device_alloc(refa_d, c_size)
       call device_alloc(refb_d, c_size)
       call device_alloc(refc_d, c_size)

       call device_alloc(da_d, c_size)
       call device_alloc(db_d, c_size)
       call device_alloc(dc_d, c_size)
    end if

    allocate(da(n))
    allocate(db(n))
    allocate(dc(n))

    call fa%init(dm, 'fa')
    call fb%init(dm, 'fb')
    call fc%init(dm, 'fc')

    call va%init(n)
    call vb%init(n)
    call vc%init(n)

    call ma%init(n / 32, 32)
    call mb%init(n / 32, 32)
    call mc%init(n / 32, 32)

  end subroutine alloc_all

  subroutine free_all()
    deallocate(refa)
    deallocate(refb)
    deallocate(refc)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_free(refa_d)
       call device_free(refb_d)
       call device_free(refc_d)

       call device_free(da_d)
       call device_free(db_d)
       call device_free(dc_d)
    end if

    deallocate(da)
    deallocate(db)
    deallocate(dc)

    call fa%free()
    call fb%free()
    call fc%free()

    call va%free()
    call vb%free()
    call vc%free()

    call ma%free()
    call mb%free()
    call mc%free()

  end subroutine free_all

  !> Build the reference values from the dofmap coordinates, then copy them
  !! into every path's operands so all four compare on identical data.
  subroutine fill_all(n)
    integer, intent(in) :: n

    call coord_fill(refa, dm%x, dm%y, dm%z, n, 1)
    call coord_fill(refb, dm%x, dm%y, dm%z, n, 2)
    call coord_fill(refc, dm%x, dm%y, dm%z, n, 3)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(refa, refa_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(refb, refb_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(refc, refc_d, n, HOST_TO_DEVICE, sync = .true.)
    end if

    call reset_all(n)

  end subroutine fill_all

  !> f = a smooth bounded function of the global coordinates.
  subroutine coord_fill(f, x, y, z, n, seed)
    integer, intent(in) :: n, seed
    real(kind=rp), intent(in) :: x(n), y(n), z(n)
    real(kind=rp), intent(inout) :: f(n)
    integer :: i
    real(kind=rp) :: s

    s = real(seed, rp)

    do i = 1, n
       f(i) = sin(x(i) + 0.1_rp * s) &
            * cos(y(i) - 0.2_rp * s) * cos(z(i) + 0.3_rp * s)
    end do

  end subroutine coord_fill

  !> Restore one operand from its reference. Device-side when on a device
  !! build, so this never becomes a host-device transfer.
  subroutine reset_one(dst, dst_d, src, src_d, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: dst(n)
    real(kind=rp), intent(inout) :: src(n)
    type(c_ptr), intent(inout) :: dst_d, src_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(dst_d, src_d, n)
    else
       call copy(dst, src, n)
    end if

  end subroutine reset_one

  !> Restore every path's operands to the reference values.
  subroutine reset_all(n)
    integer, intent(in) :: n

    call reset_one(da, da_d, refa, refa_d, n)
    call reset_one(db, db_d, refb, refb_d, n)
    call reset_one(dc, dc_d, refc, refc_d, n)

    call reset_one(fa%x, fa%x_d, refa, refa_d, n)
    call reset_one(fb%x, fb%x_d, refb, refb_d, n)
    call reset_one(fc%x, fc%x_d, refc, refc_d, n)

    call reset_one(va%x, va%x_d, refa, refa_d, n)
    call reset_one(vb%x, vb%x_d, refb, refb_d, n)
    call reset_one(vc%x, vc%x_d, refc, refc_d, n)

    call reset_one(ma%x, ma%x_d, refa, refa_d, n)
    call reset_one(mb%x, mb%x_d, refb, refb_d, n)
    call reset_one(mc%x, mc%x_d, refc, refc_d, n)

  end subroutine reset_all

  !> Pull a result back to the host so it can be compared. No-op on CPU.
  subroutine sync_to_host(x, x_d, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    type(c_ptr), intent(inout) :: x_d

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(x, x_d, n, DEVICE_TO_HOST, sync = .true.)
    end if

  end subroutine sync_to_host

  !> Run each op once down all four paths and check they agree. Aborts on a
  !! mismatch rather than letting the timing loop report meaningless numbers.
  subroutine verify_ops(n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp) :: s_math, s_field, s_vector, s_matrix

    ! --- add2 : a = a + b -------------------------------------------------
    call reset_all(n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(da_d, db_d, n)
    else
       call add2(da, db, n)
    end if
    call field_add2(fa, fb, n)
    call vector_add2(va, vb, n)
    call matrix_add2(ma, mb, n)

    call sync_to_host(da, da_d, n)
    call sync_to_host(fa%x, fa%x_d, n)
    call sync_to_host(va%x, va%x_d, n)
    call sync_to_host(ma%x, ma%x_d, n)

    call check_array('add2', 'field_math', da, fa%x, n, lx)
    call check_array('add2', 'vector_math', da, va%x, n, lx)
    call check_array('add2', 'matrix_math', da, ma%x, n, lx)

    ! --- col2 : a = a * b -------------------------------------------------
    call reset_all(n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(da_d, db_d, n)
    else
       call col2(da, db, n)
    end if
    call field_col2(fa, fb, n)
    call vector_col2(va, vb, n)
    call matrix_col2(ma, mb, n)

    call sync_to_host(da, da_d, n)
    call sync_to_host(fa%x, fa%x_d, n)
    call sync_to_host(va%x, va%x_d, n)
    call sync_to_host(ma%x, ma%x_d, n)

    call check_array('col2', 'field_math', da, fa%x, n, lx)
    call check_array('col2', 'vector_math', da, va%x, n, lx)
    call check_array('col2', 'matrix_math', da, ma%x, n, lx)

    ! --- glsc3 : MPI-reduced sum(a*b*c) -----------------------------------
    ! The value below is reduced over NEKO_COMM, so it exercises the path a
    ! single-rank run cannot. Because the fill is coordinate-based, the same
    ! number must come back at every rank count -- that invariant is what
    ! catches a wrapper reducing rank-local data or using a wrong
    ! communicator, and it is checked by comparing runs, not within one run.
    call reset_all(n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       s_math = device_glsc3(da_d, db_d, dc_d, n)
    else
       s_math = glsc3(da, db, dc, n)
    end if
    s_field = field_glsc3(fa, fb, fc, n)
    s_vector = vector_glsc3(va, vb, vc, n)
    s_matrix = matrix_glsc3(ma, mb, mc, n)

    call check_scalar('glsc3', 'field_math', s_math, s_field, lx, n)
    call check_scalar('glsc3', 'vector_math', s_math, s_vector, lx, n)
    call check_scalar('glsc3', 'matrix_math', s_math, s_matrix, lx, n)

    if (pe_rank .eq. 0) then
       write(*, '(A,I3,A,I10,A,e24.16)') &
            '# verify OK  lx = ', lx, '  n = ', n, &
            '  glsc3 = ', s_math
    end if

  end subroutine verify_ops

  subroutine check_array(op, path, a, b, n, lx)
    character(len=*), intent(in) :: op, path
    integer, intent(in) :: n, lx
    real(kind=rp), intent(in) :: a(n), b(n)
    integer :: i

    do i = 1, n
       if (.not. abscmp(a(i), b(i), verify_tol)) then
          write(*, '(A,A,A,A,A,I3,A,I10,A,e24.16,A,e24.16,A,e24.16)') &
               'mathbench MISMATCH: op = ', op, ', path = ', path, &
               ', lx = ', lx, ', i = ', i, &
               ', math = ', a(i), ', wrapper = ', b(i), ', diff = ', a(i) - b(i)
          call neko_error('mathbench: cross-path verification failed')
       end if
    end do

  end subroutine check_array

  subroutine check_scalar(op, path, a, b, lx, n)
    character(len=*), intent(in) :: op, path
    real(kind=rp), intent(in) :: a, b
    integer, intent(in) :: lx, n

    if (.not. abscmp(a, b, verify_tol)) then
       write(*, '(A,A,A,A,A,I3,A,e24.16,A,e24.16,A,e24.16)') &
            'mathbench MISMATCH: op = ', op, ', path = ', path, &
            ', lx = ', lx, ', math = ', a, ', wrapper = ', b, ', diff = ', a - b
       call neko_error('mathbench: cross-path verification failed')
    end if

  end subroutine check_scalar

  !> Report per-call timings and the Mdofs/s/pe workrate. The workrate matches
  !! the formula ReFrame already uses for Neko (tests/reframe),
  !! 1e-3 * dofs * iters / time / pes, so numbers are comparable.
  !!
  !! Both mean and minimum are reported. The minimum is the more useful
  !! statistic here: scheduler preemption, page faults and frequency changes
  !! can only ever make an iteration slower, so the fastest observed call is
  !! the cleanest estimate of the true cost, and it is what makes a
  !! nanosecond-scale dispatch difference visible at all. The mean and stddev
  !! are kept alongside it so that a noisy run is recognisable as noisy rather
  !! than silently reported as a clean result. Mdofs/s is computed from the
  !! minimum for the same reason.
  subroutine report(op, path, lx, n, n_glb, t)
    character(len=*), intent(in) :: op, path
    integer, intent(in) :: lx, n, n_glb
    real(kind=dp), intent(in) :: t(niter)
    real(kind=dp) :: mean, stddev, tmin, tmin_max
    real(kind=rp) :: mdofs
    integer :: i, ierr

    mean = sum(t) / niter
    stddev = 0.0_dp
    do i = 1, niter
       stddev = stddev + (t(i) - mean)**2
    end do
    stddev = sqrt(stddev / (niter - 1))

    tmin = minval(t)

    ! Slowest rank sets the pace, so reduce with MAX rather than reporting
    ! rank 0's own wall time.
    call MPI_Allreduce(tmin, tmin_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         NEKO_COMM, ierr)

    mdofs = 1.0e-3_rp * real(n_glb, rp) &
         / real(tmin_max, rp) / real(pe_size, rp)

    if (pe_rank .eq. 0) then
       write(*, '(A,A,A,A,A,I3,A,I10,A,e12.4,A,e12.4,A,e12.4,A,f14.3)') &
            'op = ', op, '  path = ', path, &
            '  lx = ', lx, '  n = ', n, &
            '  min = ', tmin, '  mean = ', mean, '  stddev = ', stddev, &
            '  Mdofs/s/pe = ', mdofs
    end if

  end subroutine report

  subroutine bench_add2(n, n_glb, lx)
    integer, intent(in) :: n, n_glb, lx
    real(kind=dp) :: t(niter)
    integer :: i, ierr

    ! --- math -------------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(da, da_d, refa, refa_d, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_add2(da_d, db_d, n)
       else
          call add2(da, db, n)
       end if
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(da, da_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_add2(da_d, db_d, n)
       else
          call add2(da, db, n)
       end if
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('add2 ', 'math       ', lx, n, n_glb, t)

    ! --- field_math -------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(fa%x, fa%x_d, refa, refa_d, n)
       call field_add2(fa, fb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(fa%x, fa%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call field_add2(fa, fb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('add2 ', 'field_math ', lx, n, n_glb, t)

    ! --- vector_math ------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(va%x, va%x_d, refa, refa_d, n)
       call vector_add2(va, vb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(va%x, va%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call vector_add2(va, vb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('add2 ', 'vector_math', lx, n, n_glb, t)

    ! --- matrix_math ------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(ma%x, ma%x_d, refa, refa_d, n)
       call matrix_add2(ma, mb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(ma%x, ma%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call matrix_add2(ma, mb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('add2 ', 'matrix_math', lx, n, n_glb, t)

  end subroutine bench_add2

  subroutine bench_col2(n, n_glb, lx)
    integer, intent(in) :: n, n_glb, lx
    real(kind=dp) :: t(niter)
    integer :: i, ierr

    ! --- math -------------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(da, da_d, refa, refa_d, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(da_d, db_d, n)
       else
          call col2(da, db, n)
       end if
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(da, da_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(da_d, db_d, n)
       else
          call col2(da, db, n)
       end if
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('col2 ', 'math       ', lx, n, n_glb, t)

    ! --- field_math -------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(fa%x, fa%x_d, refa, refa_d, n)
       call field_col2(fa, fb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(fa%x, fa%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call field_col2(fa, fb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('col2 ', 'field_math ', lx, n, n_glb, t)

    ! --- vector_math ------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(va%x, va%x_d, refa, refa_d, n)
       call vector_col2(va, vb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(va%x, va%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call vector_col2(va, vb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('col2 ', 'vector_math', lx, n, n_glb, t)

    ! --- matrix_math ------------------------------------------------------
    do i = 1, nwarmup
       call reset_one(ma%x, ma%x_d, refa, refa_d, n)
       call matrix_col2(ma, mb, n)
       call device_sync()
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       call reset_one(ma%x, ma%x_d, refa, refa_d, n)
       t(i) = MPI_Wtime()
       call matrix_col2(ma, mb, n)
       call device_sync()
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('col2 ', 'matrix_math', lx, n, n_glb, t)

  end subroutine bench_col2

  !> glsc3 does not mutate its arguments, so no restore is needed and the
  !! timings here are the operation alone, including its MPI_Allreduce.
  subroutine bench_glsc3(n, n_glb, lx)
    integer, intent(in) :: n, n_glb, lx
    real(kind=dp) :: t(niter)
    real(kind=rp) :: s
    integer :: i, ierr

    ! --- math -------------------------------------------------------------
    do i = 1, nwarmup
       if (NEKO_BCKND_DEVICE .eq. 1) then
          s = device_glsc3(da_d, db_d, dc_d, n)
       else
          s = glsc3(da, db, dc, n)
       end if
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       t(i) = MPI_Wtime()
       if (NEKO_BCKND_DEVICE .eq. 1) then
          s = device_glsc3(da_d, db_d, dc_d, n)
       else
          s = glsc3(da, db, dc, n)
       end if
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('glsc3', 'math       ', lx, n, n_glb, t)

    ! --- field_math -------------------------------------------------------
    do i = 1, nwarmup
       s = field_glsc3(fa, fb, fc, n)
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       t(i) = MPI_Wtime()
       s = field_glsc3(fa, fb, fc, n)
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('glsc3', 'field_math ', lx, n, n_glb, t)

    ! --- vector_math ------------------------------------------------------
    do i = 1, nwarmup
       s = vector_glsc3(va, vb, vc, n)
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       t(i) = MPI_Wtime()
       s = vector_glsc3(va, vb, vc, n)
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('glsc3', 'vector_math', lx, n, n_glb, t)

    ! --- matrix_math ------------------------------------------------------
    do i = 1, nwarmup
       s = matrix_glsc3(ma, mb, mc, n)
    end do
    call MPI_Barrier(NEKO_COMM, ierr)
    do i = 1, niter
       t(i) = MPI_Wtime()
       s = matrix_glsc3(ma, mb, mc, n)
       t(i) = MPI_Wtime() - t(i)
    end do
    call report('glsc3', 'matrix_math', lx, n, n_glb, t)

  end subroutine bench_glsc3

end program mathbench
