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
  integer :: argc, lx, n, niter, ierr
  character(len=80) :: suffix
  real(kind=dp), allocatable :: f(:), c(:), r(:), p(:), z(:)
  real(kind=dp), allocatable :: g(:, :, :, :, :)
  real(kind=sp), allocatable :: wk(:)
  
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
  allocate(wk(lx*lx))
  call semhat(wk, Xh%dx, Xh%dxt, Xh%zg, lx - 1)

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

  call MPI_Barrier(NEKO_COMM, ierr)

  !! Add cg loop here
  
  deallocate(f, c, g, r, p, z, wk)
  call space_free(Xh)
  call field_free(x)
  call field_free(msk)
  call mesh_free(msh)
  
  call neko_finalize

end program nekobone



