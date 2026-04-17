module hypre
  use num_types, only : rp, i8
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_free, device_memcpy, HOST_TO_DEVICE
  use profiler, only : profiler_start_region, profiler_end_region
  use hypre_boomeramg
  use hypre_ij_interface
  use coefs, only : coef_t
  use comm, only: pe_rank
  use, intrinsic :: iso_c_binding
  implicit none
  private

  interface
      integer (c_int) function HYPRE_Initialize() &
           bind(c, name='HYPRE_Initialize')
        use, intrinsic :: iso_c_binding
        implicit none
      end function HYPRE_Initialize

      integer (c_int) function HYPRE_Finalize() &
           bind(c, name='HYPRE_Finalize')
        use, intrinsic :: iso_c_binding
        implicit none
      end function HYPRE_Finalize
  end interface

  !TODO: this is defined in hypre_ij_interface_wrapper.c
  !      Eventually it should be moved to a more appropriate place.
  interface
      integer (c_int) function HYPRE_init_wrapper() &
           bind(c, name='HYPRE_init_wrapper')
        use, intrinsic :: iso_c_binding
        implicit none
      end function HYPRE_init_wrapper
  end interface

  type, public :: hypre_solver_t
     private
     ! pointer for HYPRE_Solver
     type(c_ptr) :: solver = C_NULL_PTR
     ! Matrix A
     ! pointer for HYPRE_IJMatrix
     type(c_ptr) :: A = C_NULL_PTR
     type(c_ptr) :: parcsr_A = C_NULL_PTR
     ! Right hand side b
     ! pointer for HYPRE_IJVector
     type(c_ptr) :: b = C_NULL_PTR
     type(c_ptr) :: par_b = C_NULL_PTR
     ! Approximated solution x
     ! pointer for HYPRE_IJVector
     type(c_ptr) :: x = C_NULL_PTR
     type(c_ptr) :: par_x = C_NULL_PTR
     ! things that should not be here but are stuck here for convenience
     ! or lazyness
     integer, allocatable :: dofs(:)! dof list here to have i4 instead of i8
     type(c_ptr) :: dofs_d = C_NULL_PTR! dof list on device
   contains
     procedure, pass(this) :: init => hypre_solver_init
     procedure, pass(this) :: free => hypre_solver_free
     procedure, pass(this) :: setup => hypre_solver_setup
     procedure, pass(this) :: solve => hypre_solve
     procedure, pass(this) :: device_solve => hypre_device_solve
     procedure, pass(this) :: set_matrix
     procedure, pass(this) :: set_vector
  end type hypre_solver_t

  public :: hypre_init, hypre_fin

contains

  subroutine hypre_init()
     integer :: ierr
     ! Initialize hypre
     ierr = HYPRE_Initialize()
     ! Set other hypre config parameters
     !call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE, ierr)
     !call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE, ierr)
     !call HYPRE_SetSpGemmUseVendor(0, ierr)
     if (NEKO_BCKND_DEVICE .eq. 1) then
       ierr = HYPRE_init_wrapper()
     end if
  end subroutine hypre_init

  subroutine hypre_fin()
     integer :: ierr
     ierr = HYPRE_Finalize()
  end subroutine hypre_fin

  !> Initialize a hypre solver object (BoomerAMG)
  subroutine hypre_solver_init(this, solver_params)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: solver_params ! Replace with some parameter object
    ! Pass and process parameters here.
    call boomeramg_init(this%solver)
  end subroutine hypre_solver_init

  !> Free the hypre solver object and related things
  subroutine hypre_solver_free(this)
    class(hypre_solver_t), intent(inout) :: this
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_free(this%dofs_d)
    end if
    deallocate(this%dofs)
    ! destroy the hypre solver object
    call boomeramg_destroy(this%solver)
  end subroutine hypre_solver_free

  subroutine set_matrix(this, A)
    class(hypre_solver_t), intent(inout) :: this
    type(c_ptr), intent(in) :: A
    this%A = A
  end subroutine set_matrix

  subroutine set_vector(this, x, b)
    class(hypre_solver_t), intent(inout) :: this
    type(c_ptr), intent(in) :: x, b
    this%b = b
    this%x = x
  end subroutine set_vector

  !> Setup the hypre solver object
  subroutine hypre_solver_setup(this, coef)
    class(hypre_solver_t), intent(inout) :: this
    type(coef_t), intent(in), target :: coef

    ! workaround to have dofs as i4 instead of i8
    call hypre_dofs_workaround(this, coef)

    ! Asseble the linear system
    call hypre_matrix_assemble_from_neko(coef, this%A, this%b, this%x)

    call hypre_matrix_get_object(this%A, this%parcsr_A)
    call hypre_vector_get_object(this%b, this%par_b)
    call hypre_vector_get_object(this%x, this%par_x)

    ! Setup BoomerAMG
    call boomeramg_setup(this%solver, this%parcsr_A, this%par_b, this%par_x)
  end subroutine hypre_solver_setup

  !> Solve system for right hand side f
  !! @param x Approximated solution. Fortran array on host.
  !! @param f Right hand side. Fortran array on host.
  !! @param n Size of vectors
  subroutine hypre_solve(this, x, f, n)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(in) :: f
    ! Copy to hypre vector
    ! (copy to x may be unneeded if zero initial guess is always used)
    call hypre_copy_to_vector(this%x, n, this%dofs, x)
    call hypre_vector_assemble(this%x)
    call hypre_copy_to_vector(this%b, n, this%dofs, f)
    call hypre_vector_assemble(this%b)
    ! Solve
    call boomeramg_solve(this%solver, this%parcsr_A, this%par_b, this%par_x)
    ! Copy from hypre vector
    call hypre_copy_from_vector(this%x, n, this%dofs, x)
  end subroutine hypre_solve

  !> Solve system for right hand side f
  !! @param x Approximated solution. c_ptr on device.
  !! @param f Right hand side. c_ptr on device.
  !! @param n Size of vectors
  subroutine hypre_device_solve(this, x, f, n)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    type(c_ptr), intent(in) :: x
    type(c_ptr), intent(in) :: f
    ! Copy to hypre vector
    ! (copy to x may be unneeded if zero initial guess is always used)
    call profiler_start_region("neko_to_hypre_X")
    call hypre_device_copy_to_vector(this%x, n, this%dofs_d, x)
    call hypre_vector_assemble(this%x)
    call profiler_end_region()
    call profiler_start_region("neko_to_hypre_B")
    call hypre_device_copy_to_vector(this%b, n, this%dofs_d, f)
    call hypre_vector_assemble(this%b)
    call profiler_end_region()
    ! Solve
    call boomeramg_solve(this%solver, this%parcsr_A, this%par_b, this%par_x)
    ! Copy from hypre vector
    call profiler_start_region("hypre_to_neko_X")
    call hypre_device_copy_from_vector(this%x, n, this%dofs_d, x)
    call profiler_end_region()
  end subroutine hypre_device_solve

  subroutine hypre_dofs_workaround(hs, coeff)
    !DIR$ INLINENEVER hypre_dofs_workaround
    type(hypre_solver_t), intent(inout) :: hs
    type(coef_t), intent(in), target :: coeff
    integer :: i, n
    n = coeff%dof%size()
    allocate(hs%dofs(n))
    do i = 1, n
       hs%dofs(i) = coeff%dof%dof(i,1,1,1)
    end do
    print *, "DOF workaround mapping", n
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(hs%dofs, hs%dofs_d, n)
       call device_memcpy(hs%dofs, hs%dofs_d, n, HOST_TO_DEVICE, .true.)
    end if
  end subroutine hypre_dofs_workaround

  subroutine hypre_matrix_assemble_from_neko(coef, A, b, x)
    !DIR$ INLINENEVER hypre_matrix_assemble_from_neko
    type(coef_t), intent(in), target :: coef
    type(c_ptr), intent(inout) :: A, b, x
    real(kind=rp), allocatable :: A_vals(:)
    integer, allocatable :: A_rows(:), A_cols(:)
    integer :: nrows
    real(kind=rp), target, allocatable :: vals(:)
    integer, target, allocatable :: ncols(:), rows(:), cols(:)
    real(kind=rp) :: tmp2
    real(kind=rp) :: tmp_vals(1)
    integer :: tmp_rows(1), tmp_cols(1), ncol(1)
    integer :: ilower, iupper, jlower, jupper
    integer :: e, k, j, i, s, l
    integer :: nnz, nelv, lx, idof
    type(c_ptr) :: rows_d, cols_d, vals_d, ncols_d
    print *, "mtx_assemble", coef%dof%size()
    lx = coef%dof%Xh%lx
    nelv = coef%dof%msh%nelv
    ! storing the matrix in (i,j,val) format needs one entry per dof contribution
    ! this means we will need elems
    nnz = nelv*lx*lx*lx*lx*3
    allocate(A_vals(nnz))
    A_vals = 0.0_rp
    allocate(A_rows(nnz))
    allocate(A_cols(nnz))
    A_rows = 0
    A_cols = 0

    !note: we cast from i8 to integer when we fill A_rows from coef%dof%dof
    associate( D => coef%dof%Xh%dx, Dt => coef%dof%Xh%dxt, &
         G11 => coef%G11, G22 => coef%G22, G33 => coef%G33, &
         G12 => coef%G12, G13 => coef%G13, G23 => coef%G23)
      idof = 1 ! fortran indexing
      ! Loop over mesh elements
      do e = 1, nelv
         ! Compute the action of the derivative operator (D u)
         ! TODO: this assumes the mesh is structured
         !       so it's missing Gij with i .ne. j
         do k = 1, lx
            do j = 1, lx
               do i = 1, lx

                 ! Hacky if statement to support 2 mpi ranks (more not supported).
                 ! NOTE: Periodic BC not supported on 2 mpi ranks.
                 if((pe_rank .eq. 0).or.(.not. coef%dof%shared_dof(i,j,k,e))) then
                  do s = 1, lx
                     tmp2 = 0.0_rp
                     do l = 1, lx
                        tmp2 = tmp2 + Dt(i,l) * D(l,s)
                     end do

                     A_vals(idof) = tmp2 * G11(s,j,k,e)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(s,j,k,e)
                     idof = idof + 1
                  end do

                  do s = 1, lx
                     tmp2 = 0.0_rp
                     do l = 1, lx
                        tmp2 = tmp2 + Dt(j,l) * D(l,s)
                     end do

                     A_vals(idof) = tmp2 * G22(i,s,k,e)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(i,s,k,e)
                     idof = idof + 1
                  end do

                  do s = 1, lx
                     tmp2 = 0.0_rp
                     do l = 1, lx
                        tmp2 = tmp2 + Dt(k,l) * D(l,s)
                     end do

                     A_vals(idof) = tmp2 * G33(i,j,s,e)
                     A_rows(idof) = coef%dof%dof(i,j,k,e)
                     A_cols(idof) = coef%dof%dof(i,j,s,e)
                     idof = idof + 1
                  end do
                 end if

               end do
            end do
         end do
      end do ! e = 1, n
    end associate

    ! Update nnz based on how many nonzeros were filled (counting duplicates)
    nnz = idof - 1

    ! Get index range of dofs on current rank
    ! (row partition) = (col partition) for square linear systems.
    ilower = minval(A_rows(1:nnz))
    iupper = maxval(A_rows(1:nnz))
    jlower = ilower
    jupper = iupper
    ! TODO: For MPI parallelism, duplicated dofs need to be assigned
    !       to a single "owning" rank, which is currently not supported
    !       by the current neko dofmap
    print *, "RANK", pe_rank, "ilower", ilower, "iupper", iupper, "jlower", jlower, "jupper", jupper

    ! Initialize matrix
    call hypre_matrix_init(1, ilower, iupper, jlower, jupper, A)

    ! Initialize vectors
    call hypre_vector_init(1, ilower, iupper, b)
    call hypre_vector_init(1, jlower, jupper, x)

    ! Overallocate here.
    allocate(rows(nnz))
    allocate(cols(nnz))
    allocate(vals(nnz))
    call count_and_fill_rows(A_rows, rows, nrows, nnz)
    allocate(ncols(nrows))
    call count_and_fill_cols_vals(A_rows, A_cols, A_vals, rows, cols, &
         nrows, ncols, vals, nnz)

    ! Fill matrix
    ! We do this the lazy and expensive way, one dof at a time.
    ! The interface expects array, so we fill some dummy arrays of size 1
    ! and pass these to the hypre interface
    if (NEKO_BCKND_DEVICE .eq. 1) then
       rows_d = C_NULL_PTR
       call device_map(rows, rows_d, nrows)
       call device_memcpy( rows, rows_d, nrows, HOST_TO_DEVICE, .true.)
       ncols_d = C_NULL_PTR
       call device_map(ncols,  ncols_d, nrows)
       call device_memcpy( ncols, ncols_d, nrows, HOST_TO_DEVICE, .true.)
       cols_d = C_NULL_PTR
       call device_map(cols, cols_d, nnz)
       call device_memcpy( cols, cols_d, nnz, HOST_TO_DEVICE, .true.)
       vals_d = C_NULL_PTR
       call device_map(vals, vals_d, nnz)
       call device_memcpy( vals, vals_d, nnz, HOST_TO_DEVICE, .true.)

       call hypre_device_matrix_update(A, nrows, ncols_d, rows_d, cols_d, vals_d)
    else
       !call hypre_matrix_update(A, nrows, ncols, rows, cols, vals)
       do i = 1, nnz
          ncol(1) = 1
          tmp_rows(1) = A_rows(i)
          tmp_cols(1) = A_cols(i)
          tmp_vals(1) = A_vals(i)
          call hypre_matrix_update(A, 1, ncol, tmp_rows, tmp_cols, tmp_vals)
       end do
    end if

    call hypre_matrix_assemble(A)
  end subroutine hypre_matrix_assemble_from_neko

  subroutine count_and_fill_rows(A_rows, rows, nrows, nnz)
    !DIR$ INLINENEVER count_and_fill_rows
    integer, intent(in) :: nnz
    integer, intent(out) :: nrows
    integer, intent(in) :: A_rows(:)
    integer, intent(inout) :: rows(:)
    integer :: e, s, i, idof
    logical :: in_rows
    idof = 0
    nrows = 0
    rows(:) = 0
    do e = 1, nnz
       in_rows = .false.
       i = A_rows(e)
       do s = 1, idof
          if (i .eq. rows(s)) then
             in_rows = .true.
          end if
       end do
       if (.not. in_rows) then
          idof = idof + 1
          nrows = nrows + 1
          rows(nrows) = i
       end if
    end do
  end subroutine count_and_fill_rows

  subroutine count_and_fill_cols_vals(A_rows, A_cols, A_vals, rows, cols, nrows, ncols, vals, nnz)
    !DIR$ INLINENEVER count_and_fill_cols_vals
    integer, intent(in) :: nnz, nrows
    integer, intent(in) :: A_rows(:), A_cols(:), rows(:)
    real(kind=rp), intent(in) :: A_vals(:)
    integer, intent(inout) :: cols(:), ncols(:)
    real(kind=rp), intent(inout) :: vals(:)
    integer :: k, e, idof
    idof = 0
    ncols(:) = 0
    do k = 1, nrows
       do e = 1, nnz
          if (A_rows(e) .eq. rows(k)) then
             ncols(k) = ncols(k) + 1
             idof = idof + 1
             cols(idof) = A_cols(e)
             vals(idof) = A_vals(e)
          end if
       end do
    end do
  end subroutine count_and_fill_cols_vals

end module hypre
