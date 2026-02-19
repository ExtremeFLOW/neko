module hypre
  use num_types, only : rp, i8
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_free, device_memcpy, HOST_TO_DEVICE
  use hypre_boomeramg
  use hypre_ij_interface
  use coefs, only : coef_t
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

  type, public :: hypre_solver_t
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
     procedure, pass(this) :: setup => hypre_solver_setup
     procedure, pass(this) :: solve => hypre_solve
     procedure, pass(this) :: device_solve => hypre_device_solve
     procedure, pass(this) :: set_matrix
     procedure, pass(this) :: set_vector
     !procedure, pass(this) :: free
     procedure, pass(this) :: set_dofs => hypre_dofs_workaround
  end type hypre_solver_t

  public :: hypre_init, hypre_fin, hypre_matrix_assemble_from_neko

contains

  subroutine hypre_init()
     integer :: ierr
     ! Initialize hypre
     ierr = HYPRE_Initialize()
     ! Set other hypre config parameters
     !call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE, ierr)
     !call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE, ierr)
     !call HYPRE_SetSpGemmUseVendor(0, ierr)
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
  subroutine hypre_solver_setup(this)
    class(hypre_solver_t), intent(inout) :: this
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
  !! @param dof_i8 Array of global DoFs on rank. Fortran array on host.
  subroutine hypre_solve(this, x, f, n, dof_i8)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(in) :: f
    integer(kind=i8), dimension(n), intent(in) :: dof_i8
    integer, dimension(n) :: dof
    ! convert to i4 as hypre takes int*
    dof = dof_i8
    ! Copy to hypre vector
    ! (copy to x may be unneeded if zero initial guess is always used)
    call hypre_copy_to_vector(this%x, n, dof, x)
    call hypre_copy_to_vector(this%b, n, dof, f)
    ! Solve
    call boomeramg_solve(this%solver, this%parcsr_A, this%par_b, this%par_x)
    ! Copy from hypre vector
    call hypre_copy_from_vector(this%x, n, dof, x)
  end subroutine hypre_solve

  !> Solve system for right hand side f
  !! @param x Approximated solution. c_ptr on device.
  !! @param f Right hand side. c_ptr on device.
  !! @param n Size of vectors
  !! @param dof Array of global DoFs on rank. c_ptr on device.
  subroutine hypre_device_solve(this, x, f, n, dof)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    type(c_ptr), intent(in) :: x
    type(c_ptr), intent(in) :: f
    type(c_ptr), intent(in) :: dof
    ! Copy to hypre vector
    ! (copy to x may be unneeded if zero initial guess is always used)
    call hypre_device_copy_to_vector(this%x, n, dof, x)
    call hypre_device_copy_to_vector(this%b, n, dof, f)
    ! Solve
    call boomeramg_solve(this%solver, this%parcsr_A, this%par_b, this%par_x)
    ! Copy from hypre vector
    call hypre_device_copy_from_vector(this%x, n, dof, x)
  end subroutine hypre_device_solve

  subroutine hypre_dofs_workaround(this, coef)
    class(hypre_solver_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: i
    allocate(this%dofs(size(coef%dof%dof)))
    do i = 1, size(coef%dof%dof)
       this%dofs(i) = coef%dof%dof(i,1,1,1)
    end do
    print *, "DOF workaround mapping"
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%dofs, this%dofs_d, size(coef%dof%dof))
       call device_memcpy(this%dofs, this%dofs_d, size(coef%dof%dof), HOST_TO_DEVICE, .true.)
    end if
  end subroutine hypre_dofs_workaround

  subroutine hypre_matrix_assemble_from_neko(coef, A, b, x)
    type(coef_t), intent(in), target :: coef
    type(c_ptr), intent(inout) :: A, b, x
    real(kind=rp), allocatable :: A_vals(:)
    real(kind=rp) :: tmp_vals(1)
    integer, allocatable :: A_rows(:), A_cols(:)
    integer :: tmp_rows(1), tmp_cols(1), ncol(1)
    real(kind=rp) :: tmp2
    integer :: ilower, iupper, jlower, jupper
    integer :: e, k, j, i, s, l
    integer :: nnz, nelv, lx, idof
    type(c_ptr) :: tmp_rows_d, tmp_cols_d, tmp_vals_d, ncol_d
    tmp_rows_d = C_NULL_PTR
    tmp_cols_d = C_NULL_PTR
    tmp_vals_d = C_NULL_PTR
    ncol_d = C_NULL_PTR
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
               end do
            end do
         end do
      end do ! e = 1, n
    end associate
    ! Get index range of dofs on current rank
    ! TODO: For MPI parallelism, duplicated dofs need to be assigned
    !       to a single "owning" rank, which is currently not supported
    !       by the current neko dofmap
    ilower = minval(A_rows)
    iupper = maxval(A_rows)
    jlower = minval(A_cols)
    jupper = maxval(A_cols)

    ! Initialize matrix
    call hypre_matrix_init(1, ilower, iupper, jlower, jupper, A)

    ! Initialize vectors
    call hypre_vector_init(1, ilower, iupper, b)
    call hypre_vector_init(1, jlower, jupper, x)

    ! Fill matrix
    ! We do this the lazy and expensive way, one dof at a time.
    ! The interface expects array, so we fill some dummy arrays of size 1
    ! and pass these to the hypre interface
    print *, "MATRIX MAPPING"
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(tmp_rows, tmp_rows_d, 1)
       call device_map(tmp_cols, tmp_cols_d, 1)
       call device_map(tmp_vals, tmp_vals_d, 1)
       call device_map(ncol, ncol_d, 1)
       print *, "MAPS DONE"
       do i = 1, nnz
         ncol(1) = 1
         tmp_rows(1) = A_rows(i)
         tmp_cols(1) = A_cols(i)
         tmp_vals(1) = A_vals(i)
         call device_memcpy( tmp_rows, tmp_rows_d, 1, HOST_TO_DEVICE, .true.)
         call device_memcpy( tmp_cols, tmp_cols_d, 1, HOST_TO_DEVICE, .true.)
         call device_memcpy( tmp_vals, tmp_vals_d, 1, HOST_TO_DEVICE, .true.)
         call device_memcpy( ncol, ncol_d, 1, HOST_TO_DEVICE, .true.)
         call hypre_device_matrix_update(A, 1, ncol_d, tmp_rows_d, tmp_cols_d, tmp_vals_d)
       end do
    else
       do i = 1, nnz
         ncol(1) = 1
         tmp_rows(1) = A_rows(i)
         tmp_cols(1) = A_cols(i)
         tmp_vals(1) = A_vals(i)
         call hypre_matrix_update(A, 1, ncol, tmp_rows, tmp_cols, tmp_vals)
       end do
    end if

    print *, "MATRIX ASSEMBLE"
    call hypre_matrix_assemble(A)
  end subroutine hypre_matrix_assemble_from_neko

end module hypre
