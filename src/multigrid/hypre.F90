module hypre
  use num_types, only : rp, i8
  use hypre_boomeramg
  use hypre_ij_interface
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
   contains
     procedure, pass(this) :: init => hypre_solver_init
     procedure, pass(this) :: setup => hypre_solver_setup
     procedure, pass(this) :: solve => hypre_solve
     procedure, pass(this) :: set_matrix
     procedure, pass(this) :: set_vector
     !procedure, pass(this) :: free
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

  subroutine hypre_solve(this, x, f, n, dof_i8)
    class(hypre_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), dimension(n), intent(in) :: f
    integer(kind=i8), dimension(n), intent(in) :: dof_i8
    integer, dimension(n) :: dof
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

end module hypre
