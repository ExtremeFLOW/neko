module hypre_ij_interface
  use, intrinsic :: iso_c_binding
  implicit none
  private

  include 'HYPRE.h'

  ! Matrix creation and assembly
  interface
     integer (c_int) function HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, matrix) &
          bind(c, name='HYPRE_IJMatrixCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: comm
       integer(c_int), value :: ilower, iupper, jlower, jupper
       type(c_ptr) :: matrix
     end function HYPRE_IJMatrixCreate
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixDestroy(matrix) &
          bind(c, name='HYPRE_IJMatrixDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
     end function HYPRE_IJMatrixDestroy
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixSetObjectType(matrix, obj_type) &
          bind(c, name='HYPRE_IJMatrixSetObjectType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
       integer(c_int), value :: obj_type
     end function HYPRE_IJMatrixSetObjectType
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixInitialize(matrix) &
          bind(c, name='HYPRE_IJMatrixInitialize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
     end function HYPRE_IJMatrixInitialize
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixSetValues(matrix, nrows, ncols, rows, cols, values) &
          bind(c, name='HYPRE_IJMatrixSetValues')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
       integer(c_int), value :: nrows
       type(c_ptr), value :: ncols
       type(c_ptr), value :: rows
       type(c_ptr), value :: cols
       type(c_ptr), value :: values
     end function HYPRE_IJMatrixSetValues
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixAssemble(matrix) &
          bind(c, name='HYPRE_IJMatrixAssemble')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
     end function HYPRE_IJMatrixAssemble
  end interface

  interface
     integer (c_int) function HYPRE_IJMatrixGetObject(matrix, object) &
          bind(c, name='HYPRE_IJMatrixGetObject')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
       type(c_ptr) :: object
     end function HYPRE_IJMatrixGetObject
  end interface

  ! Vector creation and assembly
  interface
     integer (c_int) function HYPRE_IJVectorCreate(comm, jlower, jupper, vector) &
          bind(c, name='HYPRE_IJVectorCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: comm
       integer(c_int), value :: jlower, jupper
       type(c_ptr) :: vector
     end function HYPRE_IJVectorCreate
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorDestroy(vector) &
          bind(c, name='HYPRE_IJVectorDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
     end function HYPRE_IJVectorDestroy
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorSetObjectType(vector, obj_type) &
          bind(c, name='HYPRE_IJVectorSetObjectType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
       integer(c_int), value :: obj_type
     end function HYPRE_IJVectorSetObjectType
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorInitialize(vector) &
          bind(c, name='HYPRE_IJVectorInitialize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
     end function HYPRE_IJVectorInitialize
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorSetValues(vector, nvalues, indices, values) &
          bind(c, name='HYPRE_IJVectorSetValues')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
       integer(c_int), value :: nvalues
       type(c_ptr), value :: indices
       type(c_ptr), value :: values
     end function HYPRE_IJVectorSetValues
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorGetValues(vector, nvalues, indices, values) &
          bind(c, name='HYPRE_IJVectorGetValues')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
       integer(c_int), value :: nvalues
       type(c_ptr), value :: indices
       type(c_ptr), value :: values
     end function HYPRE_IJVectorGetValues
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorAssemble(vector) &
          bind(c, name='HYPRE_IJVectorAssemble')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
     end function HYPRE_IJVectorAssemble
  end interface

  interface
     integer (c_int) function HYPRE_IJVectorGetObject(vector, object) &
          bind(c, name='HYPRE_IJVectorGetObject')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
       type(c_ptr) :: object
     end function HYPRE_IJVectorGetObject
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Matrix
  !-----------------------------------------------------------------------------

  !> Initialize a hypre matrix
  !! Each processor knows only of its own rows
  !! - the range is denoted by ijlower and ijupper.
  subroutine hypre_matrix_init(mpi_comm, ilower, iupper, jlower, jupper, A)
     integer, intent(in) :: mpi_comm, ilower, iupper, jlower, jupper
     type(c_ptr), intent(inout) :: A
     integer :: ierr
     !TODO: may need 1->0 shift on indexing
     ! Create the matrix object
     ierr = HYPRE_IJMatrixCreate(mpi_comm, ilower, iupper, jlower, jupper, A)
     ! Choose a parallel csr format storage (see hypre's user manual)
     ierr = HYPRE_IJ_MatrixSetObjectType(A, HYPRE_PARCSR)
     ! Initialize before setting coefficients
     ierr = HYPRE_IJMatrixInitialize(A)
  end subroutine hypre_matrix_init

  subroutine hypre_matrix_fill(A, nrows, ncols, rows, cols, values)
     type(c_ptr), intent(in) :: A
     integer, intent(in) :: nrows, ncols(:)
     integer, intent(in) :: rows(:), cols(:)
     double precision, intent(in) :: values(:)
     integer :: ierr
     ! Fill the matrix with values.
     ierr = HYPRE_IJMatrixSetValues(A, nrows, c_loc(ncols), &
         c_loc(rows), c_loc(cols), c_loc(values))
  end subroutine hypre_matrix_fill

  subroutine hypre_matrix_assemble(A, parcsr_A)
     type(c_ptr), intent(in) :: A
     type(c_ptr), intent(inout) :: parcsr_A
     integer :: ierr
     ! Assemble after setting the coefficients
     ierr = HYPRE_IJMatrixAssemble(A)
     ! Get parcsr matrix pointer (for solver use)
     ierr = HYPRE_IJMatrixGetObject(A, parcsr_A)
  end subroutine hypre_matrix_assemble

  subroutine hypre_matrix_destroy(A)
     type(c_ptr), intent(inout) :: A
     integer :: ierr
     ierr = HYPRE_IJMatrixDestroy(A)
  end subroutine hypre_matrix_destroy

  !-----------------------------------------------------------------------------
  ! Vector
  !-----------------------------------------------------------------------------

  !> Interface to initialize a vector for hypre
  !! mpi_comm
  !! jlower
  !! jupper
  !! v
  subroutine hypre_vector_init(mpi_comm, jlower, jupper, v)
    integer, intent(in) :: mpi_comm, jlower, jupper
    type(c_ptr), intent(inout) :: v
    integer :: ierr
    !TODO: indexing may need 1->0 shift
    ! Create the vector
    ierr = HYPRE_IJVectorCreate(mpi_comm, jlower, jupper, v)
    ierr = HYPRE_IJVectorSetObjectType(v, HYPRE_PARCSR)
    ierr = HYPRE_IJVectorInitialize(v)
  end subroutine hypre_vector_init

  subroutine hypre_vector_fill(v, nvalues, indices, values)
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues
    integer, intent(in) :: indices(:)
    double precision, intent(in) :: values(:)
    integer :: ierr
    ! Set vector values
    ierr = HYPRE_IJVectorSetValues(v, nvalues, c_loc(indices), c_loc(values))
  end subroutine hypre_vector_fill

  subroutine hypre_vector_assemble(v, par_v)
    type(c_ptr), intent(in) :: v
    type(c_ptr), intent(inout) :: par_v
    integer :: ierr
    ! Finalize and assemble vector object
    ierr = HYPRE_IJVectorAssemble(v)
    ierr = HYPRE_IJVectorGetObject(v, par_v)
  end subroutine hypre_vector_assemble

  subroutine copy_to_vector(v, par_v, nvalues, indices, values)
    type(c_ptr), intent(in) :: v
    type(c_ptr), intent(inout) :: par_v
    integer, intent(in) :: nvalues
    integer, intent(in) :: indices(:)
    double precision, intent(in) :: values(:)
    integer :: ierr
    ! re-initialize vector
    ierr = HYPRE_IJVectorInitialize(v)
    ! Set vector values
    ierr = HYPRE_IJVectorSetValues(v, nvalues, c_loc(indices), c_loc(values))
    ! Finalize and assemble vector object
    ierr = HYPRE_IJVectorAssemble(v)
    ierr = HYPRE_IJVectorGetObject(v, par_v)
  end subroutine copy_to_vector

  subroutine hypre_vector_destroy(v)
     type(c_ptr), intent(inout) :: v
     integer :: ierr
     ierr = HYPRE_IJVectorDestroy(A)
  end subroutine hypre_vector_destroy

end module hypre_ij_interface
