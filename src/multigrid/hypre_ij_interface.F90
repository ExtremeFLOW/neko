module hypre_ij_interface
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! HYPRE types
  integer(c_int), parameter :: HYPRE_PARCSR = 5555

  ! Matrix creation and assembly
  interface
     integer (c_int) function HYPRE_IJMatrixCreate(ilower, iupper, jlower, jupper, matrix) &
          bind(c, name='HYPRE_IJMatrixCreate_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: ilower, iupper, jlower, jupper
       type(c_ptr) :: matrix
     end function HYPRE_IJMatrixCreate
  end interface
  !!interface
  !!   integer (c_int) function HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, matrix) &
  !!        bind(c, name='HYPRE_IJMatrixCreate')
  !!     use, intrinsic :: iso_c_binding
  !!     implicit none
  !!     integer(c_int), value :: comm
  !!     integer(c_int), value :: ilower, iupper, jlower, jupper
  !!     type(c_ptr) :: matrix
  !!   end function HYPRE_IJMatrixCreate
  !!end interface

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
     integer (c_int) function HYPRE_IJMatrixAddToValues(matrix, nrows, ncols, rows, cols, values) &
          bind(c, name='HYPRE_IJMatrixAddToValues')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: matrix
       integer(c_int), value :: nrows
       type(c_ptr), value :: ncols
       type(c_ptr), value :: rows
       type(c_ptr), value :: cols
       type(c_ptr), value :: values
     end function HYPRE_IJMatrixAddToValues
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
     integer (c_int) function HYPRE_IJVectorCreate(jlower, jupper, vector) &
          bind(c, name='HYPRE_IJVectorCreate_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: jlower, jupper
       type(c_ptr) :: vector
     end function HYPRE_IJVectorCreate
  end interface
  !!interface
  !!   integer (c_int) function HYPRE_IJVectorCreate(comm, jlower, jupper, vector) &
  !!        bind(c, name='HYPRE_IJVectorCreate')
  !!     use, intrinsic :: iso_c_binding
  !!     implicit none
  !!     integer(c_int), value :: comm
  !!     integer(c_int), value :: jlower, jupper
  !!     type(c_ptr) :: vector
  !!   end function HYPRE_IJVectorCreate
  !!end interface

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

  interface
     integer (c_int) function HYPRE_IJVectorPrint(vector, filename) &
          bind(c, name='HYPRE_IJVectorPrint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vector
       type(c_ptr), value :: filename
     end function HYPRE_IJVectorPrint
  end interface

  public :: hypre_matrix_init, hypre_matrix_fill, hypre_matrix_assemble, hypre_matrix_update, hypre_matrix_destroy
  public :: hypre_device_matrix_update
  public :: hypre_vector_init, hypre_vector_fill, hypre_vector_assemble, hypre_vector_destroy
  public :: hypre_copy_from_vector, hypre_copy_to_vector
  public :: hypre_device_copy_from_vector, hypre_device_copy_to_vector
  public :: hypre_matrix_get_object, hypre_vector_get_object

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
     !ierr = HYPRE_IJMatrixCreate(mpi_comm, ilower, iupper, jlower, jupper, A)
     ierr = HYPRE_IJMatrixCreate(ilower, iupper, jlower, jupper, A)
     ! Choose a parallel csr format storage (see hypre's user manual)
     ierr = HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR)
     ! Initialize before setting coefficients
     ierr = HYPRE_IJMatrixInitialize(A)
  end subroutine hypre_matrix_init

  !> Fill hypre matrix from neko. Overwrites existing matrix entries.
  !! nrows may be larger than nuber of rows on rank
  !! due to duplicated DoFs in values.
  !! @param A c_ptr for HYPRE_IJMatrix
  !! @param nrows Number of rows to be set. Counts duplicates.
  !! @param ncols Number of cols on each row. Fortran array on host.
  !! @param rows Array of row indices. global DoFs. Fortran array on host.
  !! @param cols Array of col indices. global DoFs. Fortran array on host.
  !! @param values Array of values to be set. Fortran array on host.
  subroutine hypre_matrix_fill(A, nrows, ncols, rows, cols, values)
     type(c_ptr), intent(in) :: A
     integer, target, intent(in) :: nrows, ncols(:)
     integer, target, intent(in) :: rows(:), cols(:)
     double precision, target, intent(in) :: values(:)
     integer :: ierr
     ! Fill the matrix with values.
     ierr = HYPRE_IJMatrixSetValues(A, nrows, c_loc(ncols), &
         c_loc(rows), c_loc(cols), c_loc(values))
  end subroutine hypre_matrix_fill

  !> Add to hypre matrix from neko. Creates new entries or adds to existing.
  !! nrows may be larger than number of rows on rank
  !! due to duplicated DoFs in values.
  !! @param A c_ptr for HYPRE_IJMatrix
  !! @param nrows Number of rows to be set. Counts duplicates.
  !! @param ncols Number of cols on each row. Fortran array on host.
  !! @param rows Array of row indices. global DoFs. Fortran array on host.
  !! @param cols Array of col indices. global DoFs. Fortran array on host.
  !! @param values Array of values to be set. Fortran array on host.
  subroutine hypre_matrix_update(A, nrows, ncols, rows, cols, values)
     !DIR$ INLINENEVER hypre_matrix_update
     type(c_ptr), intent(in) :: A
     integer, intent(in) :: nrows
     integer, target, intent(in) :: ncols(:)
     integer, target, intent(in) :: rows(:)
     integer, target, intent(in) :: cols(:)
     double precision, target, intent(in) :: values(:)
     integer :: ierr
     ! Fill the matrix with values.
     ierr = HYPRE_IJMatrixAddToValues(A, nrows, c_loc(ncols), &
         c_loc(rows), c_loc(cols), c_loc(values))
  end subroutine hypre_matrix_update

  !> Add to hypre matrix from neko. Creates new entries or adds to existing.
  !! nrows may be larger than number of rows on rank
  !! due to duplicated DoFs in values.
  !! @param A c_ptr for HYPRE_IJMatrix
  !! @param nrows Number of rows to be set. Counts duplicates.
  !! @param ncols Number of cols on each row. c_ptr on device
  !! @param rows Array of row indices. global DoFs. c_ptr on device.
  !! @param cols Array of col indices. global DoFs. c_ptr on device.
  !! @param values Array of values to be set. c_ptr on device.
  subroutine hypre_device_matrix_update(A, nrows, ncols, rows, cols, values)
     type(c_ptr), intent(in) :: A
     integer, intent(in) :: nrows
     type(c_ptr), intent(in) :: ncols, rows, cols
     type(c_ptr), intent(in) :: values
     integer :: ierr
     ! Fill the matrix with values.
     ierr = HYPRE_IJMatrixAddToValues(A, nrows, ncols, &
         rows, cols, values)
  end subroutine hypre_device_matrix_update

  !> Finish and assemble matrix
  !! @param A c_ptr for HYPRE_IJMatrix
  subroutine hypre_matrix_assemble(A)
     type(c_ptr), intent(in) :: A
     integer :: ierr
     ! Assemble after setting the coefficients
     ierr = HYPRE_IJMatrixAssemble(A)
  end subroutine hypre_matrix_assemble

  subroutine hypre_matrix_get_object(A, parcsr_A)
     type(c_ptr), intent(in) :: A
     type(c_ptr), intent(inout) :: parcsr_A
     integer :: ierr
     ! Get parcsr matrix pointer (for solver use)
     ierr = HYPRE_IJMatrixGetObject(A, parcsr_A)
  end subroutine hypre_matrix_get_object

  subroutine hypre_matrix_destroy(A)
     type(c_ptr), intent(inout) :: A
     integer :: ierr
     ierr = HYPRE_IJMatrixDestroy(A)
  end subroutine hypre_matrix_destroy

  !-----------------------------------------------------------------------------
  ! Vector
  !-----------------------------------------------------------------------------

  !> Interface to initialize a vector for hypre
  !! @param mpi_comm Not used.
  !! @param jlower Smallest global DoF on rank.
  !! @param jupper Largest global DoF on rank.
  !! @param v HYPRE_IJVector c_ptr
  subroutine hypre_vector_init(mpi_comm, jlower, jupper, v)
    integer, intent(in) :: mpi_comm, jlower, jupper
    type(c_ptr), intent(inout) :: v
    integer :: ierr
    ! Create the vector
    !ierr = HYPRE_IJVectorCreate(mpi_comm, jlower, jupper, v)
    ierr = HYPRE_IJVectorCreate(jlower, jupper, v)
    ierr = HYPRE_IJVectorSetObjectType(v, HYPRE_PARCSR)
    ierr = HYPRE_IJVectorInitialize(v)
  end subroutine hypre_vector_init

  !> Fill hypre vector from neko.
  !! nvalues may be larger than vector length
  !! due to duplicated DoFs in indices/values.
  !! @param v c_ptr for HYPRE_IJVector
  !! @param nvalues Length of indicies/values.
  !! @param indices Array of global DoFs on rank. Fortran array on host.
  !! @param values Array of values to be set. Fortran array on host.
  subroutine hypre_vector_fill(v, nvalues, indices, values)
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues
    integer, target, intent(in) :: indices(:)
    double precision, target, intent(in) :: values(:)
    integer :: ierr
    ! Set vector values
    ierr = HYPRE_IJVectorSetValues(v, nvalues, c_loc(indices), c_loc(values))
  end subroutine hypre_vector_fill

  !> Finalze and assemble the HYPRE_IJVector object
  !! @param v c_ptr HYPRE_IJVector
  subroutine hypre_vector_assemble(v)
    type(c_ptr), intent(in) :: v
    integer :: ierr
    ! Finalize and assemble vector object
    ierr = HYPRE_IJVectorAssemble(v)
  end subroutine hypre_vector_assemble

  subroutine hypre_vector_get_object(v, par_v)
    type(c_ptr), intent(in) :: v
    type(c_ptr), intent(inout) :: par_v
    integer :: ierr
    ierr = HYPRE_IJVectorGetObject(v, par_v)
  end subroutine hypre_vector_get_object

  !> Copy vector from neko to hypre.
  !! nvalues may be larger than vector length
  !! due to duplicated DoFs in indices/values.
  !! @param v c_ptr for HYPRE_IJVector
  !! @param nvalues Length of indicies/values.
  !! @param indices Array of global DoFs on rank. Fortran array on host.
  !! @param values Array of values to be set. Fortran array on host.
  subroutine hypre_copy_to_vector(v, nvalues, indices, values, ilower)
    !DIR$ INLINENEVER hypre_copy_to_vector
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues, ilower
    integer, target, intent(in) :: indices(nvalues)
    double precision, target, intent(in) :: values(nvalues)
    integer :: ierr, i
    ! re-initialize vector
    ierr = HYPRE_IJVectorInitialize(v)
    ! Set vector values
    !ierr = HYPRE_IJVectorSetValues(v, nvalues, c_loc(indices), c_loc(values))
    do i = 1, nvalues
      if (indices(i) .ge. ilower) then
        ierr = HYPRE_IJVectorSetValues(v, 1, c_loc(indices(i)), c_loc(values(i)))
      end if
    end do
  end subroutine hypre_copy_to_vector

  !> Copy vector from neko to hypre.
  !! nvalues may be larger than vector length
  !! due to duplicated DoFs in indices/values.
  !! @param v c_ptr for HYPRE_IJVector
  !! @param nvalues Length of indicies/values.
  !! @param indices Array of global DoFs on rank. c_ptr on device.
  !! @param values Array of values to be set. c_ptr on device.
  subroutine hypre_device_copy_to_vector(v, nvalues, indices, values)
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues
    type(c_ptr), intent(in) :: indices
    type(c_ptr), intent(in) :: values
    integer :: ierr
    ! re-initialize vector
    ierr = HYPRE_IJVectorInitialize(v)
    ! Set vector values
    ierr = HYPRE_IJVectorSetValues(v, nvalues, indices, values)
  end subroutine hypre_device_copy_to_vector

  !> Copy vector from hypre to neko.
  !! nvalues may be larger than vector length
  !! due to duplicated DoFs in indices/values.
  !! @param v c_ptr for HYPRE_IJVector
  !! @param nvalues Length of indicies/values.
  !! @param indices Array of global DoFs on rank. Fortran array on host.
  !! @param values Array of values to be set. Fortran array on host.
  subroutine hypre_copy_from_vector(v, nvalues, indices, values, ilower)
    !DIR$ INLINENEVER hypre_copy_from_vector
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues, ilower
    integer, target, intent(in) :: indices(nvalues)
    double precision, target, intent(inout) :: values(nvalues)
    integer :: ierr, i
    ! Get vector values
    !ierr = HYPRE_IJVectorGetValues(v, nvalues, c_loc(indices), c_loc(values))
    do i = 1, nvalues
      if (indices(i) .ge. ilower) then
        ierr = HYPRE_IJVectorGetValues(v, 1, c_loc(indices(i)), c_loc(values(i)))
      else
        values(i) = 0.0
      end if
    end do
  end subroutine hypre_copy_from_vector

  !> Copy vector from hypre to neko.
  !! nvalues may be larger than vector length
  !! due to duplicated DoFs in indices/values.
  !! @param v c_ptr for HYPRE_IJVector
  !! @param nvalues Length of indicies/values.
  !! @param indices Array of global DoFs on rank. c_ptr on device.
  !! @param values Array of values to be set. c_ptr on device.
  subroutine hypre_device_copy_from_vector(v, nvalues, indices, values)
    type(c_ptr), intent(in) :: v
    integer, intent(in) :: nvalues
    type(c_ptr), intent(in) :: indices
    type(c_ptr), intent(in) :: values
    integer :: ierr
    ! Get vector values
    ierr = HYPRE_IJVectorGetValues(v, nvalues, indices, values)
  end subroutine hypre_device_copy_from_vector

  !> Free the hypre vector object
  subroutine hypre_vector_destroy(v)
     type(c_ptr), intent(inout) :: v
     integer :: ierr
     ierr = HYPRE_IJVectorDestroy(v)
  end subroutine hypre_vector_destroy

end module hypre_ij_interface
