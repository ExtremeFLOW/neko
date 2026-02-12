module hypre_ij_interface
  use, intrinsic :: iso_c_binding
  implicit none
  private

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

end module hypre_ij_interface
