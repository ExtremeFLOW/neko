module hypre_boomeramg
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! Common parameter constants
  integer(c_int), parameter :: HYPRE_AMG_COARSEN_FALGOUT = 6
  integer(c_int), parameter :: HYPRE_AMG_COARSEN_PMIS = 8
  integer(c_int), parameter :: HYPRE_AMG_COARSEN_HMIS = 10
  integer(c_int), parameter :: HYPRE_AMG_INTERP_EXTENDED = 6
  integer(c_int), parameter :: HYPRE_AMG_INTERP_EXT_PLUS_I = 14
  integer(c_int), parameter :: HYPRE_AMG_RELAX_JACOBI = 0
  integer(c_int), parameter :: HYPRE_AMG_RELAX_GS_SEQ = 3
  integer(c_int), parameter :: HYPRE_AMG_RELAX_GS_PAR = 4
  integer(c_int), parameter :: HYPRE_AMG_RELAX_HYBRID_GS = 6

  ! BoomerAMG solver creation and destruction
  interface
     integer (c_int) function HYPRE_BoomerAMGCreate(solver) &
           bind(c, name='HYPRE_BoomerAMGCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: solver
     end function HYPRE_BoomerAMGCreate
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGDestroy(solver) &
          bind(c, name='HYPRE_BoomerAMGDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
     end function HYPRE_BoomerAMGDestroy
  end interface

  interface
     ! Setup and solve
     integer (c_int) function HYPRE_BoomerAMGSetup(solver, A, b, x) &
          bind(c, name='HYPRE_BoomerAMGSetup')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       type(c_ptr), value :: A, b, x
     end function HYPRE_BoomerAMGSetup
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSolve(solver, A, b, x) &
          bind(c, name='HYPRE_BoomerAMGSolve')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       type(c_ptr), value :: A, b, x
     end function HYPRE_BoomerAMGSolve
  end interface

  ! Parameters
  interface
     integer (c_int) function HYPRE_BoomerAMGSetTol(solver, tol) &
          bind(c, name='HYPRE_BoomerAMGSetTol')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       real(c_double), value :: tol
     end function HYPRE_BoomerAMGSetTol
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetMaxIter(solver, max_iter) &
          bind(c, name='HYPRE_BoomerAMGSetMaxIter')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: max_iter
     end function HYPRE_BoomerAMGSetMaxIter
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetPrintLevel(solver, print_level) &
          bind(c, name='HYPRE_BoomerAMGSetPrintLevel')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: print_level
     end function HYPRE_BoomerAMGSetPrintLevel
  end interface

  ! Coarsening and interpolation
  interface
     integer (c_int) function HYPRE_BoomerAMGSetCoarsenType(solver, coarsen_type) &
          bind(c, name='HYPRE_BoomerAMGSetCoarsenType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: coarsen_type
     end function HYPRE_BoomerAMGSetCoarsenType
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetInterpType(solver, interp_type) &
          bind(c, name='HYPRE_BoomerAMGSetInterpType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: interp_type
     end function HYPRE_BoomerAMGSetInterpType
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetMaxLevels(solver, max_levels) &
           bind(c, name='HYPRE_BoomerAMGSetMaxLevels')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: max_levels
     end function HYPRE_BoomerAMGSetMaxLevels
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetStrongThreshold(solver, strong_threshold) &
          bind(c, name='HYPRE_BoomerAMGSetStrongThreshold')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       real(c_double), value :: strong_threshold
     end function HYPRE_BoomerAMGSetStrongThreshold
  end interface

  ! Relaxation (smoothing)
  interface
     integer (c_int) function HYPRE_BoomerAMGSetRelaxType(solver, relax_type) &
          bind(c, name='HYPRE_BoomerAMGSetRelaxType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: relax_type
     end function HYPRE_BoomerAMGSetRelaxType
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetNumSweeps(solver, num_sweeps) &
          bind(c, name='HYPRE_BoomerAMGSetNumSweeps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: num_sweeps
     end function HYPRE_BoomerAMGSetNumSweeps
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetCycleType(solver, cycle_type) &
          bind(c, name='HYPRE_BoomerAMGSetCycleType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: cycle_type
     end function HYPRE_BoomerAMGSetCycleType
  end interface

  ! Get diagnostic information
  interface
     integer (c_int) function HYPRE_BoomerAMGGetNumIterations(solver, num_iterations) &
          bind(c, name='HYPRE_BoomerAMGGetNumIterations')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int) :: num_iterations
     end function HYPRE_BoomerAMGGetNumIterations
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, rel_resid_norm) &
          bind(c, name='HYPRE_BoomerAMGGetFinalRelativeResidualNorm')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       real(c_double) :: rel_resid_norm
     end function HYPRE_BoomerAMGGetFinalRelativeResidualNorm
  end interface

  ! Advanced options
  interface
     integer (c_int) function HYPRE_BoomerAMGSetMaxRowSum(solver, max_row_sum) &
          bind(c, name='HYPRE_BoomerAMGSetMaxRowSum')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       real(c_double), value :: max_row_sum
     end function HYPRE_BoomerAMGSetMaxRowSum
  end interface

  interface
     integer (c_int) function HYPRE_BoomerAMGSetAggNumLevels(solver, agg_num_levels) &
          bind(c, name='HYPRE_BoomerAMGSetAggNumLevels')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: solver
       integer(c_int), value :: agg_num_levels
     end function HYPRE_BoomerAMGSetAggNumLevels
  end interface

end module hypre_boomeramg
