module hypre
  use, intrinsic :: iso_c_binding
  implicit none
  private

  include 'HYPRE.h'

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

  subroutine hypre_finalize()
     integer :: ierr
     ierr = HYPRE_Finalize()
  end subroutine hypre_finalize

end module hypre
