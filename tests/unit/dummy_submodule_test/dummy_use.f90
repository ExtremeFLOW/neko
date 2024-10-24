module dummy_use
  implicit none

  public :: dummy_use_subroutine

contains

  module subroutine dummy_use_subroutine()
    print *, "Hello from dummy used subroutine"
  end subroutine dummy_use_subroutine


end module dummy_use
