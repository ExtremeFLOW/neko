module module_test
  use dummy_use, only : dummy_use_subroutine
  implicit none
  private

  public :: module_test_subroutine

  interface
     module subroutine module_test_subroutine()
     end subroutine module_test_subroutine

  end interface

end module module_test
