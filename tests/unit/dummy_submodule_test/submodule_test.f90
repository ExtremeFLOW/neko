submodule (module_test) module_test_interface
  use dummy_use, only : dummy_use_subroutine
  implicit none

contains

  module subroutine module_test_subroutine()
    print *, "Hello from submodule"
    call dummy_use_subroutine()
  end subroutine module_test_subroutine


end submodule module_test_interface
