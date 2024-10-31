program test
  use module_test, only : module_test_subroutine
  implicit none

  print *, "Hello from main program"

  call module_test_subroutine()

end program test
