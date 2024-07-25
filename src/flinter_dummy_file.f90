!> @file flinter_dummy_file.f90
!! @brief Dummy file for testing the flinter

! Attempt to test the flinter with a dummy file

! Anti-Test: Indentation
program main
  call subroutine
end program main

! Anti-Test: Spaces inside strings
print *, "This  is a test"
print *, 'This  is a test'

! Tests
call foo(a, b)
write (*, '(A,A)') "some  string"
write (*, '(A,A)') "some  string"
write (*, '(A,A)') 'some  string'
write (*, '(A,A)') 'some  string'
write (*, '(A,A)') 'some  string' "test"

