!> Simulation driver
module simulation
  use case
  implicit none
contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), intent(inout) :: C
    real(kind=dp) :: t
    integer :: i

    t = 0d0
    do i = 1, C%params%nsteps
       call C%fluid%step()
       t = t + C%params%dt
    end do
    
  end subroutine neko_solve

end module simulation
