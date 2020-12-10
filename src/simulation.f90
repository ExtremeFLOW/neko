!> Simulation driver
module simulation
  use case
  implicit none
contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), intent(inout) :: C
    real(kind=dp) :: t
    integer :: i, tstep

    t = 0d0
    do i = 1, C%params%nsteps
       tstep = i
       call C%fluid%step(t, tstep)
       !> @todo Add call to sampler
       ! t = t + C%params%dt !< @todo Re-enable once settime is fixed
    end do
    
  end subroutine neko_solve

end module simulation
