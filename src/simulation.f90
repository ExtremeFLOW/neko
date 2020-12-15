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
    type(file_t) :: fileout

    t = 0d0
    fileout = file_t("oufluid.fld")
    do i = 1, C%params%nsteps
       tstep = i
       call C%fluid%step(t, tstep)
       if (mod(i, 3) .eq. 0) then
          write(*,*) 'Save'
          call fileout%write(C%fluid)
       end if

       !> @todo Add call to sampler
!       t = t + C%params%dt !< @todo Re-enable once settime is fixed
    end do

    
  end subroutine neko_solve

end module simulation
