module global_params
  use logger
  use num_types

  implicit none

  real(kind=rp) :: fst_ti = 3.1d-2 * 6.0! turbulence intensity
  real(kind=rp) :: fst_il = 3.2d-2! integral length scale		

  integer, parameter :: nshells = 80! No of spherical shells
  integer, parameter :: Npmax = 20  ! Npmax  -  Number of points in a shell
  real(kind=rp), parameter :: kstart = 43.77d0 ! smallest wavenumber
  real(kind=rp), parameter :: kend = 1660d0 ! largest wavenumber

  integer, parameter :: fst_modes = 2 * nshells * Npmax ! No of freestream modes  
  integer :: shell_modes(nshells) ! Modes saved per shell    

  real(kind=rp) :: kmin = kstart*0.5
  real(kind=rp) :: kmax = kend*2.0

  real(kind=rp) :: bb(fst_modes, 3)  ! RENAME THESE!!! So confusing

  integer :: k_length
  integer :: shell(fst_modes)

  real(kind=rp) :: k_num(fst_modes, 3)
  real(kind=rp) :: k_num_all(fst_modes, 3)
  real(kind=rp) :: u_hat_pn(fst_modes, 3)
  real(kind=rp) :: shell_amp(nshells)

  logical :: write_files = .true.

contains

  ! Use to print one parameter
  subroutine print_param(name, value)
    implicit none

    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, *) name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_param

  ! Use to print one parameter
  subroutine print_int(name, value)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, *) name, ": ", value
    call neko_log%message(log_buf)

  end subroutine print_int

  subroutine gparams_print_params()
    character(len=LOG_SIZE) :: log_buf

    call print_param("nshells", real(nshells, kind=rp))
    call print_param("Npmax", real(Npmax, kind=rp))
    call print_param("fst_ti", fst_ti)
    call print_param("fst_il", fst_il)

  end subroutine gparams_print_params

  real(kind=rp) FUNCTION ran2(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL(kind=rp) AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
  end function ran2
end module global_params
