! * Free stream turbulence implementation for neko
! * Author: Victor Baconnet (baconnet@kth.se)
! * Based in original implementation by Elektra Kluesberg, Prabal Negi
! 	and Philipp Schlatter. 
! * For FST theory see Schlatter (2001)

module FST

  use global_params
  use turbu

  implicit none

  type, public :: FST_t

     ! Attributes
     real(kind=rp) :: xmin
     real(kind=rp) :: xmax
     real(kind=rp) :: xstart
     real(kind=rp) :: xend
     real(kind=rp) :: delta_rise
     real(kind=rp) :: delta_fall
     real(kind=rp) :: fringe_max
     real(kind=rp) :: uinf

   contains
     procedure, pass(this) :: init => FST_init
     procedure, pass(this) :: generate => FST_generate
     procedure, pass(this) :: FST_forcing
     procedure, pass(this) :: FST_forcing_zone
     procedure, pass(this) :: free => FST_free_params
     procedure, pass(this) :: print => FST_print_params
     generic :: apply_forcing => FST_forcing, FST_forcing_zone
  end type FST_t

contains

  ! This routine intializes FST parameters.
  ! TODO: Initialize parameters from case file with default user params
  ! or create new parameters
  subroutine FST_init(this, params, xstart, xend, xmin, xmax, &
       fringe_max, delta_rise, delta_fall)
    class(FST_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    real(kind=rp), intent(in) :: xstart
    real(kind=rp), intent(in) :: xend
    real(kind=rp), intent(in) :: xmin
    real(kind=rp), intent(in) :: xmax
    real(kind=rp), intent(in) :: fringe_max
    real(kind=rp), intent(in) :: delta_rise
    real(kind=rp), intent(in) :: delta_fall

    call neko_log%section('Initializing FST')

    call params%get("case.fluid.initial_condition.value(1)", this%uinf)

    this%xstart = xstart!-0.925   !-0.185
    this%xend   = xend!-0.825   !-0.165
    this%xmin   = xmin!-0.94    !-0.188
    this%xmax   = xmax!-0.815   !-0.163
    this%fringe_max = fringe_max!0.2d0/36.0d0 * 350.
    this%delta_rise = delta_rise!0.002
    this%delta_fall = delta_fall!0.002

    call this%print() ! show parameters

    call neko_log%end_section('Done --> Intializing FST')

  end subroutine FST_init

  !! Free parameters in global params
  subroutine FST_free_params(this)
    class(FST_t), intent(inout) :: this

  end subroutine FST_free_params

  subroutine FST_print_params(this)
    class(FST_t) :: this

    call print_param("uinf", this%uinf)
    call print_param("xstart", this%xstart)
    call print_param("xend", this%xend)
    call print_param("xmin", this%xmin)
    call print_param("xmax", this%xmax)
    call print_param("fringe_max", this%fringe_max)
    call print_param("delta_rise", this%delta_rise)
    call print_param("delta_fall", this%delta_fall)
!!$    call print_param("fst_ti", this%fst_ti)
!!$    call print_param("fst_il", this%fst_il
!!$    call print_param("kstart", this%kstart)
!!$    call print_int  ("nshell", this%nshells)
!!$    call print_int  ("npmax", this%npmax)
!!$    call print_param("kstart", this%kstart)
!!$    call print_param("kend", this%kend)

  end subroutine FST_print_params

  ! Step 2 in the FST generation process, this computes
  ! all the required parameters and variables such as the
  ! wavenumbers, amplitudes, etc.
  subroutine FST_generate(this, u)

    implicit none

    class(FST_t), intent(inout) :: this
    type(field_t), intent(in) :: u
    character(len=LOG_SIZE) :: log_buf

    integer :: ierr

    call neko_log%section ('Generating FST')

    call neko_log%message("calling make_turbu")
    call make_turbu(u)
    call neko_log%message("done make_turbu")

    call MPI_Bcast(k_length , 1                   , &
         MPI_INTEGER         , 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_num_all, fst_modes*u%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(k_num    , fst_modes*u%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(u_hat_pn , fst_modes*u%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(bb       , fst_modes*u%msh%gdim, &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)
    call MPI_Bcast(shell    , fst_modes           , &
         MPI_INTEGER         , 0, NEKO_COMM, ierr)
    call MPI_Bcast(shell_amp, nshells             , &
         MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)

    call neko_log%message("FST - spectrum generation done.")

    call neko_log%end_section('Done --> Generating FST')

  end subroutine FST_generate

  ! Apply the forcing on every element
  subroutine FST_forcing(this, f, u, v, w, t)
    class(FST_t), intent(inout) :: this
    class(fluid_user_source_term_t), intent(inout) :: f
    type(field_t), intent(inout) :: u, v, w
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: tempfam

    integer :: i

    do i=1, f%dm%size()
       call FST_forcing_local(this, f%dm%x(i,1,1,1), &
            f%dm%y(i,1,1,1), f%dm%z(i,1,1,1), t, &
            u%x(i,1,1,1), v%x(i,1,1,1), w%x(i,1,1,1), &
            f%u(i,1,1,1), f%v(i,1,1,1), f%w(i,1,1,1), & 
            tempfam, i)
       f%chi(i,1,1,1) = tempfam
       !f%chi(i,1,1,1) = 0.0
    end do

  end subroutine FST_forcing
  
  ! Apply the forcing on every element
  subroutine FST_forcing_zone(this, f, u, v, w, t, zone)
    class(FST_t), intent(inout) :: this
    class(fluid_user_source_term_t), intent(inout) :: f
    type(field_t), intent(inout) :: u, v, w
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: tempfam
    class(point_zone_t), intent(inout) :: zone

    integer :: i, idx

    do idx=1, zone%size
       i = zone%mask(idx)
       call FST_forcing_local(this, f%dm%x(i,1,1,1), &
            f%dm%y(i,1,1,1), f%dm%z(i,1,1,1), t, &
            u%x(i,1,1,1), v%x(i,1,1,1), w%x(i,1,1,1), &
            f%u(i,1,1,1), f%v(i,1,1,1), f%w(i,1,1,1), & 
            tempfam, i)
       f%chi(i,1,1,1) = tempfam
       !f%chi(i,1,1,1) = 0.0
    end do

  end subroutine FST_forcing_zone

  ! Forcing to be performed on entire domain, on a local element ix
  ! Final values of the forcing are to be applied
  ! on f%u, f%v and f%w
  subroutine FST_forcing_local(FST_obj, &
       x, y, z, t, &
       u, v, w, fu, fv, fw, chi, idx)
    implicit none

    class(FST_t), intent(in) :: FST_obj
    real(kind=rp), intent(in) :: x, y, z, t
    real(kind=rp), intent(inout) :: u, v, w
    real(kind=rp), intent(inout) :: fu, fv, fw, chi
    integer, intent(in) :: idx

    integer :: l, m, i, shellno
    integer, parameter :: gdim = 3
    real(kind=rp) :: phase_shft, phi, amp, pert
    real(kind=rp) :: rand_vec(gdim)

!!$    print *, x

    if (x .ge. FST_obj%xmin .and. x .le. FST_obj%xmax) then

!!$       print *, x, "DOING IT"

       !	initialize rand_vec
       do l = 1,gdim
          rand_vec(l)= 0.0
       enddo

       !! compute m sinus modes
       do m=1,k_length
          do i=1,gdim 

             phase_shft = bb(m,1)

             phi = k_num_all(m,1)*x + &
                  k_num_all(m,2)*y + &
                  k_num_all(m,3)*z + &
                  FST_obj%uinf * k_num_all(m,1)*t + &
                  phase_shft

             shellno = shell(m)
             amp = shell_amp(shellno)
             pert = u_hat_pn(m,i)*amp*sin(phi)

             rand_vec(i) = rand_vec(i) + pert
          enddo
       enddo

       ! print *, "fringe is ", fringe(x,FST_obj)
       ! call FST_obj%print_params()


       ! ok now we split between implicit and explicit parts
       ! THIS IS EXPLICIT
       !fu = fringe(x,y,FST_obj) * (FST_obj%uinf + rand_vec(1) - u)
       !fv = fringe(x,y,FST_obj) * (rand_vec(2) - v)
       !fw = fringe(x,y,FST_obj) * (rand_vec(3) - w)
       !chi = 0.0

		 ! THIS IS IMPLICIT
       fu = fringe(x,y,FST_obj) * (FST_obj%uinf + rand_vec(1))
       fv = fringe(x,y,FST_obj) * (rand_vec(2) )
       fw = fringe(x,y,FST_obj) * (rand_vec(3) )
       chi = fringe(x,y,FST_obj) 

    else
       fu = 0.0_rp
       fv = 0.0_rp
       fw = 0.0_rp
       chi = 0.0_rp
    end if

  end subroutine FST_forcing_local

  ! Fringe function as defined in original FST.
  function fringe(x,yloc, f) result(y)
    real(kind=rp), intent(in) :: x
    type(FST_t) :: f
    real(kind=rp) :: y, S_s, S_e, x_s, x_e, yloc

	 ! Victor victor victor...
	 x_s = (x-f%xstart)/f%delta_rise
	 x_e = (x-f%xend)/f%delta_fall + 1.0
    if(x_s.le.0) then
       S_s=0
    else if (x_s.ge.1) then
       S_s=1
    else
       S_s=1/(1+exp(1/(x_s-1)+1/(x_s)))
    endif

    if(x_e.le.0) then
       S_e=0
    else if (x_e.ge.1) then
       S_e=1
    else
       S_e=1/(1+exp(1/(x_e-1)+1/(x_e)))
    endif

    y = f%fringe_max * (S_s - S_e)

  end function fringe

!!$  ! Fringe function as described in Schlatter (2001)
!!$  function fringe(x, xstart, xend, delta_rise, delta_fall, fringe_max) result(y)
!!$    implicit none
!!$
!!$    real(kind=rp), intent(in) :: x, xstart, xend, delta_rise, delta_fall, fringe_max
!!$    real(kind=rp) :: y
!!$    y = fringe_max * ( S((x-xstart)/delta_rise) - S((x-xend)/delta_fall + 1d0) )
!!$  end function fringe
!!$
!!$  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/exp(1/(x-1) + 1/x) between 0 and 1
!!$  function S(x) result(y)
!!$    implicit none
!!$
!!$    real(kind=rp), intent(in) :: x
!!$    real(kind=rp)             :: y
!!$    if ( x.le.0d0 ) then
!!$       y = 0d0
!!$    else if ( x.ge.1d0 ) then
!!$       y = 1d0
!!$    else
!!$       y = 1d0 / exp( 1d0/(x-1d0) + 1d0/x)
!!$    end if
!!$    return
!!$  end function S


end module FST
