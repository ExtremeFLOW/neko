module user
  use neko
  use ieee_arithmetic, only: ieee_is_nan

  implicit none

  integer, parameter :: num_rows = 128 ! resolution of TS wave in Cheb points
  integer, parameter :: num_columns = 7 ! y, real and imag part of TS wave in Cheb points
  integer, parameter :: num_ygll = 8 * 8 ! number of GLL points in the y direction
  real (kind=rp) :: Re = 5000.0_rp
  real (kind=rp), parameter :: pi_rp = 4.0_rp * atan (1.0_rp)
contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(usr)
    type(user_t), intent(inout) :: usr

    usr%user_mesh_setup => user_mesh_scale
    usr%fluid_user_ic => user_ic
    
  end subroutine user_setup


  ! Rescale mesh
  ! Original mesh size: (2.0, 2.0, 1.0).
  ! New mesh can easily be genreated with genmeshbox.
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i

    do i = 1, size(msh%points)
       msh%points(i)%x(1) = pi_rp * msh%points(i)%x(1) / 1.12_rp
       msh%points(i)%x(3) = pi_rp * (msh%points(i)%x(3) - 0.5_rp) &
                                  * 2.0_rp / 2.1_rp
    end do

  end subroutine user_mesh_scale

  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, j, i_y
    real(kind=rp) :: uvw(3)

    real(kind=rp), dimension(num_rows, num_columns) :: data_mode_cheb_2D
    real(kind=rp), dimension(num_rows, num_columns) :: data_mode_cheb_3D
    real(kind=rp) :: y_GLC(num_rows)
    real(kind=rp) :: TS2D_GLC(num_rows, num_columns-1)
    real(kind=rp) :: TS3D_GLC(num_rows, num_columns-1)

    real(kind=rp), dimension(num_ygll) :: y_GLL
    real(kind=rp) :: TS2D_GLL(num_ygll, num_columns-1)
    real(kind=rp) :: TS3D_GLL(num_ygll, num_columns-1)

    real(kind=rp) :: ur_2D, ui_2D, vr_2D, vi_2D
    real(kind=rp) :: ur_3D, ui_3D, vr_3D, vi_3D, wr_3D, wi_3D

    ! data reading
    open(unit = 10, file = 'TSwave_cheb_2D.bin', form = 'unformatted', &
                  access = 'stream')
       read(10) data_mode_cheb_2D
    close(10)
    y_GLC = data_mode_cheb_2D(:,1)
    TS2D_GLC = data_mode_cheb_2D(:,2:num_columns)

    open(unit = 10, file = 'TSwave_cheb_3D.bin', form = 'unformatted', &
                  access = 'stream')
       read(10) data_mode_cheb_3D
    close(10)
    TS3D_GLC = data_mode_cheb_3D(:,2:num_columns)

    ! alternative of point zone
    i_y = 1
    
    ! initialize y_GLL
    do i = 1, num_ygll
       y_GLL(i) = 0.0
    end do

    do i = 1, u%dof%size()
       if (.not. in_array(u%dof%y(i,1,1,1), y_GLL)) then
          y_GLL(i_y) = u%dof%y(i,1,1,1)
          i_y = i_y +1
       end if
    end do

    do i = 1, num_columns-1
       do j = 1, num_ygll
          TS2D_GLL(j,i) = GLC_GLL_interp(TS2D_GLC(:,i), y_GLC, y_GLL(j))
          TS3D_GLL(j,i) = GLC_GLL_interp(TS3D_GLC(:,i), y_GLC, y_GLL(j))
       end do
    end do

    do i = 1, u%dof%size()
       ur_2D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS2D_GLL(:,1))
       ui_2D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS2D_GLL(:,2))
       vr_2D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS2D_GLL(:,3))
       vi_2D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS2D_GLL(:,4))
       ur_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,1))
       ui_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,2))
       vr_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,3))
       vi_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,4))
       wr_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,5))
       wi_3D = pick_pt(u%dof%y(i,1,1,1), y_GLL, TS3D_GLL(:,6))
       uvw = channel_ic(u%dof%x(i,1,1,1), u%dof%y(i,1,1,1), u%dof%z(i,1,1,1), &
                        ur_2D, ui_2D, vr_2D, vi_2D, &
                        ur_3D, ui_3D, vr_3D, vi_3D, wr_3D, wi_3D)
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do

  end subroutine user_ic

  function in_array(y, y_list) result(is_in)
    logical :: is_in
    real(kind=rp) :: y, y_list(num_ygll)
    integer :: i
    real(kind=rp) :: tol

    tol = 1e-4_rp
    is_in = .false.

    do i = 1, num_ygll
       if (relcmp( y, y_list(i), tol)) then
          is_in = .true.
          exit
       end if
    end do

  end function in_array

  function pick_pt(y_target, y_source, Pt_source) result(Pt_target)
    real(kind=rp) :: y_target, Pt_target
    real(kind=rp), dimension(num_ygll) :: y_source, Pt_source
    real(kind=rp) :: tol
    integer :: i
    logical :: found
    real(kind=rp) :: dist
    
    tol = 1e-4_rp
    found = .false.
    dist = 100
    do i = 1, num_ygll
       if (relcmp(y_target, y_source(i), tol)) then
          Pt_target = Pt_source(i)
          dist = abs(y_target-y_source(i))
       end if
    end do
    if ( dist .le. tol ) found = .true.
    
    if (.not. found) then
       write (*,*) 'tolerence too small for picking points. dist:', dist
    end if

  end function pick_pt

  function channel_ic(x, y, z, &
                      ur_2D, ui_2D, vr_2D, vi_2D, &
                      ur_3D, ui_3D, vr_3D, vi_3D, wr_3D, wi_3D) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ub

    real(kind=rp) :: ur_2D, ui_2D, vr_2D, vi_2D
    real(kind=rp) :: ur_3D, ui_3D, vr_3D, vi_3D, wr_3D, wi_3D
    real(kind=rp) :: TS_amp_2D, TS_amp_3D, alpha, beta

    complex(kind=rp) :: spa_osci_2D
    complex(kind=rp) :: u_mode_2D, v_mode_2D
    real(kind=rp) :: u_pert_TS_2D, v_pert_TS_2D

    complex(kind=rp) :: u_mode_3D, v_mode_3D, w_mode_3D
    complex(kind=rp) :: spa_osci_3D_p, spa_osci_3D_n
    real(kind=rp) :: u_pert_TS_3D, v_pert_TS_3D, w_pert_TS_3D
    
    ! amplitude for decaying TS wave:
    ! TS_amp_2D = 1e-6_rp
    ! TS_amp_3D = 0.0_rp

    ! amplitude for secondary instability TS wave:
    TS_amp_2D = 3e-2_rp
    TS_amp_3D = 5e-4_rp

    alpha = 1.12_rp
    beta = 2.1_rp

    spa_osci_2D = exp((0.0_rp,1.0_rp)*alpha*x)
    u_mode_2D = cmplx(ur_2D, ui_2D, rp)
    v_mode_2D = cmplx(vr_2D, vi_2D, rp)

    u_pert_TS_2D = TS_amp_2D * real(u_mode_2D * spa_osci_2D, rp)
    v_pert_TS_2D = TS_amp_2D * real(v_mode_2D * spa_osci_2D, rp)

    spa_osci_3D_p = exp((0.0_rp,1.0_rp)*alpha*x+(0.0_rp,1.0_rp)*beta*z)
    spa_osci_3D_n = exp((0.0_rp,1.0_rp)*alpha*x-(0.0_rp,1.0_rp)*beta*z)
    u_mode_3D = cmplx(ur_3D, ui_3D, rp)
    v_mode_3D = cmplx(vr_3D, vi_3D, rp)
    w_mode_3D = cmplx(wr_3D, wi_3D, rp)

    u_pert_TS_3D = TS_amp_3D*real(u_mode_3D*spa_osci_3D_p, rp)
    u_pert_TS_3D = u_pert_TS_3D + TS_amp_3D*real(u_mode_3D*spa_osci_3D_n, rp)
    v_pert_TS_3D = TS_amp_3D*real(v_mode_3D*spa_osci_3D_p, rp)
    v_pert_TS_3D = v_pert_TS_3D + TS_amp_3D*real(v_mode_3D*spa_osci_3D_n, rp)
    w_pert_TS_3D = TS_amp_3D*real(w_mode_3D*spa_osci_3D_p, rp)
    w_pert_TS_3D = w_pert_TS_3D + TS_amp_3D*real(-w_mode_3D*spa_osci_3D_n, rp)

    ub = 1.0_rp-y*y

    uvw(1) = ub + u_pert_TS_2D + u_pert_TS_3D
    uvw(2) = v_pert_TS_2D + v_pert_TS_3D
    uvw(3) = w_pert_TS_3D

  end function channel_ic

  ! The raw data in .csv file is on GLC points,
  ! hence a interpolation from GLC points to GLL points is performed.
  function GLC_GLL_interp(f_GLC, x_GLC, x_GLL) result(f_GLL)
    real(kind=rp) :: x_GLL, f_GLL
    real(kind=rp), dimension(num_rows) :: f_GLC, x_GLC
    real(kind=rp), dimension(num_rows) :: chi, psi
    real(kind=rp) :: xn, xs, xe, xnc, dchebyshev
    integer :: i, N, cj
    logical :: is_nan

    N = num_rows - 1
    xn = x_GLL
    xs = x_GLC(1)
    xe = x_GLC(N+1)
    xnc = ( 2.0_rp*xn - (xs+xe) )/ &
          ( xe - xs)
    dchebyshev = eval_dchebyshev(xnc,N)
    do i = 1, N+1
       chi(i) = ( 2.0_rp*x_GLC(i) - (xs+xe) )/ &
                ( xe - xs)

       if (i .eq. 1) then
          cj = 2.0_rp
       else if (i .eq. N+1) then
          cj = 2.0_rp
       else
          cj = 1.0_rp
       end if

       psi(i) = ((-1.0_rp)**(i) * (1.0_rp-xnc*xnc) * dchebyshev)/ &
                (cj * N*N * (xnc-chi(i)))

       is_nan = ieee_is_nan(psi(i))
       if (is_nan .or. psi(i) .gt. 1e+3) then
          psi(i) = 1.0_rp
       end if

    end do

    f_GLL = 0.0_rp
    do i = 1, N+1
       f_GLL = f_GLL + psi(i)*f_GLC(i)
    end do

  end function GLC_GLL_interp

  ! Evaluate the derivative of the Chebyshev polynomials of the first kind
  function eval_dchebyshev(x,N) result(dT_Nx)
    real(kind=rp) :: x, dT_Nx, tmp
    integer :: N, j, nn ! N should be num_rows - 1
    real(kind=rp), dimension(N+1) :: dT

    dT(1) = 0.0_rp
    dT(2) = 1.0_rp

    do j = 3, N+1
       nn = j - 2
       tmp = eval_chebyshev(x, nn)
       if (nn .eq. 1) then
          dT(j) = 2.0_rp*(nn+1.0_rp)*tmp
       else
          dT(j) = 2.0_rp*(nn+1.0_rp)*tmp + (nn+1.0_rp)/(nn-1.0_rp)*dT(j-2)
       end if
    end do

    dT_Nx = dT(N+1)
  end function eval_dchebyshev
  
  ! Evaluate the Chebyshev polynomials of the first kind
  function eval_chebyshev(x,N) result(T_Nx)
    real(kind=rp) :: x, T_Nx
    integer :: N, j ! N should be num_rows - 1
    real(kind=rp), dimension(N+1) :: T
    T(1) = 1.0_rp
    T(2) = x
    do j = 3, N+1
       T(j) = 2.0_rp*x*T(j-1) - T(j-2)
    end do
    T_Nx = T(N+1)
  end function eval_chebyshev

end module user
