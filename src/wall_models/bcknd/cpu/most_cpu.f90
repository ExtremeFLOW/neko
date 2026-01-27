module most_cpu
  use num_types, only : rp
  use utils, only : neko_error
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  public :: most_compute_cpu

  ! Note:
  ! I actually tried to avoid going into the neutral regime altogether (i.e. take out the if abs(Ri)<1e-3 part)
  ! AI tells me it's not safe/correct because I approach neutrality only asymptotically, and usign slaw_m/h which were 
  ! obtained by assuming z/L=0.01. However, the results look better in this way (IMO), so I have to understand if that's 
  ! or not. 
  ! By doing it without neutral, I actually get an appreciable difference between roughlog and most, and I also get the
  ! uw flux to go exactly zero above the inversion. I also get a very slightly different v-profile than Nek though

  abstract interface
     function slaw_m_interface(z,L_ob,z0) result(slaw)
      import rp
      real(kind=rp), intent(in) :: z, L_ob, z0
      real(kind=rp) :: slaw
     end function slaw_m_interface

     function slaw_h_interface(z,L_ob,z0h) result(slaw)
      import rp
      real(kind=rp), intent(in) :: z, L_ob, z0h
      real(kind=rp) :: slaw
     end function slaw_h_interface

     function corr_m_interface(z,L_ob) result(corr)
      import rp
      real(kind=rp), intent(in) :: z, L_ob
      real(kind=rp) :: corr
     end function corr_m_interface

     function corr_h_interface(z,L_ob) result(corr)
      import rp
      real(kind=rp), intent(in) :: z, L_ob
      real(kind=rp) :: corr
     end function corr_h_interface

     function f_interface(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
      import rp, slaw_m_interface, slaw_h_interface
      real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
      real(kind=rp) :: f
      procedure(slaw_m_interface) :: slaw_m
      procedure(slaw_h_interface) :: slaw_h
     end function f_interface

     function dfdl_interface(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
      import rp, slaw_m_interface, slaw_h_interface
      real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
      real(kind=rp) :: dfdl
      procedure(slaw_m_interface) :: slaw_m
      procedure(slaw_h_interface) :: slaw_h
     end function dfdl_interface
  end interface

  !===============================
  ! 3. Procedure pointers
  !===============================
  !
  ! These will point to the correct functions
  ! depending on stability regime and bc_type.
  !

  procedure(slaw_m_interface), pointer :: slaw_m_ptr => null()
  procedure(slaw_h_interface), pointer :: slaw_h_ptr => null()
  procedure(corr_m_interface), pointer :: corr_m_ptr => null()
  procedure(corr_h_interface), pointer :: corr_h_ptr => null()
  procedure(f_interface), pointer :: f_ptr => null()
  procedure(dfdl_interface), pointer :: dfdl_ptr => null()

contains

  !================================================
  ! 4. User-facing selectors for stability regime
  !================================================

  subroutine select_bc_operators(bc_type)
    character(len=*), intent(in) :: bc_type

    select case (bc_type)
    case ("neumann")
      f_ptr    => f_neumann
      dfdl_ptr => dfdl_neumann
    case ("dirichlet")
      f_ptr    => f_dirichlet
      dfdl_ptr => dfdl_dirichlet
    case default
      call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine select_bc_operators

  subroutine compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in)    :: g, hi, ti, ts, magu, kappa
    real(kind=rp), intent(inout) :: q, Ri_b

    select case (bc_type)
    case ("neumann")
      Ri_b = - g*hi / ti*q / (magu**3*kappa**2)

    case ("dirichlet")
      Ri_b = g*hi/ti*(ti - ts)/magu**2

    case default
      call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine compute_Ri_b

  subroutine init_q(bc_type, hi, ti, ts, kappa, utau, z0h, q)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in)    :: hi, ti, ts, kappa, utau, z0h
    real(kind=rp), intent(inout) :: q

    if (bc_type == "dirichlet") then
      q = kappa*utau*(ts - ti)/log(hi/z0h)
    end if
  end subroutine init_q

  subroutine set_stability_regime(Ri_b)
    real(kind=rp), intent(in) :: Ri_b

    if (Ri_b > 0.001) then
      slaw_m_ptr => slaw_m_stable
      slaw_h_ptr => slaw_h_stable
      corr_m_ptr => corr_m_stable
      corr_h_ptr => corr_h_stable
    elseif (Ri_b < -0.001) then
      slaw_m_ptr => slaw_m_convective
      slaw_h_ptr => slaw_h_convective
      corr_m_ptr => corr_m_convective
      corr_h_ptr => corr_h_convective
    else
      slaw_m_ptr => slaw_m_neutral
      slaw_h_ptr => slaw_h_neutral
    end if
  end subroutine set_stability_regime

  !================================================
  ! 5. Main wall-flux routine
  !================================================

  subroutine most_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &  
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, z0, bc_type, zone_idx, h_idx, q, tstep)  ! q in input only temporarily
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: q   ! should this be multidimensional?
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer, intent(in) :: zone_idx  ! WARNING: only supports wall model on ONE boundary atm!
    integer :: ts_idx(3)
    integer, intent(in) :: h_idx
    integer :: i, count
    integer, parameter :: max_count = 20
    real(kind=rp) :: ui, vi, ti, ts, hi
    real(kind=rp) :: magu, utau, normu, z0h
    real(kind=rp) :: L_ob, L_upper, L_lower, L_old
    real(kind=rp) :: Ri_b, f, dfdl, fd_h, L_new, L_sign    
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: tol = 0.001_rp
    real(kind=rp), parameter :: NR_step = 0.001_rp
    character(len=LOG_SIZE) :: log_buf

    ! debug only:
    ! ts  = 300.0_rp
    ! q = 0.05_rp
    ! call neko_registry%add_field(this%coef%dof, "sampling_height", &
    !      ignore_existing=.true.)
    ! h_field => neko_registry%get_field_by_name("sampling_height")

    ! Select the ts offset based on fid
    select case (zone_idx)
      case (1)
        ts_idx = [h_idx, 0, 0 ]
      case (2)
        ts_idx = [-h_idx, 0, 0]
      case (3)
        ts_idx = [0, h_idx, 0 ]
      case (4)
        ts_idx = [0, -h_idx, 0]
      case (5)
        ts_idx = [0, 0, h_idx ]
      case (6)
        ts_idx = [0, 0, -h_idx]
      case default
        call neko_error("The face index is not correct (most_cpu.f90)")
    end select

    do i=1, n_nodes
   
      ! Sample the variables
      ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      ts = temp(ind_r(i)-ts_idx(1), ind_s(i)-ts_idx(2), ind_t(i)-ts_idx(3), ind_e(i))
      hi = h(i)

      ! Project on horizontal directions
      normu = ui * n_x(i) + vi * n_y(i)
      ui = ui - normu * n_x(i) 
      vi = vi - normu * n_y(i)

      ! Compute velocity magnitude
      magu = sqrt(ui**2 + vi**2) 

      ! utau initialisation
      if (tstep < 1) then
        utau = sqrt( sqrt( tau_x(i)**2 + tau_y(i)**2 ) ) 
      else 
        utau = magu*kappa / log(hi/z0)
   
      end if
      ! Compute thermal roughness length from Zilitinkevich 1995
      z0h = z0 * exp(-0.1_rp*sqrt((utau*z0)/1.46e-5_rp))  

      ! Get q, Ri_b, f_ptr, dfdl_ptr based on bc_type 
      ! Maybe redundant, but needed to initialise Rib
      call select_bc_operators(bc_type)
      call init_q(bc_type, hi, ti, ts, kappa, utau, z0h, q)
      call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)

      ! Compute Obukhov length
      if (tstep > 0) then
        ! Get q, Ri, f, dfdl based on bc_type 
        call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)

        if (abs(Ri_b) <= 1.0e-3_rp) then
          ! Neutral (L_ob undefined)
          L_ob = 1.0e+10_rp
        else
          ! Determine target regime sign
          if (Ri_b > 0.0_rp) then
            L_ob = hi / max(Ri_b, 0.01_rp)  ! Stable guess
            L_sign = 1.0_rp
          else
            L_ob = hi / min(Ri_b, -0.01_rp) ! Convective guess
            L_sign = -1.0_rp
          end if

          L_old = 1.0e10_rp
          count = 0
          ! Set slaw and corr pointers based on stability
          call set_stability_regime(Ri_b)

          ! Find Obukhov length
          do while  ((abs(L_old - L_ob)/abs(L_ob) > tol) .and. (count < max_count))
            ! Switch between stable and convective based on bulk Richardson (Ri_b)
            L_old = L_ob
            count = count + 1

            fd_h = NR_step*L_ob
            L_upper = L_ob + fd_h
            L_lower = L_ob - fd_h

            ! Compute L_ob based on stability and bc_type
            if (.not. associated(f_ptr) .or. .not. associated(dfdl_ptr)) then
              call neko_error("Unassociated pointer for f or dfdl")
            end if

            f = f_ptr(Ri_b, hi, z0, z0h, L_ob, slaw_m_ptr, slaw_h_ptr)
            dfdl = dfdl_ptr(l_upper, l_lower, hi, z0, z0h, L_ob, slaw_m_ptr, slaw_h_ptr, fd_h)
            if (abs(dfdl) < 1.0e-12_rp) call neko_error("Division by zero in dfdl")
            L_new = L_ob - f/dfdl

            !!! Optional bounding
            ! ! Avoid regime crossing during Newton iter
            ! if (L_new*L_sign <= 0.0_rp) then
            !   ! "damp update" (stay on same side)
            !   L_new = 0.5_rp * L_ob
            ! end if

            ! ! Bound L_ob
            ! L_ob = sign(max(abs(L_new), 1.0e-4_rp), L_sign)
            ! L_ob = sign(min(abs(L_ob), 1.0e6_rp), L_sign)
          end do

          if (abs(L_ob) > 5e4_rp .or. abs(L_ob) < 1e-5_rp) then
            count = max_count
            call neko_error("Obukhov length did not converge (MOST wall model)") 
          end if
        end if

        ! Based on stability and bc_type, compute utau/q 
        call set_stability_regime(Ri_b) 
        select case (bc_type)
          case ("neumann")
            ! Compute u* with the new Obukhov length
            utau = kappa*magu/slaw_m_ptr(hi, L_ob, z0)
          case ("dirichlet")
            ! Compute u* with the new Obukhov length
            utau = kappa*magu/slaw_m_ptr(hi, L_ob, z0)
            ! and compute q from here
            q = kappa*utau*(ts - ti)/slaw_h_ptr(L_ob, hi, z0h)  ! z0h placeholder
          case default
            call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
        end select 

      end if

      ! Distribute according to the velocity vector and bound magu to avoid 0 division
      magu = max(magu, 1.0e-6_rp)
      tau_x(i) = -utau**2 * ui / magu
      tau_y(i) = -utau**2 * vi / magu
      tau_z(i) = 0
    end do

    call neko_log%section('Wall model quick look')
    write(log_buf, '(A,E15.7)') 'Ri_b: ', Ri_b
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'L_ob: ', L_ob
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'utau: ', utau
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'magu: ', magu
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'ts: ', ts
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'ti: ', ti
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'q: ', q
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'hi: ', hi
    call neko_log%message(trim(log_buf))
    call neko_log%end_section()

  end subroutine most_compute_cpu

  !================================================
  ! 6. Concrete similarity functions
  !================================================

  !--------------- Stable ----------------

  function slaw_m_stable(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z,L_ob,z0
    real(kind=rp) :: slaw

    slaw = log(z/z0)-corr_m_stable(z,L_ob)+corr_m_stable(z0,L_ob)
  end function slaw_m_stable

  function slaw_h_stable(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z,L_ob,z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h)-corr_h_stable(z,L_ob)+corr_h_stable(z0h,L_ob)
  end function slaw_h_stable

  function corr_m_stable(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z,L_ob
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d  ! parameter
    real(kind=rp) :: zeta

    zeta = z/L_ob
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = - a*zeta - b*(zeta-c/d)*exp(-d*zeta) - b*c/d
  end function corr_m_stable

  function corr_h_stable(z,L_ob) result(corr)    ! identical to corr_m_stable...?
    real(kind=rp), intent(in) :: z,L_ob
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d
    real(kind=rp) :: zeta

    zeta = z/L_ob
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = -b * (zeta-c/d)*exp(-d*zeta)-(1.0_rp+ 0.6666666_rp * a * zeta)**1.5_rp-b*c/d + 1.0_rp
  end function corr_h_stable

  !--------------- Convective ----------------

  function slaw_m_convective(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0
    real(kind=rp) :: slaw

    slaw = log(z/z0) - corr_m_convective(z, L_ob) + corr_m_convective(z0, L_ob)
  end function slaw_m_convective

  function slaw_h_convective(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h) - corr_h_convective(z, L_ob) + corr_h_convective(z0h, L_ob)
  end function slaw_h_convective

  function corr_m_convective(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z, L_ob
    real(kind=rp) :: xi, pi, zeta
    real(kind=rp) :: corr

    zeta = z/L_ob
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi)) + log(0.5_rp*(1 + xi**2)) - 2*atan(xi) + pi/2
  end function corr_m_convective

  function corr_h_convective(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z, L_ob
    real(kind=rp) :: zeta, pi, xi
    real(kind=rp) :: corr

    zeta = z/L_ob
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi**2))
  end function corr_h_convective

  !--------------- Neutral ----------------

  function slaw_m_neutral(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0
    real(kind=rp) :: slaw

    slaw = log(z/z0)
  end function slaw_m_neutral

  function slaw_h_neutral(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h)
  end function slaw_h_neutral

  !------------- Similarity laws --------------

  function f_neumann(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - z/L_ob/slaw_m(z, L_ob, z0)**3) 
  end function f_neumann

  function dfdl_neumann(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-z/l_upper/slaw_m(z, l_upper, z0)**3)  ! conv
    dfdl = dfdl + (z/l_lower/slaw_m(z, l_lower, z0)**3)  ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_neumann

  function f_dirichlet(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - z/L_ob*slaw_h(z, L_ob, z0h)/slaw_m(z, L_ob, z0)**2)  ! conv
  end function f_dirichlet

  function dfdl_dirichlet(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-z/l_upper*slaw_h(z, l_upper, z0h)/slaw_m(z, l_upper, z0)**2)  ! conv
    dfdl=dfdl + (z/l_lower*slaw_h(z, l_lower, z0h)/slaw_m(z, l_lower, z0)**2)  ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_dirichlet
  

end module most_cpu
