module most_cpu
  use num_types, only : rp
  use utils, only : neko_error
  implicit none
  private

  public :: most_compute_cpu

  ! TO DO:
  ! - where/when should I sample ts?  and z0t?
  ! - create a neutral branch if Rib==0

  abstract interface
     function slaw_m_interface(h,L_o,z0) result(slaw)
        real(kind=rp), intent(in) :: h, L_o, z0
        real(kind=rp) :: slaw
     end function slaw_m_interface

     function slaw_h_interface(h,L_o,z0t) result(slaw)
        real(kind=rp), intent(in) :: h, L_o, z0t
        real(kind=rp) :: slaw
     end function slaw_h_interface

     function corr_m_interface(h,L_o) result(corr)
        real(kind=rp), intent(in) :: h, L_o
        real(kind=rp) :: corr
     end function corr_m_interface

     function corr_h_interface(h,L_o) result(corr)
        real(kind=rp), intent(in) :: h, L_o
        real(kind=rp) :: corr
     end function corr_h_interface

     function f_interface(Ri_b, h, z0, z0t, L_ob, slaw_m, slaw_h) result(f)
        real(kind=rp), intent(in) :: Ri_b, h, z0, z0t, L_ob
        procedure(slaw_m_interface) :: slaw_m
        procedure(slaw_h_interface) :: slaw_h
     end function f_interface

     function dfdl_interface(l_upper, l_lower, h, z0, z0t, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
        real(kind=rp), intent(in) :: l_upper, l_lower, h, z0, z0t, L_ob, fd_h
        procedure(slaw_m_interface) :: slaw_m
        procedure(slaw_h_interface) :: slaw_h
        real(kind=rp) :: dfdl
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

  subroutine split_bc_type(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b, flag)
    real(kind=rp), intent(in) :: g, hi, ti, ts, magu, kappa
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: Ri_b, q
    integer :: flag

    select case (bc_type)
    case ("neumann")
      if (flag /= 0), then
        Ri_b = - g*hi / ti*q / (magu**3*kappa**2)
        f_ptr => f_neumann
        dfdl_ptr => dfdl_neumann
      end if
    case ("dirichlet")
      q =  kappa*utau*(ts - ti)/log(hi/z0t)
      if (flag /= 0), then
        Ri_b = g*hi/ti*(ti - ts)/magu**2
        f_ptr => f_dirichlet
        dfdl_ptr => dfdl_dirichlet
      end if
    case default
      call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
      ! perhaps I should split this in two "cases" so that i can have two error messages 
      ! (one in case the user forget ts or q and one in case it forgets bc_type) ?
    end select
  end subroutine split_bc_type

  subroutine set_stability_regime(Ri_b)
    real(kind=rp), intent(in) :: Ri_b

    if (Ri_b > 0.0) then
       slaw_m_ptr => slaw_m_stable
       slaw_h_ptr => slaw_h_stable
       corr_m_ptr => corr_m_stable
       corr_h_ptr => corr_h_stable
    elseif (Ri_b < 0.0) then
       slaw_m_ptr => slaw_m_convective
       slaw_h_ptr => slaw_h_convective
       corr_m_ptr => corr_m_convective
       corr_h_ptr => corr_h_convective

!!> INCLUDE NEUTRAL WITH A BETTER IMPLEMENTATION (COPY ROUGHLOG)
    ! else
    !    slaw_m_ptr => slaw_m_neutral   ! include or not?
    !    slaw_h_ptr => slaw_h_neutral
    !    corr_m_ptr => corr_m_neutral
    !    corr_h_ptr => corr_h_neutral
    end if
  end subroutine set_stability_regime

  !================================================
  ! 5. Main wall-flux routine
  !================================================

  subroutine most_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &  
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, z0, bc_type, ts, q, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0, z0t
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: ts,q   ! should this be multidimensional?
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: i, count, count_def
    integer, parameter :: max_count = 20
    real(kind=rp) :: ui, vi, wi, ti, hi
    real(kind=rp) :: magu, utau, normu
    real(kind=rp) :: L_ob, L_ob_def, L_upper, L_lower, L_backup, L_old
    real(kind=rp) :: Ri_b, Ri_b_def, f, dfdl, fd_h    
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: tol = 0.001_rp
    real(kind=rp), parameter :: NR_step = 0.001_rp

    do i=1, n_nodes
   
      ! Sample the variables
      ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      hi = h(i)

      ! Project on horizontal directions
      normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
      ui = ui - normu * n_x(i) 
      vi = vi - normu * n_y(i)
      wi = wi - normu * n_z(i)

      ! Compute velocity magnitude
      magu = sqrt(ui**2 + wi**2)  ! does every number need to be _rp?

      ! utau initialisation
      if (tstep < 1) then
        utau = sqrt( sqrt( tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2 ) ) 
      else 
        utau = magu*kappa / log(hi/z0)
      end if

      ! Get q, Ri_b, f_ptr, dfdl_ptr based on bc_type 
      ! Maybe redundant, but needed to initialise Rib
      call split_bc_type(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b, 0)

      ! Compute Obukhov length
      if (tstep > 0) then
        ! Get q, Ri, f, dfdl based on bc_type 
        call split_bc_type(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b, 1)

        ! Obukhov _ based on the previous-step utau
        L_ob= -(ti*utau**3)/(kappa*g*q)  ! this is only initialisation (ts,ti,300 doesn't matter)
        L_ob_def = L_ob
        L_backup = L_ob

        if (Ri_b == 0) then
          ! Neutral
          L_ob = 0   ! improve this once the neutral option is implemented

        else:
          L_old = 0
          count = 0
          do while  ((abs(L_old - L_ob)/abs(L_ob) .gt. tol) .and. (count .lt. max_count))
            ! Switch between stable and convective based on bulk Richardson (Ri_b)
            L_old = L_ob
            count = count + 1
            fd_h = NR_step*L_ob
            L_upper = L_ob + fd_h
            L_lower = L_ob - fd_h

            ! Set slaw and corr pointers based on stability
            call set_stability_regime(Ri_b)

            ! Compute L_ob based on stability and bc_type
            f = f_ptr(Ri_b, hi, z0, z0t, L_ob, slaw_m_ptr, slaw_h_ptr)
            dfdl = dfdl_ptr(l_upper, l_lower, h, z0, z0t, L_ob, slaw_m_ptr, slaw_h_ptr, fd_h)
            L_ob = L_ob - f/dfdl

            if (abs(L_ob) > 20000 .or. abs(L_ob) < 1e-5_rp) then
                count = max_count
            end if
          end if
        end do

        ! Error handling 
        if (count .eq. max_count) then
            L_ob_def = L_backup
            call neko_error("Obukhov length did not converge (MOST wall model)") 
        end if

        L_ob_def = L_ob
        Ri_b_def = Ri_b
        count_def = count

        ! Based on stability and bc_type, compute utau/q 
        call set_stability_regime(Ri_b) 
        select case (bc_type)
          case ("neumann")
            ! Compute u* with the new Obukhov length
            utau = kappa*magu/slaw_m_ptr(L_ob, hi, z0)
          case ("dirichlet")
            ! Compute u* with the new Obukhov length
            utau = kappa*magu/slaw_m_ptr(L_ob, hi, z0)
            ! and compute q from here
            q = kappa*utau*(ts - ti)/slaw_h_ptr(L_ob, hi, z0t)  ! z0t placeholder
          case default
            call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
        end select 

      end if

      ! Distribute according to the velocity vector
      tau_x(i) = -utau**2 * ui / magu
      tau_y(i) = -utau**2 * vi / magu
      tau_z(i) = -utau**2 * wi / magu
      ! export q as well?
    end do

  end subroutine most_compute_cpu

  !================================================
  ! 6. Concrete similarity functions
  !================================================

  !--------------- Neutral ----------------
  ! function slaw_m_neutral(zeta) result(slaw)
  !   real(kind=rp), intent(in) :: zeta
  !   real(kind=rp) :: slaw
  !   slaw = 1.0                        include neutral or not?
  ! end function slaw_m_neutral

  ! function slaw_h_neutral(zeta) result(slaw)
  !   real(kind=rp), intent(in) :: zeta
  !   real(kind=rp) :: slaw
  !   slaw = 1.0
  ! end function slaw_h_neutral

  ! function corr_m_neutral(zeta) result(corr)
  !   real(kind=rp), intent(in) :: zeta
  !   real(kind=rp) :: corr
  !   corr = 0.0
  ! end function corr_m_neutral

  ! function corr_h_neutral(zeta) result(corr)
  !   real(kind=rp), intent(in) :: zeta
  !   real(kind=rp) :: corr
  !   corr = 0.0
  ! end function corr_h_neutral


  !--------------- Stable ----------------

  function slaw_m_stable(h,L_o,z_0) result(slaw)
    real(kind=rp), intent(in) :: h,L_o,z_0
    real(kind=rp) :: slaw

    slaw = log(h/z0)-corr_m_stable(h,L_o)+corr_m_stable(z0,L_o)
  end function slaw_m_stable

  function slaw_h_stable(h,L_o,z0t) result(slaw)
    real(kind=rp), intent(in) :: h,L_o,z0t
    real(kind=rp) :: slaw

    slaw = log(h/z0t)-corr_h_stable(h,L_o)+corr_h_stable(z0t,L_o)
  end function slaw_h_stable

  function corr_m_stable(h,L_o) result(corr)
    real(kind=rp), intent(in) :: h,L_o
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d  ! parameter
    real(kind=rp) :: zeta

    zeta = h/L_o
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = - a*zeta - b*(zeta-c/d)*exp(-d*zeta) - b*c/d
  end function corr_m_stable

  function corr_h_stable(h,L_o) result(corr)    ! identical to corr_m_stable...?
    real(kind=rp), intent(in) :: h,L_o
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d
    real(kind=rp) :: zeta

    zeta = h/L_o
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = -b * (zeta-c/d)*exp(-d*zeta)-(1.0_rp+ 0.6666666_rp * a * zeta)**1.5_rp-b*c/d + 1.0_rp
  end function corr_h_stable

  !--------------- Convective ----------------

  function slaw_m_convective(h,L_o,z0) result(slaw)
    real(kind=rp), intent(in) :: h, L_o, z0
    real(kind=rp) :: slaw

    slaw = log(h/z0) - corr_m_convective(h, L_o) + corr_m_convective(z0, L_o)
  end function slaw_m_convective

  function slaw_h_convective(h,L_o,z0t) result(slaw)
    real(kind=rp), intent(in) :: h, L_o, z0t
    real(kind=rp) :: slaw

    slaw = log(h/z0t) - corr_h_conv(h, L_o) + corr_h_conv(z0t, L_o)
  end function slaw_h_convective

  function corr_m_convective(h,L_o) result(corr)
    real(kind=rp), intent(in) :: h, L_o
    real(kind=rp) :: xi, pi, zeta
    real(kind=rp) :: corr

    zeta = h/L_o
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi)) + log(0.5_rp*(1 + xi**2)) - 2*atan(xi) + pi/2
  end function corr_m_convective

  function corr_h_convective(h,L_o) result(corr)
    real(kind=rp), intent(in) :: h, L_o
    real(kind=rp) :: zeta, pi, xi
    real(kind=rp) :: corr

    zeta = h/L_o
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi**2))
  end function corr_h_convective

!!!_-----------------------------------------_!!!

  function f_neumann(Ri_b, h, z0, z0t, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, h, z0, z0t, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - h/L_ob/slaw_m(h, L_ob, z0)**3) !! ChatGPT reports that this is wrong in terms of sign?
  end function f_neumann

  function dfdl_neumann(l_upper, l_lower, h, z0, z0t, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, h, z0, z0t, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-h/l_upper/slaw_m(h, l_upper, z0)**3)  ! conv
    dfdl = dfdl + (h/l_lower/slaw_m(h, l_lower, z0)**3)  ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_neumann

  function f_dirichlet(Ri_b, h, z0, z0t, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, h, z0, z0t, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - h/L_ob*slaw_h(h, L_ob, z0t)/slaw_m(h, L_ob, z0)**2)  ! conv
  end function f_dirichlet

  function dfdl_dirichlet(l_upper, l_lower, h, z0, z0t, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, h, z0, z0t, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-h/l_upper*slaw_h(h, l_upper, z0t)/slaw_m(h, l_upper, z0)**2)  ! conv
    dfdl=dfdl + (h/l_lower*slaw_h(h, l_lower, z0t)/slaw_m(h, l_lower, z0)**2)  ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_dirichlet
  

end module most_cpu
