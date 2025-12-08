module most_cpu
  use num_types, only : rp
  use utils, only : neko_error
  implicit none
  private

  public :: most_compute_cpu

  ! REWIRTE STUFF ACCORDING OT NEKO CONVENTION


  !===============================
  ! 2. Abstract interfaces
  !===============================
  !
  ! These specify the *shape* of the functions
  ! without tying them to a specific expression.
  !

  abstract interface
     function slaw_m_interface(z,L_o,z0) result(slaw)
       real(kind=rp), intent(in) :: z, L_o, z1
       real(kind=rp) :: slaw
     end function slaw_m_interface

     function slaw_h_interface(z,L_o,z1) result(slaw)
       real(kind=rp), intent(in) :: z, L_o, z1
       real(kind=rp) :: slaw
     end function slaw_h_interface

     function corr_m_interface(z,L_o) result(corr)
       real(kind=rp), intent(in) :: z, L_o
       real(kind=rp) :: corr
     end function corr_m_interface

     function corr_h_interface(z,L_o) result(corr)
       real(kind=rp), intent(in) :: z, L_o
       real(kind=rp) :: corr
     end function corr_h_interface
  end interface

  !===============================
  ! 3. Procedure pointers
  !===============================
  !
  ! These will point to the correct similarity
  ! functions depending on stability regime.
  !

  procedure(slaw_m_interface), pointer :: slaw_m_ptr => null()
  procedure(slaw_h_interface), pointer :: slaw_h_ptr => null()
  procedure(corr_m_interface), pointer :: corr_m_ptr => null()
  procedure(corr_h_interface), pointer :: corr_h_ptr => null()

contains

  !================================================
  ! 4. User-facing selectors for stability regime
  !================================================

  subroutine get_Ri(bc_type, g, hi, ti, ts, magu, kappa, q, rib)
    real(kind=rp), intent(in) :: g, hi, ti, ts, , magu, kappa
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: rib, q

    select case (bc_type)
    case ("neumann")
      ! Compute Bulk Richardson number with q form bc
      rib = - g*hi / ti*q / (magu**3*kappa**2)
    case ("dirichlet")
      ! Obtain q from similarity
      q =  kappa*utau*(ts - ti)/log(h/z1)
      ! and compute Ri from there
      rib = g*hi/ti*(ti - ts)/magu**2
    case default
      call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
      ! perhaps I should split this in two "cases" so that i can have two error messages 
      ! (one in case the user forget ts or q and one in case it forgets bc_type) ?
    end select
  end subroutine get_Ri

  subroutine set_stability_regime(rib)
    real(kind=rp), intent(in) :: rib

    if (rib > 0.0) then
       slaw_m_ptr => slaw_m_stable
       slaw_h_ptr => slaw_h_stable
       corr_m_ptr => corr_m_stable
       corr_h_ptr => corr_h_stable
    elseif (rib < 0.0) then
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
  ! 5. Main wall-flux routine calling the pointers
  !================================================

  subroutine most_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &  
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, z0, bc_type, ts, q, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: ts,q
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: i, count
    real(kind=rp) :: ui, vi, wi, ti, magu, utau, normu
    real(kind=rp) :: L_ob, L_upper, L_lower, L_backup, L_old
    real(kind=rp) :: rib, f, dfdl, fd_h    
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: tol = 0.001_rp
    real(kind=rp), parameter :: NR_step = 0.001_rp
    integer, parameter :: max_count = 20

    do i=1, n_nodes
   
      ! Sample the variables
      ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
      ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

      ! Project on horizontal directions
      normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
      ui = ui - normu * n_x(i) 
      vi = vi - normu * n_y(i)
      wi = wi - normu * n_z(i)

      ! Compute velocity magnitude
      magu = sqrt(ui**2 + vi**2 + wi**2)  ! does every number need to be _rp?

      ! Get uncorrected utau for a first guess other wise from last  (e.g. if timestep < 3)
      if (tstep < 1) then
        utau = sqrt( sqrt( tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2 ) ) 
      else 
        utau = magu*kappa / log(h(i)/z0)
      end if

      ! Get heat flux from the bc    
      call get_Ri(bc_type, g, h(i), ti, ts, magu, kappa, q, rib)

      ! Compute Obukhov length
      if (tstep > 0) then
        ! Get Ri/q/ts from this timestep
        call get_Ri(bc_type, g, h(i), ti, ts, magu, kappa, q, rib)
        !!!! here is the first case-select (can we join them?)

        ! Obukhov l based on the previous-step utau
        L_ob = -(ts*utau**3)/(kappa*g*q) 
        L_backup = L_ob

        L_old = 0
        count = 0
        do while  ((abs(L_old - L_ob)/abs(L_ob) .gt. tol) .and. (count .lt. max_count))
          if rib == 0:
            ! Neutral
            L_ob = 0   ! imrpove this once the neutral option is implemented
          else:
            ! Switch between neutral and convective based on bulk Richardson (rib)
            L_old = L_ob
            count = count + 1

            fd_h = NR_step*L_ob
            L_upper = L_ob + fd_h
            L_lower = L_ob - fd_h


            !> Compute f, dfdl, l_ob based on stability
            call set_stability_regime(rib)
            
            select case (bc_type)
              case ("neumann")
                L_ob = L_ob - f_neumann/dfdl_neumann
              case ("dirichlet")
                L_ob = L_ob - f_dirichlet/dfdl_dirichlet
              case default
                call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
            end select

            if (abs(L_ob) .gt. 20000 .or. abs(L_ob) .lt. 1e-5) then
                count = 20
            end if

        end do

        ! Error handling 
        if (count .eq. 20) then
            call neko_error("Obukhov length did not converge (convective MOST wall model)") 
        end if

        ! Based on stability and bc_type, compute utau/q   DOES THIS MAKE SENSE?
        call set_stability_regime(rib) !necessary to call it again or what?
        select case (bc_type)
          case ("neumann")
            ! Compute u* with the new Obukhov length
            utau = kappa*magvh/slaw_m_ptr(l_obukhov, h, z0)
          case ("dirichlet")
            ! Compute u* with the new Obukhov length
            utau = kappa*magvh/slaw_m_ptr(l_obukhov, h, z0)
            ! and compute q from here
            q = kappa*utau*(ts - th)/slaw_h_ptr(l_obukhov, h, z1)
          case default
            call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
        end select 

      end if















        ! Compute utau (should I update Obkuhov length?)
        utau = kappa*magu/slaw(L_ob, h(i), z0)    

      end if

      ! Distribute according to the velocity vector
      tau_x(i) = -utau**2 * ui / magu
      tau_y(i) = -utau**2 * vi / magu
      tau_z(i) = -utau**2 * wi / magu
      ! export q as well?
    end do

  end subroutine most_convective_compute_cpu

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

  function slaw_m_stable(z,L_o,z_0) result(slaw)
    real(kind=rp), intent(in) :: z,L_o,z_0
    real(kind=rp) :: slaw
    real(kind=rp) :: corr_m_stable

    slaw = log(z/z0)-corr_m_stable(z,L_o)+corr_m_stable(z0,L_o)
  end function slaw_m_stable

  function slaw_h_stable(z,L_o,z1) result(slaw)
    real(kind=rp), intent(in) :: z,L_o,z1
    real(kind=rp) :: slaw
    real(kind=rp) :: corr_h_stable

    slaw = log(z/z1)-corr_q_stable(z,L_o)+corr_q_stable(z1,L_o)
  end function slaw_h_stable

  function corr_m_stable(z,L_o) result(corr)
    real(kind=rp), intent(in) :: z,L_o
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d  ! parameter

    a = 0.7
    b = 0.75
    c = 5
    d = 0.35
    corr = a*z/L_o + b*(z/L_o-c/d)*exp(-d*z/L_o) + b*c/d
  end function corr_m_stable

  function corr_h_stable(z,L_o) result(corr)    ! identical to corr_m_stable...?
    real(kind=rp), intent(in) :: z,L_o
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d

    a = 0.7
    b = 0.75
    c = 5
    d = 0.35
    corr = a*z/L_o + b*(z/L_o-c/d)*exp(-d*z/L_o) + b*c/d
  end function corr_h_stable

  !--------------- Convective ----------------

  function slaw_m_convective(z,L_o,z0) result(slaw)
    real(kind=rp), intent(in) :: z, L_o, z0
    real(kind=rp) :: slaw

    slaw = log(z/z0) - corr_u_conv(z, L_o) + corr_u_conv(z0, L_o)
  end function slaw_m_convective

  function slaw_h_convective(z,L_o,z1) result(slaw)
    real(kind=rp), intent(in) :: z, L_o, z1
    real(kind=rp) :: slaw

    slaw = log(z/z1) - corr_h_conv(z, L_o) + corr_h_conv(z1, L_o)
  end function slaw_h_convective

  function corr_m_convective(z,L_o) result(corr)
    real(kind=rp), intent(in) :: z, L_o
    real(kind=rp) :: xi, pi, zeta
    real(kind=rp) :: corr

    ztea = z/L_o
    pi = 4*atan(1.0)
    xi = (1.0 - 16.0*zeta)**0.25
    corr_m_conv = 2*log(0.5*(1 + xi)) + log(0.5*(1 + xi**2)) - 2*atan(xi) + pi/2
  end function corr_m_convective

  function corr_h_convective(z,L_o) result(corr)
    real(kind=rp), intent(in) :: z, L_o
    real(kind=rp) :: zeta, pi, xi
    real(kind=rp) :: corr

    zeta = z/L_o
    pi = 4*atan(1.0)
    xi = (1.0 - 16.0*zeta)**0.25
    corr_h_conv = 2*log(0.5*(1 + xi**2)))
  end function corr_h_convective

!!!_-----------------------------------------_!!!

!!!! TO BE COMPLETED !!!!!

  function f_neumann()
   f = (rib - h/l_obukhov/
     $        similarity_law_u_conv(l_obukhov, h, z0)**3)

  end function f_neumann()

  function dfdl_neumann()
              dfdl = (-h/l_upper
     $               /similarity_law_u_conv(l_upper, h, z0)**3)
              dfdl = dfdl + (h/l_lower/
     $        similarity_law_u_conv(l_lower, h, z0)**3)
              dfdl = dfdl/(2*fd_h)
  end function dfdl_neumann()

  function f_dirichlet()
f=(rib - h/l_obukhov
     $          *similarity_law_q_conv(l_obukhov, h, z1)
     $          /similarity_law_u_conv(l_obukhov, h, z0)**2)

  end function f_dirichlet

  function dfdl_dirichlet()
              dfdl = (-h/l_upper*similarity_law_q_conv(l_upper, h, z1)
     $             /similarity_law_u_conv(l_upper, h, z0)**2)
              dfdl=dfdl + (h/l_lower
     $             *similarity_law_q_conv(l_lower, h, z1)
     $             /similarity_law_u_conv(l_lower, h, z0)**2)
              dfdl = dfdl/(2*fd_h)
  end function dfdl_dirchlet
  

end module most_wall_model
