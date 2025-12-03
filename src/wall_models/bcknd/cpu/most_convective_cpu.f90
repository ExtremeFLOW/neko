! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the CPU kernel for the `most_convective_t` type.



!>===============================<!!
!!! WIP !!!
! Things to do (prio order):
! 1) read ts
! 2) read q
! 3) improve error handling
! 4) change arguments to the subroutine 
!    in the rest of neko
!>===============================<!!


module most_convective_cpu
  use num_types, only : rp
  use utils, only : neko_error
  implicit none
  private

  public :: most_convective_compute_cpu

contains
  !> Compute the wall shear stress on CPU using rough convective MOST model.
  !! @param tstep The current time-step.
  subroutine most_convective_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &  
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, z0, tstep)!  I added s   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: i, count
    real(kind=rp) :: ui, vi, wi, magu, utau, normu
    real(kind=rp) :: l_obukhov, l_upper, l_lower, l_backup, l_old
    real(kind=rp) :: f, dfdl, fd_h    
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: tol = 0.001_rp
    integer, parameter :: max_count = 20
    real(kind=rp) :: rib

    real(kind=rp) :: ti !, intent(in) q   ! will need to deal with surface temperature sooner or later
    real(kind=rp), parameter :: q = 0.05_rp  ! placeholder for heat flux, to be passed as argument later  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rp), parameter :: ts = 300_rp! placeholder for heat flux, to be passed as argument later   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1, n_nodes

       ! Sample the variables
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Get heat flux from the bc     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Project on horizontal directions
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i) 
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       ! Compute velocity magnitude
       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Compute Bulk Richardson number
       rib = - g*h(i) / ti*q / (magu**3*kappa**2)

       ! Get uncorrected utau for a first guess other wise from last  (e.g. if timestep < 3)
        if (tstep < 3) then
          utau = sqrt( sqrt( tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2 ) ) 
        else 
          utau = magu*kappa / log(h(i)/z0)
        end if
   
       ! Compute Obukhov length  (when timestep > e.g. 3)
        if (tstep > 2) then
          l_obukhov = -(ts*utau**3)/(kappa*g*q)
          l_backup = l_obukhov
 
          l_old = 0
          count = 0
          do while  ((abs(l_old - l_obukhov)/abs(l_obukhov) .gt. tol) .and. (count .lt. max_count))
            l_old = l_obukhov
            count = count + 1

            fd_h = 1e-3*l_obukhov
            l_upper = l_obukhov + fd_h
            l_lower = l_obukhov - fd_h

            f = (rib - h(i)/l_obukhov/similarity_law(l_obukhov, h(i), z0)**3)
            dfdl = (-h(i)/l_upper/similarity_law(l_upper, h(i), z0)**3)
            dfdl = dfdl + (h(i)/l_lower/similarity_law(l_lower, h(i), z0)**3)
            dfdl = dfdl/(2*fd_h)

            l_obukhov = l_obukhov - f/dfdl
            
            if (abs(l_obukhov) .gt. 20000 .or. abs(l_obukhov) .lt. 1e-5) then
                count = 20
            end if

          end do

       ! Error handling (if convergence not achieved)
       ! Temporary: I think it's better to just make the sim crash (?)
        if (count .eq. 20) then
           call neko_error("Obukhov length did not converge (convective MOST wall model)")   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  l_obukhov = l_backup
        end if

       ! Compute utau (should I update Obkuhov length?)
        utau = kappa*magu/similarity_law(l_obukhov, h(i), z0)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if

       ! Distribute according to the velocity vector
       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine most_convective_compute_cpu

!> Compute correction for the u log law in the convective case
  real(kind=rp) function correction(z,l) result(corr)
   implicit none
   real(kind=rp), intent(in) :: z, l
   real(kind=rp) :: xi, pi

   pi = 4.0_rp*atan(1.0_rp)
   xi = (1.0_rp - 16.0_rp*z/l)**0.25_rp
   corr = 2.0_rp*log(0.5_rp*(1.0_rp + xi)) + log(0.5_rp*(1.0_rp + xi**2.0_rp)) - 2.0_rp*atan(xi) + pi/2.0_rp
  end function correction

!> Compute the similarity law for velocity
  real(kind=rp) function similarity_law(l_obukhov, z, z0) result(slaw)
    implicit none
    real(kind=rp), intent(in) :: l_obukhov, z, z0

    slaw = log(z/z0) - correction(z, l_obukhov) + correction(z0, l_obukhov)
  end function

end module most_convective_cpu