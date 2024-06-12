! Copyright (c) 2024, The Neko Authors
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
!
!> Implements `spalding_t`.
module spalding
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get_or_default
  use logger, only : neko_log, NEKO_LOG_DEBUG
  use utils, only : neko_error

  implicit none
  private

  !> Wall model based on Spalding's law of the wall.
  !! Reference: http://dx.doi.org/10.1115/1.3641728
  type, public, extends(wall_model_t) :: spalding_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => spalding_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
       spalding_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => spalding_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => spalding_compute
     !> Solve for the friction velocity
     procedure, private, pass(this) :: solve
  end type spalding_t

contains
  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  subroutine spalding_init(this, coef, msk, facet, nu, h_index, json)
    class(spalding_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B

    call json_get_or_default(json, "kappa", kappa, 0.41_rp)
    call json_get_or_default(json, "B", B, 5.2_rp)

    call this%init_from_components(coef, msk, facet, nu, h_index, kappa, B)
  end subroutine spalding_init

  !> Constructor from components.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  subroutine spalding_init_from_components(this, coef, msk, facet, nu, h_index,&
                                           kappa, B)
    class(spalding_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: nu
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Spalding's law is only available on the CPU backend.")
    end if

    call this%init_base(coef, msk, facet, nu, h_index)

    this%kappa = kappa
    this%B = B
  end subroutine spalding_init_from_components


  !> Destructor for the spalding_t (base) class.
  subroutine spalding_free(this)
    class(spalding_t), intent(inout) :: this

    call this%free_base()

  end subroutine spalding_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine spalding_compute(this, t, tstep)
    class(spalding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    do i=1, this%n_nodes
      ! Sample the velocity
      ui = u%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      vi = v%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
      wi = w%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))

      ! Project on tangential direction
      normu = ui * this%n_x%x(i) + vi * this%n_y%x(i) + wi * this%n_z%x(i)

      ui = ui - normu * this%n_x%x(i)
      vi = vi - normu * this%n_y%x(i)
      wi = wi - normu * this%n_z%x(i)

      magu = sqrt(ui**2 + vi**2 + wi**2)

      ! Get initial guess for Newton solver
      if (tstep .eq. 1) then
         guess = sqrt(magu * this%nu / this%h%x(i))
      else
         guess = this%tau_x(i)**2 + this%tau_y(i)**2 + this%tau_z(i)**2
         guess = sqrt(sqrt(guess))
      end if

      utau =  this%solve(magu, this%h%x(i), guess)

      ! Distribute according to the velocity vector
      this%tau_x(i) = -utau**2 * ui / magu
      this%tau_y(i) = -utau**2 * vi / magu
      this%tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine spalding_compute

  !> Newton solver for the algebraic equation defined by the law.
  !! @param u The velocity value.
  !! @param y The wall-normal distance.
  !! @param guess Initial guess.
  function solve(this, u,  y, guess) result(utau)
    class(spalding_t), intent(inout) :: this
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: guess
    real(kind=rp) :: yp, up, kappa, B, utau
    real(kind=rp) :: error, f, df, old
    integer :: niter, k, maxiter

    utau = guess
    kappa = this%kappa
    B = this%B

    maxiter = 100

    do k=1, maxiter
      up = u / utau
      yp = y * utau / this%nu
      niter = k
      old = utau

      ! Evaluate function and its derivative
      f = (up + exp(-kappa*B)* &
          (exp(kappa*up) - 1.0_rp - kappa*up - 0.5_rp*(kappa*up)**2 - &
           1.0_rp/6*(kappa*up)**3) - yp)

      df = (-y / this%nu - u/utau**2 - kappa*up/utau*exp(-kappa*B) * &
           (exp(kappa*up) - 1 - kappa*up - 0.5*(kappa*up)**2))

      ! Update solution
      utau = utau - f / df

      error = abs((old - utau)/old)

      if (error < 1e-3) then
        exit
      endif

    enddo

    if ((niter .eq. maxiter) .and. (neko_log%level_ .eq. NEKO_LOG_DEBUG)) then
       write(*,*) "Newton not converged", error, f, utau, old, guess
    end if
end function solve


end module spalding