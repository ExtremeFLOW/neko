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
!> Subroutines to add advection terms to the RHS of a transport equation.
module adv_oifs
  use advection, only: advection_t
  use num_types, only: rp
  use space, only: space_t, GL
  use field, only: field_t
  use coefs, only: coef_t
  use math, only: copy, rzero
  use operators, only: runge_kutta, set_convect_rst
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM
  use interpolation, only: interpolator_t
  use time_interpolator, only: time_interpolator_t
  use field_series, only: field_series_t
  use field_list, only: field_list_t
  use time_scheme_controller, only: time_scheme_controller_t
  use device, only: device_map, device_free
  use device_math, only: device_addcol3s2, device_rzero
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  implicit none
  private

  !! Type encapsulating operator-integration-factor splitting advection
  !! routines with dealiasing applied
  !! Literature:
  !! https://www.mcs.anl.gov/~fischer/nek5000/oifs.pdf
  !! https://publications.anl.gov/anlpubs/2017/12/140626.pdf
  !! https://dl.acm.org/doi/abs/10.1007/BF01063118
  type, public, extends(advection_t) :: adv_oifs_t
     !> Number of RK4 sub-steps
     integer :: ntaubd
     !> Coeffs of the higher-order space
     type(coef_t) :: coef_GL
     !> Coeffs of the original space in the simulation
     type(coef_t), pointer :: coef_GLL => null()
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     !> The additional higher-order space used in dealiasing
     type(space_t) :: Xh_GL
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL => null()
     !> The time interpolator scheme
     type(time_interpolator_t) :: dtime
     !> The lagged velocity and scalar fields
     type(field_series_t), pointer :: ulag, vlag, wlag, slag => null()
     !> The times corresponding to the lagged fields
     real(kind=rp), pointer :: ctlag(:) => null()
     !> The time-steps corresponding to the lagged fields
     real(kind=rp), pointer :: dctlag(:) => null()
     !> The time scheme controller for the oifs scheme
     type(time_scheme_controller_t), pointer :: oifs_scheme => null()
     !> The current convecting field in GL space and rst format
     type(field_t) :: cr_GL, cs_GL, ct_GL
     !> The convecting field series in GL space and rst format
     type(field_series_t) :: convr_GL, convs_GL, convt_GL
     !> The time interpolated convecting field used in Runge_Kutta method
     type(field_t), pointer :: cr_k1, cs_k1, ct_k1
     type(field_t), pointer :: cr_k23, cs_k23, ct_k23
     type(field_t), pointer :: cr_k4, cs_k4, ct_k4
     !> The field_list containing the time interpolated convecting field
     type(field_list_t) :: conv_k1, conv_k23, conv_k4
     !> The convecting velocity field in GL space
     real(kind=rp), allocatable :: cx(:), cy(:), cz(:)
     !> Device pointers for cx, cy, cz
     type(c_ptr) :: cx_d = C_NULL_PTR, cy_d = C_NULL_PTR, cz_d = C_NULL_PTR

   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS
     procedure, pass(this) :: compute => adv_oifs_compute
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS
     procedure, pass(this) :: compute_scalar => adv_oifs_compute_scalar
     !> Constructor
     procedure, pass(this) :: init => adv_oifs_init
     !> setting characteristic convecting field
     procedure, pass(this) :: set_conv_velocity_fst
     !> Destructor
     procedure, pass(this) :: free => adv_oifs_free
  end type adv_oifs_t

contains

  !> Constructor
  !! @param lxd The polynomial order of the space used in the dealiasing.
  !! @param coef The coefficients of the (space, mesh) pair.
  !! @param ctarget Target CFL number.
  !! @param ulag The x component of lagged velocity.
  !! @param vlag The y component of lagged velocity.
  !! @param wlag The z component of lagged velocity.
  !! @param dtlag Lagged time steps.
  !! @param tlag Lagged simulation times.
  !! @param time_scheme The bdf-ext time scheme used in the method.
  !! @param slag The lagged scalar field.
  subroutine adv_oifs_init(this, lxd, coef, ctarget, ulag, vlag, wlag, &
       dtlag, tlag, time_scheme, slag)
    implicit none
    class(adv_oifs_t) :: this
    integer, intent(in) :: lxd
    type(coef_t), target :: coef
    real(kind=rp), intent(in) :: ctarget
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    real(kind=rp), target, intent(in) :: dtlag(10)
    real(kind=rp), target, intent(in) :: tlag(10)
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_series_t), target, optional :: slag
    integer :: nel, n_GL, n, idx, idy, idz
    real(kind=rp) :: max_cfl_rk4

    ! stability limit for RK4 including safety factor
    max_cfl_rk4 = 2.0
    this%ntaubd = max(int(ctarget/max_cfl_rk4),1)

    call this%Xh_GL%init(GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)

    call this%coef_GL%init(this%Xh_GL, coef%msh)

    call this%cr_GL%init(coef%msh, this%Xh_GL)
    call this%cs_GL%init(coef%msh, this%Xh_GL)
    call this%ct_GL%init(coef%msh, this%Xh_GL)

    nel = coef%msh%nelv
    n_GL = nel*this%Xh_GL%lxyz
    n = nel*coef%Xh%lxyz

    call this%GLL_to_GL%map(this%coef_GL%drdx, coef%drdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdx, coef%dsdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdx, coef%dtdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdy, coef%drdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdy, coef%dsdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdy, coef%dtdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdz, coef%drdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdz, coef%dsdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdz, coef%dtdz, nel, this%Xh_GL)


    allocate(this%cx(n_GL))
    allocate(this%cy(n_GL))
    allocate(this%cz(n_GL))

    allocate(this%cr_k1)
    allocate(this%cs_k1)
    allocate(this%ct_k1)
    allocate(this%cr_k23)
    allocate(this%cs_k23)
    allocate(this%ct_k23)
    allocate(this%cr_k4)
    allocate(this%cs_k4)
    allocate(this%ct_k4)


    call this%cr_k1%init(coef%msh, this%Xh_GL)
    call this%cs_k1%init(coef%msh, this%Xh_GL)
    call this%ct_k1%init(coef%msh, this%Xh_GL)

    call this%cr_k23%init(coef%msh, this%Xh_GL)
    call this%cs_k23%init(coef%msh, this%Xh_GL)
    call this%ct_k23%init(coef%msh, this%Xh_GL)

    call this%cr_k4%init(coef%msh, this%Xh_GL)
    call this%cs_k4%init(coef%msh, this%Xh_GL)
    call this%ct_k4%init(coef%msh, this%Xh_GL)

    call this%conv_k1%init(3)
    call this%conv_k23%init(3)
    call this%conv_k4%init(3)

    call this%conv_k1%assign(1, this%cr_k1)
    call this%conv_k1%assign(2, this%cs_k1)
    call this%conv_k1%assign(3, this%ct_k1)

    call this%conv_k23%assign(1, this%cr_k23)
    call this%conv_k23%assign(2, this%cs_k23)
    call this%conv_k23%assign(3, this%ct_k23)

    call this%conv_k4%assign(1, this%cr_k4)
    call this%conv_k4%assign(2, this%cs_k4)
    call this%conv_k4%assign(3, this%ct_k4)

    call this%dtime%init(1)
    this%ulag => ulag
    this%vlag => vlag
    this%wlag => wlag
    this%ctlag => tlag
    this%dctlag => dtlag
    this%oifs_scheme => time_scheme

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%cx, this%cx_d, n_GL)
       call device_map(this%cy, this%cy_d, n_GL)
       call device_map(this%cz, this%cz_d, n_GL)
    end if

    ! Initializing the convecting fields
    ! Map the velocity fields from GLL space to GL space
    call this%GLL_to_GL%map(this%cx, this%ulag%f%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cy, this%vlag%f%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cz, this%wlag%f%x, nel, this%Xh_GL)

    ! Set the convecting field in the rst format
    call set_convect_rst(this%cr_GL, this%cs_GL, this%ct_GL, &
         this%cx, this%cy, this%cz, this%Xh_GL, this%coef_GL)

    ! Set the convecting field series
    call this%convr_GL%init(this%cr_GL, 3)
    call this%convs_GL%init(this%cs_GL, 3)
    call this%convt_GL%init(this%ct_GL, 3)

    ! Repeat for previous time-steps
    call this%GLL_to_GL%map(this%cx, this%ulag%lf(1)%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cy, this%vlag%lf(1)%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cz, this%wlag%lf(1)%x, nel, this%Xh_GL)

    call set_convect_rst(this%cr_GL, this%cs_GL, this%ct_GL, &
         this%cx, this%cy, this%cz, this%Xh_GL, this%coef_GL)

    this%convr_GL%lf(1) = this%cr_GL
    this%convs_GL%lf(1) = this%cs_GL
    this%convt_GL%lf(1) = this%ct_GL

    call this%GLL_to_GL%map(this%cx, this%ulag%lf(2)%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cy, this%vlag%lf(2)%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cz, this%wlag%lf(2)%x, nel, this%Xh_GL)

    call set_convect_rst(this%cr_GL, this%cs_GL, this%ct_GL, &
         this%cx, this%cy, this%cz, this%Xh_GL, this%coef_GL)

    this%convr_GL%lf(2) = this%cr_GL
    this%convs_GL%lf(2) = this%cs_GL
    this%convt_GL%lf(2) = this%ct_GL

    ! Initilize the lagged scalar field, if present.
    if (present(slag)) then
       this%slag => slag
    end if

  end subroutine adv_oifs_init

  !> Destructor
  subroutine adv_oifs_free(this)
    class(adv_oifs_t), intent(inout) :: this

    call this%coef_GL%free()

    nullify(this%coef_GLL)

    call this%GLL_to_GL%free()

    call this%Xh_GL%free()

    nullify(this%Xh_GLL)

    call this%dtime%free()

    call this%cr_GL%free()
    call this%cs_GL%free()
    call this%ct_GL%free()

    call this%convr_GL%free()
    call this%convs_GL%free()
    call this%convt_GL%free()

    call this%conv_k1%free()
    call this%conv_k23%free()
    call this%conv_k4%free()

    if (associated(this%cr_k1)) then
       deallocate(this%cr_k1)
    end if
    if (associated(this%cs_k1)) then
       deallocate(this%cs_k1)
    end if
    if (associated(this%ct_k1)) then
       deallocate(this%ct_k1)
    end if
    if (associated(this%cr_k23)) then
       deallocate(this%cr_k23)
    end if
    if (associated(this%cs_k23)) then
       deallocate(this%cs_k23)
    end if
    if (associated(this%ct_k23)) then
       deallocate(this%ct_k23)
    end if
    if (associated(this%cr_k4)) then
       deallocate(this%cr_k4)
    end if
    if (associated(this%cs_k4)) then
       deallocate(this%cs_k4)
    end if
    if (associated(this%ct_k4)) then
       deallocate(this%ct_k4)
    end if

    nullify(this%ulag)
    nullify(this%vlag)
    nullify(this%wlag)
    nullify(this%slag)
    nullify(this%ctlag)
    nullify(this%dctlag)
    nullify(this%oifs_scheme)
    nullify(this%cr_k1)
    nullify(this%cs_k1)
    nullify(this%ct_k1)
    nullify(this%cr_k23)
    nullify(this%cs_k23)
    nullify(this%ct_k23)
    nullify(this%cr_k4)
    nullify(this%cs_k4)
    nullify(this%ct_k4)

    if (allocated(this%cx)) then
       deallocate(this%cx)
    end if
    if (allocated(this%cy)) then
       deallocate(this%cy)
    end if
    if (allocated(this%cz)) then
       deallocate(this%cz)
    end if
    if (c_associated(this%cx_d)) then
       call device_free(this%cx_d)
    end if
    if (c_associated(this%cy_d)) then
       call device_free(this%cy_d)
    end if
    if (c_associated(this%cz_d)) then
       call device_free(this%cz_d)
    end if

  end subroutine adv_oifs_free

  !> Mapping the velocity fields to GL space and transforming them to the rst format
  !! @param u Velocity component in x-direction
  !! @param v Velocity component in y-direction
  !! @param w Velocity component in z-direction
  !! @note similar to set_ct_cvx in NEK5000
  subroutine set_conv_velocity_fst(this, u, v, w)
    implicit none
    class(adv_oifs_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w
    integer :: i, nel, n_GL, idx, idy, idz

    nel = this%coef_GLL%msh%nelv
    n_GL = nel*this%Xh_GL%lxyz

    call this%convr_GL%update()
    call this%convs_GL%update()
    call this%convt_GL%update()

    call this%GLL_to_GL%map(this%cx, u%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cy, v%x, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%cz, w%x, nel, this%Xh_GL)

    call set_convect_rst(this%cr_GL, this%cs_GL, this%ct_GL, &
         this%cx, this%cy, this%cz, this%Xh_GL, this%coef_GL)

    this%convr_GL%f = this%cr_GL
    this%convs_GL%f = this%cs_GL
    this%convt_GL%f = this%ct_GL

  end subroutine set_conv_velocity_fst


  !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
  !! the RHS using the OIFS method.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step.
  subroutine adv_oifs_compute(this, vx, vy, vz, fx, fy, fz, Xh, coef, n, dt)
    implicit none
    class(adv_oifs_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    real(kind=rp) :: tau, tau1, th, dtau
    integer :: i, ilag, itau, nel, n_GL

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    associate(ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         ctlag => this%ctlag, dctlag => this%dctlag, dtime => this%dtime, &
         Xh_GL => this%Xh_GL, coef_GL => this%coef_GL, ntaubd => this%ntaubd, &
         GLL_to_GL => this%GLL_to_GL, oifs_scheme => this%oifs_scheme, &
         cr_k1 => this%cr_K1, cs_k1 => this%cs_K1, ct_k1 => this%ct_K1, &
         cr_k23 => this%cr_K23, cs_k23 => this%cs_K23, ct_k23 => this%ct_K23, &
         cr_k4 => this%cr_K4, cs_k4 => this%cs_K4, ct_k4 => this%ct_K4, &
         convr_GL => this%convr_GL, convs_GL => this%convs_GL, &
         convt_GL => this%convt_GL, conv_k1 => this%conv_k1, &
         conv_k23 => this%conv_k23, conv_k4 => this%conv_k4)

      call dtime%init(oifs_scheme%ndiff)

      tau = ctlag(oifs_scheme%ndiff)

      call this%set_conv_velocity_fst(vx, vy, vz)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_rzero(fx%x_d,n)
         call device_rzero(fy%x_d,n)
         call device_rzero(fz%x_d,n)
      else
         call rzero(fx%x,n)
         call rzero(fy%x,n)
         call rzero(fz%x,n)
      end if

      do ilag = oifs_scheme%ndiff, 1, -1
         if (NEKO_BCKND_DEVICE .eq. 1) then
            if (ilag .eq. 1) then
               call device_addcol3s2(fx%x_d, vx%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(2), n)
               call device_addcol3s2(fy%x_d, vy%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(2), n)
               call device_addcol3s2(fz%x_d, vz%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(2), n)
            else
               call device_addcol3s2(fx%x_d, ulag%lf(ilag-1)%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(ilag+1), n)
               call device_addcol3s2(fy%x_d, vlag%lf(ilag-1)%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(ilag+1), n)
               call device_addcol3s2(fz%x_d, wlag%lf(ilag-1)%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(ilag+1), n)
            end if
         else
            if (ilag .eq. 1) then
               do i = 1, n
                  fx%x(i,1,1,1) = fx%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(2) &
                       * vx%x(i,1,1,1) * coef%B(i,1,1,1)
                  fy%x(i,1,1,1) = fy%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(2) &
                       * vy%x(i,1,1,1) * coef%B(i,1,1,1)
                  fz%x(i,1,1,1) = fz%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(2) &
                       * vz%x(i,1,1,1) * coef%B(i,1,1,1)
               end do
            else
               do i = 1, n
                  fx%x(i,1,1,1) = fx%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(ilag+1) &
                       * ulag%lf(ilag-1)%x(i,1,1,1) &
                       * coef%B(i,1,1,1)
                  fy%x(i,1,1,1) = fy%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(ilag+1) &
                       * vlag%lf(ilag-1)%x(i,1,1,1) &
                       * coef%B(i,1,1,1)
                  fz%x(i,1,1,1) = fz%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(ilag+1) &
                       * wlag%lf(ilag-1)%x(i,1,1,1) &
                       * coef%B(i,1,1,1)
               end do
            end if
         end if
         dtau = dctlag(ilag)/real(ntaubd)
         do itau = 1, ntaubd
            th = tau + dtau/2.
            tau1 = tau + dtau
            call dtime%interpolate_scalar(tau, cr_k1, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau, cs_k1, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau, ct_k1, convt_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, cr_k23, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, cs_k23, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, ct_k23, convt_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, cr_k4, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, cs_k4, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, ct_k4, convt_GL, ctlag, n_GL)
            call runge_kutta(fx, conv_k1, conv_k23, conv_k4, Xh, Xh_GL, &
                 coef, coef_GL, GLL_to_GL, tau, dtau, &
                 n, nel, n_GL)
            call runge_kutta(fy, conv_k1, conv_k23, conv_k4, Xh, Xh_GL, &
                 coef, coef_GL, GLL_to_GL, tau, dtau, &
                 n, nel, n_GL)
            call runge_kutta(fz, conv_k1, conv_k23, conv_k4, Xh, Xh_GL, &
                 coef, coef_GL, GLL_to_GL, tau, dtau, &
                 n, nel, n_GL)
            tau = tau1
         end do
      end do

    end associate

  end subroutine adv_oifs_compute
  !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to the
  !! RHS.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param s The scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step.
  subroutine adv_oifs_compute_scalar(this, vx, vy, vz, s, fs, Xh, coef, n, dt)
    implicit none
    class(adv_oifs_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fs
    type(field_t), intent(inout) :: s
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    real(kind=rp) :: tau, tau1, th, dtau
    integer :: i, ilag, itau, nel, n_GL
    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    associate(slag => this%slag, ctlag => this%ctlag, dctlag => this%dctlag, &
         dtime => this%dtime, Xh_GL => this%Xh_GL, coef_GL => this%coef_GL, &
         ntaubd => this%ntaubd, GLL_to_GL => this%GLL_to_GL, &
         oifs_scheme => this%oifs_scheme, cr_k1 => this%cr_K1, &
         cs_k1 => this%cs_K1, ct_k1 => this%ct_K1, cr_k23 => this%cr_K23, &
         cs_k23 => this%cs_K23, ct_k23 => this%ct_K23, cr_k4 => this%cr_K4, &
         cs_k4 => this%cs_K4, ct_k4 => this%ct_K4, &
         convr_GL => this%convr_GL, convs_GL => this%convs_GL, &
         convt_GL => this%convt_GL, conv_k1 => this%conv_k1, &
         conv_k23 => this%conv_k23, conv_k4 => this%conv_k4)

      call dtime%init(oifs_scheme%ndiff)

      tau = ctlag(oifs_scheme%ndiff)

      call this%set_conv_velocity_fst(vx, vy, vz)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_rzero(fs%x_d,n)
      else
         call rzero(fs%x,n)
      end if

      do ilag = oifs_scheme%ndiff, 1, -1
         if (NEKO_BCKND_DEVICE .eq. 1) then
            if (ilag .eq. 1) then
               call device_addcol3s2(fs%x_d, s%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(2), n)
            else
               call device_addcol3s2(fs%x_d, slag%lf(ilag-1)%x_d, coef%B_d, &
                    oifs_scheme%diffusion_coeffs(ilag+1), n)
            end if
         else
            if (ilag .eq. 1) then
               do i = 1, n
                  fs%x(i,1,1,1) = fs%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(2) &
                       * s%x(i,1,1,1) * coef%B(i,1,1,1)
               end do
            else
               do i = 1, n
                  fs%x(i,1,1,1) = fs%x(i,1,1,1) + &
                       oifs_scheme%diffusion_coeffs(ilag+1) &
                       * slag%lf(ilag-1)%x(i,1,1,1) * coef%B(i,1,1,1)
               end do
            end if
         end if
         dtau = dctlag(ilag)/real(ntaubd)
         do itau = 1, ntaubd
            th = tau + dtau/2.
            tau1 = tau + dtau
            call dtime%interpolate_scalar(tau, cr_k1, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau, cs_k1, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau, ct_k1, convt_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, cr_k23, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, cs_k23, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(th, ct_k23, convt_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, cr_k4, convr_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, cs_k4, convs_GL, ctlag, n_GL)
            call dtime%interpolate_scalar(tau1, ct_k4, convt_GL, ctlag, n_GL)
            call runge_kutta(fs, conv_k1, conv_k23, conv_k4, Xh, Xh_GL, &
                 coef, coef_GL, GLL_to_GL, tau, dtau, &
                 n, nel, n_GL)
            tau = tau1
         end do
      end do

    end associate

  end subroutine adv_oifs_compute_scalar

end module adv_oifs
