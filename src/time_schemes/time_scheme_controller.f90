! Copyright (c) 2023, The Neko Authors
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
!> Compound scheme for the advection and diffusion operators in a transport
!! equation. 
module time_scheme_controller
  use neko_config
  use num_types, only : rp
  use bdf_time_scheme, only: bdf_time_scheme_t
  use ext_time_scheme, only: ext_time_scheme_t
  use ab_time_scheme, only: ab_time_scheme_t
  use device, only : device_free, device_map, device_memcpy, HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !> Implements the logic to compute the time coefficients for the advection and
  !! diffusion operators in a transport equation.
  !! @details
  !! Uses the BDF scheme for the diffusion, where as the term for advection 
  !! the scheme depends on the orders of the BDF and advection schemes.
  !! - Order 1 advection
  !!   - BDF1 for diffusion -> Adam-Bashforth scheme.
  !! - Order 2 advection
  !!   - BDF1 for diffusion -> Adam-Bashforth scheme.
  !!   - BDF2 for diffusion -> Explcit extrapolation scheme.
  !! - Order 3 advection
  !!   - BDF1 for diffusion -> Adam-Bashforth scheme.
  !!   - BDF2 for diffusion -> Modified explcit extrapolation scheme.
  !!   - BDF3 for diffusion -> Explciit extrapolation scheme.
  !! The order of the BDF scheme in the above logic is set by the user, whereas
  !! the advection scheme is set to forward Euler when BDF is order 1, 
  !! and otherwise to a 3rd order scheme (excluding the first 2 timesteps).
  !! This means that some of the options in the above list never get realized,
  !! particularly order 2 and 3 advection for 1st order BDF. They remain in the
  !! code so as to have the orinigal Nek5000 logic in place for possible
  !! adoption in the future.
  !! An important detail here is the handling of the first timesteps where a 
  !! high-order scheme cannot be constructed. The parameters `nadv` and `ndiff`,
  !! which are initialized to 0, hold the current order of the respective
  !! scheme.
  !! @note the advection scheme also applies to source terms. 
  type, public :: time_scheme_controller_t
     type(ext_time_scheme_t) :: ext
     type(ab_time_scheme_t) :: ab
     type(bdf_time_scheme_t) :: bdf
  
     !> Time coefficients for the advection operator
     real(kind=rp) advection_coeffs(4)
     !> Time coefficients for the diffusion operator
     real(kind=rp) diffusion_coeffs(4)
     !> Controls the actual order of the diffusion scheme,
     !!  e.g. 1 at the first time-step
     integer :: ndiff = 0
     !> Controls the actual order of the advection scheme,
     !!  e.g. 1 at the first time-step
     integer :: nadv = 0
     !> Order of the advection scheme
     integer :: advection_time_order = 3
     !> Order of the diffusion scheme
     integer :: diffusion_time_order
     !> Device pointer for `advection_coeffs`
     type(c_ptr) :: advection_coeffs_d = C_NULL_PTR 
     !> Device pointer for `diffusion_coeffs`
     type(c_ptr) :: diffusion_coeffs_d = C_NULL_PTR 
     
   contains
     !> Constructor
     procedure, pass(this) :: init => time_scheme_controller_init
     !> Destructor
     procedure, pass(this) :: free => time_scheme_controller_free
     !> Set the time coefficients
     procedure, pass(this) :: set_coeffs => &
       time_scheme_controller_set_coeffs
  end type time_scheme_controller_t

  contains

  !> Contructor
  !! @param torder Desired order of the scheme: 1, 2, 3.
  !! This sets the order for the diffusion scheme only.
  subroutine time_scheme_controller_init(this, torder)
    implicit none
    class(time_scheme_controller_t) :: this
    integer :: torder 
  
    this%diffusion_time_order = torder
    
    ! Force 1st order advection when diffusion is 1st order
    if (torder .eq. 1) then
       this%advection_time_order = 1
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%advection_coeffs, this%advection_coeffs_d, 4)
       call device_map(this%diffusion_coeffs, this%diffusion_coeffs_d, 4)
    end if
  end subroutine time_scheme_controller_init

  !> Destructor
  subroutine time_scheme_controller_free(this)
    implicit none
    class(time_scheme_controller_t) :: this

    if (c_associated(this%advection_coeffs_d)) then
       call device_free(this%advection_coeffs_d)
    end if
    if (c_associated(this%diffusion_coeffs_d)) then
       call device_free(this%diffusion_coeffs_d)
    end if
  end subroutine time_scheme_controller_free

  !> Set the time coefficients
  !! @details Implements all necessary logic to handle  
  !! @param t Timestep values, first element is the current timestep.
  subroutine time_scheme_controller_set_coeffs(this, dt)
    implicit none
    class(time_scheme_controller_t) :: this
    real(kind=rp), intent(inout), dimension(10) :: dt
    real(kind=rp), dimension(4) :: adv_coeffs_old
    real(kind=rp), dimension(4) :: diff_coeffs_old

    associate( &
      nadv          => this%nadv, &
      ndiff         => this%ndiff, &
      adv_coeffs    => this%advection_coeffs, &
      adv_coeffs_d  => this%advection_coeffs_d, &
      diff_coeffs   => this%diffusion_coeffs, &
      diff_coeffs_d => this%diffusion_coeffs_d)
    
      ! Increment the order of the scheme if below time_order
      ndiff = ndiff + 1
      ndiff = min(ndiff, this%diffusion_time_order)
      nadv = nadv + 1
      nadv = min(nadv, this%advection_time_order)
      
      call this%bdf%compute_coeffs(diff_coeffs, dt, ndiff)
      
      if (nadv .eq. 1) then
         ! Forward euler
         call this%ext%compute_coeffs(adv_coeffs, dt, nadv)
      else if (nadv .eq. 2) then
         if (ndiff .eq. 1) then
            ! 2nd order Adam-Bashforth, currently never used
            call this%ab%compute_coeffs(adv_coeffs, dt, nadv)
         else
            ! Linear extrapolation
            call this%ext%compute_coeffs(adv_coeffs, dt, nadv)
         end if
      else if (nadv .eq. 3) then
         if (ndiff .eq. 1) then
            ! 3rd order Adam-Bashforth, currently never used
            call this%ab%compute_coeffs(adv_coeffs, dt, nadv)
         else if (ndiff .eq. 2) then
            ! The modified EXT scheme
            call this%ext%compute_modified_coeffs(adv_coeffs, dt)
         else
            ! Quadratic extrapolation
            call this%ext%compute_coeffs(adv_coeffs, dt, nadv)
         end if
      end if
      

      if (c_associated(adv_coeffs_d)) then
         if (maxval(abs(adv_coeffs - adv_coeffs_old)) .gt. 1e-10_rp) then
            call device_memcpy(adv_coeffs, adv_coeffs_d, 4, HOST_TO_DEVICE)
         end if
      end if

      if (c_associated(diff_coeffs_d)) then
         if (maxval(abs(diff_coeffs - diff_coeffs_old)) .gt. 1e-10_rp) then
            call device_memcpy(diff_coeffs, diff_coeffs_d, 4, HOST_TO_DEVICE)
         end if
      end if
    end associate

  end subroutine time_scheme_controller_set_coeffs
end module time_scheme_controller