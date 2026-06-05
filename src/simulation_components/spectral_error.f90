! Copyright (c) 2022, The Neko Authors
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
!> Implements type spectral_error_t.
module spectral_error
  use num_types, only : rp, dp
  use field, only : field_t
  use coefs, only : coef_t
  use field_list, only : field_list_t
  use math, only : rzero, copy, add3s2
  use file, only : file_t, file_free
  use time_state, only : time_state_t
  use tensor, only : tnsr3d
  use device_math, only : device_copy
  use neko_config, only : NEKO_BCKND_HIP, NEKO_BCKND_CUDA, NEKO_BCKND_OPENCL, &
       NEKO_BCKND_DEVICE
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use device, only : DEVICE_TO_HOST, HOST_TO_DEVICE, device_memcpy
  use comm, only : pe_rank
  use utils, only : NEKO_FNAME_LEN, NEKO_VARNAME_LEN, neko_error, neko_warning
  use field_writer, only : field_writer_t
  use simulation_component, only : simulation_component_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use case, only : case_t
  use registry, only : neko_registry
  use vector_list, only : vector_list_t
  use scratch_registry, only : neko_scratch_registry
  use amr_reconstruct, only : amr_reconstruct_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Provides tools to calculate the spectral error indicator
  !! @details
  !! This is a posteriori error measure, based on the local properties of
  !! the spectral solution, which was developed by Mavriplis. This method
  !! formally only gives an indication of the error.
  type, public, extends(simulation_component_t) :: spectral_error_t
     !> Pointers to main fields
     type(field_list_t) :: field
     !> Transformed fields
     type(field_list_t) :: field_hat
     !> Spectral error indicator per element
     type(vector_list_t) :: eind
     !> Fit coefficients per element
     type(vector_list_t) :: sig
     !> Averaged spectral error indicator per element
     type(vector_list_t) :: eind_av
     !> Averaged fit coefficients per element
     type(vector_list_t) :: sig_av
     !> Number of fields
     integer :: nfld
     !> Field names
     character(NEKO_VARNAME_LEN), dimension(:), allocatable :: field_names
     !> Vector length
     integer :: nelv
     !> Simulation time of the first averaging step
     real(dp) :: time_start
     !> Simulation time of previous averaging step
     real(dp) :: time_previous
     !> Restart file name
     character(:), allocatable :: restart_file
     !> Configuration of spectral error calculation
     real(kind=rp) :: SERI_SMALL = 1.e-14
     !> used for ratios
     real(kind=rp) :: SERI_SMALLR = 1.e-10
     !> used for gradients
     real(kind=rp) :: SERI_SMALLG = 1.e-5
     !> used for sigma and rtmp in error calculations
     real(kind=rp) :: SERI_SMALLS = 0.2
     !> number of points in fitting
     integer :: SERI_NP = 4
     integer :: SERI_NP_MAX = 4
     !> last modes skipped
     integer :: SERI_ELR = 0
   contains
     !> Constructor.
     procedure, pass(this) :: init => spectral_error_init
     !> Destructor.
     procedure, pass(this) :: free => spectral_error_free
     !> Compute the indicator.
     procedure, pass(this) :: compute_ => spectral_error_compute
     !> Restart the simcomp.
     procedure, pass(this) :: restart_ => spectral_error_restart
     !> Set start simulation time for averaging
     procedure, pass(this) :: set_time_start => spectral_error_set_time_start
     !> Set previous simulation time for averaging
     procedure, pass(this) :: set_time_previous => &
          spectral_error_set_time_previous
     !> Get averaging time
     procedure, pass(this) :: get_collect_time => &
          spectral_error_get_collect_time
     !> Reset averaged variables to zero
     procedure, pass(this) :: average_reset => spectral_error_average_reset
     !> AMR restart
     procedure, pass(this) :: amr_restart => spectral_error_amr_restart
  end type spectral_error_t

contains

  !> Constructor.
  subroutine spectral_error_init(this, json, case)
    class(spectral_error_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=NEKO_VARNAME_LEN), allocatable, dimension(:) :: field_names
    character(len=:), allocatable :: name, restart_file
    real(rp) :: start_time

    call this%free()

    call json_get_or_default(json, "name", name, "spectral_error")
    this%name = trim(name)
    call json_get_or_default(json, "restart_file", restart_file, "no restart")
    this%restart_file = trim(restart_file)

    call json_get(json, "fields", field_names)
    if (.not. allocated(field_names)) &
         call neko_error('Spectral error: missing fields')

    call json_get_or_default(json, 'start_time', start_time, 0.0_rp)

    call this%init_base(json, case)

    call spectral_error_init_from_components(this, case%fluid%c_Xh, &
         field_names, start_time)

    deallocate(field_names)

  end subroutine spectral_error_init

  !> Actual constructor.
  !! @param coef type with all geometrical variables.
  subroutine spectral_error_init_from_components(this, coef, field_names, &
       start_time)
    class(spectral_error_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=NEKO_VARNAME_LEN), dimension(:), intent(in) :: field_names
    real(rp), intent(in) :: start_time
    integer :: il, jl, nfld

    ! array sizes
    this%nfld = size(field_names)
    this%nelv = coef%msh%nelv

    ! save field names
    allocate(this%field_names(this%nfld))
    do il = 1, this%nfld
       this%field_names(il) = trim(field_names(il))
    end do

    ! allocate pointer space
    call this%field%init(this%nfld)
    call this%field_hat%init(this%nfld)
    call this%eind%init(this%nfld)
    call this%sig%init(this%nfld)
    call this%eind_av%init(this%nfld)
    call this%sig_av%init(this%nfld)

    ! add missing fields and vectors to registry
    do il = 1, this%nfld
       call neko_registry%add_field(coef%dof, trim(field_names(il))//'_hat', &
            ignore_existing = .false.)
       call neko_registry%add_vector(this%nelv, &
            trim(field_names(il))//'_eind', ignore_existing = .false.)
       call neko_registry%add_vector(this%nelv, &
            trim(field_names(il))//'_sig', ignore_existing = .false.)
       call neko_registry%add_vector(this%nelv, &
            trim(field_names(il))//'_eind_av', ignore_existing = .false.)
       call neko_registry%add_vector(this%nelv, &
            trim(field_names(il))//'_sig_av', ignore_existing = .false.)
    end do

    ! get pointers
    do il = 1, this%nfld
       this%field%items(il)%ptr => &
            neko_registry%get_field_by_name(trim(field_names(il)))
       this%field_hat%items(il)%ptr => &
            neko_registry%get_field_by_name(trim(field_names(il))//'_hat')
       this%eind%items(il)%ptr => &
            neko_registry%get_vector_by_name(trim(field_names(il))//'_eind')
       this%sig%items(il)%ptr => &
            neko_registry%get_vector_by_name(trim(field_names(il))//'_sig')
       this%eind_av%items(il)%ptr => &
            neko_registry%get_vector_by_name(trim(field_names(il))//'_eind_av')
       this%sig_av%items(il)%ptr => &
            neko_registry%get_vector_by_name(trim(field_names(il))//'_sig_av')
    end do

    ! zero averaged variables
    call this%average_reset()

    ! set start and previous time
    this%time_start = start_time
    this%time_previous = this%time_start

    !> The following code has been lifted from Adam's implementation
    associate(LX1 => coef%Xh%lx, LY1 => coef%Xh%ly, &
         LZ1 => coef%Xh%lz, &
         SERI_SMALL => this%SERI_SMALL, &
         SERI_SMALLR => this%SERI_SMALLR, &
         SERI_SMALLG => this%SERI_SMALLG, &
         SERI_SMALLS => this%SERI_SMALLS, &
         SERI_NP => this%SERI_NP, &
         SERI_NP_MAX => this%SERI_NP_MAX, &
         SERI_ELR => this%SERI_ELR &
         )
      ! correctness check
      if (SERI_NP .gt. SERI_NP_MAX) then
         call neko_log%message('SETI_NP greater than SERI_NP_MAX')
      end if
      il = SERI_NP + SERI_ELR
      jl = min(LX1, LY1)
      jl = min(jl, LZ1)
      if (il .gt. jl) then
         call neko_log%message('SERI_NP+SERI_ELR greater than L?1')
      end if
    end associate

  end subroutine spectral_error_init_from_components

  !> Destructor
  subroutine spectral_error_free(this)
    class(spectral_error_t), intent(inout) :: this

    this%nfld = 0
    this%nelv = 0
    this%time_start = 0.0_dp
    this%time_previous = 0.0_dp

    if (allocated(this%field_names)) deallocate(this%field_names)
    if (allocated(this%restart_file)) deallocate(this%restart_file)

    call this%field%free()
    call this%field_hat%free()
    call this%eind%free()
    call this%sig%free()
    call this%eind_av%free()
    call this%sig_av%free()

    call this%free_base()

    call this%free_amr_base()

  end subroutine spectral_error_free

  !> Compute the spectral error indicator.
  subroutine spectral_error_compute(this, time)
    class(spectral_error_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: wk
    integer :: il, idx
    real(rp) :: alpha, beta
    real(dp) :: collect, increment

    ! Check consistency
    if (this%time_start .gt. this%time_previous) &
         call neko_error('Spectral error; time_start grater than &
         &time_previous')

    ! Skip steps wit simulation steps smaller than time_previous
    if (time%t .le. this%time_previous) return

    call neko_scratch_registry%request_field(wk, idx, .false.)

    associate(coef => this%case%fluid%c_Xh)
      ! Generate the field_hat (legendre coeff)
      do il = 1, this%nfld
         call transform_to_spec_or_phys(this%field_hat%items(il)%ptr, &
              this%field%items(il)%ptr, wk, coef, 'spec')
      end do

      ! Get the spectral error indicator
      do il = 1, this%nfld
         call calculate_indicators(this, coef, this%eind%items(il)%ptr%x, &
              this%sig%items(il)%ptr%x, this%nelv, coef%Xh%lx, coef%Xh%ly, &
              coef%Xh%lz, this%field_hat%items(il)%ptr%x)
      end do

      collect = time%t - this%time_start
      increment = time%t - this%time_previous
      this%time_previous = time%t
      if (collect .gt. 0.0_dp) then
         beta = real(increment / collect, rp)
      else
         beta = 1.0_rp
      end if
      alpha = 1.0_rp - beta

      ! Get time averages
      do il = 1, this%nfld
         call add3s2(this%eind_av%items(il)%ptr%x, &
              this%eind_av%items(il)%ptr%x, this%eind%items(il)%ptr%x, &
              alpha, beta, this%nelv)
      end do
      do il = 1, this%nfld
         call add3s2(this%sig_av%items(il)%ptr%x, &
              this%sig_av%items(il)%ptr%x, this%sig%items(il)%ptr%x, &
              alpha, beta, this%nelv)
      end do

    end associate

    call neko_scratch_registry%relinquish_field(idx)

  end subroutine spectral_error_compute

  !> Restart the simcomp.
  subroutine spectral_error_restart(this, time)
    class(spectral_error_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (trim(this%restart_file) .eq. "no restart") then
       if (time%t .gt. this%time_start) then
          this%time_start = time%t
          this%time_previous = this%time_start
       end if
    else
       call neko_error('Spectral error restart not done yet')
       this%time_previous = time%t
    end if
  end subroutine spectral_error_restart

  !> Set start simulation time for averaging
  subroutine spectral_error_set_time_start(this, time)
    class(spectral_error_t), intent(inout) :: this
    real(dp), intent(in) :: time

    this%time_start = time
    if (this%time_start .gt. this%time_previous) &
         call neko_warning('Spectral error; time_start grater than &
         &time_previous')

  end subroutine spectral_error_set_time_start

  !> Set previous simulation time for averaging
  subroutine spectral_error_set_time_previous(this, time)
    class(spectral_error_t), intent(inout) :: this
    real(dp), intent(in) :: time

    this%time_previous = time
    if (this%time_start .gt. this%time_previous) &
         call neko_warning('Spectral error; time_start grater than &
         &time_previous')

  end subroutine spectral_error_set_time_previous

  !> Get averaging time
  pure function spectral_error_get_collect_time(this) result(time)
    class(spectral_error_t), intent(in) :: this
    real(dp) :: time

    time = this%time_previous - this%time_start

  end function spectral_error_get_collect_time

  !> Reset averaged variables to zero
  subroutine spectral_error_average_reset(this)
    class(spectral_error_t), intent(inout) :: this
    integer :: il

    ! Reset averaged variables
    do il = 1, this%nfld
       if (associated(this%eind_av%items(il)%ptr)) then
          if (allocated(this%eind_av%items(il)%ptr%x)) then
             this%eind_av%items(il)%ptr%x(:) = 0.0_rp
          end if
       end if
    end do
    do il = 1, this%nfld
       if (associated(this%sig_av%items(il)%ptr)) then
          if (allocated(this%sig_av%items(il)%ptr%x)) then
             this%sig_av%items(il)%ptr%x(:) = 0.0_rp
          end if
       end if
    end do

  end subroutine spectral_error_average_reset

  !> Transform a field u to u_hat into physical or spectral space
  !! the result of the transformation is in u_hat.
  !! @param u_hat Transformed field (output).
  !! @param u Field to transform (input).
  !! @param wk Working field.
  !! @param coef Type coef for mesh parameters.
  !! @param space String that indicates which space to transform, "spec" or
  !! "phys".
  subroutine transform_to_spec_or_phys(u_hat, u, wk, coef, space)
    type(field_t), intent(inout) :: u_hat
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: wk
    type(coef_t), intent(inout) :: coef
    character(len=4), intent(in) :: space
    integer :: i, j, k, e, nxyz, nelv, n

    !> Define some constants
    nxyz = coef%Xh%lx * coef%Xh%lx * coef%Xh%lx
    nelv = coef%msh%nelv
    n = nxyz*nelv

    !> Copy field to working array
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_copy(wk%x_d, u%x_d, n)
    else
       call copy(wk%x, u%x, n)
    end if

    select case (space)
    case ('spec')
       call tnsr3d(u_hat%x, coef%Xh%lx, wk%x, &
            coef%Xh%lx, coef%Xh%vinv, &
            coef%Xh%vinvt, coef%Xh%vinvt, nelv)
    case ('phys')
       call tnsr3d(u_hat%x, coef%Xh%lx, wk%x, &
            coef%Xh%lx, coef%Xh%v, &
            coef%Xh%vt, coef%Xh%vt, nelv)
    end select

    ! Synchronize
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then

       call device_memcpy(u_hat%x, u_hat%x_d, n, &
            DEVICE_TO_HOST, sync = .true.)
    end if

  end subroutine transform_to_spec_or_phys

  !> Wrapper for old fortran 77 subroutines
  !! @param coef coef type
  !! @param eind spectral indicator
  !! @param sig coefficient of the exponential fit
  !! @param lnelt number of elements
  !! @param LX1 gll points in x
  !! @param LY1 gll points in y
  !! @param LZ1 gll points in z
  !! @paran var variable to calculate indicator
  subroutine calculate_indicators(this, coef, eind, sig, lnelt, LX1, LY1, LZ1, &
       var)
    type(spectral_error_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: lnelt
    integer :: LX1
    integer :: LY1
    integer :: LZ1
    real(kind=rp) :: eind(lnelt)
    real(kind=rp) :: sig(lnelt)
    real(kind=rp) :: var(LX1, LY1, LZ1, lnelt)

    real(kind=rp) :: xa(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz)
    real(kind=rp) :: xb(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz)
    integer :: i, e

    ! zero arrays
    call rzero(eind, lnelt)
    call rzero(sig, lnelt)

    ! Get the indicator
    call speri_var(this, eind, sig, var, lnelt, xa, xb, LX1, LY1, LZ1)

  end subroutine calculate_indicators

  !> Calculate the indicator in a specified variable
  !! @param est spectral indicator
  !! @param sig coefficient of the exponential fit
  !! @param nell number of elements
  !! @paran var variable to calculate indicator
  !! @paran xa work array
  !! @paran xb work array
  !! @param LX1 gll points in x
  !! @param LY1 gll points in y
  !! @param LZ1 gll points in z
  subroutine speri_var(this, est, sig, var, nell, xa, xb, LX1, LY1, LZ1)
    type(spectral_error_t), intent(inout) :: this
    integer :: nell
    integer :: LX1
    integer :: LY1
    integer :: LZ1
    real(kind=rp) :: est(nell)
    real(kind=rp) :: sig(nell)
    real(kind=rp) :: var(LX1, LY1, LZ1, nell)
    real(kind=rp) :: xa(LX1, LY1, LZ1)
    real(kind=rp) :: xb(LX1, LY1, LZ1)

    !> local variables
    integer :: il, jl, kl, ll, j_st, j_en, ii
    !> polynomial coefficients
    real(kind=rp) :: coeff(LX1, LY1, LZ1)
    !> Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) :: coef11
    !> copy of last SERI_NP columns of coefficients
    real(kind=rp) :: coefx(this%SERI_NP_MAX, LY1, LZ1), &
         coefy(this%SERI_NP_MAX, LX1, LZ1), &
         coefz(this%SERI_NP_MAX, LX1, LY1)
    !> estimated error
    real(kind=rp) :: estx, esty, estz
    !> estimated decay rate
    real(kind=rp) :: sigx, sigy, sigz
    real(kind=rp) :: third
    parameter (third = 1.0/3.0)

    !> loop over elements
    do il = 1, nell
       !> go to Legendre space (done in two operations)
       !! and square the coefficient
       do ii = 1, LX1*LY1*LZ1
          coeff(ii, 1, 1) = var(ii, 1, 1, il) * var(ii, 1, 1, il)
       end do

       !> lower left corner
       coef11 = coeff(1, 1, 1)

       !> small value; nothing to od
       if (coef11 .ge. this%SERI_SMALL) then
          !> extrapolate coefficients
          !> X - direction
          !> copy last SERI_NP collumns (or less if NX1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1, LX1 - this%SERI_NP + 1 - this%SERI_ELR)
          j_en = max(1, LX1 - this%SERI_ELR)
          do ll = 1, LZ1
             do kl = 1, LY1
                do jl = j_st, j_en
                   coefx(j_en - jl + 1, kl, ll) = coeff(jl, kl, ll)
                end do
             end do
          end do
          !> get extrapolated values
          call speri_extrap(this, estx, sigx, coef11, coefx, &
               j_st, j_en, LY1, LZ1)

          !> Y - direction
          !> copy last SERI_NP collumns (or less if NY1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1, LY1 - this%SERI_NP + 1 - this%SERI_ELR)
          j_en = max(1, LY1 - this%SERI_ELR)
          do ll = 1, LZ1
             do kl = j_st, j_en
                do jl = 1, LX1
                   coefy(j_en - kl + 1, jl, ll) = coeff(jl, kl, ll)
                end do
             end do
          end do
          !> get extrapolated values
          call speri_extrap(this, esty, sigy, coef11, coefy, &
               j_st, j_en, LX1, LZ1)

          !> Z - direction
          !> copy last SERI_NP collumns (or less if NZ1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1, LZ1 - this%SERI_NP + 1 - this%SERI_ELR)
          j_en = max(1, LZ1 - this%SERI_ELR)
          do ll = j_st, j_en
             do kl = 1, LY1
                do jl = 1, LX1
                   coefz(j_en - ll + 1, jl, kl) = coeff(jl, kl, ll)
                end do
             end do
          end do
          !> get extrapolated values
          call speri_extrap(this, estz, sigz, coef11, coefz, &
               j_st, j_en, LX1, LY1)

          !> average
          est(il) = sqrt(estx + esty + estz)
          sig(il) = third*(sigx + sigy + sigz)

       else
          !> for testing
          estx = 0.0
          esty = 0.0
          estz = 0.0
          sigx = -1.0
          sigy = -1.0
          sigz = -1.0
          !> for testing; end

          est(il) = 0.0
          sig(il) = -1.0
       end if

    end do

  end subroutine speri_var

  !> Extrapolate the Legendre spectrum from the last points
  !! @param estx spectral indicator
  !! @param sigx coefficient of the exponential fit
  !! @param coef11 legendre coefficients
  !! @paran coef legendre coefficients
  !! @paran ix_st argument list
  !! @paran ix_en argument list
  !! @param nyl argument list
  !! @param nzl argument list
  subroutine speri_extrap(this, estx, sigx, coef11, coef, &
       ix_st, ix_en, nyl, nzl)
    implicit none
    type(spectral_error_t), intent(inout) :: this
    !> argument list
    integer :: ix_st, ix_en, nyl, nzl
    !> Legendre coefficients; last SERI_NP columns
    real(kind=rp) :: coef(this%SERI_NP_MAX, nyl, nzl)
    !> Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) :: coef11
    !> estimated error and decay rate
    real(kind=rp) :: estx, sigx

    !> local variables
    integer :: il, jl, kl, ll ! loop index
    integer :: nsigt, pnr, nzlt
    real(kind=rp) :: sigt, smallr, cmin, cmax, cnm, rtmp, rtmp2, rtmp3
    real(kind=rp) :: sumtmp(4), cffl(this%SERI_NP_MAX)
    real(kind=rp) :: stmp, estt, clog, ctmp, cave, erlog
    logical :: cuse(this%SERI_NP_MAX)

    associate(SERI_SMALL => this%SERI_SMALL, &
         SERI_SMALLR => this%SERI_SMALLR, &
         SERI_SMALLG => this%SERI_SMALLG, &
         SERI_SMALLS => this%SERI_SMALLS, &
         SERI_NP => this%SERI_NP, &
         SERI_NP_MAX => this%SERI_NP_MAX, &
         SERI_ELR => this%SERI_ELR &
         )
      ! initial values
      estx = 0.0
      sigx = -1.0

      ! relative cutoff
      smallr = coef11*SERI_SMALLR

      ! number of points
      pnr = ix_en - ix_st + 1

      ! to few points to interpolate
      !if ((ix_en - ix_st).le.1) return

      ! for averaging, initial values
      sigt = 0.0
      nsigt = 0

      ! loop over all face points
      nzlt = max(1, nzl - SERI_ELR) !  for 2D runs
      do il = 1, nzlt
         ! weight
         rtmp3 = 1.0/(2.0*(il - 1) + 1.0)
         do jl = 1, nyl - SERI_ELR

            ! find min and max coef along single row
            cffl(1) = coef(1, jl, il)
            cmin = cffl(1)
            cmax = cmin
            do kl = 2, pnr
               cffl(kl) = coef(kl, jl, il)
               cmin = min(cmin, cffl(kl))
               cmax = max(cmax, cffl(kl))
            end do

            ! are coefficients sufficiently big
            if ((cmin .gt. 0.0) .and. (cmax .gt. smallr)) then
               ! mark array position we use in iterpolation
               do kl = 1, pnr
                  cuse(kl) = .TRUE.
               end do
               ! max n for polynomial order
               cnm = real(ix_en)

               ! check if all the points should be taken into account
               ! in original code by Catherine Mavriplis this part is written
               ! for 4 points, so I place if statement first
               if (pnr .eq. 4) then
                  ! should we neglect last values
                  if ((cffl(1) .lt. smallr) .and. &
                       (cffl(2) .lt. smallr)) then
                     if (cffl(3) .lt. smallr) then
                        cuse(1) = .FALSE.
                        cuse(2) = .FALSE.
                        cnm = real(ix_en - 2)
                     else
                        cuse(1) = .FALSE.
                        cnm = real(ix_en - 1)
                     end if
                  else
                     ! should we take stronger gradient
                     if ((cffl(1)/cffl(2) .lt. SERI_SMALLG) .and. &
                          (cffl(3)/cffl(4) .lt. SERI_SMALLG)) then
                        cuse(1) = .FALSE.
                        cuse(3) = .FALSE.
                        cnm = real(ix_en - 1)
                     elseif ((cffl(2)/cffl(1) .lt. SERI_SMALLG) .and. &
                          (cffl(4)/cffl(3) .lt. SERI_SMALLG)) then
                        cuse(2) = .FALSE.
                        cuse(4) = .FALSE.
                     end if
                  end if
               end if

               ! get sigma for given face point
               do kl = 1, 4
                  sumtmp(kl) = 0.0
               end do
               ! find new min and count number of points
               cmin = cmax
               cmax = 0.0
               do kl = 1, pnr
                  if (cuse(kl)) then
                     rtmp = real(ix_en - kl)
                     rtmp2 = log(cffl(kl))
                     sumtmp(1) = sumtmp(1) + rtmp2
                     sumtmp(2) = sumtmp(2) + rtmp
                     sumtmp(3) = sumtmp(3) + rtmp*rtmp
                     sumtmp(4) = sumtmp(4) + rtmp2*rtmp
                     ! find new min and count used points
                     cmin = min(cffl(kl), cmin)
                     cmax = cmax + 1.0
                  end if
               end do
               ! decay rate along single row
               stmp = (sumtmp(1)*sumtmp(2) - sumtmp(4)*cmax)/ &
                    (sumtmp(3)*cmax - sumtmp(2)*sumtmp(2))
               ! for averaging
               sigt = sigt + stmp
               nsigt = nsigt + 1

               ! get error estimator depending on calculated decay rate
               estt = 0.0
               if (stmp .lt. SERI_SMALLS) then
                  estt = cmin
               else
                  ! get averaged constant in front of c*exp(-sig*n)
                  clog = (sumtmp(1)+stmp*sumtmp(2))/cmax
                  ctmp = exp(clog)
                  ! average exponent
                  cave = sumtmp(1)/cmax
                  ! check quality of approximation comparing is to the
                  ! constant cave
                  do kl = 1, 2
                     sumtmp(kl) = 0.0
                  end do
                  do kl = 1, pnr
                     if (cuse(kl)) then
                        erlog = clog - stmp*real(ix_en - kl)
                        sumtmp(1) = sumtmp(1) + &
                             (erlog - log(cffl(kl)))**2
                        sumtmp(2) = sumtmp(2) + &
                             (erlog - cave)**2
                     end if
                  end do
                  rtmp = 1.0 - sumtmp(1)/sumtmp(2)
                  if (rtmp .lt. SERI_SMALLS) then
                     estt = cmin
                  else
                     ! last coefficient is not included in error estimator
                     estt = ctmp/stmp*exp(-stmp*cnm)
                  end if
               end if
               ! add contribution to error estimator; variable weight
               estx = estx + estt/(2.0*(jl - 1) + 1.0)*rtmp3
            end if ! if((cmin.gt.0.0).and.(cmax.gt.smallr))
         end do
      end do
      ! constant weight
      ! Multiplication by 4 in 2D / 8 in 3D
      ! Normalization of the error by the volume of the reference element
      ! which is equal to 4 in 2D / 8 in 3D
      ! ==> Both operations cancel each other
      estx = estx/(2.0*(ix_en - 1) + 1.0)

      ! final everaging
      ! sigt = 2*sigma so we divide by 2
      if (nsigt .gt. 0) then
         sigx = 0.5*sigt/nsigt
      end if

    end associate

  end subroutine speri_extrap

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine spectral_error_amr_restart(this, reconstruct, counter, tstep)
    class(spectral_error_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    character(len=LOG_SIZE) :: log_buf
    integer :: il

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    log_buf = trim(this%name)
    call neko_log%section(log_buf, NEKO_LOG_VERBOSE)

    ! These should be already restarted, but AMR restart prevents
    ! recursive restarting, so it is safe to call it here
    call this%field%amr_restart(reconstruct, counter, tstep)

    ! These I reallocate here assuming former values do not matter
    call this%field_hat%amr_reallocate(reconstruct, counter, tstep)

    ! Reallocate arrays
    if (reconstruct%nold .ne. reconstruct%nnew) then

       this%nelv = reconstruct%nnew

       do il = 1, this%nfld
          if (associated(this%eind%items(il)%ptr)) then
             if (allocated(this%eind%items(il)%ptr%x)) then
                deallocate(this%eind%items(il)%ptr%x)
                allocate(this%eind%items(il)%ptr%x(reconstruct%nnew))
             end if
          end if
       end do
       do il = 1, this%nfld
          if (associated(this%sig%items(il)%ptr)) then
             if (allocated(this%sig%items(il)%ptr%x)) then
                deallocate(this%sig%items(il)%ptr%x)
                allocate(this%sig%items(il)%ptr%x(reconstruct%nnew))
             end if
          end if
       end do
       do il = 1, this%nfld
          if (associated(this%eind_av%items(il)%ptr)) then
             if (allocated(this%eind_av%items(il)%ptr%x)) then
                deallocate(this%eind_av%items(il)%ptr%x)
                allocate(this%eind_av%items(il)%ptr%x(reconstruct%nnew))
             end if
          end if
       end do
       do il = 1, this%nfld
          if (associated(this%sig_av%items(il)%ptr)) then
             if (allocated(this%sig_av%items(il)%ptr%x)) then
                deallocate(this%sig_av%items(il)%ptr%x)
                allocate(this%sig_av%items(il)%ptr%x(reconstruct%nnew))
             end if
          end if
       end do
    end if

    ! Reset averaging time
    this%time_start = this%time_previous

    ! Reset averaged variables
    call this%average_reset()

    call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)

  end subroutine spectral_error_amr_restart

end module spectral_error
