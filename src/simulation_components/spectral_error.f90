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
  use num_types, only: rp
  use field, only: field_t
  use coefs, only: coef_t
  use field_list, only: field_list_t
  use math, only: rzero, copy
  use file, only: file_t, file_free
  use tensor, only: tnsr3d
  use device_math, only: device_copy
  use gather_scatter
  use neko_config
  use logger, only: neko_log
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE, device_memcpy
  use comm, only: pe_rank
  use utils, only: NEKO_FNAME_LEN, neko_error
  use field_writer, only: field_writer_t
  use simulation_component, only: simulation_component_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use case, only: case_t
  use field_registry, only: neko_field_registry

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
     type(field_t), pointer :: u  => null()
     type(field_t), pointer :: v  => null()
     type(field_t), pointer :: w  => null()
     !> Transformed fields
     type(field_t), pointer :: u_hat => null()
     type(field_t), pointer :: v_hat => null()
     type(field_t), pointer :: w_hat => null()
     !> Working field - Consider making this a simple array
     type(field_t) :: wk
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
     !> spectral error indicator per element
     real(kind=rp), allocatable :: eind_u(:), eind_v(:), eind_w(:)
     !> fit coeficients per element
     real(kind=rp), allocatable :: sig_u(:), sig_v(:), sig_w(:)
     !> List to write the spectral error indicator as a field
     type(field_list_t) :: speri_l
     !> File to write
     type(file_t) :: mf_speri
     !> Field writer controller for the output
     type(field_writer_t) :: writer
   contains
     !> Constructor.
     procedure, pass(this) :: init => spectral_error_init
     !> Destructor.
     procedure, pass(this) :: free => spectral_error_free
     !> Compute the indicator (called according to the simcomp controller).
     procedure, pass(this) :: compute_ => spectral_error_compute
     !> Calculate the indicator.
     procedure, pass(this) :: get_indicators => spectral_error_get_indicators

  end type spectral_error_t

contains

  !> Constructor.
  subroutine spectral_error_init(this, json, case)
    class(spectral_error_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    character(len=20) :: fields(3)

    real(kind=rp) :: fluid_output_value, val
    character(len=:), allocatable :: fluid_output_control, str

    !> Add keyword "fields" to the json so that the field writer
    ! picks it up. Will also add those fields to the registry.
    fields(1) = "u_hat"
    fields(2) = "v_hat"
    fields(3) = "w_hat"
    call json%add("fields", fields)

    call json_get(case%params, "case.fluid.output_control", &
         fluid_output_control)
    call json_get(case%params, "case.fluid.output_value", &
         fluid_output_value)

    ! See if the user has set a compute control, otherwise
    ! set it to the fluid output control.
    call json_get_or_default(json, "compute_control", str, fluid_output_control)
    call json_get_or_default(json, "compute_value", val, fluid_output_value)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call spectral_error_init_from_attributes(this, case%fluid%c_Xh)

  end subroutine spectral_error_init

  !> Actual constructor.
  !! @param coef type with all geometrical variables.
  subroutine spectral_error_init_from_attributes(this, coef)
    class(spectral_error_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: il, jl, aa
    character(len=NEKO_FNAME_LEN) :: fname_speri

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    this%u_hat => neko_field_registry%get_field("u_hat")
    this%v_hat => neko_field_registry%get_field("v_hat")
    this%w_hat => neko_field_registry%get_field("w_hat")

    !> Initialize fields and copy data from proper one
    this%wk = this%u

    !> Allocate arrays (Consider moving some to coef)
    allocate(this%eind_u(coef%msh%nelv))
    allocate(this%eind_v(coef%msh%nelv))
    allocate(this%eind_w(coef%msh%nelv))
    allocate(this%sig_u(coef%msh%nelv))
    allocate(this%sig_v(coef%msh%nelv))
    allocate(this%sig_w(coef%msh%nelv))

    !> The following code has been lifted from Adam's implementation
    associate(LX1 => coef%Xh%lx, LY1 => coef%Xh%ly, &
      LZ1 => coef%Xh%lz, &
      SERI_SMALL  => this%SERI_SMALL,  &
      SERI_SMALLR => this%SERI_SMALLR, &
      SERI_SMALLG => this%SERI_SMALLG, &
      SERI_SMALLS => this%SERI_SMALLS, &
      SERI_NP     => this%SERI_NP,     &
      SERI_NP_MAX => this%SERI_NP_MAX, &
      SERI_ELR    => this%SERI_ELR     &
      )
      ! correctness check
      if (SERI_NP.gt.SERI_NP_MAX) then
         call neko_log%message('SETI_NP greater than SERI_NP_MAX')
      endif
      il = SERI_NP+SERI_ELR
      jl = min(LX1,LY1)
      jl = min(jl,LZ1)
      if (il.gt.jl) then
         call neko_log%message('SERI_NP+SERI_ELR greater than L?1')
      endif
    end associate

  end subroutine spectral_error_init_from_attributes

  !> Destructor
  subroutine spectral_error_free(this)
    class(spectral_error_t), intent(inout) :: this

    if(allocated(this%eind_u)) then
       deallocate(this%eind_u)
    end if

    if(allocated(this%eind_v)) then
       deallocate(this%eind_v)
    end if

    if(allocated(this%eind_w)) then
       deallocate(this%eind_w)
    end if

    if(allocated(this%sig_u)) then
       deallocate(this%sig_u)
    end if

    if(allocated(this%sig_w)) then
       deallocate(this%sig_w)
    end if

    if(allocated(this%sig_w)) then
       deallocate(this%sig_w)
    end if

    call this%wk%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%u_hat)
    nullify(this%v_hat)
    nullify(this%w_hat)

    call this%writer%free()
    call this%free_base()

  end subroutine spectral_error_free

  !> Compute the spectral error indicator.
  subroutine spectral_error_compute(this, t, tstep)
    class(spectral_error_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    integer :: e, i, lx, ly, lz, nelv, n

    print *, "COMPUTEEEEEE"

    call this%get_indicators(this%case%fluid%c_Xh)

    lx = this%u_hat%Xh%lx
    ly = this%u_hat%Xh%ly
    lz = this%u_hat%Xh%lz
    nelv = this%u_hat%msh%nelv

    !> Copy the element indicator into all points of the field
    do e = 1,nelv
       do i = 1,lx*ly*lx
          this%u_hat%x(i,1,1,e) = this%eind_u(e)
          this%v_hat%x(i,1,1,e) = this%eind_v(e)
          this%w_hat%x(i,1,1,e) = this%eind_w(e)
       end do
    end do

    ! We need this copy to GPU since the field writer does an internal copy
    ! GPU -> CPU before writing the field files. This overwrites the *_hat
    ! host arrays that contain the values.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%u_hat%x, this%u_hat%x_d, lx*ly*lz*nelv, &
                         HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%v_hat%x, this%v_hat%x_d, lx*ly*lz*nelv, &
                         HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%w_hat%x, this%w_hat%x_d, lx*ly*lz*nelv, &
                         HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine spectral_error_compute

  !> Transform a field u to u_hat into physical or spectral space
  !! the result of the transformation is in u_hat.
  !! @param u_hat Transformed field (output).
  !! @param u Field to transform (input).
  !! @param wk Working field.
  !! @param coef Type coef for mesh parameters.
  !! @param space String that indicates which space to transform, "spec" or "phys".
  subroutine transform_to_spec_or_phys(u_hat, u, wk, coef, space)
    type(field_t), intent(inout) :: u_hat
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: wk
    type(coef_t), intent(inout) :: coef
    character(len=4), intent(in) :: space
    integer :: i, j, k, e, nxyz, nelv, n

    !> Define some constants
    nxyz = coef%Xh%lx*coef%Xh%lx*coef%Xh%lx
    nelv = coef%msh%nelv
    n    = nxyz*nelv

    !> Copy field to working array
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
      (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_copy(wk%x_d, u%x_d, n)
    else
       call copy(wk%x,u%x,n)
    end if

    select case(space)
    case('spec')
       call tnsr3d(u_hat%x, coef%Xh%lx, wk%x, &
                    coef%Xh%lx,coef%Xh%vinv, &
                    coef%Xh%vinvt, coef%Xh%vinvt, nelv)
    case('phys')
       call tnsr3d(u_hat%x, coef%Xh%lx, wk%x, &
                    coef%Xh%lx,coef%Xh%v, &
                    coef%Xh%vt, coef%Xh%vt, nelv)
    end select

    ! Synchronize
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then

       call device_memcpy(u_hat%x,u_hat%x_d, n, &
                         DEVICE_TO_HOST, sync=.true.)
    end if

  end subroutine transform_to_spec_or_phys

  !> Transform and get the spectral error indicators
  !! @param coef type coef for mesh parameters and space
  subroutine spectral_error_get_indicators(this,coef)
    class(spectral_error_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: i

    ! Generate the uvwhat field (legendre coeff)
    call transform_to_spec_or_phys(this%u_hat, this%u, this%wk, coef, 'spec')
    call transform_to_spec_or_phys(this%v_hat, this%v, this%wk, coef, 'spec')
    call transform_to_spec_or_phys(this%w_hat, this%w, this%wk, coef, 'spec')

    ! Get the spectral error indicator
    call calculate_indicators(this, coef, this%eind_u, this%sig_u, &
         coef%msh%nelv, coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
         this%u_hat%x)
    call calculate_indicators(this, coef, this%eind_v, this%sig_v, &
         coef%msh%nelv, coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
         this%v_hat%x)
    call calculate_indicators(this, coef, this%eind_w, this%sig_w, &
         coef%msh%nelv, coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
         this%w_hat%x)

  end subroutine spectral_error_get_indicators

  !> Write error indicators in a field file.
  !! @param t Current simulation time.
  subroutine spectral_error_write(this, t)
    class(spectral_error_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    integer i, e
    integer lx, ly, lz, nelv

    !> Write the file
    !! Remember that the list is already ponting to the fields
    !! that were just modified.
    call this%mf_speri%write(this%speri_l,t)

  end subroutine spectral_error_write

  !> Wrapper for old fortran 77 subroutines
  !! @param coef coef type
  !! @param eind spectral indicator
  !! @param sig coefficient of the exponential fit
  !! @param lnelt number of elements
  !! @param LX1 gll points in x
  !! @param LY1 gll points in y
  !! @param LZ1 gll points in z
  !! @paran var variable to calculate indicator
  subroutine calculate_indicators(this, coef, eind, sig, lnelt, LX1, LY1, LZ1, var)
    type(spectral_error_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: lnelt
    integer :: LX1
    integer :: LY1
    integer :: LZ1
    real(kind=rp) :: eind(lnelt)
    real(kind=rp) :: sig(lnelt)
    real(kind=rp) :: var(LX1,LY1,LZ1,lnelt)

    real(kind=rp) :: xa(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz)
    real(kind=rp) :: xb(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz)
    integer :: i, e

    ! zero arrays
    call rzero(eind,lnelt)
    call rzero(sig,lnelt)

    ! Get the indicator
    call speri_var(this, eind,sig,var,lnelt,xa,xb, LX1, LY1, LZ1)

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
  subroutine speri_var(this, est,sig,var,nell,xa,xb,LX1,LY1,LZ1)
    type(spectral_error_t), intent(inout) :: this
    integer :: nell
    integer :: LX1
    integer :: LY1
    integer :: LZ1
    real(kind=rp) :: est(nell)
    real(kind=rp) :: sig(nell)
    real(kind=rp) :: var(LX1,LY1,LZ1,nell)
    real(kind=rp) :: xa(LX1,LY1,LZ1)
    real(kind=rp) :: xb(LX1,LY1,LZ1)

    !> local variables
    integer :: il, jl, kl, ll, j_st, j_en, ii
    !> polynomial coefficients
    real(kind=rp) :: coeff(LX1,LY1,LZ1)
    !> Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) ::  coef11
    !> copy of last SERI_NP columns of coefficients
    real(kind=rp) ::  coefx(this%SERI_NP_MAX,LY1,LZ1), &
                      coefy(this%SERI_NP_MAX,LX1,LZ1), &
                      coefz(this%SERI_NP_MAX,LX1,LY1)
    !> estimated error
    real(kind=rp) ::  estx, esty, estz
    !> estimated decay rate
    real(kind=rp) ::  sigx, sigy, sigz
    real(kind=rp) ::  third
    parameter (third = 1.0/3.0)

    !> loop over elements
    do il = 1,nell
       !> go to Legendre space (done in two operations)
       !! and square the coefficient
       do ii = 1, LX1*LY1*LZ1
          coeff(ii,1,1) = var(ii,1,1,il) * var(ii,1,1,il)
       end do

       !> lower left corner
       coef11 = coeff(1,1,1)

       !> small value; nothing to od
       if (coef11.ge.this%SERI_SMALL) then
          !> extrapolate coefficients
          !> X - direction
          !> copy last SERI_NP collumns (or less if NX1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1,LX1-this%SERI_NP+1-this%SERI_ELR)
          j_en = max(1,LX1-this%SERI_ELR)
          do ll = 1,LZ1
             do kl = 1,LY1
                do jl = j_st,j_en
                   coefx(j_en-jl+1,kl,ll) = coeff(jl,kl,ll)
                enddo
             enddo
          enddo
          !> get extrapolated values
          call speri_extrap(this,estx,sigx,coef11,coefx, &
            j_st,j_en,LY1,LZ1)

          !> Y - direction
          !> copy last SERI_NP collumns (or less if NY1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1,LY1-this%SERI_NP+1-this%SERI_ELR)
          j_en = max(1,LY1-this%SERI_ELR)
          do ll = 1,LZ1
             do kl = j_st,j_en
                do jl = 1,LX1
                   coefy(j_en-kl+1,jl,ll) = coeff(jl,kl,ll)
                enddo
             enddo
          enddo
          !> get extrapolated values
          call speri_extrap(this, esty,sigy,coef11,coefy, &
             j_st,j_en,LX1,LZ1)

          !> Z - direction
          !> copy last SERI_NP collumns (or less if NZ1 is smaller)
          !> SERI_ELR allows to exclude last row
          j_st = max(1,LZ1-this%SERI_NP+1-this%SERI_ELR)
          j_en = max(1,LZ1-this%SERI_ELR)
          do ll = j_st,j_en
             do kl = 1,LY1
                do jl = 1,LX1
                   coefz(j_en-ll+1,jl,kl) = coeff(jl,kl,ll)
                enddo
             enddo
          enddo
          !> get extrapolated values
          call speri_extrap(this, estz,sigz,coef11,coefz, &
        j_st,j_en,LX1,LY1)

          !> average
          est(il) =  sqrt(estx + esty + estz)
          sig(il) =  third*(sigx + sigy + sigz)

       else
          !> for testing
          estx = 0.0
          esty = 0.0
          estz = 0.0
          sigx = -1.0
          sigy = -1.0
          sigz = -1.0
          !> for testing; end

          est(il) =  0.0
          sig(il) = -1.0
       endif

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
  subroutine speri_extrap(this,estx,sigx,coef11,coef, &
                ix_st,ix_en,nyl,nzl)
    implicit none
    type(spectral_error_t), intent(inout) :: this
    !> argument list
    integer :: ix_st,ix_en,nyl,nzl
    !> Legendre coefficients; last SERI_NP columns
    real(kind=rp) :: coef(this%SERI_NP_MAX,nyl,nzl)
    !> Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) :: coef11
    !> estimated error and decay rate
    real(kind=rp) :: estx, sigx

    !> local variables
    integer :: il, jl, kl, ll  ! loop index
    integer :: nsigt, pnr, nzlt
    real(kind=rp) :: sigt, smallr, cmin, cmax, cnm, rtmp, rtmp2, rtmp3
    real(kind=rp) :: sumtmp(4), cffl(this%SERI_NP_MAX)
    real(kind=rp) :: stmp, estt, clog, ctmp, cave, erlog
    logical :: cuse(this%SERI_NP_MAX)

    associate(SERI_SMALL  => this%SERI_SMALL,  &
             SERI_SMALLR => this%SERI_SMALLR, &
             SERI_SMALLG => this%SERI_SMALLG, &
             SERI_SMALLS => this%SERI_SMALLS, &
             SERI_NP     => this%SERI_NP,     &
             SERI_NP_MAX => this%SERI_NP_MAX, &
             SERI_ELR    => this%SERI_ELR     &
            )
      ! initial values
      estx =  0.0
      sigx = -1.0

      ! relative cutoff
      smallr = coef11*SERI_SMALLR

      ! number of points
      pnr = ix_en - ix_st +1

      ! to few points to interpolate
      !if ((ix_en - ix_st).le.1) return

      ! for averaging, initial values
      sigt = 0.0
      nsigt = 0

      ! loop over all face points
      nzlt = max(1,nzl - SERI_ELR) !  for 2D runs
      do il=1,nzlt
         ! weight
         rtmp3 = 1.0/(2.0*(il-1)+1.0)
         do jl=1,nyl - SERI_ELR

            ! find min and max coef along single row
            cffl(1) = coef(1,jl,il)
            cmin = cffl(1)
            cmax = cmin
            do kl =2,pnr
               cffl(kl) = coef(kl,jl,il)
               cmin = min(cmin,cffl(kl))
               cmax = max(cmax,cffl(kl))
            enddo

            ! are coefficients sufficiently big
            if((cmin.gt.0.0).and.(cmax.gt.smallr)) then
               ! mark array position we use in iterpolation
               do kl =1,pnr
                  cuse(kl) = .TRUE.
               enddo
               ! max n for polynomial order
               cnm = real(ix_en)

               ! check if all the points should be taken into account
               ! in original code by Catherine Mavriplis this part is written
               ! for 4 points, so I place if statement first
               if (pnr.eq.4) then
                  ! should we neglect last values
                  if ((cffl(1).lt.smallr).and. &
                 (cffl(2).lt.smallr)) then
                     if (cffl(3).lt.smallr) then
                        cuse(1) = .FALSE.
                        cuse(2) = .FALSE.
                        cnm = real(ix_en-2)
                     else
                        cuse(1) = .FALSE.
                        cnm = real(ix_en-1)
                     endif
                  else
                     ! should we take stronger gradient
                     if ((cffl(1)/cffl(2).lt.SERI_SMALLG).and. &
                   (cffl(3)/cffl(4).lt.SERI_SMALLG)) then
                        cuse(1) = .FALSE.
                        cuse(3) = .FALSE.
                        cnm = real(ix_en-1)
                     elseif ((cffl(2)/cffl(1).lt.SERI_SMALLG).and. &
                       (cffl(4)/cffl(3).lt.SERI_SMALLG)) then
                        cuse(2) = .FALSE.
                        cuse(4) = .FALSE.
                     endif
                  endif
               endif

               ! get sigma for given face point
               do kl =1,4
                  sumtmp(kl) = 0.0
               enddo
               ! find new min and count number of points
               cmin = cmax
               cmax = 0.0
               do kl =1,pnr
                  if(cuse(kl)) then
                     rtmp  = real(ix_en-kl)
                     rtmp2 = log(cffl(kl))
                     sumtmp(1) = sumtmp(1) +rtmp2
                     sumtmp(2) = sumtmp(2) +rtmp
                     sumtmp(3) = sumtmp(3) +rtmp*rtmp
                     sumtmp(4) = sumtmp(4) +rtmp2*rtmp
                     ! find new min and count used points
                     cmin = min(cffl(kl),cmin)
                     cmax = cmax + 1.0
                  endif
               enddo
               ! decay rate along single row
               stmp = (sumtmp(1)*sumtmp(2) - sumtmp(4)*cmax)/ &
                   (sumtmp(3)*cmax - sumtmp(2)*sumtmp(2))
               ! for averaging
               sigt = sigt + stmp
               nsigt = nsigt + 1

               ! get error estimator depending on calculated decay rate
               estt = 0.0
               if (stmp.lt.SERI_SMALLS) then
                  estt = cmin
               else
                  ! get averaged constant in front of c*exp(-sig*n)
                  clog = (sumtmp(1)+stmp*sumtmp(2))/cmax
                  ctmp = exp(clog)
                  ! average exponent
                  cave = sumtmp(1)/cmax
                  ! check quality of approximation comparing is to the constant cave
                  do kl =1,2
                     sumtmp(kl) = 0.0
                  enddo
                  do kl =1,pnr
                     if(cuse(kl)) then
                        erlog = clog - stmp*real(ix_en-kl)
                        sumtmp(1) = sumtmp(1)+ &
                            (erlog-log(cffl(kl)))**2
                        sumtmp(2) = sumtmp(2)+ &
                  (erlog-cave)**2
                     endif
                  enddo
                  rtmp = 1.0 - sumtmp(1)/sumtmp(2)
                  if (rtmp.lt.SERI_SMALLS) then
                     estt = cmin
                  else
                     ! last coefficient is not included in error estimator
                     estt = ctmp/stmp*exp(-stmp*cnm)
                  endif
               endif
               ! add contribution to error estimator; variable weight
               estx = estx + estt/(2.0*(jl-1)+1.0)*rtmp3
            endif  ! if((cmin.gt.0.0).and.(cmax.gt.smallr))
         enddo
      enddo
      ! constant weight
      ! Multiplication by 4 in 2D / 8 in 3D
      ! Normalization of the error by the volume of the reference element
      ! which is equal to 4 in 2D / 8 in 3D
      ! ==> Both operations cancel each other
      estx = estx/(2.0*(ix_en-1)+1.0)

      ! final everaging
      ! sigt = 2*sigma so we divide by 2
      if (nsigt.gt.0) then
         sigx = 0.5*sigt/nsigt
      endif

    end associate

  end subroutine speri_extrap

end module spectral_error
