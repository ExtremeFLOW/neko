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
!> spectral_error_indicator
module spectral_error_indicator
  use num_types, only: rp
  use logger, only: neko_log, LOG_SIZE
  use field, only: field_t
  use coefs, only: coef_t
  use field_list, only: field_list_t
  use math, only: rzero
  use file, only: file_t, file_free
  use tensor, only: tnsr3d
  use device_math, only: device_copy
  use gather_scatter
  use neko_config
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> include information needed for compressing fields
  type, public :: spectral_error_indicator_t
     !> Pointers to main fields 
     type(field_t), pointer :: u  => null()
     type(field_t), pointer :: v  => null()
     type(field_t), pointer :: w  => null()
     !> Transformed fields
     type(field_t) :: u_hat
     type(field_t) :: v_hat
     type(field_t) :: w_hat
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
     !> Device pointers

     contains
       !> Initialize object.
       procedure, pass(this) :: init => spec_err_ind_init
       !> Destructor
       procedure, pass(this) :: free => spec_err_ind_free
       !> Calculate the indicator
       procedure, pass(this) :: get_indicators => spec_err_ind_get
       !> Calculate the indicator
       procedure, pass(this) :: write_as_field => spec_err_ind_write

  end type spectral_error_indicator_t

contains

  !> Initialize the object
  !! @param u u velocity field
  !! @param v v velocity field
  !! @param w w velocity field
  !! @param coef type with all geometrical variables
  subroutine spec_err_ind_init(this, u,v,w,coef)
    class(spectral_error_indicator_t), intent(inout) :: this
    type(field_t), intent(in), target :: u
    type(field_t), intent(in), target :: v
    type(field_t), intent(in), target :: w
    type(coef_t), intent(in) :: coef
    integer :: il, jl, aa
    character(len=NEKO_FNAME_LEN) :: fname_speri

    !> Assign the pointers
    this%u => u
    this%v => v
    this%w => w
    !> Initialize fields and copy data from proper one
    this%u_hat = u
    this%v_hat = v
    this%w_hat = w
    this%wk = u
    !> Allocate arrays (Consider moving some to coef)
    allocate(this%eind_u(coef%msh%nelv))
    allocate(this%eind_v(coef%msh%nelv))
    allocate(this%eind_w(coef%msh%nelv))
    allocate(this%sig_u(coef%msh%nelv))
    allocate(this%sig_v(coef%msh%nelv))
    allocate(this%sig_w(coef%msh%nelv))

    !> The following code has been lifted from Adams implementation

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
         if (pe_rank.eq.0) write(*,*) 'SETI_NP greater than SERI_NP_MAX' 
         endif
         il = SERI_NP+SERI_ELR
         jl = min(LX1,LY1)
         jl = min(jl,LZ1)
         if (il.gt.jl) then
            if (pe_rank.eq.0) write(*,*) 'SERI_NP+SERI_ELR greater than L?1'
         endif
    end associate

    !> Initialize the list that holds the fields to write
    call list_init3(this%speri_l,this%u_hat,this%v_hat, &
                    this%w_hat)
    ! Initialize the file
    fname_speri = 'speri.fld'
    this%mf_speri =  file_t(fname_speri)

  end subroutine spec_err_ind_init


  !> Detructor
  subroutine spec_err_ind_free(this)
    class(spectral_error_indicator_t), intent(inout) :: this

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

    call this%u_hat%free()
    call this%v_hat%free()
    call this%w_hat%free()
    call this%wk%free()


    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
   
    !> finalize data related to writing 
    call list_final3(this%speri_l)
    call file_free(this%mf_speri)
                
    ! Cleanup the device (if present)
                

  end subroutine spec_err_ind_free

  !> Transform a field u > u_hat into physical or spectral space
  ! Currentlly I need to pass the speri type but it should be in coed
  !the result of the transform is given in fldhat
  subroutine transform_to_spec_or_phys(u_hat, u, wk, coef, space)
    type(field_t), intent(inout) :: u_hat
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: wk
    type(coef_t), intent(inout) :: coef
    character(len=4), intent(in) :: space             
    integer :: i, j, k, e, nxyz, nelv, n
    character(len=LOG_SIZE) :: log_buf 

    ! Define some constants
    nxyz = coef%Xh%lx*coef%Xh%lx*coef%Xh%lx
    nelv = coef%msh%nelv
    n    = nxyz*nelv

    ! Copy field to working array
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

  !> Transform and get the indicators
  subroutine spec_err_ind_get(this,coef)
    class(spectral_error_indicator_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: i

    ! Generate the uvwhat field (legendre coeff)
    call transform_to_spec_or_phys(this%u_hat, this%u, this%wk, coef, 'spec')
    call transform_to_spec_or_phys(this%v_hat, this%v, this%wk, coef, 'spec')
    call transform_to_spec_or_phys(this%w_hat, this%w, this%wk, coef, 'spec')
   
    
    ! Get the spectral error indicator
    call calculate_indicators(this, coef, this%eind_u, this%sig_u, coef%msh%nelv, &
                              coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
                              this%u_hat%x)
    call calculate_indicators(this, coef, this%eind_v, this%sig_v, coef%msh%nelv, &
                              coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
                              this%v_hat%x)
    call calculate_indicators(this, coef, this%eind_w, this%sig_w, coef%msh%nelv, &
                              coef%Xh%lx,  coef%Xh%ly,  coef%Xh%lz, &
                              this%w_hat%x)
                             
  end subroutine spec_err_ind_get

  !> Initialize user defined variables.
  !! @param t Current simulation time.
  subroutine spec_err_ind_write(this, t)
    class(spectral_error_indicator_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    
    integer i, e
    integer lx, ly, lz, nelv
    
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

    !> Write the file
    !! Remember that the list is already ponting to the fields
    !! that were just modified.
    call this%mf_speri%write(this%speri_l,t)      

  end subroutine spec_err_ind_write

  !> Wrapper for old fortran 77 subroutines
  subroutine calculate_indicators(this, coef, eind, sig, lnelt, LX1, LY1, LZ1, var)
    type(spectral_error_indicator_t), intent(inout) :: this
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

  subroutine speri_var(this, est,sig,var,nell,xa,xb,LX1,LY1,LZ1)
    type(spectral_error_indicator_t), intent(inout) :: this
    integer :: nell
    integer :: LX1
    integer :: LY1
    integer :: LZ1
    real(kind=rp) :: est(nell)
    real(kind=rp) :: sig(nell)
    real(kind=rp) :: var(LX1,LY1,LZ1,nell)
    real(kind=rp) :: xa(LX1,LY1,LZ1)
    real(kind=rp) :: xb(LX1,LY1,LZ1)
    
    ! local variables
    integer :: il, jl, kl, ll, j_st, j_en, ii
    ! polynomial coefficients
    real(kind=rp) :: coeff(LX1,LY1,LZ1)
    ! Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) ::  coef11
    ! copy of last SERI_NP columns of coefficients
    real(kind=rp) ::  coefx(this%SERI_NP_MAX,LY1,LZ1), & 
                      coefy(this%SERI_NP_MAX,LX1,LZ1), &
                      coefz(this%SERI_NP_MAX,LX1,LY1)
    ! estimated error
    real(kind=rp) ::  estx, esty, estz
    ! estimated decay rate
    real(kind=rp) ::  sigx, sigy, sigz
    real(kind=rp) ::  third
    parameter (third = 1.0/3.0)

    ! loop over elements
    do il = 1,nell
        ! go to Legendre space (done in two operations)
        ! and square the coefficient
        do ii = 1, LX1*LY1*LZ1
           coeff(ii,1,1) = var(ii,1,1,il) * var(ii,1,1,il) 
        end do

        ! lower left corner
        coef11 = coeff(1,1,1)

        ! small value; nothing to od
        if (coef11.ge.this%SERI_SMALL) then
           ! extrapolate coefficients
           ! X - direction
           ! copy last SERI_NP collumns (or less if NX1 is smaller)
           ! SERI_ELR allows to exclude last row
            j_st = max(1,LX1-this%SERI_NP+1-this%SERI_ELR)
            j_en = max(1,LX1-this%SERI_ELR)
            do ll = 1,LZ1
                do kl = 1,LY1
                    do jl = j_st,j_en
                        coefx(j_en-jl+1,kl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo
            ! get extrapolated values
            call speri_extrap(this,estx,sigx,coef11,coefx, &
                 j_st,j_en,LY1,LZ1)
         
            ! Y - direction
            ! copy last SERI_NP collumns (or less if NY1 is smaller)
            ! SERI_ELR allows to exclude last row
            j_st = max(1,LY1-this%SERI_NP+1-this%SERI_ELR)
            j_en = max(1,LY1-this%SERI_ELR)
            do ll = 1,LZ1
                do kl = j_st,j_en
                    do jl = 1,LX1
                        coefy(j_en-kl+1,jl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

            ! get extrapolated values
            call speri_extrap(this, esty,sigy,coef11,coefy, &
                j_st,j_en,LX1,LZ1)
   
            ! Z - direction
            ! copy last SERI_NP collumns (or less if NZ1 is smaller)
            ! SERI_ELR allows to exclude last row
            j_st = max(1,LZ1-this%SERI_NP+1-this%SERI_ELR)
            j_en = max(1,LZ1-this%SERI_ELR)
            do ll = j_st,j_en
                do kl = 1,LY1
                    do jl = 1,LX1
                        coefz(j_en-ll+1,jl,kl) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

            ! get extrapolated values
            call speri_extrap(this, estz,sigz,coef11,coefz, &
               j_st,j_en,LX1,LY1)

            ! average
            est(il) =  sqrt(estx + esty + estz)
            sig(il) =  third*(sigx + sigy + sigz)

        else
            ! for testing
            estx = 0.0
            esty = 0.0
            estz = 0.0
            sigx = -1.0
            sigy = -1.0
            sigz = -1.0
            ! for testing; end

            est(il) =  0.0
            sig(il) = -1.0
        endif

    end do

  end subroutine speri_var


  subroutine speri_extrap(this,estx,sigx,coef11,coef, &
                ix_st,ix_en,nyl,nzl)
    implicit none
    type(spectral_error_indicator_t), intent(inout) :: this
    ! argument list
    integer :: ix_st,ix_en,nyl,nzl
    ! Legendre coefficients; last SERI_NP columns
    real(kind=rp) :: coef(this%SERI_NP_MAX,nyl,nzl)
    ! Legendre coefficients; first value coeff(1,1,1)
    real(kind=rp) :: coef11
    ! estimated error and decay rate
    real(kind=rp) :: estx, sigx

    ! local variables
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


  subroutine list_init3(list,uu,vv,ww)
    type(field_list_t), intent(inout) :: list
    type(field_t) , target:: uu
    type(field_t) , target:: vv
    type(field_t) , target:: ww
    !> Initialize field lists
    allocate(list%fields(3))
    list%fields(1)%f => uu
    list%fields(2)%f => vv
    list%fields(3)%f => ww
  end subroutine list_init3


  subroutine list_final3(list)
    type(field_list_t), intent(inout) :: list
    !> Deallocate field lists
    deallocate(list%fields)
  end subroutine list_final3


end module spectral_error_indicator



