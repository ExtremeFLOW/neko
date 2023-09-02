module spectral_error_indicator

  use neko
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> include information needed for compressing fields
  type, public :: spec_err_ind_t
     !> Transformation matrices
     real(kind=rp), allocatable :: v(:,:) 
     real(kind=rp), allocatable :: vt(:,:) 
     real(kind=rp), allocatable :: vinv(:,:)  
     real(kind=rp), allocatable :: vinvt(:,:)  
     !> Legendre weights in matrix form
     real(kind=rp), allocatable :: w(:,:)
     !> Place holder matrices for easy transformation 
     real(kind=rp), allocatable :: specmat(:,:) !< Transformation matrix
     real(kind=rp), allocatable :: specmatt(:,:) !< Transformation matrix
     !> Pointers to main fields 
     type(field_t), pointer :: u_  => null()
     type(field_t), pointer :: v_  => null()
     type(field_t), pointer :: w_  => null()
     !> Transformed fields
     type(field_t) :: u_hat
     type(field_t) :: v_hat
     type(field_t) :: w_hat
     !> Working field - Consider making this a simple array
     type(field_t) :: wk
     !> Configuration of spectral error calculation
     real(kind=rp) :: SERI_SMALL
     !!> used for ratios
     real(kind=rp) :: SERI_SMALLR
     !!> used for gradients
     real(kind=rp) :: SERI_SMALLG
     !!> used for sigma and rtmp in error calculations
     real(kind=rp) :: SERI_SMALLS
     !!> number of points in fitting
     integer :: SERI_NP
     integer :: SERI_NP_MAX
     !!> last modes skipped
     integer :: SERI_ELR
     !!> spectral error indicator per element
     real(kind=rp), allocatable :: eind_u(:), eind_v(:), eind_w(:)
     !!> fit coeficients per element
     real(kind=rp), allocatable :: sig_u(:), sig_v(:), sig_w(:)
     !> Device pointers
     type(c_ptr) :: v_d = C_NULL_PTR
     type(c_ptr) :: vt_d = C_NULL_PTR
     type(c_ptr) :: vinv_d = C_NULL_PTR
     type(c_ptr) :: vinvt_d = C_NULL_PTR
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: specmat_d = C_NULL_PTR
     type(c_ptr) :: specmatt_d = C_NULL_PTR

     contains
       !> Initialize object.
       procedure, pass(this) :: init => spec_err_ind_init
       !> Destructor
       procedure, pass(this) :: free => spec_err_ind_free
       !> Calculate the indicator
       procedure, pass(this) :: get_indicators => spec_err_ind_get

  end type spec_err_ind_t

contains

  !> Initialize the object
  !! @param u_ u velocity field
  !! @param v_ v velocity field
  !! @param w_ w velocity field
  !! @param coef type with all geometrical variables
  subroutine spec_err_ind_init(this, u_,v_,w_,coef)
    class(spec_err_ind_t), intent(inout) :: this
    type(field_t), intent(in), target :: u_
    type(field_t), intent(in), target :: v_
    type(field_t), intent(in), target :: w_
    type(coef_t), intent(in) :: coef
    integer :: il, jl, aa

    !> Assign the pointers
    this%u_ => u_
    this%v_ => v_
    this%w_ => w_
    !> Initialize fields and copy data from proper one
    this%u_hat = u_
    this%v_hat = v_
    this%w_hat = w_
    this%wk = u_
    !> Allocate arrays (Consider moving some to coef)
    allocate(this%v(coef%Xh%lx, coef%Xh%lx))
    allocate(this%vt(coef%Xh%lx, coef%Xh%lx))
    allocate(this%vinv(coef%Xh%lx, coef%Xh%lx))
    allocate(this%vinvt(coef%Xh%lx, coef%Xh%lx))
    allocate(this%w(coef%Xh%lx, coef%Xh%lx))
    allocate(this%specmat(coef%Xh%lx, coef%Xh%lx))
    allocate(this%specmatt(coef%Xh%lx, coef%Xh%lx))
    allocate(this%eind_u(coef%msh%nelv))
    allocate(this%eind_v(coef%msh%nelv))
    allocate(this%eind_w(coef%msh%nelv))
    allocate(this%sig_u(coef%msh%nelv))
    allocate(this%sig_v(coef%msh%nelv))
    allocate(this%sig_w(coef%msh%nelv))

    ! Initialize all the transformation matrices
    call generate_transformation_matrices(this,coef)

    !> The following code has been lifted from Adams implementation
    ! set cutoff parameters
    ! used for values
    this%SERI_SMALL = 1.e-14
    ! used for ratios
    this%SERI_SMALLR = 1.e-10
    ! used for gradients
    this%SERI_SMALLG = 1.e-5
    ! used for sigma and rtmp in error calculations
    this%SERI_SMALLS = 0.2
    ! number of points in fitting
    this%SERI_NP = 4
    this%SERI_NP_MAX = 4
    ! last modes skipped
    this%SERI_ELR = 0

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

  end subroutine spec_err_ind_init


  !> Detructor
  subroutine spec_err_ind_free(this)
    class(spec_err_ind_t), intent(inout) :: this

    if(allocated(this%v)) then
       deallocate(this%v)
    end if

    if(allocated(this%vt)) then
       deallocate(this%vt)
    end if

    if(allocated(this%vinv)) then
       deallocate(this%vinv)
    end if

    if(allocated(this%vinvt)) then
       deallocate(this%vinvt)
    end if

    if(allocated(this%w)) then
       deallocate(this%w)
    end if

    if(allocated(this%specmat)) then
       deallocate(this%specmat)
    end if

    if(allocated(this%specmatt)) then
       deallocate(this%specmatt)
    end if

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

    call field_free(this%u_hat)
    call field_free(this%v_hat)
    call field_free(this%w_hat)
    call field_free(this%wk)

    nullify(this%u_)
    nullify(this%v_)
    nullify(this%w_)
                
    ! Cleanup the device (if present)
                
    if (c_associated(this%v_d)) then
       call device_free(this%v_d)
    end if

    if (c_associated(this%vt_d)) then
       call device_free(this%vt_d)
    end if

    if (c_associated(this%vinv_d)) then
       call device_free(this%vinv_d)
    end if

    if (c_associated(this%vinvt_d)) then
       call device_free(this%vinvt_d)
    end if

    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if

    if (c_associated(this%specmat_d)) then
       call device_free(this%specmat_d)
    end if

    if (c_associated(this%specmatt_d)) then
       call device_free(this%specmatt_d)
    end if

  end subroutine spec_err_ind_free


  !> Generate spectral tranform matrices
  !! @param coef type with all geometrical variables
  subroutine generate_transformation_matrices(this,coef)
    type(spec_err_ind_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    
    real(kind=rp) :: L(0:coef%Xh%lx-1)
    real(kind=rp) :: delta(coef%Xh%lx)
    integer :: i, kj, j, j2, kk
    character(len=LOG_SIZE) :: log_buf 

    associate(Xh => coef%Xh, v=> this%v, vt => this%vt, &
      vinv => this%vinv, vinvt => this%vinvt, w => this%w)
      ! Get the Legendre polynomials for each point
      ! Then proceed to compose the transform matrix
      kj = 0
      do j = 1, Xh%lx
         L(0) = 1.
         L(1) = Xh%zg(j,1)
         do j2 = 2, Xh%lx-1
            L(j2) = ( (2*j2-1) * Xh%zg(j,1) * L(j2-1) &
                  - (j2-1) * L(j2-2) ) / j2 
         end do
         do kk = 1, Xh%lx
            kj = kj+1
            v(kj,1) = L(KK-1)
         end do
      end do

      ! transpose the matrix
      call trsp1(v, Xh%lx) !< non orthogonal wrt weights

      ! Calculate the nominal scaling factors
      do i = 1, Xh%lx
         delta(i) = 2.0_rp / (2*(i-1)+1)
      end do
      ! modify last entry  
      delta(Xh%lx) = 2.0_rp / (Xh%lx-1)

      ! calculate the inverse to multiply the matrix
      do i = 1, Xh%lx
         delta(i) = sqrt(1.0_rp / delta(i))
      end do
      ! scale the matrix      
      do i = 1, Xh%lx
         do j = 1, Xh%lx
            v(i,j) = v(i,j) * delta(j) ! orthogonal wrt weights
         end do
      end do

      ! get the trasposed
      call copy(vt, v, Xh%lx * Xh%lx)
      call trsp1(vt, Xh%lx)

      !populate the mass matrix
      kk = 1
      do i = 1, Xh%lx
         do j = 1, Xh%lx
            if (i .eq. j) then
               w(i,j) = Xh%wx(kk)
               kk = kk+1
            else
               w(i,j) = 0
            end if
         end do
      end do

      !Get the inverse of the transform matrix
      call mxm(vt, Xh%lx, w, Xh%lx, vinv, Xh%lx)

      !get the transposed of the inverse
      call copy(vinvt, vinv, Xh%lx * Xh%lx)
      call trsp1(vinvt, Xh%lx)
    end associate

    ! Copy the data to the GPU
    ! Move all this to space.f90 to for next version 
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
    (NEKO_BCKND_OPENCL .eq. 1)) then 
       call device_map(this%v,     this%v_d,     coef%Xh%lxy)
       call device_map(this%vt,    this%vt_d,    coef%Xh%lxy)
       call device_map(this%vinv,  this%vinv_d,  coef%Xh%lxy)
       call device_map(this%vinvt, this%vinvt_d, coef%Xh%lxy)
       call device_map(this%w,     this%w_d,     coef%Xh%lxy)
       !Map the following pointers but do not copy data for them
       call device_map(this%specmat,  this%specmat_d,  coef%Xh%lxy)
       call device_map(this%specmatt, this%specmatt_d, coef%Xh%lxy)

       call device_memcpy(this%v,     this%v_d,     coef%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(this%vt,    this%vt_d,    coef%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(this%vinv,  this%vinv_d,  coef%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(this%vinvt, this%vinvt_d, coef%Xh%lxy, &
                          HOST_TO_DEVICE)
       call device_memcpy(this%w,     this%w_d,     coef%Xh%lxy, &
                          HOST_TO_DEVICE)

    end if

  end subroutine generate_transformation_matrices


  !> Transform a field u > u_hat into physical or spectral space
  ! Currentlly I need to pass the speri type but it should be in coed
  !the result of the transform is given in fldhat
  subroutine transform_to_spec_or_phys(u_hat, u, wk, coef, sei, space)
    type(field_t), intent(inout) :: u_hat
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: wk
    type(spec_err_ind_t), intent(inout) :: sei
    type(coef_t), intent(in) :: coef
                
    integer :: i, j, k, e, nxyz, nelv, n
    character(len=LOG_SIZE) :: log_buf 
    character(len=4) :: space 

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
                      coef%Xh%lx,sei%vinv, &
                      sei%vinvt, sei%vinvt, nelv)
       case('phys') 
          call tnsr3d(u_hat%x, coef%Xh%lx, wk%x, &
                      coef%Xh%lx,sei%v, &
                      sei%vt, sei%vt, nelv)
    end select

    ! Synchronize
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then 

       call device_memcpy(u_hat%x,u_hat%x_d, n, &
                          DEVICE_TO_HOST)
    end if

  end subroutine transform_to_spec_or_phys

  !> Trnasform and get the indicators
  subroutine spec_err_ind_get(this,coef)
    class(spec_err_ind_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    
    ! Generate the uvwhat field (legendre coeff)
    call transform_to_spec_or_phys(this%u_hat, this%u_, this%wk, coef, this, 'spec')
    call transform_to_spec_or_phys(this%v_hat, this%v_, this%wk, coef, this, 'spec')
    call transform_to_spec_or_phys(this%w_hat, this%w_, this%wk, coef, this, 'spec')
    
  
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


  !> Wrapper for old fortran 77 subroutines
  subroutine calculate_indicators(this, coef, eind, sig, lnelt, LX1, LY1, LZ1, var)
    type(spec_err_ind_t), intent(inout) :: this
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

    ! Put the indicator in the fieldhat to write it
    do e = 1,lnelt
       do i = 1,LX1*LY1*LZ1
          var(i,1,1,e) = eind(e)
       end do
    end do

  end subroutine calculate_indicators 

  subroutine speri_var(this, est,sig,var,nell,xa,xb,LX1,LY1,LZ1)
    type(spec_err_ind_t), intent(inout) :: this
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
      type(spec_err_ind_t), intent(inout) :: this
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
!-----------------------------------------------------------------------
      
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
!      if ((ix_en - ix_st).le.1) return

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


end module spectral_error_indicator

module user
  use neko
  use spectral_error_indicator
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0


  !> ========== Needed for Probes =================
  
  !> Probe type
  type(probes_t) :: pb

  !> Output variables
  type(file_t) :: fout
  type(matrix_t) :: mat_out

  !> Case IO parameters  
  integer            :: n_fields
  character(len=:), allocatable  :: output_file

  !> Output control
  logical :: write_output = .false.

  !> =============================================

  !> Spectral error indicator
  type(spec_err_ind_t) :: speri
  type(field_list_t) :: speri_l
  type(file_t) :: mf_speri
  character(len=NEKO_FNAME_LEN) :: fname_speri

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => check
  end subroutine user_setup
 
  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    !> Variables for bias
    real(kind=rp) :: arg, bias
    
    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z

  end subroutine set_scalar_boundary_conditions

  subroutine set_initial_conditions_for_u_and_s(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand, r,z
    s => neko_field_registry%get_field('s')

    !> Initialize with zero velocity
    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())

    !> Initialize with rand perturbations on temperature
    call rzero(s%x,w%dof%size())
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) 
    end do
    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e+s%msh%offset_el,rp)*real(i*j*k,rp))
                r = sqrt(s%dof%x(i,j,k,e)**2+s%dof%y(i,j,k,e)**2)
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001*rand*s%dof%x(i,j,k,e)*&
                                                    sin(3*pi*r/0.05_rp)*sin(10*pi*z)
             end do
          end do
       end do
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine set_initial_conditions_for_u_and_s

  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    !> Log variable
    character(len=LOG_SIZE) :: log_buf ! For logging status

    !> Support variables for probes 
    integer :: i
    type(matrix_t) :: mat_coords

    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    write(log_buf,*) 'Rayleigh Number is Ra=', Ra
    call neko_log%message(log_buf)
    
!    if (pe_rank.eq.0) write(*,*) 

    !> ========== Needed for Probes =================
    
    !> Read the output information
    call json_get(params, 'case.probes.output_file', output_file) 

    !> Probe set up
    !! Read probe info and initialize the controller, arrays, etc.
    call pb%init(t, params)
    !! Perform the set up of gslib_findpts
    call pb%setup(coef)
    !! Find the probes in the mesh. Map from xyz -> rst
    call pb%map(coef)
    !> Write a summary of the probe info
    call pb%show()

    !> Initialize the output
    fout = file_t(trim(output_file))
    call mat_out%init(pb%n_probes, pb%n_fields)

    !> Write coordinates in output file (imitate nek5000)
    !! Initialize the arrays
    call mat_coords%init(pb%n_probes,3)
    !! Array them as rows
    call transpose(mat_coords%x, pb%n_probes, pb%xyz, 3)
    !! Write the data to the file
    call fout%write(mat_coords)
    !! Free the memory 
    call mat_coords%free

    !> ==============================================
    
    !> Initialize the list of pointers to write out files
    call speri%init(u,v,w,coef)
    call list_init3(speri_l,speri%u_hat,speri%v_hat,speri%w_hat)
    fname_speri = 'speri.fld'
    mf_speri =  file_t(fname_speri)


  end subroutine user_initialize

  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param

    !> ========== Needed for Probes =================
    
    call pb%free
    call mat_out%free
    call file_free(fout)
    
    !> ==============================================
    
    call speri%free()
    call list_final3(speri_l)
    call file_free(mf_speri)
     
  end subroutine user_finalize

  subroutine set_bousinesq_forcing_term(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    !> ========== Needed for Probes =================

    !> Interpolate the desired fields
    call pb%interpolate(t,tstep, write_output)
    !! Write if the interpolate function returs write_output=.true.
    if (write_output) then
       call transpose(mat_out%x, pb%n_probes, pb%out_fields, pb%n_fields)
       call fout%write(mat_out, t)
       write_output = .false.
    
       call speri%get_indicators(coef)
       call mf_speri%write(speri_l,t)      
    end if

    !> ==============================================

  end subroutine check
  
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

end module user
