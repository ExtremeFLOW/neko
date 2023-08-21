module speri

  use neko
  use, intrinsic :: iso_c_binding
  implicit none
  
  !> include information needed for compressing fields
  type :: speri_t
     real(kind=rp), allocatable :: v(:,:) !< Transformation matrix
     real(kind=rp), allocatable :: vt(:,:) !< Transformation matrix transposed
     real(kind=rp), allocatable :: vinv(:,:) !< Transformation matrix inversed 
     real(kind=rp), allocatable :: vinvt(:,:) !< Transformation matrix
     !! inversed and transposed 
     real(kind=rp), allocatable :: w(:,:) !< Diagonal matrix with weights
     real(kind=rp), allocatable :: specmat(:,:) !< Transformation matrix
     real(kind=rp), allocatable :: specmatt(:,:) !< Transformation matrix
     type(field_t), pointer :: fld  => null()
     type(field_t) :: fldhat
     type(field_t) :: wk
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dof => null()
     !> From adam
     real(kind=rp) :: SERI_SMALL
     ! used for ratios
     real(kind=rp) :: SERI_SMALLR
     ! used for gradients
     real(kind=rp) :: SERI_SMALLG
     ! used for sigma and rtmp in error calculations
     real(kind=rp) :: SERI_SMALLS
     ! number of points in fitting
     integer :: SERI_NP
     integer :: SERI_NP_MAX
     ! last modes skipped
     integer :: SERI_ELR
     real(kind=rp), allocatable :: eind(:) !<
     real(kind=rp), allocatable :: sig(:) !<

     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: v_d = C_NULL_PTR
     type(c_ptr) :: vt_d = C_NULL_PTR
     type(c_ptr) :: vinv_d = C_NULL_PTR
     type(c_ptr) :: vinvt_d = C_NULL_PTR
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: specmat_d = C_NULL_PTR
     type(c_ptr) :: specmatt_d = C_NULL_PTR

  end type speri_t

  interface speri_init
          module procedure speri_init_all
  end interface speri_init

  public :: speri_init, speri_free

contains

        !> Initialize
        subroutine speri_init_all(speri, u)
                type(speri_t), intent(inout) :: speri
                type(field_t), intent(in), target :: u
                integer :: il, jl, aa

                call speri_free(speri)

                speri%fld => u
                speri%fldhat = u
                speri%wk = u
                speri%msh => u%msh
                speri%Xh => u%Xh
                speri%dof => u%dof

                ! Allocate arrays
                allocate(speri%v(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%vt(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%vinv(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%vinvt(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%w(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%specmat(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%specmatt(speri%Xh%lx, speri%Xh%lx))
                allocate(speri%eind(speri%msh%nelv))
                allocate(speri%sig(speri%msh%nelv))

                ! Initialize all the matrices
                call speri_generate_specmat(speri)

                ! Generate the uhat field (legendre coeff)

                call speri_goto_space(speri,'spec') !< 'spec' / 'phys'

                !> From ADAM
                ! set cutoff parameters
                ! used for values
                speri%SERI_SMALL = 1.e-14
                ! used for ratios
                speri%SERI_SMALLR = 1.e-10
                ! used for gradients
                speri%SERI_SMALLG = 1.e-5
                ! used for sigma and rtmp in error calculations
                speri%SERI_SMALLS = 0.2
                ! number of points in fitting
                speri%SERI_NP = 4
                speri%SERI_NP_MAX = 4
                ! last modes skipped
                speri%SERI_ELR = 0

                associate(LX1 => speri%Xh%lx, LY1 => speri%Xh%ly, &
                                LZ1 => speri%Xh%lz, &
                                SERI_SMALL  => speri%SERI_SMALL,  &
                                SERI_SMALLR => speri%SERI_SMALLR, &
                                SERI_SMALLG => speri%SERI_SMALLG, & 
                                SERI_SMALLS => speri%SERI_SMALLS, & 
                                SERI_NP     => speri%SERI_NP,     &
                                SERI_NP_MAX => speri%SERI_NP_MAX, &
                                SERI_ELR    => speri%SERI_ELR     &
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

        end subroutine speri_init_all


        !> Deallocate coefficients
        subroutine speri_free(speri)
                type(speri_t), intent(inout) :: speri

                if(allocated(speri%v)) then
                        deallocate(speri%v)
                end if

                if(allocated(speri%vt)) then
                        deallocate(speri%vt)
                end if

                if(allocated(speri%vinv)) then
                        deallocate(speri%vinv)
                end if

                if(allocated(speri%vinvt)) then
                        deallocate(speri%vinvt)
                end if

                if(allocated(speri%w)) then
                        deallocate(speri%w)
                end if

                if(allocated(speri%specmat)) then
                        deallocate(speri%specmat)
                end if

                if(allocated(speri%specmatt)) then
                        deallocate(speri%specmatt)
                end if

                if(allocated(speri%eind)) then
                        deallocate(speri%eind)
                end if

                if(allocated(speri%sig)) then
                        deallocate(speri%sig)
                end if

                call field_free(speri%fldhat)

                call field_free(speri%wk)

                nullify(speri%fld)
                nullify(speri%msh)
                nullify(speri%Xh)
                nullify(speri%dof)


                !
                ! Cleanup the device (if present)
                !

                if (c_associated(speri%v_d)) then
                        call device_free(speri%v_d)
                end if

                if (c_associated(speri%vt_d)) then
                        call device_free(speri%vt_d)
                end if

                if (c_associated(speri%vinv_d)) then
                        call device_free(speri%vinv_d)
                end if

                if (c_associated(speri%vinvt_d)) then
                        call device_free(speri%vinvt_d)
                end if

                if (c_associated(speri%w_d)) then
                        call device_free(speri%w_d)
                end if

                if (c_associated(speri%specmat_d)) then
                        call device_free(speri%specmat_d)
                end if

                if (c_associated(speri%specmatt_d)) then
                        call device_free(speri%specmatt_d)
                end if

        end subroutine speri_free


        !> Generate spectral tranform matrices
        subroutine speri_generate_specmat(speri)
                type(speri_t), intent(inout) :: speri
                real(kind=rp) :: L(0:speri%Xh%lx-1)
                real(kind=rp) :: delta(speri%Xh%lx)
                integer :: i, kj, j, j2, kk
                character(len=LOG_SIZE) :: log_buf 

                associate(Xh => speri%Xh, v=> speri%v, vt => speri%vt, &
                                vinv => speri%vinv, vinvt => speri%vinvt, w => speri%w)
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
                                speri%w(i,j) = 0
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
                        call device_map(speri%v,     speri%v_d,     speri%Xh%lxy)
                        call device_map(speri%vt,    speri%vt_d,    speri%Xh%lxy)
                        call device_map(speri%vinv,  speri%vinv_d,  speri%Xh%lxy)
                        call device_map(speri%vinvt, speri%vinvt_d, speri%Xh%lxy)
                        call device_map(speri%w,     speri%w_d,     speri%Xh%lxy)
                        !Map the following pointers but do not copy data for them
                        call device_map(speri%specmat,  speri%specmat_d,  speri%Xh%lxy)
                        call device_map(speri%specmatt, speri%specmatt_d, speri%Xh%lxy)


                        call device_memcpy(speri%v,     speri%v_d,     speri%Xh%lxy, &
                                HOST_TO_DEVICE)
                        call device_memcpy(speri%vt,    speri%vt_d,    speri%Xh%lxy, &
                                HOST_TO_DEVICE)
                        call device_memcpy(speri%vinv,  speri%vinv_d,  speri%Xh%lxy, &
                                HOST_TO_DEVICE)
                        call device_memcpy(speri%vinvt, speri%vinvt_d, speri%Xh%lxy, &
                                HOST_TO_DEVICE)
                        call device_memcpy(speri%w,     speri%w_d,     speri%Xh%lxy, &
                                HOST_TO_DEVICE)

                end if

        end subroutine speri_generate_specmat


        !> Tranform to spectral space (using tensor product)
        !the result of the transform is given in fldhat
        subroutine speri_goto_space(speri, space)
                type(speri_t), intent(inout) :: speri
                integer :: i, j, k, e, nxyz, nelv, n
                character(len=LOG_SIZE) :: log_buf 
                character(len=4) :: space 

                ! define some constants
                nxyz = speri%Xh%lx*speri%Xh%lx*speri%Xh%lx
                nelv = speri%msh%nelv
                n    = nxyz*nelv

                if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
                        (NEKO_BCKND_OPENCL .eq. 1)) then 

                        if (pe_rank.eq.0) write(*,*) 'SpErInd: Transform in the GPU'

                        ! Define the matrix according to which transform to do 
                        if (space .eq. 'spec') then
                                call device_copy(speri%specmat_d,  speri%vinv_d,  speri%Xh%lxy)
                                call device_copy(speri%specmatt_d, speri%vinvt_d, speri%Xh%lxy)
                                call device_copy(speri%wk%x_d, speri%fld%x_d, n)
                        endif
                        if (space .eq. 'phys') then
                                call device_copy(speri%specmat_d,  speri%v_d,  speri%Xh%lxy)
                                call device_copy(speri%specmatt_d, speri%vt_d, speri%Xh%lxy)
                                call device_copy(speri%wk%x_d, speri%fldhat%x_d, n)
                        endif

                else

                        if (pe_rank.eq.0) write(*,*) 'SpErInd: Transform in the CPU'

       ! Define the matrix according to which transform to do 
       if (space .eq. 'spec') then
          call copy(speri%specmat, speri%vinv, speri%Xh%lx*speri%Xh%lx)
          call copy(speri%specmatt, speri%vinvt, speri%Xh%lx*speri%Xh%lx)
          call copy(speri%wk%x,speri%fld%x,n)
       endif
       if (space .eq. 'phys') then
          call copy(speri%specmat, speri%v, speri%Xh%lx*speri%Xh%lx)
          call copy(speri%specmatt, speri%vt, speri%Xh%lx*speri%Xh%lx)
          call copy(speri%wk%x,speri%fldhat%x,n)
       endif

    end if

    call tnsr3d(speri%fldhat%x, speri%Xh%lx, speri%wk%x, &
                speri%Xh%lx,speri%specmat, &
                speri%specmatt, speri%specmatt, nelv)

    !! Synchronize
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1)) then 

       call device_memcpy(speri%fldhat%x,speri%fldhat%x_d, n, &
                          DEVICE_TO_HOST)
    end if

  end subroutine speri_goto_space

  subroutine speri_get(speri)
    type(speri_t), intent(inout) :: speri
    real(kind=rp) :: xa(speri%Xh%lx,speri%Xh%ly,speri%Xh%lz)
    real(kind=rp) :: xb(speri%Xh%lx,speri%Xh%ly,speri%Xh%lz)
    integer :: i, e

    associate(eind        => speri%eind, &
              sig         => speri%sig , &
              lnelt       => speri%msh%nelv, &
              LX1         => speri%Xh%lx, &
              LY1         => speri%Xh%ly, &
              LZ1         => speri%Xh%lz, &
              var         => speri%fldhat%x &
             )   
      
       ! zero arrays
       call rzero(eind,lnelt)
       call rzero(sig,lnelt)
       
       ! Get the indicator
       call speri_var(speri, eind,sig,var,lnelt,xa,xb, LX1, LY1, LZ1)

       ! Put the indicator in the fieldhat to write it
       do e = 1,lnelt
          do i = 1,LX1*LY1*LZ1
             var(i,1,1,e) = eind(e)
          end do
       end do

     end associate

  end subroutine

  subroutine speri_var(speri, est,sig,var,nell,xa,xb,LX1,LY1,LZ1)
    type(speri_t), intent(inout) :: speri
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
    real(kind=rp) ::  coefx(speri%SERI_NP_MAX,LY1,LZ1), & 
                      coefy(speri%SERI_NP_MAX,LX1,LZ1), &
                      coefz(speri%SERI_NP_MAX,LX1,LY1)
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
        if (coef11.ge.speri%SERI_SMALL) then
           ! extrapolate coefficients
           ! X - direction
           ! copy last SERI_NP collumns (or less if NX1 is smaller)
           ! SERI_ELR allows to exclude last row
            j_st = max(1,LX1-speri%SERI_NP+1-speri%SERI_ELR)
            j_en = max(1,LX1-speri%SERI_ELR)
            do ll = 1,LZ1
                do kl = 1,LY1
                    do jl = j_st,j_en
                        coefx(j_en-jl+1,kl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo
            ! get extrapolated values
            call speri_extrap(speri,estx,sigx,coef11,coefx, &
                 j_st,j_en,LY1,LZ1)
         
            ! Y - direction
            ! copy last SERI_NP collumns (or less if NY1 is smaller)
            ! SERI_ELR allows to exclude last row
            j_st = max(1,LY1-speri%SERI_NP+1-speri%SERI_ELR)
            j_en = max(1,LY1-speri%SERI_ELR)
            do ll = 1,LZ1
                do kl = j_st,j_en
                    do jl = 1,LX1
                        coefy(j_en-kl+1,jl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

            ! get extrapolated values
            call speri_extrap(speri, esty,sigy,coef11,coefy, &
                j_st,j_en,LX1,LZ1)
   
            ! Z - direction
            ! copy last SERI_NP collumns (or less if NZ1 is smaller)
            ! SERI_ELR allows to exclude last row
            j_st = max(1,LZ1-speri%SERI_NP+1-speri%SERI_ELR)
            j_en = max(1,LZ1-speri%SERI_ELR)
            do ll = j_st,j_en
                do kl = 1,LY1
                    do jl = 1,LX1
                        coefz(j_en-ll+1,jl,kl) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

            ! get extrapolated values
            call speri_extrap(speri, estz,sigz,coef11,coefz, &
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

  end subroutine


  subroutine speri_extrap(speri,estx,sigx,coef11,coef, &
                ix_st,ix_en,nyl,nzl)
      implicit none
      type(speri_t), intent(inout) :: speri
      ! argument list
      integer :: ix_st,ix_en,nyl,nzl
      ! Legendre coefficients; last SERI_NP columns
      real(kind=rp) :: coef(speri%SERI_NP_MAX,nyl,nzl)
      ! Legendre coefficients; first value coeff(1,1,1)
      real(kind=rp) :: coef11
      ! estimated error and decay rate
      real(kind=rp) :: estx, sigx

      ! local variables
      integer :: il, jl, kl, ll  ! loop index
      integer :: nsigt, pnr, nzlt
      real(kind=rp) :: sigt, smallr, cmin, cmax, cnm, rtmp, rtmp2, rtmp3
      real(kind=rp) :: sumtmp(4), cffl(speri%SERI_NP_MAX)
      real(kind=rp) :: stmp, estt, clog, ctmp, cave, erlog
      logical :: cuse(speri%SERI_NP_MAX)
!-----------------------------------------------------------------------
      
     associate(LX1 => speri%Xh%lx, LY1 => speri%Xh%ly, &
              LZ1 => speri%Xh%lz, &
              SERI_SMALL  => speri%SERI_SMALL,  &
              SERI_SMALLR => speri%SERI_SMALLR, &
              SERI_SMALLG => speri%SERI_SMALLG, & 
              SERI_SMALLS => speri%SERI_SMALLS, & 
              SERI_NP     => speri%SERI_NP,     &
              SERI_NP_MAX => speri%SERI_NP_MAX, &
              SERI_ELR    => speri%SERI_ELR     &
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

      end subroutine



end module speri



module user
  use neko
  use speri
  implicit none

  !> Variables to store the Rayleigh and Prandlt numbers
  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Re = 0
  real(kind=rp) :: Pr = 0

  !> Arrays asociated with Method#1 for nusselt calculation
  type(field_t) :: work_field ! Field to perform operations
  type(field_t) :: uzt ! u_z * T
  type(field_t) :: dtdx ! Derivative of scalar wrt x
  type(field_t) :: dtdy ! Detivative of scalar wrt y
  type(field_t) :: dtdz ! Derivative of scalar wrt z
  type(field_t) :: dtdn ! Derivative of scalar wrt normal
  type(field_t) :: mass_area_top ! mass matrix for area on top wall
  type(field_t) :: mass_area_bot ! mass matrix for area on bottom wall
  type(field_t) :: mass_area_side ! mass matrix for area on top wall
  type(field_t) :: bm1 ! mass matrix for area on bottom wall
  real(kind=rp) :: bar_uzt ! Volume integral
  real(kind=rp) :: top_wall_bar_dtdz ! area integral
  real(kind=rp) :: bot_wall_bar_dtdz ! area integral
  real(kind=rp) :: top_wall_area ! area of top and bottom wall
  real(kind=rp) :: bot_wall_area ! area of top and bottom wall
  real(kind=rp) :: side_wall_area ! area of top and bottom wall
  
  !> Boundary conditions  
  integer :: istep = 1

  !> Variables to write extra files
  character(len=NEKO_FNAME_LEN) :: fname_dtdx
  character(len=NEKO_FNAME_LEN) :: fname_dtdy
  character(len=NEKO_FNAME_LEN) :: fname_dtdz
  character(len=NEKO_FNAME_LEN) :: fname_top_area
  character(len=NEKO_FNAME_LEN) :: fname_bot_area
  character(len=NEKO_FNAME_LEN) :: fname_side_area
  character(len=NEKO_FNAME_LEN) :: fname_bm1
  character(len=NEKO_FNAME_LEN) :: fname_speri_x
  character(len=NEKO_FNAME_LEN) :: fname_speri_y
  character(len=NEKO_FNAME_LEN) :: fname_speri_z
  type(file_t) :: mf_dtdx
  type(file_t) :: mf_dtdy
  type(file_t) :: mf_dtdz
  type(file_t) :: mf_top_area
  type(file_t) :: mf_bot_area
  type(file_t) :: mf_side_area
  type(file_t) :: mf_bm1
  type(file_t) :: mf_speri_x
  type(file_t) :: mf_speri_y
  type(file_t) :: mf_speri_z

  !> List of elements and facets in uper and lower boundary
  type(stack_i4t4_t) :: sidewall_facet
  type(stack_i4t4_t) :: wall_facet
  type(tuple4_i4_t)  :: facet_el_type_0

  !> Spectral error indicator
  type(speri_t) :: speri_u
  type(speri_t) :: speri_v
  type(speri_t) :: speri_w

  !> Mesh deformation
  real(kind=rp) :: delta_mesh = 0.0_rp !How much to deform to the top and bottom plate. =0 means no deformation

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_mesh_setup => deform_mesh
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
    u%fluid_user_ic => set_initial_conditions_for_u_and_s
    u%scalar_user_bc => set_scalar_boundary_conditions
    u%fluid_user_f_vector => set_bousinesq_forcing_term
    u%user_check => calculate_nusselt
  end subroutine user_setup

  subroutine deform_mesh(msh)
    type(mesh_t), intent(inout) :: msh
    msh%apply_deform => redistribute_elements
  end subroutine deform_mesh
  
  subroutine redistribute_elements(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    type(tuple_i4_t) :: el_and_facet
    real(kind=rp) :: th
    integer :: e, i, j ,k, l,  facet, nxyze
    real(kind=rp) :: Betaz, x_pt, y_pt, z_pt, z_pt_def
    
    nxyze = lx*ly*lz*msh%nelv
    if (delta_mesh .gt. 1e-8_rp) then
       Betaz = delta_mesh 
       do i = 1, nxyze
          z_pt = z(i,1,1,1)
          z_pt_def = 0.5*(tanh(Betaz*(2*z_pt-1.0))/tanh(Betaz) + 1.0)
          z(i,1,1,1) = z_pt_def
       end do
    end if
  end subroutine redistribute_elements
  
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
    
    arg  = -istep*0.0005
    bias = x*0.2*exp(arg)
    bias = 0

    ! This will be used on all zones without labels
    ! e.g. the ones hardcoded to 'v', 'w', etcetc
    s = 1.0_rp - z + bias


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
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
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

    !> Parameters to populate the list of elements and facets
    real(kind=rp) :: stack_size(1) !Do it like this so it is cast
    real(kind=rp) :: stack_size_global
    integer :: e, facet, typ, e_counter, facet_counter, i,j,k,n
    real(kind=rp) :: normal(3)
    real(kind=rp) :: area_mass, area_rank(1),area_global(1)
    type(tuple4_i4_t), pointer :: wall_facet_list(:)
    type(tuple4_i4_t), pointer :: sidewall_facet_list(:)
    real(kind=rp) :: voll(1), voll_temp(1)
    integer :: lx, ly, lz
    
    !> Recalculate the non dimensional parameters
    call json_get(params, 'case.scalar.Pr', Pr)
    call json_get(params, 'case.fluid.Re', Re)
    Ra = (Re**2)*Pr
    if (pe_rank.eq.0) write(*,*) 'Rayleigh Number is Ra=', Ra
    
    !> Initialize variables related to nusselt calculation
    call field_init(work_field, u%dof, 'work_field')
    call field_init(uzt, u%dof, 'uzt')
    call field_init(dtdx, u%dof, 'dtdx')
    call field_init(dtdy, u%dof, 'dtdy')
    call field_init(dtdz, u%dof, 'dtdz')
    call field_init(dtdn, u%dof, 'dtdn')
    call field_init(mass_area_top, u%dof, 'mat')
    call field_init(mass_area_bot, u%dof, 'mab')
    call field_init(mass_area_side, u%dof, 'masd')
    call field_init(bm1, u%dof, 'mass_mat')

    !> Initialize the file
    fname_dtdx = 'dtdx.fld'
    fname_dtdy = 'dtdy.fld'
    fname_dtdz = 'dtdz.fld'
    fname_top_area = 'top_area.fld'
    fname_bot_area = 'bot_area.fld'
    fname_side_area = 'side_area.fld'
    fname_bm1 = 'bm1.fld'
    fname_speri_x = 'speri_x.fld'
    fname_speri_y = 'speri_y.fld'
    fname_speri_z = 'speri_z.fld'
    mf_dtdx =  file_t(fname_dtdx)
    mf_dtdy =  file_t(fname_dtdy)
    mf_dtdz =  file_t(fname_dtdz)
    mf_top_area =  file_t(fname_top_area)
    mf_bot_area =  file_t(fname_bot_area)
    mf_side_area =  file_t(fname_side_area)
    mf_bm1 =  file_t(fname_bm1)
    mf_speri_x =  file_t(fname_speri_x)
    mf_speri_y =  file_t(fname_speri_y)
    mf_speri_z =  file_t(fname_speri_z)


    !> Initialize list with upper and lower wall facets
    call wall_facet%init()
    call sidewall_facet%init()

    !> Populate the list with upper and lower and wall facets 
    do e = 1, u%msh%nelv !Go over all elements
      do facet = 1, 6 ! Go over all facets of hex element
        normal = coef_get_normal(coef,1,1,1,e,facet) ! Get facet normal
        typ =  u%msh%facet_type(facet,e)             ! Get facet type
        if (typ.ne.0) then !then it is a boundary facet
          if (abs(normal(3)).ge.1-1e-5) then !then it is on the plates
            !write(*,*) 'In boundary: facet=', facet, 'and element=', e
            !write(*,*) 'BC facet type ', typ 
            !write(*,*) 'normal to this facet is: ', normal   
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call wall_facet%push(facet_el_type_0)
          else !then it is on the sidewall
            !write(*,*) 'In boundary: facet=', facet, 'and element=', e
            !write(*,*) 'BC facet type ', typ 
            !write(*,*) 'normal to this facet is: ', normal   
            facet_el_type_0%x = (/facet, e, typ, 0/)
            call wall_facet%push(facet_el_type_0)
          end if
        end if
      end do
    end do
    !! Determine the number of facets in the walls
    stack_size = real(wall_facet%top_,kind=rp)
    stack_size_global = glsum(stack_size,1)
    if (pe_rank .eq. 0) then
       write(*,*) 'Facets at the wall in this rank: ', stack_size
       write(*,*) 'Facets at the wall global: ', stack_size_global
    end if

    !> Fill up arrays to serve as area mass matrix for integration
    !! we want <dt/dz>_A,t at z=0,1.
    !! Average in Area is: (dt/dz * normal * Area_mass) / Area_total
    wall_facet_list => wall_facet%array()
    sidewall_facet_list => sidewall_facet%array()
    !! Initialize the mass as 0. This serves as a mask for nodes that are not integrated
    call rzero(mass_area_top%x,u%dof%size())
    call rzero(mass_area_bot%x,u%dof%size())
    call rzero(mass_area_side%x,u%dof%size())
    !! Fill up the top and bottom mass matrices
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz
    do e_counter = 1, wall_facet%top_
      e = wall_facet_list(e_counter)%x(2)
      facet = wall_facet_list(e_counter)%x(1)
      select case(facet)
         case(1)
            do j = 1 , ly
               do k = 1 , lz
                  mass_area_side%x(1,j,k,e) = coef%area(j,k,facet,e) 
               end do
            end do
         case(2)
            do j = 1 , ly
               do k = 1 , lz
                  mass_area_side%x(lx,j,k,e) = coef%area(j,k,facet,e) 
               end do
            end do
         case(3)
            do i = 1 , lx
               do k = 1 , lz
                  mass_area_side%x(i,1,k,e) = coef%area(i,k,facet,e)    
               end do
            end do
         case(4)
            do i = 1 , lx
               do k = 1 , lz
                  mass_area_side%x(i,ly,k,e) = coef%area(i,k,facet,e)       
               end do
            end do
         case(5)
            normal = coef_get_normal(coef,1,1,1,e,facet)
            do i = 1 , lx
               do j = 1 , ly
                  mass_area_bot%x(i,j,1,e) = coef%area(i,j,facet,e)*normal(3) 
               end do
            end do
         case(6)
            normal = coef_get_normal(coef,1,1,1,e,facet)
            do i = 1 , lx
               do j = 1 , ly
                  mass_area_top%x(i,j,lz,e) = coef%area(i,j,facet,e)*normal(3)       
               end do
            end do
      end select
    end do

    n = size(coef%B)
    top_wall_area = glsum(mass_area_top%x,n)
    bot_wall_area = glsum(mass_area_bot%x,n)
    side_wall_area = glsum(mass_area_side%x,n)
    !! Write it out
    if (pe_rank .eq. 0) then
       write(*,*) 'Area of top wall * normal= ', top_wall_area
       write(*,*) 'Area of bottom wall * normal= ', bot_wall_area
       write(*,*) 'Area of side wall= ', side_wall_area
    end if 
    !! Put it also on device
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_memcpy(mass_area_top%x,mass_area_top%x_d, &
                          n,HOST_TO_DEVICE)
       call device_memcpy(mass_area_bot%x,mass_area_bot%x_d, &
                          n,HOST_TO_DEVICE)
       call device_memcpy(mass_area_side%x,mass_area_side%x_d, &
                          n,HOST_TO_DEVICE)
    end if

    !> Store volume mass matrix for post processing
    n = size(coef%B) 
    call rzero(bm1%x,u%dof%size())
    call copy(bm1%x,coef%B,n)
    voll = glsum(bm1%x,n)
    if (pe_rank .eq. 0) then
       write(*,*) 'volume=', voll
    end if 
    !! rescale mass if volume is too small
    do i = 1,n
       bm1%x(i,1,1,1)=bm1%x(i,1,1,1) * 1e3
    end do
    voll = glsum(bm1%x,n)
    if (pe_rank .eq. 0) then
       write(*,*) 'rescaled volume=', voll
    end if 

    !> Perform IO
    call mf_top_area%write(mass_area_top,t)
    call mf_bot_area%write(mass_area_bot,t)
    call mf_side_area%write(mass_area_side,t)
    call mf_bm1%write(bm1,t)      

  end subroutine user_initialize

  pure function get_area_mass(coef, i, j, k, e, facet) result(area_mass)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: i, j, k, e, facet
    real(kind=rp) :: area_mass
      
    select case (facet)               
      case(1,2)
        area_mass = coef%area(j, k, facet, e)
      case(3,4)
        area_mass = coef%area(i, k, facet, e)
      case(5,6)
        area_mass = coef%area(i, j, facet, e)
      end select
  end function get_area_mass


  subroutine user_finalize(t, param)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: param

    ! Finalize variables related to nusselt calculation
    call field_free(work_field)
    call field_free(uzt)
    call field_free(dtdx)
    call field_free(dtdy)
    call field_free(dtdz)
    call field_free(dtdn)
    call field_free(mass_area_top)
    call field_free(mass_area_bot)
    call field_free(mass_area_side)
    call field_free(bm1)
    
    ! Finilize list that contains uper and lower wall facets
    call wall_facet%free()
    call sidewall_facet%free()
  
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

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term

  subroutine calculate_nusselt(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(field_t), pointer :: s
    integer :: n,ntot, i,j,k, facet, e, lx,ly,lz
    integer :: index(4)
    type (tuple_i4_t) :: facet_el
    real(kind=rp) :: normal(3)
    real(kind=rp) :: voll(1), voll_temp(1)
    ! output control
    integer :: calc_freq = 0 ! How frequently should we calculate Nu
    integer :: calc_frequency = 0 ! How frequently should we calculate Nu
    logical :: verify_bc = .false. ! write boundary conditions
    logical :: calculate_now = .false. ! write boundary conditions
    character(len=:), allocatable :: oc
    real(kind=rp) :: dt
    
    !> This value is used for breaking symtetries in bc
    istep = istep + 1

    !> Get the control parameters
    call json_get(params, 'case.timestep', dt)
    call json_get(params, 'case.monitor.output_control', oc)
    call json_get(params, 'case.monitor.calc_frequency', calc_freq)
    call json_get(params, 'case.monitor.verify_bc', verify_bc)
 
    if (oc(1:14).eq.'simulationtime') then
        calc_frequency = int(calc_freq / dt)
    else
        calc_frequency = int(calc_freq)
    end if

    if (mod(tstep,calc_frequency).ne.0) return

    s => neko_field_registry%get_field('s')
    n = size(coef%B)
    ntot = coef%dof%size()
    lx = u%Xh%lx
    ly = u%Xh%ly
    lz = u%Xh%lz

    !> ------ Method #1 for nusselt calculation -----
    !> Nu_v = 1 + sqrt(Ra) * <u_z * T>_{v,t}
    !> Steps:
       !1.    Calculate the convective current T*u_z
       !2.    Get <u_z * T>_{v}(t) (Average in volume)
       !!2.1. Multiply field with mass matrix
       !!2.2. Perform a global sum to get the integral
       !!2.3. Normalize with total volume to get the average
       !3.    Get <u_z * T>_{v,t} by averaging time signal
       !4.    Multiply average by sqrt(Ra) and sum 1
    !> Steps 1. and 2. are done here. Do 3. and 4. in post  
    if (NEKO_BCKND_DEVICE .eq. 1) then 
       call device_col3(uzt%x_d,w%x_d,s%x_d,n)             !1.  
       call device_col3(work_field%x_d,uzt%x_d,coef%B_d,n) !2.1.   
       bar_uzt = device_glsum(work_field%x_d,n)            !2.2.
       bar_uzt = bar_uzt / coef%volume                     !2.3. 
    else
       call col3(uzt%x,w%x,s%x,n)                          !1.
       call col3(work_field%x,uzt%x,coef%B,n)              !2.1.
       bar_uzt = glsum(work_field%x,n)                     !2.2.
       bar_uzt = bar_uzt / coef%volume                     !2.3.
    end if

    !> ------ Method #2 for nusselt calculation -----
    ! Calculate derivatives. Automatically in device with opgrad
    call dudxyz (dtdx%x, s%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (dtdy%x, s%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (dtdz%x, s%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then 
       ! Calculate for top wall
       call device_col3(work_field%x_d,dtdz%x_d,mass_area_top%x_d,n)      
       top_wall_bar_dtdz = device_glsum(work_field%x_d,n)                  
       top_wall_bar_dtdz = top_wall_bar_dtdz / abs(top_wall_area)            
       ! Calculate for bot wall
       call device_col3(work_field%x_d,dtdz%x_d,mass_area_bot%x_d,n)      
       bot_wall_bar_dtdz = device_glsum(work_field%x_d,n)                  
       bot_wall_bar_dtdz = bot_wall_bar_dtdz / abs(bot_wall_area)            
    else
       ! Calculate for top wall
       call col3(work_field%x,dtdz%x,mass_area_top%x,n)      
       top_wall_bar_dtdz = glsum(work_field%x,n)                  
       top_wall_bar_dtdz = top_wall_bar_dtdz / abs(top_wall_area)            
       ! Calculate for bot wall
       call col3(work_field%x,dtdz%x,mass_area_bot%x,n)      
       bot_wall_bar_dtdz = glsum(work_field%x,n)                  
       bot_wall_bar_dtdz = bot_wall_bar_dtdz / abs(bot_wall_area)            
    end if
    
    !> write variables to monitor
    !! Integral quantities
    if (pe_rank .eq. 0) then
       open(10,file="nusselt.txt",position="append")
       write(10,*) t,'', bar_uzt, '', top_wall_bar_dtdz, '', &
                   bot_wall_bar_dtdz
       close(10)
    end if
    !! Fields
    call device_memcpy(dtdx%x,dtdx%x_d, n,DEVICE_TO_HOST)
    call device_memcpy(dtdy%x,dtdy%x_d, n,DEVICE_TO_HOST)
    call device_memcpy(dtdz%x,dtdz%x_d, n,DEVICE_TO_HOST)
    call mf_dtdz%write(dtdz,t)


    if (verify_bc.eqv..true.) then
    !> Calculate temp flux on each element boundary and chech boundary to verify BC
     !! Zero out the vector
     call rzero(dtdn%x,dtdn%dof%size())
     !! Calculate the normal component for each facet in the domain
     !! based on functions coef_get_normal and index_is_on_facet
     do e = 1, u%msh%nelv !Go over all elements
      do facet = 1, 6 ! Go over all facets of hex element
         select case(facet)
            case(1)
               do j = 1 , dtdn%Xh%ly
                  do k = 1 , dtdn%Xh%lz
                     dtdn%x(1,j,k,e) = coef%nx(j, k, facet, e) &
                                    *dtdx%x(1, j, k, e)  &
                                    +coef%ny(j, k, facet, e)  &
                                    *dtdy%x(1, j, k, e)  &
                                    +coef%nz(j, k, facet, e)  &
                                    *dtdz%x(1, j, k, e) 
                  end do
               end do
            case(2)
               do j = 1 , dtdn%Xh%ly
                  do k = 1 , dtdn%Xh%lz
                     dtdn%x(dtdn%Xh%lx,j,k,e) = coef%nx(j, k, facet, e) &
                                    *dtdx%x(lx, j, k, e)  &
                                    +coef%ny(j, k, facet, e)  &
                                    *dtdy%x(lx, j, k, e)  &
                                    +coef%nz(j, k, facet, e)  &
                                    *dtdz%x(lx, j, k, e) 
                  end do
               end do
            case(3)
               do i = 1 , dtdn%Xh%lx
                  do k = 1 , dtdn%Xh%lz
                     dtdn%x(i,1,k,e) = coef%nx(i, k, facet, e) &
                                    *dtdx%x(i, 1, k, e)  &
                                    +coef%ny(i, k, facet, e)  &
                                    *dtdy%x(i, 1, k, e)  &
                                    +coef%nz(i, k, facet, e)  &
                                    *dtdz%x(i, 1, k, e) 
                  end do
               end do
            case(4)
               do i = 1 , dtdn%Xh%lx
                  do k = 1 , dtdn%Xh%lz
                     dtdn%x(i,dtdn%Xh%ly,k,e) = coef%nx(i, k, facet, e) &
                                    *dtdx%x(i, ly, k, e)  &
                                    +coef%ny(i, k, facet, e)  &
                                    *dtdy%x(i, ly, k, e)  &
                                    +coef%nz(i, k, facet, e)  &
                                    *dtdz%x(i, ly, k, e) 
                  end do
               end do
            
            case(5)
               do i = 1 , dtdn%Xh%lx
                  do j = 1 , dtdn%Xh%ly
                     dtdn%x(i,j,1,e) = coef%nx(i, j, facet, e) &
                                    *dtdx%x(i, j, 1, e)  &
                                    +coef%ny(i, j, facet, e)  &
                                    *dtdy%x(i, j, 1, e)  &
                                    +coef%nz(i, j, facet, e)  &
                                    *dtdz%x(i, j, 1, e) 
                  end do
               end do
            case(6)
               do i = 1 , dtdn%Xh%lx
                  do j = 1 , dtdn%Xh%ly
                     dtdn%x(i,j,dtdn%Xh%lz,e) = coef%nx(i, j, facet, e) &
                                    *dtdx%x(i, j, lz, e)  &
                                    +coef%ny(i, j, facet, e)  &
                                    *dtdy%x(i, j, lz, e)  &
                                    +coef%nz(i, j, facet, e)  &
                                    *dtdz%x(i, j, lz, e) 
                  end do
               end do
         end select
      end do
     end do
     call mf_dtdx%write(dtdn,t)
    end if

    !> Get the spectral error indicators for the mesh
    !! Initialize spectral error indicator
    call speri_init(speri_u,u)
    call speri_init(speri_v,v)
    call speri_init(speri_w,w)
    !! calculate spectral error indicator
    call speri_get(speri_u)
    call speri_get(speri_v)
    call speri_get(speri_w)
    !! Write the data into a field
    call mf_speri_x%write(speri_u%fldhat,t)      
    call mf_speri_y%write(speri_v%fldhat,t)      
    call mf_speri_z%write(speri_w%fldhat,t)      
    !! Finalize specral error indicator
    call speri_free(speri_u)
    call speri_free(speri_v)
    call speri_free(speri_w)

  end subroutine calculate_nusselt

end module user
