!--------------------------------------------------------------------------!
! from KTH framework which was the backbone for this
!--------------------------------------------------------------------------!
!> @file trip.f
!! @ingroup trip_line
!! @brief Tripping function for AMR version of nek5000
!! @note  This version uses developed framework parts. This is because
!!   I'm in a hurry and I want to save some time writing the code. So
!!   I reuse already tested code and focuse important parts. For the
!!   same reason for now only lines parallel to z axis are considered. 
!!   The tripping is based on a similar implementation in the SIMSON code
!!   (Chevalier et al. 2007, KTH Mechanics), and is described in detail 
!!   in the paper Schlatter & Örlü, JFM 2012, DOI 10.1017/jfm.2012.324.
!! @author Adam Peplinski
!! @date May 03, 2018
!--------------------------------------------------------------------------!

! General note:
! trip was originally converted from Nek5000 to Neko by Martin Karp
! it was later separated into a separate module by Victor Baconnet
! at every stage this was done with a specific usage in mind, therefore
! the implementation is a bit hacky and inconsistent with the more modern parts of the code

!==========================================================================!
!! NOTE
!!
!! The tripping only works for 1 line parallel to the y or z axis (can be
!! extended to x-axis as well if needed). See above explanation for details.
!!
!! Usage in user file:
!!  - Depending on if you need along y or z axis, change trip_1dprj_y or
!!    trip_1dprj_z at lines 241 and 355.
!!  - add `use trip` under `use neko`
!!  - declare global object of type(trip_t) (e.g. `tripper`)
!!  - call tripper%init in `user_init_modules`
!!  - call tripper%apply_device() in fluid_user_f_vector
!==========================================================================!


module trip
  use neko
  use device_inhom_dirichlet
  implicit none

! max number of lines and Fourier modes
  integer, parameter :: trip_nline_max=2
  integer, parameter :: trip_nmode_max=500
!!  real(kind=rp), parameter :: pi = 4.0_rp*atan(1.0_rp)
! max number of random phase sets stored; 1- time independent, 2, 3 and 4 - time dependent
! I keep two old random pahase sets to get correct restart after AMR refinement
  integer, parameter :: trip_nset_max = 4
  integer, parameter :: ldim = 3 !lets hardcode this for now gosh


  type, public :: trip_t 

     type(dofmap_t), pointer :: dof
     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh
     integer :: id
   
     ! timer id
     integer :: tmr_id
   
   ! initialisation flag
     logical :: ifinit = .false.
   
   ! runtime parameter part
   ! section id
     integer :: sec_id    
   ! parameter section
     integer :: nline                  !< @var number of tripping lines
     integer :: nline_id
     real(kind=rp) :: tiamp                     !< @var time independent amplitude
     integer :: tiamp_id
     real(kind=rp) :: tdamp                     !< @var time dependent amplitude
     integer :: tdamp_id
     real(kind=rp) :: spos(3,trip_nline_max) !< @var coordinates of starting point of tripping line
     integer :: spos_id(3,trip_nline_max)
     real(kind=rp) :: epos(3,trip_nline_max) !< @var coordinates of ending point of tripping line
     integer :: epos_id(3,trip_nline_max)
     real(kind=rp) :: smth(3,trip_nline_max) !< @var smoothing radius
     integer :: smth_id(3,trip_nline_max)
     logical :: lext(trip_nline_max)    !< @var do we extend a line beyond starting and endig points
     integer :: lext_id(trip_nline_max)
     real(kind=rp) :: rota(trip_nline_max)      !< @var elipse rotation angle
     integer :: rota_id(trip_nline_max)
     integer :: nmode(trip_nline_max)  !< @var number of Fourier modes
     integer :: nmode_id(trip_nline_max)
     real(kind=rp) :: tdt(trip_nline_max)       !< @var time step for tripping
     integer :: tdt_id(trip_nline_max)
   
   ! inverse line length
     real(kind=rp) :: ilngt(trip_nline_max)
   
   ! inverse smoothing radius
     real(kind=rp) :: ismth(3,trip_nline_max)
     
   ! projection of 3D pionts on 1D line
     real(kind=rp), allocatable :: prj(:,:)
   
   ! number of points in 1D projection
     integer npoint(trip_nline_max)
     
   ! mapping of 3D array to 1D projection array
     integer, allocatable :: map(:,:,:,:,:)
   
   ! function for smoothing of the forcing
     real(kind=rp), allocatable :: fsmth(:,:,:,:,:)
   
   ! mask for tripping
     integer, allocatable :: mask(:)
     type(c_ptr) :: mask_d = C_NULL_PTR
     real(kind=rp), allocatable :: fsmth_mask(:)
     type(c_ptr) :: fsmth_mask_d = C_NULL_PTR
   ! forces for trippping
     type(c_ptr) :: ftripx_d = C_NULL_PTR
     type(c_ptr) :: ftripy_d = C_NULL_PTR
     type(c_ptr) :: ftripz_d = C_NULL_PTR
   ! force for interpolation
     type(c_ptr) :: f_interpolate_d(trip_nset_max)
   ! seed for random number generator; different for each line
     integer :: seed(trip_nline_max)
   
   ! number of tripping time intervals
     integer :: ntdt(trip_nline_max), ntdt_old(trip_nline_max)
     
   ! set of random phases (static, current and prevoious)
     real(kind=rp) :: rphs(trip_nmode_max,trip_nset_max,trip_nline_max)
   
   ! set of forcing arrays (static, current and prevoious)
     real(kind=rp), allocatable :: frcs(:,:,:)
   
   ! tripping array; interpolated value to set in 3D arrays
     real(kind=rp), allocatable :: ftrp(:,:)
    
     integer :: iff(trip_nline_max), iy(trip_nline_max)
     integer :: ir(97,trip_nline_max)
 contains
      procedure, pass(this) :: apply => trip_forcing
      procedure, pass(this) :: apply_device => trip_forcing_device
!!!      procedure, pass(this) :: apply => test_forcing
!!!      procedure, pass(this) :: apply_device => test_forcing_device
      procedure, pass(this) :: init => trip_init
      procedure, pass(this) :: update => trip_update
      procedure, pass(this) :: reset => trip_reset
      procedure, pass(this) :: ran2 => trip_ran2
   end type trip_t

!=======================================================================
!> @brief Register tripping module
!! @ingroup trip_line
!! @note This routine should be called in frame_usr_register
contains
!=======================================================================
!> @brief Initilise tripping module
!! @ingroup trip_line
!! @note This routine should be called in frame_usr_init
  subroutine trip_init(this, dof, nline, nmode, tiamp, tdamp, spos, epos, smth, lext, rota, tdt, time)
    class(trip_t) :: this
    integer, intent(in) :: nline, nmode(trip_nline_max)
    real(kind=rp), intent(inout), dimension(3,trip_nline_max) :: spos, epos, smth
    real(kind=rp), intent(inout) ::  rota(trip_nline_max), tdt(trip_nline_max)
    logical, intent(inout) :: lext(trip_nline_max)
    real(kind=rp), intent(in) :: time, tiamp, tdamp
    type(dofmap_t), target :: dof
    integer :: itmp
    real(kind=rp) :: rtmp
    logical :: ltmp

    integer :: il, jl


    this%dof => dof
    this%msh => dof%msh
    this%Xh => dof%Xh

    ! get runtime parameters
    this%nline = nline
    this%nmode = nmode
    this%tiamp = tiamp
    this%tdamp = tdamp
    this%spos = spos
    this%epos = epos
    this%smth = smth
    this%lext = lext
    this%rota = rota
    this%tdt = tdt
    allocate(this%prj(dof%size(),trip_nline_max))
    allocate(this%map(this%Xh%lz,dof%Xh%ly,dof%Xh%lz,dof%msh%nelv,trip_nline_max))
    allocate(this%fsmth(this%Xh%lz,this%Xh%ly,this%Xh%lz,dof%msh%nelv,trip_nline_max))
    allocate(this%frcs(dof%size(),trip_nset_max,trip_nline_max))
    allocate(this%ftrp(dof%size(),trip_nline_max))


    ! get sure z position of stating point is lower than ending point position
    do il=1,this%nline
       if (this%spos(ldim,il).gt.this%epos(ldim,il)) then
          do jl=1,LDIM
             rtmp = this%spos(jl,il)
             this%spos(jl,il) = this%epos(jl,il)
             this%epos(jl,il) = rtmp
          enddo
       endif
    enddo


    ! get inverse line lengths and smoothing radius
    do il=1,this%nline
       this%ilngt(il) = 0.0
       do jl=1,LDIM
          this%ilngt(il) = this%ilngt(il) + (this%epos(jl,il)-&
                           this%spos(jl,il))**2
       enddo
       if (this%ilngt(il).gt.0.0) then
          this%ilngt(il) = 1.0/sqrt(this%ilngt(il))
       else
          this%ilngt(il) = 1.0
       endif
       do jl=1,LDIM
          if (this%smth(jl,il).gt.0.0) then
             this%ismth(jl,il) = 1.0/this%smth(jl,il)
          else
             this%ismth(jl,il) = 1.0
          endif
       enddo
    enddo


    ! get 1D projection and array mapping
    call trip_1dprj(this)


    ! initialise random generator seed and number of time intervals
    do il=1,this%nline
       this%seed(il) = -32*il
       this%ntdt(il) = 1 - trip_nset_max
       this%ntdt_old(il) = this%ntdt(il)
       this%iff(il) = 0.0
    enddo


    ! generate random phases (time independent and time dependent)
    call trip_rphs_get(this, time)
!    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
!       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call trip_init_device(this)
    end if
    
    ! get forcing
    call trip_frcs_get(this, time, .true.)


    
    ! everything is initialised
    this%ifinit=.true.

  end subroutine
!=======================================================================
!> @brief Update tripping
!! @ingroup trip_line
  subroutine trip_update(this, time)
    class(trip_t), intent(inout) :: this
    ! local variables
    real(kind=rp) :: time

!-----------------------------------------------------------------------
    ! update random phases (time independent and time dependent)
    call trip_rphs_get(this,time)

    ! update forcing
    call trip_frcs_get(this,time,.false.)



  end subroutine trip_update 
!=======================================================================
!> @brief Compute tripping forcing
!! @ingroup trip_line
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    iel             local element number
  subroutine trip_forcing(this, ffx,ffy,ffz,ix,iy,iz,iel)
!!!    class(trip_t), intent(inout) :: this
    class(trip_t) :: this
    real(kind=rp), intent(inout) :: ffx, ffy, ffz
    integer, intent(in) :: ix,iy,iz,iel
    integer :: ipos,il
    real(kind=rp) :: ffn


    do il= 1, this%nline
       ffn = this%fsmth(ix,iy,iz,iel,il)
       if (ffn.gt.0.0) then

          ipos = this%map(ix,iy,iz,iel,il)
          ffn = this%ftrp(ipos,il)*ffn

          ffx =ffx - ffn*sin(this%rota(il))
          ffy =ffy  +ffn*cos(this%rota(il))
       endif
    enddo
  end subroutine trip_forcing

!=======================================================================
!> @brief Compute tripping forcing
!! @ingroup trip_line
!! @param[inout] ffx,ffy,ffz     forcing; x,y,z component
!! @param[in]    ix,iy,iz        GLL point index
!! @param[in]    iel             local element number
  subroutine trip_forcing_device( this, fx_d,fy_d,fz_d)
    class(trip_t) :: this
    type(c_ptr), intent(inout) :: fx_d, fy_d, fz_d


    call device_rzero(fx_d,this%dof%size())
    call device_rzero(fy_d,this%dof%size())
    call device_rzero(fz_d,this%dof%size())
    
    if (this%mask(0) .gt. 0) &
       call device_inhom_dirichlet_apply_vector(this%mask_d,fx_d,fy_d,fz_d,&
            this%ftripx_d,this%ftripy_d,this%ftripz_d,this%mask(0), glb_cmd_queue)
	

  end subroutine trip_forcing_device
!=======================================================================
!> @brief Reset tripping
!! @ingroup trip_line
  subroutine trip_reset(this, time)
    ! local variables
    class(trip_t) :: this
    real(kind=rp) :: time
      
    ! get 1D projection and array mapping
    call trip_1dprj(this)
    
    ! update forcing
    call trip_frcs_get(this,time,.true.)

  end subroutine trip_reset
!=======================================================================
!> @brief Get 1D projection, array mapping and forcing smoothing
!! @ingroup trip_line
!! @details This routine is just a simple version supporting only lines
!!   paralles to y axis. In future it can be generalised.
!! @remark This routine uses global scratch space \a CTMP0 and \a CTMP1
  subroutine trip_1dprj_y(this)
    class(trip_t) :: this
    real(kind=rp), allocatable :: lcoord(:)
    integer, allocatable :: lmap(:)
    integer :: npxy, npel, nptot, itmp, jtmp, ktmp, eltmp, istart
    integer :: il, jl, nx1
    real(kind=rp) :: xl, yl, zl, xr, yr, rota, rtmp, ptmp
    real(kind=rp), parameter :: epsl = 1.0d-10

    npxy = this%Xh%lxy
    npel = this%Xh%lxyz
    nx1 = this%Xh%lx
    nptot = this%dof%size()
    allocate(lcoord(nptot), lmap(nptot))

    ! for each line
    do il=1,this%nline
       ! reset mapping array
       !call ifill(this%map(1,1,1,1,il),-1,nptot)
       do jl = 1, nptot
          this%map(jl,1,1,1,il) = -1
       end do

       ! Get coordinates and sort them
       call copy(lcoord,this%dof%y,nptot)
       call sort(lcoord,lmap,nptot)
       ! get smoothing profile
       rota = this%rota(il)
       ! initialize smoothing factor
       call rzero(this%fsmth(1,1,1,1,il),nptot)
       this%npoint(il) = 0

       if (.not.this%lext(il) .and.&
          (lcoord(nptot).lt. (this%spos(2,il)-3.0*this%smth(2,il)) .or. &
          (lcoord(1).gt.&
          (this%epos(2,il)+3.0*this%smth(2,il))))) then
          exit
       end if

       ! if we do not extend a line exclude points below line start (y coordinate matters only)
       ! this cannot be mixed with Gauss profile
       istart = 1
       if (.not.this%lext(il)) then
          do jl=1,nptot
             if (lcoord(jl).lt. &
                 (this%spos(2,il)-3.0_rp*this%smth(2,il))) then
                istart = istart+1
             else
                exit
             endif
          enddo
       endif


       ! find unique entrances and provide mapping
       this%npoint(il) = 1
       this%prj(this%npoint(il),il) = lcoord(istart)
       itmp = lmap(istart)-1
       eltmp = itmp/npel + 1
       itmp = itmp - npel*(eltmp-1)
       ktmp = itmp/npxy + 1
       itmp = itmp - npxy*(ktmp-1)
       jtmp = itmp/nx1 + 1
       itmp = itmp - nx1*(jtmp-1) + 1
       this%map(itmp,jtmp,ktmp,eltmp,il) = this%npoint(il)
       do jl=istart+1,nptot
          ! if line is not extended finish at proper position
          if (.not.this%lext(il).and.(lcoord(jl).gt. &
              (this%epos(2,il)+3.0_rp*this%smth(2,il)))) exit

          if((lcoord(jl)-this%prj(this%npoint(il),il)).gt. &
              max(epsl,abs(epsl*lcoord(jl)))) then
             this%npoint(il) = this%npoint(il) + 1
             this%prj(this%npoint(il),il) = lcoord(jl)
          endif

          itmp = lmap(jl)-1
          eltmp = itmp/npel + 1
          itmp = itmp - npel*(eltmp-1)
          ktmp = itmp/npxy + 1
          itmp = itmp - npxy*(ktmp-1)
          jtmp = itmp/nx1 + 1
          itmp = itmp - nx1*(jtmp-1) + 1
          this%map(itmp,jtmp,ktmp,eltmp,il) = this%npoint(il)
       enddo

       ! rescale 1D array
       do jl=1,this%npoint(il)
          this%prj(jl,il) = (this%prj(jl,il) - this%spos(2,il))&
              *this%ilngt(il)
       enddo

       ! get smoothing profile
       rota = this%rota(il)
       ! initialize smoothing factor
       call rzero(this%fsmth(1,1,1,1,il),nptot)

       do jl=1,nptot
          itmp = jl-1
          eltmp = itmp/npel + 1
          itmp = itmp - npel*(eltmp-1)
          ktmp = itmp/npxy + 1
          itmp = itmp - npxy*(ktmp-1)
          jtmp = itmp/nx1 + 1
          itmp = itmp - nx1*(jtmp-1) + 1

          ! take only mapped points
          istart = this%map(itmp,jtmp,ktmp,eltmp,il)
          if (istart.gt.0) then

             ! rotation
             xl = this%dof%x(itmp,jtmp,ktmp,eltmp)-this%spos(1,il)
             yl = this%dof%z(itmp,jtmp,ktmp,eltmp)-this%spos(3,il)

             xr = xl*cos(rota)+yl*sin(rota)
             yr = -xl*sin(rota)+yl*cos(rota)

             rtmp = (xr*this%ismth(1,il))**2+(yr*this%ismth(3,il))**2
             ! do we extend a line beyond its ends
             if (.not.this%lext(il)) then
                if (this%prj(istart,il).lt.0.0_rp) then
                    zl = this%dof%y(itmp,jtmp,ktmp,eltmp)-this%spos(2,il)
                   rtmp = rtmp+(zl*this%ismth(2,il))**2
                elseif(this%prj(istart,il).gt.1.0_rp) then
                    zl = this%dof%y(itmp,jtmp,ktmp,eltmp)-this%epos(2,il)
                   rtmp = rtmp+(zl*this%ismth(2,il))**2
                endif
             endif
             ! Gauss; cannot be used with lines not extended beyond their ending points
             !trip_fsmth(itmp,jtmp,ktmp,eltmp,il) = exp(-4.0*rtmp)
             ! limited support
             if (rtmp.lt.1.0_rp) then
                this%fsmth(itmp,jtmp,ktmp,eltmp,il) = &
                    exp(-rtmp)*(1.0_rp-rtmp)**2.0_rp
             else
                this%fsmth(itmp,jtmp,ktmp,eltmp,il) = 0.0_rp
             endif
          endif

       enddo
    enddo
    deallocate(lcoord, lmap)
  end subroutine trip_1dprj_y

!=======================================================================
!> @brief Get 1D projection, array mapping and forcing smoothing
!! @ingroup trip_line
!! @details This routine is just a simple version supporting only lines
!!   paralles to z axis. In future it can be generalised.
!! @remark This routine uses global scratch space \a CTMP0 and \a CTMP1
  subroutine trip_1dprj(this)
    class(trip_t) :: this
    real(kind=rp), allocatable :: lcoord(:)
    integer, allocatable :: lmap(:)
    integer :: npxy, npel, nptot, itmp, jtmp, ktmp, eltmp, istart
    integer :: il, jl, nx1
    real(kind=rp) :: xl, yl, zl, xr, yr, rota, rtmp, ptmp
    real(kind=rp), parameter :: epsl = 1.0d-10

    npxy = this%Xh%lxy
    npel = this%Xh%lxyz
    nx1 = this%Xh%lx
    nptot = this%dof%size()
    allocate(lcoord(nptot), lmap(nptot))
    
    ! for each line
    do il=1,this%nline
       ! reset mapping array
       !call ifill(this%map(1,1,1,1,il),-1,nptot)
       do jl = 1, nptot
          this%map(jl,1,1,1,il) = -1
       end do
  
       ! Get coordinates and sort them
       call copy(lcoord,this%dof%z,nptot)
       call sort(lcoord,lmap,nptot)
       ! get smoothing profile
       rota = this%rota(il)
       ! initialize smoothing factor
       call rzero(this%fsmth(1,1,1,1,il),nptot)
       this%npoint(il) = 0
       if (.not.this%lext(il) .and.&
          (lcoord(nptot).lt. (this%spos(ldim,il)-3.0*this%smth(ldim,il)) .or. &
          (lcoord(1).gt.&
          (this%epos(ldim,il)+3.0*this%smth(ldim,il))))) then
          exit
       end if
  
       ! if we do not extend a line exclude points below line start (z coordinate matters only)
       ! this cannot be mixed with Gauss profile
       istart = 1
       if (.not.this%lext(il)) then
          do jl=1,nptot
             if (lcoord(jl).lt. &
                 (this%spos(ldim,il)-3.0_rp*this%smth(ldim,il))) then
                istart = istart+1
             else
                exit
             endif
          enddo
       endif

       
       ! find unique entrances and provide mapping
       this%npoint(il) = 1
       this%prj(this%npoint(il),il) = lcoord(istart)
       itmp = lmap(istart)-1
       eltmp = itmp/npel + 1
       itmp = itmp - npel*(eltmp-1)
       ktmp = itmp/npxy + 1
       itmp = itmp - npxy*(ktmp-1)
       jtmp = itmp/nx1 + 1
       itmp = itmp - nx1*(jtmp-1) + 1
       this%map(itmp,jtmp,ktmp,eltmp,il) = this%npoint(il)
       do jl=istart+1,nptot
          ! if line is not extended finish at proper position
          if (.not.this%lext(il).and.(lcoord(jl).gt. &
              (this%epos(ldim,il)+3.0_rp*this%smth(ldim,il)))) exit
  
          if((lcoord(jl)-this%prj(this%npoint(il),il)).gt. &
              max(epsl,abs(epsl*lcoord(jl)))) then
             this%npoint(il) = this%npoint(il) + 1
             this%prj(this%npoint(il),il) = lcoord(jl)
          endif
  
          itmp = lmap(jl)-1
          eltmp = itmp/npel + 1
          itmp = itmp - npel*(eltmp-1)
          ktmp = itmp/npxy + 1
          itmp = itmp - npxy*(ktmp-1)
          jtmp = itmp/nx1 + 1
          itmp = itmp - nx1*(jtmp-1) + 1
          this%map(itmp,jtmp,ktmp,eltmp,il) = this%npoint(il)
       enddo
           
       ! rescale 1D array
       do jl=1,this%npoint(il)
          this%prj(jl,il) = (this%prj(jl,il) - this%spos(ldim,il))&
              *this%ilngt(il)
       enddo
       
       ! get smoothing profile
       rota = this%rota(il)
       ! initialize smoothing factor
       call rzero(this%fsmth(1,1,1,1,il),nptot)
       
       do jl=1,nptot
          itmp = jl-1
          eltmp = itmp/npel + 1
          itmp = itmp - npel*(eltmp-1)
          ktmp = itmp/npxy + 1
          itmp = itmp - npxy*(ktmp-1)
          jtmp = itmp/nx1 + 1
          itmp = itmp - nx1*(jtmp-1) + 1
  
          ! take only mapped points
          istart = this%map(itmp,jtmp,ktmp,eltmp,il)
          if (istart.gt.0) then
  
             ! rotation
             xl = this%dof%x(itmp,jtmp,ktmp,eltmp)-this%spos(1,il)
             yl = this%dof%y(itmp,jtmp,ktmp,eltmp)-this%spos(2,il)
  
             xr = xl*cos(rota)+yl*sin(rota)
             yr = -xl*sin(rota)+yl*cos(rota)
  
             rtmp = (xr*this%ismth(1,il))**2+(yr*this%ismth(2,il))**2
             ! do we extend a line beyond its ends
             if (.not.this%lext(il)) then
                if (this%prj(istart,il).lt.0.0_rp) then
                    zl = this%dof%z(itmp,jtmp,ktmp,eltmp)-this%spos(ldim,il)
                   rtmp = rtmp+(zl*this%ismth(ldim,il))**2
                elseif(this%prj(istart,il).gt.1.0_rp) then
                    zl = this%dof%z(itmp,jtmp,ktmp,eltmp)-this%epos(ldim,il)
                   rtmp = rtmp+(zl*this%ismth(ldim,il))**2
                endif
             endif
             ! Gauss; cannot be used with lines not extended beyond their ending points
             !trip_fsmth(itmp,jtmp,ktmp,eltmp,il) = exp(-4.0*rtmp)
             ! limited support
             if (rtmp.lt.1.0_rp) then
                this%fsmth(itmp,jtmp,ktmp,eltmp,il) = &
                    exp(-rtmp)*(1.0_rp-rtmp)**2.0_rp
             else
                this%fsmth(itmp,jtmp,ktmp,eltmp,il) = 0.0_rp
             endif
          endif
  
       enddo
    enddo
    deallocate(lcoord, lmap)
  end subroutine trip_1dprj
  !=======================================================================
  !> @brief Generate set of random phases
  !! @ingroup trip_line
  subroutine trip_rphs_get(this, time)
    class(trip_t), intent(inout) :: this
    real(kind=rp), intent(in) :: time
    integer :: il, jl, kl
    integer :: itmp
  
!#ifdef DEBUG
!    character*3 str1, str2
!    integer iunit, ierr
!    ! call number
!    integer icalldl
!    save icalldl
!    data icalldl /0/
!#endif
    ! time independent part
    if (this%tiamp.gt.0.0.and..not.this%ifinit) then
       do il = 1, this%nline
          do jl=1, this%nmode(il)
             this%rphs(jl,1,il) = 2.0*pi*this%ran2(il)
          enddo
       enddo
    endif
  
    ! time dependent part
    do il = 1, this%nline
       itmp = int(time/this%tdt(il))
       !call bcast(itmp,ISIZE) ! just for safety
       do kl= this%ntdt(il)+1, itmp
          do jl= trip_nset_max,3,-1
             call copy(this%rphs(1,jl,il),this%rphs(1,jl-1,il), &
                  this%nmode(il))
          enddo
          do jl=1, this%nmode(il)
          this%rphs(jl,2,il) = 2.0_rp*pi*this%ran2(il)
          enddo
       enddo
       ! update time interval
       this%ntdt_old(il) = this%ntdt(il)
       this%ntdt(il) = itmp
    enddo

!#ifdef DEBUG
!    ! for testing
!    ! to output refinement
!    icalldl = icalldl+1
!    call io_file_freeid(iunit, ierr)
!    write(str1,'(i3.3)') NID
!    write(str2,'(i3.3)') icalldl
!    open(unit=iunit,file='trp_rps.txt'//str1//'i'//str2)
!  
!    do il=1,trip_nmode(1)
!       write(iunit,*) il,trip_rphs(il,1:4,1)
!    enddo
!  
!    close(iunit)
!#endif

  end subroutine trip_rphs_get
!=======================================================================
!> @brief A simple portable random number generator
!! @ingroup trip_line
!! @details  Requires 32-bit integer arithmetic. Taken from Numerical
!!   Recipes, William Press et al. Gives correlation free random
!!   numbers but does not have a very large dynamic range, i.e only
!!   generates 714025 different numbers. Set seed negative for
!!   initialization
!! @param[in]   il      line number
!! @return      ran
real(kind=rp) function trip_ran2(this,il)
    class(trip_t) :: this
    integer, intent(in) :: il
    ! local variables
    integer, parameter :: m=714025
    integer, parameter :: ia=1366
    integer, parameter :: ic=150889
    real, parameter :: rm=1./m
    integer :: j
    associate( seed => this%seed, iff => this%iff, iy => this%iy, ir => this%ir) 
    ! initialise
    if (seed(il).lt.0.or.iff(il).eq.0) then
       iff(il)=1
       seed(il)=mod(ic-seed(il),m)
       do j=1,97
          seed(il)=mod(ia*seed(il)+ic,m)
          ir(j,il)=seed(il)
       end do
       seed(il)=mod(ia*seed(il)+ic,m)
       iy(il)=seed(il)
    end if
    
    ! generate random number
    j=1+(97*iy(il))/m
    iy(il)=ir(j,il)
    trip_ran2=iy(il)*rm
    seed(il)=mod(ia*seed(il)+ic,m)
    ir(j,il)=seed(il)
    end associate

  end function
!=======================================================================
!> @brief Generate forcing along 1D line
!! @ingroup trip_line
!! @param[in] ifreset    reset flag
  subroutine trip_frcs_get(this, time, ifreset)
    ! argument list
    class(trip_t), intent(inout) :: this
    logical, intent(in) ::  ifreset
    real(kind=rp), intent(in) :: time
    integer :: il, jl, kl, ll
    integer :: istart, m
    real(kind=rp) :: theta0, theta
    logical :: ifntdt_dif
!#ifdef TRIP_PR_RST
!    ! variables necessary to reset pressure projection for P_n-P_n-2
!    integer nprv(2)
!    common /orthbi/ nprv
!
!    ! variables necessary to reset velocity projection for P_n-P_n-2
!    include 'VPROJ'
!#endif      
    ! local variables

!#ifdef DEBUG
!    character*3 str1, str2
!    integer iunit, ierr
!    ! call number
!    integer icalldl
!    save icalldl
!    data icalldl /0/
!#endif
    ! reset all
    if (ifreset) then
       if (this%tiamp.gt.0.0) then
          istart = 1
       else
          istart = 2
       endif
       do il= 1, this%nline
          do jl = istart, trip_nset_max
             call rzero(this%frcs(1,jl,il),this%npoint(il))
             do kl= 1, this%npoint(il)
                theta0 = 2*pi*this%prj(kl,il)
                do ll= 1, this%nmode(il)
                   theta = theta0*ll
                   this%frcs(kl,jl,il) = this%frcs(kl,jl,il) + &
                       sin(theta+this%rphs(ll,jl,il))
                enddo
             enddo
          enddo
       enddo
       ! rescale time independent part
       if (this%tiamp.gt.0.0) then
          do il= 1, this%nline
             call cmult(this%frcs(1,1,il),this%tiamp,this%npoint(il))
          enddo
       endif
       if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
          call trip_update_forces_device(this)
       end if

    else
       ! reset only time dependent part if needed
       ifntdt_dif = .FALSE.
       do il= 1, this%nline
          if (this%ntdt(il).ne.this%ntdt_old(il)) then
             ifntdt_dif = .TRUE.
             do jl= trip_nset_max,3,-1
                call copy(this%frcs(1,jl,il),this%frcs(1,jl-1,il), &
                    this%npoint(il))
             enddo
             call rzero(this%frcs(1,2,il),this%npoint(il))
             do jl= 1, this%npoint(il)
                theta0 = 2*pi*this%prj(jl,il)
                do kl= 1, this%nmode(il)
                   theta = theta0*kl
                   this%frcs(jl,2,il) = this%frcs(jl,2,il) + &
                       sin(theta+this%rphs(kl,2,il))
                enddo
             enddo
             if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
             .or. (NEKO_BCKND_OPENCL .eq. 1)) then
                call trip_update_forces_device(this)
             end if

          endif
       enddo
!         if (ifntdt_dif) then
!#ifdef TRIP_PR_RST
!            ! reset projection space
!            ! pressure
!            if (int(PARAM(95)).gt.0) then
!               PARAM(95) = ISTEP
!               nprv(1) = 0      ! veloctiy field only
!            endif
!            ! velocity
!            if (int(PARAM(94)).gt.0) then
!               PARAM(94) = ISTEP!+2
!               ivproj(2,1) = 0
!               ivproj(2,2) = 0
!               if (IF3D) ivproj(2,3) = 0
!            endif
!#endif
!         endif
    endif
      
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
    .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       ! get tripping for current time stepa
       m = this%mask(0)
       if ( m.gt. 0) then
          if (this%tiamp.gt.0.0) then
               call device_copy(this%ftripx_d,this%f_interpolate_d(1),m)
          else
               call device_rzero(this%ftripx_d,m)
          endif
          !> We only support one line for now!
          il = 1
          ! interpolation in time
          theta0= time/this%tdt(il)-real(this%ntdt(il))
          if (theta0.gt.0.0) then
             theta0=theta0*theta0*(3.0-2.0*theta0)
             theta = (1.0-theta0)*this%tdamp
             call device_add2s2(this%ftripx_d,this%f_interpolate_d(3),theta,m)
             theta = theta0*this%tdamp
             call device_add2s2(this%ftripx_d,this%f_interpolate_d(2),theta,m)
          else
             theta0=theta0+1.0
             theta0=theta0*theta0*(3.0-2.0*theta0)
             theta = (1.0-theta0)*this%tdamp
             call device_add2s2(this%ftripx_d,this%f_interpolate_d(4),theta,m)
             theta = theta0*this%tdamp
             call device_add2s2(this%ftripx_d,this%f_interpolate_d(3),theta,m)
          endif
          call device_col2(this%ftripx_d,this%fsmth_mask_d,m)
          call device_cmult2(this%ftripy_d,this%ftripx_d,cos(this%rota(1)),m)
          call device_cmult(this%ftripx_d,-sin(this%rota(1)),m)
          call device_rzero(this%ftripz_d,m)
       end if
     else
       ! get tripping for current time step
       if (this%tiamp.gt.0.0) then
          do il= 1, this%nline
            call copy(this%ftrp(1,il),this%frcs(1,1,il),this%npoint(il))
          enddo
       else
          do il= 1, this%nline
             call rzero(this%ftrp(1,il),this%npoint(il))
          enddo
       endif
       ! interpolation in time
       do il = 1, this%nline
          theta0= time/this%tdt(il)-real(this%ntdt(il))
          if (theta0.gt.0.0) then
             theta0=theta0*theta0*(3.0-2.0*theta0)
             !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
             do jl= 1, this%npoint(il)
                this%ftrp(jl,il) = this%ftrp(jl,il) + &
                    this%tdamp*((1.0-theta0)*this%frcs(jl,3,il) + &
                    theta0*this%frcs(jl,2,il))
             enddo
          else
             theta0=theta0+1.0
             theta0=theta0*theta0*(3.0-2.0*theta0)
             !theta0=theta0*theta0*theta0*(10.0+(6.0*theta0-15.0)*theta0)
             do jl= 1, this%npoint(il)
                this%ftrp(jl,il) = this%ftrp(jl,il) + &
                     this%tdamp*((1.0-theta0)*this%frcs(jl,4,il) + &
                     theta0*this%frcs(jl,3,il))
             enddo
          endif
       enddo
    end if

!#efdef DEBUG
!      ! for testing
!      ! to output refinement
!      icalldl = icalldl+1
!      call io_file_freeid(iunit, ierr)
!      write(str1,'(i3.3)') NID
!      write(str2,'(i3.3)') icalldl
!      open(unit=iunit,file='trp_fcr.txt'//str1//'i'//str2)
!
!      do il=1,trip_npoint(1)
!         write(iunit,*) il,trip_prj(il,1),trip_ftrp(il,1),
!     $        trip_frcs(il,1:4,1)
!      enddo
!
!      close(iunit)
!#endif
  end subroutine
  subroutine trip_init_device(this)
    type(trip_t) :: this
    integer, allocatable :: mask_temp(:)
    real(kind=rp), allocatable :: fsmth_mask_temp(:)
    integer :: i, j,n, m = 0
    integer(c_size_t) :: array_size
    n = this%dof%size()
    allocate(mask_temp(n), fsmth_mask_temp(n))
    
    do i = 1, n
       if (this%fsmth(i,1,1,1,1) .gt. 0.0) then
          m = m + 1
          mask_temp(m) = i
          fsmth_mask_temp(m) = this%fsmth(i,1,1,1,1)
       end if
    end do

    allocate(this%mask(0:m))
    allocate(this%fsmth_mask(m))
    this%mask(0) = m
    print *, m
    do i = 1, m
       this%mask(i) = mask_temp(i)
       this%fsmth_mask(i) = fsmth_mask_temp(i)
    end do

    deallocate(mask_temp, fsmth_mask_temp)
    call device_map(this%mask, this%mask_d, m+1)
    call device_map(this%fsmth_mask, this%fsmth_mask_d, m)
    call device_memcpy(this%mask, this%mask_d, m+1,HOST_TO_DEVICE,sync=.true.)
    call device_memcpy(this%fsmth_mask, this%fsmth_mask_d, m,HOST_TO_DEVICE,sync=.true.)
    array_size = rp*m
    do i = 1, trip_nset_max
       call device_alloc(this%f_interpolate_d(i),array_size)
       if (m .gt. 0) call device_rzero(this%f_interpolate_d(i),m)
    end do
    call device_alloc(this%ftripx_d,array_size)
    call device_alloc(this%ftripy_d,array_size)
    call device_alloc(this%ftripz_d,array_size)

  end subroutine trip_init_device


  subroutine trip_update_forces_device(this)
    type(trip_t) :: this
    real(kind=rp), allocatable :: frc_mask_temp(:,:)
    integer :: m, ipos, il, i, j
    m = this%mask(0)
    il = 1
    allocate(frc_mask_temp(m,trip_nset_max))
    call rzero(frc_mask_temp,m*trip_nset_max)
    do i = 1, m
       ipos = this%map(this%mask(i),1,1,1,1) 
       do j = 1, trip_nset_max
          frc_mask_temp(i, j) = this%frcs(ipos,j,il)
       end do
    end do
    do j = 1, trip_nset_max
       call device_memcpy(frc_mask_temp(:,j),this%f_interpolate_d(j),m,HOST_TO_DEVICE, sync=.true.)
    end do

  end subroutine trip_update_forces_device


end module trip
