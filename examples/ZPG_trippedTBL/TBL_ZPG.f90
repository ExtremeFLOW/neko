module user
  use neko
  use trip
  implicit none

  type(trip_t) :: tripper
  logical :: enable_tripping

  !-----------------------------------------------------------------------------------------!
  ! trip stuff
  real(kind=rp), parameter :: Lz = 78.0_rp
  real(kind=rp), parameter :: delta_star = 0.325_rp ! displacement thickness of our inflow BL
                                                    ! this is based on the inflow d99 of 1 and d*/d99 ratio of around 0.325 for Blasius

  ! amplitude of tripping force
  real(kind=rp), parameter :: TIAMP   = 0.00_rp         ! time independent amplitude. set to zero
  real(kind=rp), parameter :: TDAMP   = 0.05_rp         ! amplitude of time dependent term. 
                                                        ! there is no guideline for this, but this was chosen for this specific example after trial and error
  ! location of the trip
  real(kind=rp), parameter :: SPOSX01 = 10.0_rp         ! 10*d99 downstream of inflow
  real(kind=rp), parameter :: SPOSY01 = 0.0_rp          ! wall
  real(kind=rp), parameter :: SPOSZ01 = 0.0_rp          ! z=0
  real(kind=rp), parameter :: EPOSX01 = SPOSX01         
  real(kind=rp), parameter :: EPOSY01 = SPOSY01          
  real(kind=rp), parameter :: EPOSZ01 = Lz         ! z_max=78*d99

  ! smoothing and spread of the tripping force
  ! Parameters SMTH*01 taken from Schlatter & Orlu 2012 doi:10.1017/jfm.2012.324
  real(kind=rp), parameter :: SMTHX01 = 4.0_rp * delta_star
  real(kind=rp), parameter :: SMTHY01 = delta_star      
  real(kind=rp), parameter :: SMTHZ01 = 0.0_rp          ! not needed, so set to zero

  ! rotation
  real(kind=rp), parameter :: ROTA01  = 0.0_rp          ! wall normal, so 0

  ! number of Fourier modes
  integer, parameter :: NMODE01 = Lz / (1.7_rp*delta_star)              ! Lz/(1.7* d*)
  
  ! time scale of the variations
  real(kind=rp), parameter :: TDT01   = 4.0_rp * delta_star / 1.0_rp    ! 4 * d* / Uinf

  ! additional parameters
  integer, parameter :: nline = 1               ! number of tripping lines
  logical, parameter :: LEXT01 = .false.        ! whether to extend the line
  !-----------------------------------------------------------------------------------------!

contains

!=============================================================================!
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => user_dirichlet_conditions
    user%source_term => userf
    user%initialize => initialize
    user%finalize => finalize
    user%compute => compute
    user%mesh_setup => user_mesh_scale
  end subroutine user_setup

!=============================================================================!
  ! User-defined finalization routine called at the end of the simulation
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time


  end subroutine finalize

 
!=============================================================================!
  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    integer :: nmode(trip_nline_max), ierr
    real(kind=rp), dimension(3,trip_nline_max) :: spos, epos, smth
    real(kind=rp) :: rota(trip_nline_max), tdt(trip_nline_max)
    logical :: lext(trip_nline_max)


    type(field_t), pointer :: u
    real(kind=rp) :: t
    integer :: ntot

    !
    ! Tripping
    !
    spos(1,1) = SPOSX01
    spos(2,1) = SPOSY01
    spos(3,1) = SPOSZ01
    epos(1,1)= EPOSX01
    epos(2,1)= EPOSY01
    epos(3,1)= EPOSZ01
    smth(1,1)= SMTHX01
    smth(2,1)= SMTHY01
    smth(3,1)= SMTHZ01
    rota(1) = ROTA01
    nmode(1) = NMODE01
    tdt(1) = TDT01
    lext(1) = LEXT01

    u => neko_registry%get_field("u")
    t = time%t

    ! initializes the tripping parameters
    ! this includes finding a map of the GLL points to the tripping line, etc.
    call tripper%init( u%dof , nline, nmode, tiamp, tdamp, &
                       spos, epos, smth, lext, rota, tdt, t)

  end subroutine initialize

!=============================================================================!
   subroutine userf(scheme_name, rhs, time)
     character(len=*), intent(in) :: scheme_name
     type(field_list_t), intent(inout) :: rhs
     type(time_state_t), intent(in) :: time

     type(field_t), pointer :: rhs_u, rhs_v, rhs_w
     real(kind=rp) :: y
     integer :: ix, iy, iz, iel
     real(kind=rp) :: ffx,ffy,ffz

     if (scheme_name .eq. 'fluid') then

       rhs_u => rhs%get_by_index(1)
       rhs_v => rhs%get_by_index(2)
       rhs_w => rhs%get_by_index(3)

       call field_rzero(rhs_u,tripper%dof%size())
       call field_rzero(rhs_v,tripper%dof%size())
       call field_rzero(rhs_w,tripper%dof%size())


       ! a hacky implementation
       if (NEKO_BCKND_DEVICE .eq. 1) then
           ! tripper%apply_device uses c pointers
           call tripper%apply_device(rhs_u%x_d, rhs_v%x_d, rhs_w%x_d)
       else
          ! tripline_org only sets the temp values to zero and calls tripper%apply
          ! tripper%apply uses arrays
          do iel = 1, tripper%msh%nelv
             do iz = 1, tripper%Xh%lz
                do iy = 1, tripper%Xh%ly
                   do ix = 1, tripper%Xh%lx
                      call tripline_org(ffx,ffy,ffz,ix,iy,iz,iel)
                      rhs_u%x(ix,iy,iz,iel) = ffx
                      rhs_v%x(ix,iy,iz,iel) = ffy
                      rhs_w%x(ix,iy,iz,iel) = ffz
                   end do
                end do
             end do
          end do
       end if

     end if
 
   end subroutine userf

!=============================================================================!
   !> Tripline forcing
   ! subroutine that only sets the values to zero and calls tripper%apply
   subroutine tripline_org(u, v, w, j, k, l, e )
     real(kind=rp), intent(inout) :: u
     real(kind=rp), intent(inout) :: v
     real(kind=rp), intent(inout) :: w
     integer, intent(in) :: j
     integer, intent(in) :: k
     integer, intent(in) :: l
     integer, intent(in) :: e

     u = 0.0_rp
     v = 0.0_rp
     w = 0.0_rp

     call tripper%apply(u,v,w,j,k,l,e)

   end subroutine tripline_org
!=============================================================================!
  ! New mesh can easily be genreated with genmeshbox
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time

    integer :: i, p, nvert
    real(kind=rp) :: y
    integer :: el_in_y
    integer :: ymax
    real(kind=rp) :: y_shift
    real(kind=rp) :: gamma
    
    
    ! mesh for Schlatter and Orlu, JFM, 2010 can be generated using: 
    !                   genmeshbox 0.0 1950.0 0.0 65.0 0.0 78.0 1152 64 112 .false. .false. .true.
    ! need to be checked
    ymax = 65.0_rp 
    el_in_y = 64

    ! manually adjusted parameters for the stretching function
    y_shift = -0.9_rp   ! -1 means half channel. can use -.5 < values < -1 if want finer grid on top boundary. gamma=-.5 becomes channel.
    gamma = 3.0_rp
  
    ! tanh stretching 
    nvert = size(msh%points)
    do i = 1, nvert
       y = msh%points(i)%x(2) 
       y = y / ymax + y_shift
       y = ymax * &
        (tanh(gamma*y)-tanh(gamma*y_shift)) / &
        (tanh(gamma*(1+y_shift))-tanh(gamma*y_shift))
       msh%points(i)%x(2) = y
    end do

  end subroutine user_mesh_scale

!=============================================================================!
! this makes a "fake" sine-shaped profile for the inflow as a proxy for the Blasius profile
   subroutine user_dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(dofmap_t ), pointer :: dof
    type(field_t), pointer :: u, v, w
    real(kind=rp) :: x, y, z
    integer :: i, msk_ind

    real(kind=rp) ::  fluid_nu,dist,th, yy
    real(kind=rp) ::  Uinf,delta99

    ! Only do this at the first time step since our BCs are constants.
    if (time%tstep .ne. 1) return

    ! Grab the dofmap from the first field in the list.
    dof => fields%dof(1)

    ! We only have the velocity solver here, so we know the contents of the
    ! field list, i.e. that it holds u, v, w.
    u => fields%get("u")
    v => fields%get("v")
    w => fields%get("w")

    Uinf = 1.0_rp

    do i = 1, bc%msk(0)
      msk_ind = bc%msk(i)
      x = dof%x(msk_ind, 1, 1, 1)
      y = dof%y(msk_ind, 1, 1, 1)
      z = dof%z(msk_ind, 1, 1, 1)

      u%x(msk_ind,1,1,1) = Uinf
      v%x(msk_ind,1,1,1) = 0.0_rp
      w%x(msk_ind,1,1,1) = 0.0_rp
      if (y.lt.1.0_rp) then
        u%x(msk_ind,1,1,1) = sin(y*PI/2)
      end if
    enddo

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.false.)
    end if

   end subroutine user_dirichlet_conditions

!=============================================================================!
! This is an extension of the sine-shaped profile inside the domain
  ! User-defined initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(dofmap_t), pointer :: dof
    type (field_t), pointer :: u, v, w

    integer :: i, ntot 
    real(kind=rp) :: uvw(3)

    dof => fields%dof(1)
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")

    ntot = dof%size()
    do i = 1, ntot
       uvw = blasius_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do

  end subroutine initial_conditions

!=============================================================================!
  ! application of a sine-shaped profile as an approximation for the Blasius profile
  function blasius_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) ::  fluid_nu,dist,th, yy
    real(kind=rp) ::  Uinf,delta99

    fluid_nu = 0.00070
    Uinf = 1.0_rp

    dist = (1.0/5.29)**2 * Uinf / fluid_nu          ! find the offset from the inlet
    delta99 = 5.29 * sqrt(fluid_nu*(x+dist)/Uinf)   ! calculate local delta99

    uvw(1) = Uinf
    if (y.lt.delta99) then
      uvw(1) = sin(y/delta99*PI/2)
    end if

    uvw(2) = 0.0_rp     ! v is not zero in reality, but for simplicity
    uvw(3) = 0.0_rp

  end function blasius_ic


!=============================================================================!
  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    ! updates the time dependent value of the trip
    ! IMPORTANT: if you restart from a field file (instead of checkpoint), 
    !            you will need to use the correct time here for tripping.
    call tripper%update( time%t )

  end subroutine compute
!=============================================================================!

end module user
