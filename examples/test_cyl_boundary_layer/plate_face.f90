!>User file to set up at rotating cylinder. Martin Karp 30/5-2022

module user
  use neko
  implicit none
  ! Case parameters
  real(kind=rp), parameter :: h = 1.0
  real(kind=rp), parameter :: gam = 20.0
  real(kind=rp), parameter :: rad = h/gam
  real(kind=rp), parameter :: n = 7
  real(kind=rp), parameter :: pw = 1/n
  real(kind=rp), parameter :: ucl = 1
  real(kind=rp), parameter :: alpha = 1
  real(kind=rp), parameter :: u_th2 = ucl*alpha
  real(kind=rp), parameter :: u_rho = 0.0
  real(kind=rp), parameter :: u_axial = 0.0
  real(kind=rp), parameter :: y0 = 0.0
  real(kind=rp), parameter :: y1 = 0.0
  real(kind=rp), parameter :: delta = 0.005*h

!  real(kind=rp), allocatable, dimension(:,:) :: yz_ele_plane
!  real(kind=rp), allocatable, dimension(:) :: yz_map
  integer :: tcounter=0  ! Temp. counter to substitute tstep in id_face 
contains
  
  
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
!    u%fluid_user_ic => user_ic
    u%fluid_user_if => user_inflow_eval
!    !u%user_check => user_do_stuff
!    u%user_mesh_setup => cylinder_deform
    u%user_mesh_setup => check_face
  end subroutine user_setup
 
  subroutine cylinder_deform(msh)
    type(mesh_t), intent(inout) :: msh
    msh%apply_deform => cylinder_gen_curve
  end subroutine cylinder_deform

  subroutine check_face(msh)
    type(mesh_t), intent(inout) :: msh
    msh%apply_deform => id_face
  end subroutine check_face

  !> Create a map between elements on the inflow velocity plane(s) and the inlet elements in Neko
  subroutine id_face(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    type(tuple_i4_t) :: el_and_facet
    real(kind=rp) :: th
    integer :: e, i, j ,k, l,  facet
    integer :: nely, nelz, yz_ele
!    real(kind=rp) :: yz_ele_plane(2,100), yz_map(100)
    real(kind=rp), allocatable :: yz_ele_plane(:,:)
    real(kind=rp), allocatable  :: yz_map(:) 

    nely = 10 ! Total no. elements in y
    nelz = 10 ! Total no. elements in z
    !if (tstep==0) then 
    ! NOTE: Time stepper is only in simulation.f90 so cannot use this check now. Later I should move this out of case.f90
    ! and perform this mapping only at tstep==0
    ! Also allocating the below once like this using tcounter is wrong if I use multiple ranks. I should use tstep==0 for this as well
    ! For now using tcounter==0
    ! I defined tcounter inside the module but outside all these subrountines, so now this works with multi-cpus 
    if (tcounter==0) then
     allocate(yz_ele_plane(2,nely*nelz)) ! y,z gll coordinates of lower left corner of all element
     allocate(yz_map(1:nely*nelz))
     open(unit=10,file='vel_plane_ele',status='old',form='unformatted',access="stream")
     read(10) yz_ele_plane ! size(2,nely*nelz). (1, nely*nelz) is y coord and (2,nely*nelz) is z coord
     close(10)
     ! allocate(yz_map(1:nely*nelz))
     yz_map(1:nely*nelz) = 0.0 ! Shows which glbal element no. corr. to each element in inflow vel. plane
     ! Create map to find global element nos. that corr. to elements on vel. plane
     ! Check on all cpus, and bcast to all, and add face element global number in an array ordered correctly as the vel. plane.
     ! Later we can use findloc to get the index on vel slice corr. to a given global ele. no.
     do yz_ele=1,nely*nelz ! Loop over elements in vel. plane
      do l = 1,msh%labeled_zones(1)%size ! Loop over elements in Neko inflow plane 
       el_and_facet = msh%labeled_zones(1)%facet_el(l)
       facet = el_and_facet%x(1)
       e = el_and_facet%x(2)
       if (index_is_on_facet(1,1,1,lx,ly,lz, facet)) then
!        if ((y(1,1,1,e).eq.yz_ele_plane(1,yz_ele)).and.(z(1,1,1,e).eq.yz_ele_plane(2,yz_ele))) then
        if ((abs(y(1,1,1,e)-yz_ele_plane(1,yz_ele))<1e-8).and.(abs(z(1,1,1,e)-yz_ele_plane(2,yz_ele))<1e-8)) then
            write(*,001) 'Plane:',yz_ele,yz_ele_plane(1,yz_ele),yz_ele_plane(2,yz_ele),&
                & ' Neko:',e+msh%offset_el,y(1,1,1,e),z(1,1,1,e)
            yz_map(yz_ele) = e+msh%offset_el ! Assigning the global element number in Neko inlet corresponding to this element in the velocity plane
        end if      
       end if
      end do  
     end do

    !The inlet zone number is 1
    do l = 1,msh%labeled_zones(1)%size
       el_and_facet = msh%labeled_zones(1)%facet_el(l)
       facet = el_and_facet%x(1)
       e = el_and_facet%x(2)
       do k = 1,1!lz
          do j = 1,1!ly
              do i = 1, 1!lx
                 if (index_is_on_facet(i,j,k,lx,ly,lz, facet)) then
                    write(*,*) 'x=',x(i,j,k,e),', y=',y(i,j,k,e),', z=',z(i,j,k,e),', e=',e,', offset=',msh%offset_el, &
                        & 'gl=',e+msh%offset_el,',pid=', pe_rank
          !          th = atan2(z(i,j,k,e), x(i,j,k,e))
          !          x(i,j,k,e) = rad * cos(th)
          !          z(i,j,k,e) = rad * sin(th) 
                 end if
              end do
          end do
       end do
    end do
   end if !
   tcounter = 1
001 format(a6,i4,f10.4,f10.4,a5,i6,f10.4,f10.4)    
  end subroutine id_face

  !> Make the cylinder more circular!
  subroutine cylinder_gen_curve(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    type(tuple_i4_t) :: el_and_facet
    real(kind=rp) :: th
    integer :: e, i, j ,k, l,  facet

    !The cylinders zone number is 7
    do l = 1,msh%labeled_zones(7)%size
       el_and_facet = msh%labeled_zones(7)%facet_el(l)
       facet = el_and_facet%x(1)
       e = el_and_facet%x(2)
       do k = 1, lz
          do j = 1, ly
              do i = 1, lx
                 if (index_is_on_facet(i,j,k,lx,ly,lz, facet)) then
                    th = atan2(z(i,j,k,e), x(i,j,k,e))
                    x(i,j,k,e) = rad * cos(th)
                    z(i,j,k,e) = rad * sin(th) 
                 end if
              end do
          end do
       end do
    end do
  end subroutine cylinder_gen_curve
  
  subroutine user_inflow_eval(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
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
    real(kind=rp) ::  u_th,dist,th, yy
    real(kind=rp) ::  arg

!   Two different regions (inflow & cyl) have the label 'v  '
!   Let compute the distance from the (0,0) in the x-y plane
!   to identify the proper one
    dist = sqrt(x**2 + z**2)

! --- INFLOW
    if (dist .gt. 1.1*rad) then
       u =  ucl*y**pw
    end if
! --- 

    w = 0.0
    v = 0.0
! --- SPINNING CYLINDER

    if (dist.lt.1.5*rad .and. y.gt. 0.1) then                      
       th = atan2(z,x)
       u = cos(th)*u_rho - sin(th)*u_th2
       w = sin(th)*u_rho + cos(th)*u_th2   
    end if     
                    
! --- 


!     Smoothing function for the velocity u_th on the spinning cylinder
!     to avoid gap in the at the bottom wall

!     u_th is smoothed if z0 < z < delta
!     u_th=1 if z >= delta
 

    yy = y + abs(y0) ! coordinate shift 

    if (dist .lt. 1.5*rad) then 
       if (yy.lt.delta) then
          arg  = yy/delta
          u_th = u_th2/(1.0_rp+exp(1.0_rp/(arg-1.0_rp)+1.0_rp/arg))
       else
          u_th = u_th2
       endif

       th = atan2(z,x)

       u = cos(th)*u_rho - sin(th)*u_th
       w = sin(th)*u_rho + cos(th)*u_th  
    end if
  end subroutine user_inflow_eval

  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    !type(param_t), intent(inout) :: params
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: y

    do i = 1, u%dof%size()
       y = u%dof%y(i,1,1,1)
       u%x(i,1,1,1) =  ucl*y**pw
       v%x(i,1,1,1) = 0.0
       w%x(i,1,1,1) = 0.0
    end do
  end subroutine user_ic

end module user
