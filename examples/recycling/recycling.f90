! Martin Karp 13/3-2023
! updated initial condition Philipp Schlatter 09/07/2024
module user
  use neko
  implicit none

  type(global_interpolation_t) :: interpolate
  type(matrix_t) :: xyz
  type(vector_t) :: res, B
  real(kind=rp) :: vol
  integer :: n_pts = 0
  logical :: init = .false.

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%initial_conditions => user_ic
    u%mesh_setup => user_mesh_scale
    u%dirichlet_conditions => dirichlet_update
  end subroutine user_setup

  !Observe this might be called multiple times per time step by
  !different solvers
  subroutine dirichlet_update(field_bc_list, bc, time)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    integer :: i, pt, msk_idx
    real(kind=rp) :: y, z, scale
    type(coef_t), pointer :: coef
    type(field_t), pointer :: field

    coef => neko_user_access%case%fluid%c_Xh

    ! Only do this at the first time step since our BCs are constants.
    if (time%tstep .eq. 1 .and. .not. init) then
       n_pts = bc%msk(0)

       call xyz%init(3,n_pts)
       call B%init(n_pts)
       pt = 0
       do i = 1, n_pts
          ! Get idx of point and its coords and store in contigous xyz array
          msk_idx = bc%msk(i)
          ! We want to recycle the flow from 10 units upstream
          xyz%x(1,i) = coef%dof%x(msk_idx,1,1,1) + 10.0
          xyz%x(2,i) = coef%dof%y(msk_idx,1,1,1)
          xyz%x(3,i) = coef%dof%z(msk_idx,1,1,1)
          B%x(i) = coef%B(msk_idx,1,1,1)
       end do
       vol = glsum(B%x,n_pts)
       call xyz%copy_from(HOST_TO_DEVICE,.false.)
       call B%copy_from(HOST_TO_DEVICE,.false.)
       !Initialize interpolator
       call interpolate%init(coef%dof)
       ! Find the the points we want to interpolate on the inflow
       call interpolate%find_points(xyz%x, n_pts)
       ! Initialize a temporary array
       call res%init(n_pts)
       init = .true.
    end if

    ! Check that we are being called by `fluid`
    if (field_bc_list%items(1)%ptr%name .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, &
            w => field_bc_list%items(3)%ptr)

         field => neko_field_registry%get_field('u')
         ! get the x-velocity (u) values 10 units upstream into res
         call interpolate%evaluate(res%x, field%x, .false.)
         !Enforce that bulk velocity is 1
         !If we run on GPU, do this on the GPU directly
         !Moving data to and from the GPU kills performance
         if (NEKO_BCKND_DEVICE .eq. 1) then
            scale = 1.0_rp/(device_glsc2(res%x_d,B%x_d,n_pts)/vol)
            call device_cmult(res%x_d,scale, n_pts)
         else
            scale = 1.0_rp/(glsc2(res%x,B%x,n_pts)/vol)
            call cmult(res%x,scale, n_pts)
         end if
         ! scatter the values in res into the correct spots in u
         call field_masked_scatter_copy_0(u, res%x, bc%msk, u%size(), n_pts)
         ! repeat for v and w
         field => neko_field_registry%get_field('v')
         call interpolate%evaluate(res%x, field%x, .false.)
         call field_masked_scatter_copy_0(v, res%x, bc%msk, u%size(), n_pts)

         field => neko_field_registry%get_field('w')
         call interpolate%evaluate(res%x, field%x, .false.)
         call field_masked_scatter_copy_0(w, res%x, bc%msk, u%size(), n_pts)

       end associate
    end if

  end subroutine dirichlet_update


  ! Rescale mesh, we create a mesh with some refinement close to the wall.
  ! initial mesh: 0..8, -1..1, 0..1.5
  ! mesh size (8*pi,2*delta,4/3*pi)
  ! New mesh can easily be genreated with genmeshbox
  ! OBS refinement is not smooth and the constant values are a bit ad hoc.
  ! OBS only moves element vertices
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i, p, nvert

    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h
    real(kind=rp) :: center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y
    real(kind=rp) :: llx, llz

    !constants to scale mesh size
    llx = 4.*pi
    llz = 4./3.*pi

    ! rescale mesh
    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = llx/4.*msh%points(i)%x(1)
       y = msh%points(i)%x(2)
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)- &
               el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = 2./3.*llz*msh%points(i)%x(3)
    end do

  end subroutine user_mesh_scale

  ! User defined initial condition
  subroutine user_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(dofmap_t), pointer :: dof
    type (field_t), pointer :: u, v, w, p
    real(kind=rp) :: uvw(3)
    integer :: i

    dof => fields%dof(1)
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    do i = 1, u%dof%size()
       uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic

  ! Kind of brute force with rather large initial disturbances
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran

    real(kind=rp) :: llx, llz

    llx = 4.*pi
    llz = 4./3.*pi

    Re_tau = 180
    C = 5.17
    k = 0.41
    Re_b = 2800

    yp = (1-y)*Re_tau
    if (y .lt. 0) yp = (1+y)*Re_tau

    ! Reichardt function
    ux = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
         (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
    ux = ux * Re_tau/Re_b

    ! actually, sometimes one may not use the turbulent profile, but
    ! rather the parabolic lamianr one
    ! ux = 1.5*(1-y**2)

    ! add perturbations to trigger turbulence
    ! base flow
    uvw(1) = ux
    uvw(2) = 0
    uvw(3) = 0

    ! first, large scale perturbation
    eps = 0.05
    kx = 3
    kz = 4
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! second, small scale perturbation
    eps = 0.005
    kx = 17
    kz = 13
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! finally, random perturbations only in y
    eps1 = 0.001
    ran = sin(-20*x*z+y**3*tan(x*z**2)+100*z*y-20*sin(x*y*z)**5)
    uvw(2) = uvw(2) + eps1*ran

  end function channel_ic

end module user
