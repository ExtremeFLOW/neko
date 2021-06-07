!> Coefficients 
module coefs
  use gather_scatter
  use num_types
  use space  
  use math
  use mesh
  use mxm_wrapper
  implicit none
  private
  
  !> Coefficients defined on a given (mesh, \f$ X_h \f$) tuple
  type, public :: coef_t     
     real(kind=rp), allocatable :: G11(:,:,:,:) !< Geometric data at index 1,1
     real(kind=rp), allocatable :: G22(:,:,:,:) !< Geometric data at index 2,2
     real(kind=rp), allocatable :: G33(:,:,:,:) !< Geometric data at index 3,3
     real(kind=rp), allocatable :: G12(:,:,:,:) !< Geometric data at index 1,2
     real(kind=rp), allocatable :: G13(:,:,:,:) !< Geometric data at index 1,3
     real(kind=rp), allocatable :: G23(:,:,:,:) !< Geometric data at index 2,3

     real(kind=rp), allocatable :: mult(:,:,:,:) !< Multiplicity
     ! generate mapping data between element and reference element 
     !! \f$ dx/dr, dy/dr, dz/dr \f$
     !! \f$ dx/ds, dy/ds, dz/ds \f$
     !! \f$ dx/dt, dy/dt, dz/dt \f$
     real(kind=rp), allocatable :: dxdr(:,:,:,:), dydr(:,:,:,:), dzdr(:,:,:,:) 
     real(kind=rp), allocatable :: dxds(:,:,:,:), dyds(:,:,:,:), dzds(:,:,:,:)
     real(kind=rp), allocatable :: dxdt(:,:,:,:), dydt(:,:,:,:), dzdt(:,:,:,:) 
     !< \f$ dr/dx, dr/dy, dr/dz \f$
     !! \f$ ds/dx, ds/dy, ds/dz \f$
     !! \f$ dt/dx, dt/dy, dt/dz \f$
     real(kind=rp), allocatable :: drdx(:,:,:,:), drdy(:,:,:,:), drdz(:,:,:,:) 
     real(kind=rp), allocatable :: dsdx(:,:,:,:), dsdy(:,:,:,:), dsdz(:,:,:,:)
     real(kind=rp), allocatable :: dtdx(:,:,:,:), dtdy(:,:,:,:), dtdz(:,:,:,:) 
     
     real(kind=rp), allocatable :: h1(:,:,:,:) 
     real(kind=rp), allocatable :: h2(:,:,:,:)
     logical :: ifh2
     
     real(kind=rp), allocatable :: jac(:,:,:,:) !< Jacobian
     real(kind=rp), allocatable :: jacinv(:,:,:,:) !< Inverted Jacobian
     real(kind=rp), allocatable :: B(:,:,:,:) !< Mass matrix/volume matrix
     real(kind=rp), allocatable :: Binv(:,:,:,:) !< Inverted Mass matrix/volume matrix

     real(kind=rp), allocatable :: area(:,:,:,:) !< Facet area
     real(kind=rp), allocatable :: nx(:,:,:,:)   !< x-direction of facet normal
     real(kind=rp), allocatable :: ny(:,:,:,:)   !< y-direction of facet normal
     real(kind=rp), allocatable :: nz(:,:,:,:)   !< z-direction of facet normal
     
     real(kind=rp) :: volume
     
     type(space_t), pointer :: Xh => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dof => null()
     type(gs_t), pointer :: gs_h=> null()
  end type coef_t

  public :: coef_init, coef_free
  
contains

  !> Initialize coefficients
  subroutine coef_init(coef, gs_h)
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout), target :: gs_h
    integer :: n

    call coef_free(coef)
    
    coef%msh => gs_h%dofmap%msh
    coef%Xh => gs_h%dofmap%Xh
    coef%dof => gs_h%dofmap
    coef%gs_h => gs_h

    !
    ! Allocate arrays for geometric data
    !
    !>@todo Be clever and try to avoid allocating zeroed geom. factors
    allocate(coef%G11(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G22(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G33(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G12(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G13(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G23(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dxdr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dxds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dxdt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dydr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dyds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dydt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%dzdr(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dzds(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dzdt(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdx(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdy(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%drdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dsdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%dtdz(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
    allocate(coef%jac(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%jacinv(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    

    allocate(coef%area(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%nx(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%ny(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    allocate(coef%nz(coef%Xh%lx, coef%Xh%ly, 6, coef%msh%nelv))
    
    call coef_generate_dxyzdrst(coef)
    
    call coef_generate_geo(coef)

    call coef_generate_area_and_normal(coef)

    !
    ! Set up multiplicity
    !
    n = coef%Xh%lx * coef%Xh%ly * coef%Xh%lz * coef%msh%nelv
    allocate(coef%mult(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    call rone(coef%mult, n)   
    call gs_op_vector(gs_h, coef%mult, n, GS_OP_ADD)
    call invcol1(coef%mult, n)
    
    allocate(coef%B(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%Binv(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    call coef_generate_mass(coef)
    
    allocate(coef%h1(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%h2(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    ! This is a placeholder, just for now
    ! We can probably find a prettier solution
    call rone(coef%h1,n)
    call rone(coef%h2,n)
    coef%ifh2 = .false.

  end subroutine coef_init

  !> Deallocate coefficients
  subroutine coef_free(coef)
    type(coef_t), intent(inout) :: coef

    if (allocated(coef%G11)) then
       deallocate(coef%G11)
    end if

    if (allocated(coef%G22)) then
       deallocate(coef%G22)
    end if

    if (allocated(coef%G33)) then
       deallocate(coef%G33)
    end if

    if (allocated(coef%G12)) then
       deallocate(coef%G12)
    end if

    if (allocated(coef%G13)) then
       deallocate(coef%G13)
    end if
    
    if (allocated(coef%G23)) then
       deallocate(coef%G23)
    end if

    if (allocated(coef%mult)) then
       deallocate(coef%mult)
    end if
    
    if (allocated(coef%B)) then
       deallocate(coef%B)
    end if
    
    if (allocated(coef%Binv)) then
       deallocate(coef%Binv)
    end if
    
    if(allocated(coef%dxdr)) then
       deallocate(coef%dxdr)
    end if
    
    if(allocated(coef%dxds)) then
       deallocate(coef%dxds)
    end if
    
    if(allocated(coef%dxdt)) then
       deallocate(coef%dxdt)
    end if
    
    if(allocated(coef%dydr)) then
       deallocate(coef%dydr)
    end if
    
    if(allocated(coef%dyds)) then
       deallocate(coef%dyds)
    end if
    
    if(allocated(coef%dydt)) then
       deallocate(coef%dydt)
    end if
    
    if(allocated(coef%dzdr)) then
       deallocate(coef%dzdr)
    end if
    
    if(allocated(coef%dzds)) then
       deallocate(coef%dzds)
    end if
    
    if(allocated(coef%dzdt)) then
       deallocate(coef%dzdt)
    end if
    
    if(allocated(coef%drdx)) then
       deallocate(coef%drdx)
    end if
    
    if(allocated(coef%dsdx)) then
       deallocate(coef%dsdx)
    end if
    
    if(allocated(coef%dtdx)) then
       deallocate(coef%dtdx)
    end if
    
    if(allocated(coef%drdy)) then
       deallocate(coef%drdy)
    end if
    
    if(allocated(coef%dsdy)) then
       deallocate(coef%dsdy)
    end if
    
    if(allocated(coef%dtdy)) then
       deallocate(coef%dtdy)
    end if
    
    if(allocated(coef%drdz)) then
       deallocate(coef%drdz)
    end if
    
    if(allocated(coef%dsdz)) then
       deallocate(coef%dsdz)
    end if
    
    if(allocated(coef%dtdz)) then
       deallocate(coef%dtdz)
    end if
    
    if(allocated(coef%jac)) then
       deallocate(coef%jac)
    end if
    
    if(allocated(coef%jacinv)) then
       deallocate(coef%jacinv)
    end if
    
    if(allocated(coef%h1)) then
       deallocate(coef%h1)
    end if
    
    if(allocated(coef%h2)) then
       deallocate(coef%h2)
    end if

    if (allocated(coef%area)) then
       deallocate(coef%area)
    end if

    if (allocated(coef%nx)) then
       deallocate(coef%nx)
    end if

    if (allocated(coef%ny)) then
       deallocate(coef%ny)
    end if

    if (allocated(coef%nz)) then
       deallocate(coef%nz)
    end if
    
    nullify(coef%msh)
    nullify(coef%Xh)
    nullify(coef%dof)
    
  end subroutine coef_free

  subroutine coef_generate_dxyzdrst(c)
    type(coef_t), intent(inout) :: c
    integer :: e,i,lxy,lyz
    
    lxy=c%Xh%lx*c%Xh%ly
    lyz=c%Xh%ly*c%Xh%lz
    do e =1,c%msh%nelv
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%x(1,1,1,e),c%Xh%lx,c%dxdr(1,1,1,e),lyz)
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%y(1,1,1,e),c%Xh%lx,c%dydr(1,1,1,e),lyz)
       call mxm(c%Xh%dx,c%Xh%lx,c%dof%z(1,1,1,e),c%Xh%lx,c%dzdr(1,1,1,e),lyz)

       DO i=1,c%Xh%lz
          call mxm(c%dof%x(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dxds(1,1,i,e),c%Xh%ly)
          call mxm(c%dof%y(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dyds(1,1,i,e),c%Xh%ly)
          call mxm(c%dof%z(1,1,i,e),c%Xh%lx,c%Xh%dyt,c%Xh%ly,c%dzds(1,1,i,e),c%Xh%ly)
       end do
       !> We actually take 2d into account, wow, need to do that for the rest.
       if(c%msh%gdim .eq. 3) then
          call mxm(c%dof%x(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dxdt(1,1,1,e),c%Xh%lz)
          call mxm(c%dof%y(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dydt(1,1,1,e),c%Xh%lz)
          call mxm(c%dof%z(1,1,1,e),lxy,c%Xh%dzt,c%Xh%lz,c%dzdt(1,1,1,e),c%Xh%lz)
       else
          call rzero(c%dxdt(1,1,1,e),lxy)
          call rzero(c%dydt(1,1,1,e),lxy)
          call rone(c%dzdt(1,1,1,e),lxy)
       endif
    end do
    
    if (c%msh%gdim .eq. 2) then
       call rzero   (c%jac,c%dof%n_dofs)
       call addcol3 (c%jac,c%dxdr,c%dyds,c%dof%n_dofs)
       call subcol3 (c%jac,c%dxds,c%dydr,c%dof%n_dofs)
       call copy    (c%drdx,c%dyds,c%dof%n_dofs)
       call copy    (c%drdy,c%dxds,c%dof%n_dofs)
       call chsign  (c%drdy,c%dof%n_dofs)
       call copy    (c%dsdx,c%dydr,c%dof%n_dofs)
       call chsign  (c%dsdx,c%dof%n_dofs)
       call copy    (c%dsdy,c%dxdr,c%dof%n_dofs)
       call rzero   (c%drdz,c%dof%n_dofs)
       call rzero   (c%dsdz,c%dof%n_dofs)
       call RONE    (c%dtdz,c%dof%n_dofs)
   else
       call rzero   (c%jac,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxdr,c%dyds,c%dzdt,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxdt,c%dydr,c%dzds,c%dof%n_dofs)
       call addcol4 (c%jac,c%dxds,c%dydt,c%dzdr,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxdr,c%dydt,c%dzds,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxds,c%dydr,c%dzdt,c%dof%n_dofs)
       call subcol4 (c%jac,c%dxdt,c%dyds,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%drdx,c%dyds,c%dzdt,c%dydt,c%dzds,c%dof%n_dofs)
       call ascol5  (c%drdy,c%dxdt,c%dzds,c%dxds,c%dzdt,c%dof%n_dofs)
       call ascol5  (c%drdz,c%dxds,c%dydt,c%dxdt,c%dyds,c%dof%n_dofs)
       call ascol5  (c%dsdx,c%dydt,c%dzdr,c%dydr,c%dzdt,c%dof%n_dofs)
       call ascol5  (c%dsdy,c%dxdr,c%dzdt,c%dxdt,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%dsdz,c%dxdt,c%dydr,c%dxdr,c%dydt,c%dof%n_dofs)
       call ascol5  (c%dtdx,c%dydr,c%dzds,c%dyds,c%dzdr,c%dof%n_dofs)
       call ascol5  (c%dtdy,c%dxds,c%dzdr,c%dxdr,c%dzds,c%dof%n_dofs)
       call ascol5  (c%dtdz,c%dxdr,c%dyds,c%dxds,c%dydr,c%dof%n_dofs)
    end if

    call invers2(c%jacinv,c%jac,c%dof%n_dofs)

  end subroutine coef_generate_dxyzdrst
  !> Generate geometric data for the given mesh
  !! @note Current implementation assumes regular shaped hex elements
  subroutine coef_generate_geo(c)
    type(coef_t), intent(inout) :: c
    integer :: e, lxyz

    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    
    if(c%msh%gdim .eq. 2) then
      call vdot2(c%G11,c%drdx,c%drdy,c%drdx,c%drdy,c%dof%n_dofs)
      call vdot2(c%G22,c%dsdx,c%dsdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call vdot2(c%G12,c%drdx,c%drdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call  col2(c%G11,c%jacinv,c%dof%n_dofs)
      call  col2(c%G22,c%jacinv,c%dof%n_dofs)
      call  col2(c%G12,c%jacinv,c%dof%n_dofs)
      call rzero(c%G33,c%dof%n_dofs)
      call rzero(c%G13,c%dof%n_dofs)
      call rzero(c%G23,c%dof%n_dofs)
    else
      call vdot3(c%G11,c%drdx,c%drdy,c%drdz,c%drdx,c%drdy,c%drdz,c%dof%n_dofs)
      call vdot3(c%G22,c%dsdx,c%dsdy,c%dsdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G33,c%dtdx,c%dtdy,c%dtdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G12,c%drdx,c%drdy,c%drdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G13,c%drdx,c%drdy,c%drdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G23,c%dsdx,c%dsdy,c%dsdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      
      call col2(c%G11,c%jacinv,c%dof%n_dofs)
      call col2(c%G22,c%jacinv,c%dof%n_dofs)
      call col2(c%G33,c%jacinv,c%dof%n_dofs)
      call col2(c%G12,c%jacinv,c%dof%n_dofs)
      call col2(c%G13,c%jacinv,c%dof%n_dofs)
      call col2(c%G23,c%jacinv,c%dof%n_dofs)
    end if
    do e=1,c%msh%nelv
       call col2(c%G11(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G22(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G12(1,1,1,e),c%Xh%w3,lxyz)
       if (c%msh%gdim .eq. 3) then
         call col2(c%G33(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G13(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G23(1,1,1,e),c%Xh%w3,lxyz)
       end if
    end do
  end subroutine coef_generate_geo
 
  !> Generate mass matrix B for the given mesh and space
  !! @note This is also a stapleholder, we need to go through the coef class properly.
  subroutine coef_generate_mass(c)
    type(coef_t), intent(inout) :: c
    integer :: e, j, k, l, lxyz
    
    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    
    call rone(c%B,c%dof%n_dofs)
    do e=1,c%msh%nelv
       ! Here we need to handle things differently for axis symmetric elements
       call col3(c%B(1,1,1,e),c%jac(1,1,1,e),c%Xh%w3,lxyz)
    end do
    
    call copy(c%Binv,c%B,c%dof%n_dofs)
    call gs_op_vector(c%gs_h,c%Binv, c%dof%n_dofs,GS_OP_ADD)
    call invcol1(c%Binv,c%dof%n_dofs)

    c%volume = glsum(c%B,c%dof%n_dofs)
  end subroutine coef_generate_mass

  !> Generate facet area and surface normals
  subroutine coef_generate_area_and_normal(coef)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), allocatable :: a(:,:,:,:)
    real(kind=rp), allocatable :: b(:,:,:,:)
    real(kind=rp), allocatable :: c(:,:,:,:)
    real(kind=rp), allocatable :: dot(:,:,:,:)
    integer :: n, e, j, k, l, lx
    real(kind=rp) :: weight, len
    n = coef%dof%n_dofs
    lx = coef%Xh%lx
    
    allocate(a(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(b(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(c(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(dot(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))

    call vcross(a,b,c, coef%dxds, coef%dyds, coef%dzds, &
         coef%dxdt, coef%dydt, coef%dzdt, n)
    call vdot3(dot, a, b, c, a, b, c, n)

    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wy(j) * coef%Xh%wz(k)
             coef%area(j, k, 2, e) = sqrt(dot(lx, j, k, e)) * weight
             coef%area(j, k, 1, e) = sqrt(dot(1, j, k, e)) * weight
             coef%nx(j,k, 1, e) = -A(1, j, k, e)
             coef%nx(j,k, 2, e) =  A(lx, j, k, e)
             coef%ny(j,k, 1, e) = -B(1, j, k, e)
             coef%ny(j,k, 2, e) =  B(lx, j, k, e)
             coef%nz(j,k, 1, e) = -C(1, j, k, e)
             coef%nz(j,k, 2, e) =  C(lx, j, k, e)
          end do
       end do
    end do

    call vcross(a,b,c, coef%dxdr, coef%dydr, coef%dzdr, &
         coef%dxdt, coef%dydt, coef%dzdt, n)
    call vdot3(dot, a, b, c, a, b, c, n)
    
    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wx(j) * coef%Xh%wz(k)
             coef%area(j, k, 3, e) = sqrt(dot(j, 1, k, e)) * weight
             coef%area(j, k, 4, e) = sqrt(dot(j, lx, k, e)) * weight
             coef%nx(j,k, 3, e) =  A(j, 1, k, e)
             coef%nx(j,k, 4, e) = -A(j, lx, k, e)
             coef%ny(j,k, 3, e) =  B(j, 1, k, e)
             coef%ny(j,k, 4, e) = -B(j, lx, k, e)
             coef%nz(j,k, 3, e) =  C(j, 1, k, e)
             coef%nz(j,k, 4, e) = -C(j, lx, k, e)             
          end do
       end do
    end do


    call vcross(a,b,c, coef%dxdr, coef%dydr, coef%dzdr, &
         coef%dxds, coef%dyds, coef%dzds, n)
    call vdot3(dot, a, b, c, a, b, c, n)
    
    do e = 1, coef%msh%nelv
       do k = 1, coef%Xh%lx
          do j = 1, coef%Xh%lx
             weight = coef%Xh%wx(j) * coef%Xh%wy(k)
             coef%area(j, k, 5, e) = sqrt(dot(j, k, 1, e)) * weight
             coef%area(j, k, 6, e) = sqrt(dot(j, j, lx, e)) * weight
             coef%nx(j,k, 5, e) = -A(j, k, 1, e)
             coef%nx(j,k, 6, e) =  A(j, k, lx, e)
             coef%ny(j,k, 5, e) = -B(j, k, 1, e)
             coef%ny(j,k, 6, e) =  B(j, k, lx, e)
             coef%nz(j,k, 5, e) = -C(j, k, 1, e)
             coef%nz(j,k, 6, e) =  C(j, k, lx, e)             
          end do
       end do
    end do

    ! Normalize
    n = size(coef%nz)
    do j = 1, n
       len = sqrt(coef%nx(j,1,1,1)**2 + &
            coef%ny(j,1,1,1)**2 + coef%nz(j,1,1,1)**2)
       if (len .gt. NEKO_EPS) then
          coef%nx(j,1,1,1) = coef%nx(j,1,1,1) / len
          coef%ny(j,1,1,1) = coef%ny(j,1,1,1) / len
          coef%nz(j,1,1,1) = coef%nz(j,1,1,1) / len
       end if
    end do

    deallocate(dot)
    deallocate(c)
    deallocate(b)
    deallocate(a)
    
  end subroutine coef_generate_area_and_normal
  
end module coefs
