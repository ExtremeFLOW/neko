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
     
     real(kind=dp), allocatable :: G1(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G2(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G3(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G4(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G5(:,:,:,:) !< Geometric data
     real(kind=dp), allocatable :: G6(:,:,:,:) !< Geometric data

     real(kind=dp), allocatable :: mult(:,:,:,:) !< Multiplicity
     ! generate mapping data between element and reference element 
     !! \f$ dx/dr, dy/dr, dz/dr \f$
     !! \f$ dx/ds, dy/ds, dz/ds \f$
     !! \f$ dx/dt, dy/dt, dz/dt \f$
     real(kind=dp), allocatable :: dxdr(:,:,:,:), dydr(:,:,:,:), dzdr(:,:,:,:) 
     real(kind=dp), allocatable :: dxds(:,:,:,:), dyds(:,:,:,:), dzds(:,:,:,:)
     real(kind=dp), allocatable :: dxdt(:,:,:,:), dydt(:,:,:,:), dzdt(:,:,:,:) 
     !< \f$ dr/dx, dr/dy, dr/dz \f$
     !! \f$ ds/dx, ds/dy, ds/dz \f$
     !! \f$ dt/dx, dt/dy, dt/dz \f$
     real(kind=dp), allocatable :: drdx(:,:,:,:), drdy(:,:,:,:), drdz(:,:,:,:) 
     real(kind=dp), allocatable :: dsdx(:,:,:,:), dsdy(:,:,:,:), dsdz(:,:,:,:)
     real(kind=dp), allocatable :: dtdx(:,:,:,:), dtdy(:,:,:,:), dtdz(:,:,:,:) 
     
     real(kind=dp), allocatable :: h1(:,:,:,:) 
     real(kind=dp), allocatable :: h2(:,:,:,:)
     logical :: ifh2
     
     real(kind=dp), allocatable :: jac(:,:,:,:) !< Jacobian
     real(kind=dp), allocatable :: jacinv(:,:,:,:) !< Inverted Jacobian
     real(kind=dp), allocatable :: B(:,:,:,:) !< Mass matrix/volume matrix
     real(kind=dp), allocatable :: Binv(:,:,:,:) !< Inverted Mass matrix/volume matrix
     
     real(kind=dp) :: volume
     
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
    allocate(coef%G1(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G2(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G3(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G4(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G5(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    allocate(coef%G6(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    
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
    
    
    call coef_generate_dxyzdrst(coef)
    
    call coef_generate_geo(coef)

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

    if (allocated(coef%G1)) then
       deallocate(coef%G1)
    end if

    if (allocated(coef%G2)) then
       deallocate(coef%G2)
    end if

    if (allocated(coef%G3)) then
       deallocate(coef%G3)
    end if

    if (allocated(coef%G4)) then
       deallocate(coef%G4)
    end if

    if (allocated(coef%G5)) then
       deallocate(coef%G5)
    end if
    
    if (allocated(coef%G6)) then
       deallocate(coef%G6)
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
    if(allocated(coef%dxdr)) deallocate(coef%dxdr)
    if(allocated(coef%dxds)) deallocate(coef%dxds)
    if(allocated(coef%dxdt)) deallocate(coef%dxdt)
    if(allocated(coef%dydr)) deallocate(coef%dydr)
    if(allocated(coef%dyds)) deallocate(coef%dyds)
    if(allocated(coef%dydt)) deallocate(coef%dydt)
    if(allocated(coef%dzdr)) deallocate(coef%dzdr)
    if(allocated(coef%dzds)) deallocate(coef%dzds)
    if(allocated(coef%dzdt)) deallocate(coef%dzdt)
    if(allocated(coef%drdx)) deallocate(coef%drdx)
    if(allocated(coef%dsdx)) deallocate(coef%dsdx)
    if(allocated(coef%dtdx)) deallocate(coef%dtdx)
    if(allocated(coef%drdy)) deallocate(coef%drdy)
    if(allocated(coef%dsdy)) deallocate(coef%dsdy)
    if(allocated(coef%dtdy)) deallocate(coef%dtdy)
    if(allocated(coef%drdz)) deallocate(coef%drdz)
    if(allocated(coef%dsdz)) deallocate(coef%dsdz)
    if(allocated(coef%dtdz)) deallocate(coef%dtdz)
    
    if(allocated(coef%jac)) deallocate(coef%jac)
    if(allocated(coef%jacinv)) deallocate(coef%jacinv)
    
    if(allocated(coef%h1)) deallocate(coef%h1)
    if(allocated(coef%h2)) deallocate(coef%h2)
    

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
      call vdot2(c%G1,c%drdx,c%drdy,c%drdx,c%drdy,c%dof%n_dofs)
      call vdot2(c%G2,c%dsdx,c%dsdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call vdot2(c%G4,c%drdx,c%drdy,c%dsdx,c%dsdy,c%dof%n_dofs)
      call  col2(c%G1,c%jacinv,c%dof%n_dofs)
      call  col2(c%G2,c%jacinv,c%dof%n_dofs)
      call  col2(c%G4,c%jacinv,c%dof%n_dofs)
      call rzero(c%G3,c%dof%n_dofs)
      call rzero(c%G5,c%dof%n_dofs)
      call rzero(c%G6,c%dof%n_dofs)
    else
      call vdot3(c%G1,c%drdx,c%drdy,c%drdz,c%drdx,c%drdy,c%drdz,c%dof%n_dofs)
      call vdot3(c%G2,c%dsdx,c%dsdy,c%dsdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G3,c%dtdx,c%dtdy,c%dtdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G4,c%drdx,c%drdy,c%drdz,c%dsdx,c%dsdy,c%dsdz,c%dof%n_dofs)
      call vdot3(c%G5,c%drdx,c%drdy,c%drdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      call vdot3(c%G6,c%dsdx,c%dsdy,c%dsdz,c%dtdx,c%dtdy,c%dtdz,c%dof%n_dofs)
      
      call col2(c%G1,c%jacinv,c%dof%n_dofs)
      call col2(c%G2,c%jacinv,c%dof%n_dofs)
      call col2(c%G3,c%jacinv,c%dof%n_dofs)
      call col2(c%G4,c%jacinv,c%dof%n_dofs)
      call col2(c%G5,c%jacinv,c%dof%n_dofs)
      call col2(c%G6,c%jacinv,c%dof%n_dofs)
    end if
    do e=1,c%msh%nelv
       call col2(c%G1(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G2(1,1,1,e),c%Xh%w3,lxyz)
       call col2(c%G4(1,1,1,e),c%Xh%w3,lxyz)
       if (c%msh%gdim .eq. 3) then
         call col2(c%G3(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G5(1,1,1,e),c%Xh%w3,lxyz)
         call col2(c%G6(1,1,1,e),c%Xh%w3,lxyz)
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
  
end module coefs
