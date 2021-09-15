!> Coefficients 
module coefs
  use gather_scatter
  use num_types
  use space  
  use math
  use mesh
  use mamba
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

     type(mmbArray) :: mba_g11
     type(mmbArray) :: mba_g22
     type(mmbArray) :: mba_g33
     type(mmbArray) :: mba_g12
     type(mmbArray) :: mba_g13
     type(mmbArray) :: mba_g23

     type(mmbLayout) :: layout

     type(mmbTileIterator) :: mba_g11_it
     type(mmbTileIterator) :: mba_g22_it
     type(mmbTileIterator) :: mba_g33_it
     type(mmbTileIterator) :: mba_g12_it
     type(mmbTileIterator) :: mba_g13_it
     type(mmbTileIterator) :: mba_g23_it

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
    
    associate(G11 => c%G11, G12 => c%G12, G13 => c%G13, &
         G22 => c%G22, G23 => c%G23, G33 => c%G33, &
         drdx => c%drdx, drdy => c%drdy, drdz => c%drdz, &
         dsdx => c%dsdx, dsdy => c%dsdy, dsdz => c%dsdz, &
         dtdx => c%dtdx, dtdy => c%dtdy, dtdz => c%dtdz, &
         dxdr => c%dxdr, dydr => c%dydr, dzdr => c%dzdr, &
         dxds => c%dxds, dyds => c%dyds, dzds => c%dzds, &
         dxdt => c%dxdt, dydt => c%dydt, dzdt => c%dzdt, &
         dx => c%Xh%dx, dy => c%Xh%dy, dz => c%Xh%dz, &
         x => c%dof%x, y => c%dof%y, z => c%dof%z, &
         lx => c%Xh%lx, ly => c%Xh%ly, lz => c%Xh%lz, &
         dyt => c%Xh%dyt, dzt => c%Xh%dzt, &
         jacinv => c%jacinv, jac => c%jac, n_dofs => c%dof%n_dofs)

      do e = 1, c%msh%nelv
         call mxm(dx, lx, x(1,1,1,e), lx, dxdr(1,1,1,e), lyz)
         call mxm(dx, lx, y(1,1,1,e), lx, dydr(1,1,1,e), lyz)
         call mxm(dx, lx, z(1,1,1,e), lx, dzdr(1,1,1,e), lyz)
         
         do i = 1, lz
            call mxm(x(1,1,i,e), lx, dyt, ly, dxds(1,1,i,e), ly)
            call mxm(y(1,1,i,e), lx, dyt, ly, dyds(1,1,i,e), ly)
            call mxm(z(1,1,i,e), lx, dyt, ly, dzds(1,1,i,e), ly)
         end do
       
         ! We actually take 2d into account, wow, need to do that for the rest.
         if(c%msh%gdim .eq. 3) then
            call mxm(x(1,1,1,e), lxy, dzt, lz, dxdt(1,1,1,e), lz)
            call mxm(y(1,1,1,e), lxy, dzt, lz, dydt(1,1,1,e), lz)
            call mxm(z(1,1,1,e), lxy, dzt, lz, dzdt(1,1,1,e), lz)
         else
            call rzero(dxdt(1,1,1,e), lxy)
            call rzero(dydt(1,1,1,e), lxy)
            call rone(dzdt(1,1,1,e), lxy)
         end if
      end do
      
      if (c%msh%gdim .eq. 2) then
         call rzero   (jac, n_dofs)
         call addcol3 (jac, dxdr, dyds, n_dofs)
         call subcol3 (jac, dxds, dydr, n_dofs)
         call copy    (drdx, dyds, n_dofs)
         call copy    (drdy, dxds, n_dofs)
         call chsign  (drdy, n_dofs)
         call copy    (dsdx, dydr, n_dofs)
         call chsign  (dsdx, n_dofs)
         call copy    (dsdy, dxdr, n_dofs)
         call rzero   (drdz, n_dofs)
         call rzero   (dsdz, n_dofs)
         call rone    (dtdz, n_dofs)
      else
         call rzero   (jac, n_dofs)
         call addcol4 (jac, dxdr, dyds, dzdt, n_dofs)
         call addcol4 (jac, dxdt, dydr, dzds, n_dofs)
         call addcol4 (jac, dxds, dydt, dzdr, n_dofs)
         call subcol4 (jac, dxdr, dydt, dzds, n_dofs)
         call subcol4 (jac, dxds, dydr, dzdt, n_dofs)
         call subcol4 (jac, dxdt, dyds, dzdr, n_dofs)
         call ascol5  (drdx, dyds, dzdt, dydt, dzds, n_dofs)
         call ascol5  (drdy, dxdt, dzds, dxds, dzdt, n_dofs)
         call ascol5  (drdz, dxds, dydt, dxdt, dyds, n_dofs)
         call ascol5  (dsdx, dydt, dzdr, dydr, dzdt, n_dofs)
         call ascol5  (dsdy, dxdr, dzdt, dxdt, dzdr, n_dofs)
         call ascol5  (dsdz, dxdt, dydr, dxdr, dydt, n_dofs)
         call ascol5  (dtdx, dydr, dzds, dyds, dzdr, n_dofs)
         call ascol5  (dtdy, dxds, dzdr, dxdr, dzds, n_dofs)
         call ascol5  (dtdz, dxdr, dyds, dxds, dydr, n_dofs)
      end if
      
      call invers2(jacinv, jac, n_dofs)
      
    end associate
    
  end subroutine coef_generate_dxyzdrst
  
  !> Generate geometric data for the given mesh
  !! @note Current implementation assumes regular shaped hex elements
  subroutine coef_generate_geo(c)
    type(coef_t), intent(inout) :: c
    integer :: e, lxyz
    integer(mmbErrorKind) :: err
    integer(mmbIndexKind), dimension(4) :: dims

    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    associate(G11 => c%G11, G12 => c%G12, G13 => c%G13, &
         G22 => c%G22, G23 => c%G23, G33 => c%G33, &
         drdx => c%drdx, drdy => c%drdy, drdz => c%drdz, &
         dsdx => c%dsdx, dsdy => c%dsdy, dsdz => c%dsdz, &
         dtdx => c%dtdx, dtdy => c%dtdy, dtdz => c%dtdz, &
         jacinv => c%jacinv, n_dofs => c%dof%n_dofs, w3 => c%Xh%w3)
    
      if(c%msh%gdim .eq. 2) then
         call vdot2(G11, drdx, drdy, drdx, drdy, n_dofs)
         call vdot2(G22, dsdx, dsdy, dsdx, dsdy, n_dofs)
         call vdot2(G12, drdx, drdy, dsdx, dsdy, n_dofs)
         call  col2(G11, jacinv, n_dofs)
         call  col2(G22, jacinv, n_dofs)
         call  col2(G12, jacinv, n_dofs)
         call rzero(G33, n_dofs)
         call rzero(G13, n_dofs)
         call rzero(G23, n_dofs)
      else
         call vdot3(G11, drdx, drdy, drdz, drdx, drdy, drdz, n_dofs)
         call vdot3(G22, dsdx, dsdy, dsdz, dsdx, dsdy, dsdz, n_dofs)
         call vdot3(G33, dtdx, dtdy, dtdz, dtdx, dtdy, dtdz, n_dofs)
         call vdot3(G12, drdx, drdy, drdz, dsdx, dsdy, dsdz, n_dofs)
         call vdot3(G13, drdx, drdy, drdz, dtdx, dtdy, dtdz, n_dofs)
         call vdot3(G23, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, n_dofs)
         
         call col2(G11, jacinv, n_dofs)
         call col2(G22, jacinv, n_dofs)
         call col2(G33, jacinv, n_dofs)
         call col2(G12, jacinv, n_dofs)
         call col2(G13, jacinv, n_dofs)
         call col2(G23, jacinv, n_dofs)
      end if
      do e = 1, c%msh%nelv
         call col2(G11(1,1,1,e), w3, lxyz)
         call col2(G22(1,1,1,e), w3, lxyz)
         call col2(G12(1,1,1,e), w3, lxyz)
         if (c%msh%gdim .eq. 3) then
            call col2(G33(1,1,1,e), w3, lxyz)
            call col2(G13(1,1,1,e), w3, lxyz)
            call col2(G23(1,1,1,e), w3, lxyz)
         end if
      end do

      dims = [c%Xh%lx, c%Xh%lx, c%Xh%lx, c%msh%nelv]
      call mmb_layout_create_regular_nd(int(storage_size(1.0_rp)/8, mmbSizeKind), &
           4_mmbSizeKind, MMB_COLMAJOR, mmb_layout_padding_create_zero(),&
           c%layout, err)

      call mmb_array_create_wrapped(G11, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g11, err)
      call mmb_array_create_wrapped(G22, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g22, err)
      call mmb_array_create_wrapped(G33, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g33, err)
      call mmb_array_create_wrapped(G12, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g12, err)
      call mmb_array_create_wrapped(G13, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g13, err)
      call mmb_array_create_wrapped(G23, dims, c%layout, &
         dram_interface, MMB_READ_WRITE, c%mba_g23, err)

      call mmb_array_tile(c%mba_g11, dims, err)
      call mmb_array_tile(c%mba_g22, dims, err)
      call mmb_array_tile(c%mba_g33, dims, err)
      call mmb_array_tile(c%mba_g12, dims, err)
      call mmb_array_tile(c%mba_g13, dims, err)
      call mmb_array_tile(c%mba_g23, dims, err)

      call mmb_tile_iterator_create(c%mba_g11, c%mba_g11_it, err)
      call mmb_tile_iterator_create(c%mba_g22, c%mba_g22_it, err)
      call mmb_tile_iterator_create(c%mba_g33, c%mba_g33_it, err)
      call mmb_tile_iterator_create(c%mba_g12, c%mba_g12_it, err)
      call mmb_tile_iterator_create(c%mba_g13, c%mba_g13_it, err)
      call mmb_tile_iterator_create(c%mba_g23, c%mba_g23_it, err)

    end associate
  end subroutine coef_generate_geo
 
  !> Generate mass matrix B for the given mesh and space
  !! @note This is also a stapleholder, we need to go through the coef class properly.
  subroutine coef_generate_mass(c)
    type(coef_t), intent(inout) :: c
    integer :: e, j, k, l, lxyz
    
    lxyz = c%Xh%lx * c%Xh%ly * c%Xh%lz
    
    call rone(c%B,c%dof%n_dofs)
    do e = 1, c%msh%nelv
       ! Here we need to handle things differently for axis symmetric elements
       call col3(c%B(1,1,1,e), c%jac(1,1,1,e), c%Xh%w3, lxyz)
    end do
    
    call copy(c%Binv, c%B, c%dof%n_dofs)
    call gs_op_vector(c%gs_h, c%Binv, c%dof%n_dofs, GS_OP_ADD)
    call invcol1(c%Binv, c%dof%n_dofs)

    c%volume = glsum(c%B, c%dof%n_dofs)
    
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
